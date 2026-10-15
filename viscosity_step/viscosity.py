# -*- coding: utf-8 -*-

"""Non-graphical part of the Viscosity step in a SEAMM flowchart

The shear viscosity from the fluctuations of the pressure tensor in equilibrium
MD, by Green-Kubo and from the Helfand moments. The five independent shear
components of the pressure tensor (Pxy, Pxz, Pyz, (Pxx-Pyy)/2 and (Pyy-Pzz)/2)
are analyzed separately and averaged.
"""

import importlib.resources
import logging
from pathlib import Path
import sys
import traceback

import numpy as np
from seamm_md_analysis import (
    SHEAR_COMPONENTS,
    autocorrelation_functions,
    component_helfand_moments,
    fit_green_kubo_integral,
    fit_helfand_moment,
    open_md,
    shear_stress_components,
)
from seamm_md_analysis.plots import plot_components
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    md_inputs,
    publish_viscosity,
    read_files,
    request_sampling,
    run_own_md,
)
import seamm_md_analysis.stats
from tabulate import tabulate

import viscosity_step
import molsystem
import seamm
import seamm_util
from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Viscosity")

# Add this module's properties to the standard properties
path = importlib.resources.files("viscosity_step") / "data"
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)

#: The longest acceptable interval between samples of the pressure tensor
MAX_INTERVAL = Q_(5.0, "fs")

#: The fewest samples in each block when a single run is divided into blocks
MIN_BLOCK_SAMPLES = 2000

#: The labels of the components and their average, for tables and results
COMPONENTS = (*SHEAR_COMPONENTS, "average")


def fmt_err(value, err, precision=2):
    """Format a value and its error, with 3 decimals when the error is unusable."""
    return seamm_md_analysis.stats.fmt_err(
        value, err, precision=precision, fallback_decimals=3
    )


class Viscosity(seamm.Node):
    """
    The non-graphical part of a Viscosity step in a flowchart.

    Attributes
    ----------
    subflowchart : seamm.Flowchart
        The MD to run, when this step runs its own MD.

    parameters : ViscosityParameters
        The control parameters for Viscosity.

    See Also
    --------
    TkViscosity,
    Viscosity, ViscosityParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Viscosity",
        namespace="org.molssi.seamm",
        extension=None,
        logger=logger,
    ):
        """A step for Viscosity in a SEAMM flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.
        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Viscosity {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Viscosity", namespace=namespace
        )

        super().__init__(
            flowchart=flowchart,
            title="Viscosity",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = viscosity_step.metadata
        self.parameters = viscosity_step.ViscosityParameters()
        self._file_handler = None
        self._reset()

    @property
    def version(self):
        """The semantic version of this module."""
        return viscosity_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return viscosity_step.__git_revision__

    def _reset(self):
        """Clear the data of earlier runs, e.g. from an earlier loop iteration."""
        self.T = []  # The mean temperature of each run, K
        self.V = []  # The mean volume of each run, Å^3
        self.timestep = []  # The interval between samples of each run, fs
        self.acfs = []  # The autocorrelation functions per run, 5 x m, MPa^2
        self.integrals = []  # The Green-Kubo running integrals per run, 5 x m
        self.Ms = []  # The Helfand moments per run, 5 x m, mPa*s*ps

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "viscosity-step"
        parser = seamm_util.getParser()

        self.subflowchart._parser = self.flowchart._parser

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        super().create_parser(name=parser_name)

        if not parser_exists:
            parser.add_argument(
                parser_name,
                "--html",
                action="store_true",
                help="whether to write out html files for graphs, etc.",
            )

        # Now need to walk through the steps in the subflowchart...
        self.subflowchart.reset_visited()
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            node = node.create_parser()

        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        short : bool = False
            Leave out the description of the sub-flowchart.

        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        source = P.get("md source", run_own_md)
        if source == run_own_md:
            text = (
                "Calculate the shear viscosity by Green-Kubo and from the Helfand "
                f"moments, averaging {P['nruns']} runs of the MD below.\n\n"
            )
        else:
            text = (
                "Calculate the shear viscosity by Green-Kubo and from the Helfand "
                "moments of "
            )
            if source == read_files:
                text += f"the MD in {P['md files']}.\n\n"
            else:
                text += f"{source} earlier in the flowchart.\n\n"
            short = True

        if not short:
            self.subflowchart.root_directory = self.flowchart.root_directory

            node = self.subflowchart.get_node("1").next()
            while node is not None:
                try:
                    text += __(node.description_text()).__str__()
                except Exception as e:
                    print(f"Error describing viscosity flowchart: {e} in {node}")
                    self.logger.critical(
                        f"Error describing viscosity flowchart: {e} in {node}"
                    )
                    raise
                except:  # noqa: E722
                    print(
                        "Unexpected error describing viscosity flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    self.logger.critical(
                        "Unexpected error describing viscosity flowchart: "
                        f"{sys.exc_info()[0]} in {str(node)}"
                    )
                    raise
                text += "\n"
                node = node.next()

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def run(self):
        """Run a Viscosity step.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        next_node = super().run(printer)
        self._reset()

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Print what we are doing
        printer.important(__(self.description_text(P, short=True), indent=self.indent))

        own_md = P["md source"] == run_own_md
        _, configuration = self.get_system_configuration()
        if own_md and configuration is not None and configuration.periodicity != 3:
            raise RuntimeError(
                "The Viscosity step needs a periodic system: the pressure tensor "
                "is divided by the volume of the periodic cell."
            )

        # Cite the methods
        self.references.cite(
            raw=self._bibliography["Daivis1994"],
            alias="Daivis1994",
            module="viscosity_step",
            level=1,
            note="Averaging the five shear components of the pressure tensor.",
        )
        self.references.cite(
            raw=self._bibliography["Helfand1960"],
            alias="Helfand1960",
            module="viscosity_step",
            level=2,
            note="The Helfand moments.",
        )

        # Find the handler for job.out and set the level up
        job_handler = None
        out_handler = None
        job_level = None
        out_level = None
        for handler in job.handlers:
            if (
                isinstance(handler, logging.FileHandler)
                and "job.out" in handler.baseFilename
            ):
                job_handler = handler
                job_level = job_handler.level
                job_handler.setLevel(printing.JOB)
            elif isinstance(handler, logging.StreamHandler):
                out_handler = handler
                out_level = out_handler.level
                out_handler.setLevel(printing.JOB)
        handlers = (job_handler, job_level, out_handler, out_level)

        if own_md:
            self._run_own_md(P, handlers)
        else:
            self._analyze_earlier_md(P, handlers)

        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)

        if len(self.Ms) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")

        self.analyze(P=P, style="full", configuration=configuration)

        return next_node

    def _run_own_md(self, P, handlers):
        """Run the MD in the sub-flowchart several times and process each run."""
        job_handler, job_level, out_handler, out_level = handlers

        # Turn on any sampling the analysis needs but the MD lacks
        changes, missing = request_sampling(self.subflowchart, self.md_requirements())
        for change in changes:
            printer.job(__(f"Note: {change}.", indent=self.indent + 4 * " "))
        if len(missing) > 0:
            raise RuntimeError(
                "The MD in the sub-flowchart cannot provide what is needed: "
                + " ".join(missing)
            )

        first_node = self.subflowchart.get_node("1").next()
        if first_node is None:
            raise RuntimeError("The sub-flowchart of the Viscosity step is empty.")

        # Ensure the nodes have their options
        node = first_node
        while node is not None:
            node.all_options = self.all_options
            node = node.next()

        # And the subflowchart has the executor
        self.subflowchart.executor = self.flowchart.executor

        nruns = P["nruns"]
        fmt = f"0{len(str(nruns))}"
        for run in range(1, nruns + 1):
            # Direct most output to Run.out
            run_dir = Path(self.directory) / f"run_{run:{fmt}}"
            run_dir.mkdir(parents=True, exist_ok=True)

            if self._file_handler is not None:
                self._file_handler.close()
                job.removeHandler(self._file_handler)
            path = run_dir / "Run.out"
            path.unlink(missing_ok=True)
            self._file_handler = logging.FileHandler(path)
            self._file_handler.setLevel(printing.NORMAL)
            formatter = logging.Formatter(fmt="{message:s}", style="{")
            self._file_handler.setFormatter(formatter)
            job.addHandler(self._file_handler)

            # Add the run to the ids so the directory structure is reasonable
            self.flowchart.reset_visited()
            self.set_subids((*self._id, f"run_{run:{fmt}}"))

            node = first_node
            try:
                while node is not None:
                    node = node.run()
            except DeprecationWarning as e:
                printer.normal("\nDeprecation warning: " + str(e))
                traceback.print_exc(file=sys.stderr)
                traceback.print_exc(file=sys.stdout)
            except Exception as e:
                printer.job(f"Caught exception in run {run}: {str(e)}")
                with open(run_dir / "stderr.out", "a") as fd:
                    traceback.print_exc(file=fd)
                if "continue" in P["errors"]:
                    continue
                elif "exit" in P["errors"]:
                    break
                else:
                    raise
            else:
                self.process_run(self._read_run(run, run_dir))
                self._print_progress(run, handlers)

            self.logger.debug(f"End of run {run}")

        # Remove any redirection of printing.
        if self._file_handler is not None:
            self._file_handler.close()
            job.removeHandler(self._file_handler)
            self._file_handler = None

    def _analyze_earlier_md(self, P, handlers):
        """Analyze MD run earlier in the flowchart, or in other jobs."""
        paths = md_inputs(self, P, name="pressure_tensor")
        runs = [self._read_md_file(path) for path in paths]

        text = f"Analyzing {len(paths)} MD run{'s' if len(paths) > 1 else ''}"
        what = "runs"
        if len(runs) == 1:
            what = "blocks"
            runs = self._split_into_blocks(runs[0], P["number of blocks"])
            text += f", divided into {len(runs)} blocks to estimate the errors"
        text += ": " + ", ".join(str(p) for p in paths)
        printer.job(__(text, indent=self.indent + 4 * " "))
        printer.job("")

        for run, data in enumerate(runs, start=1):
            self.process_run(data)
            self._print_progress(run, handlers, what=what)

    def _print_progress(self, run, handlers, what="runs"):
        """Print the running estimate of the viscosity after each run or block."""
        job_handler, job_level, out_handler, out_level = handlers
        if len(self.Ms) < 2:
            return
        if job_handler is not None:
            job_handler.setLevel(job_level)
        if out_handler is not None:
            out_handler.setLevel(out_level)
        try:
            results = self._fit(fit_components=False)
            gk = results["Green-Kubo"][-1]
            v, e = fmt_err(gk["value"], self._err95(gk))
            text = f"After {run} {what}: viscosity = {v} ± {e} mPa*s (Green-Kubo, 95%)"
        except Exception as e:
            text = f"After {run} {what} the viscosity could not be fitted: {e}"
        printer.normal(__(text, indent=self.indent + 4 * " "))
        if job_handler is not None:
            job_handler.setLevel(printing.JOB)
        if out_handler is not None:
            out_handler.setLevel(printing.JOB)

    def process_run(self, data):
        """Compute the correlation functions and Helfand moments of one run.

        Parameters
        ----------
        data : (numpy.ndarray(n, 3, 3), pint.Quantity, float, float)
            The pressure tensor in Pa, the interval between samples, the mean
            temperature in K and the mean volume in Å^3.
        """
        tensor, dt, T, V = data
        self.T.append(T)
        self.V.append(V)
        self.timestep.append(dt.m_as("fs"))

        X = shear_stress_components(tensor)
        n = X.shape[1]
        m = min(n // 10, 10000)

        T = Q_(T, "K")
        V = Q_(V, "Å^3")
        k_B = Q_("k_B")
        Psq = Q_(1.0, "Pa^2")

        # The shear stress of a fluid has a mean of zero, so it is not subtracted:
        # the sample mean of a finite run would bias the viscosity low.
        acf, integral = autocorrelation_functions(X, m=m, demean=False)
        self.acfs.append(acf * Psq.m_as("MPa^2"))
        self.integrals.append(integral * (V / (k_B * T) * dt * Psq).m_as("mPa*s"))

        constants = V * dt**2 / (2 * k_B * T) * Psq
        self.Ms.append(
            component_helfand_moments(X, m=m, demean=False) * constants.m_as("mPa*s*ps")
        )

    @staticmethod
    def _with_average(rows):
        """Add the average of the components as a last row."""
        rows = np.asarray(rows)
        return np.vstack([rows, rows.mean(axis=0)])

    @staticmethod
    def _err95(fit):
        stderr = fit["stderr_runs"]
        return None if stderr is None else 2 * stderr

    def _curves(self, per_run):
        """The mean over the runs of each component and their average, and the
        standard deviation over the runs (the error bands of the plots)."""
        rows = np.stack([self._with_average(r) for r in per_run])
        return rows.mean(axis=0), rows.std(axis=0)

    def _fit(self, fit_components=True):
        """Fit the Green-Kubo integrals and Helfand moments.

        The value of each component and of their average is fitted to the mean
        over the runs. Its error is the standard error of the values fitted to
        each run separately: the errors of the fit itself treat strongly correlated
        points as independent, so they are far too small.

        Returns
        -------
        {str: [dict]}
            For "Green-Kubo" and "Helfand", a dict for each component and their
            average with the value, the errors and the fitted curve. A component
            that could not be fitted has None.
        """
        dt = self.timestep[0] * 1.0e-3  # ps
        integral, integral_err = self._curves(self.integrals)
        M, M_err = self._curves(self.Ms)
        ts = (np.arange(integral.shape[1]) * dt).tolist()
        indices = range(len(COMPONENTS)) if fit_components else [len(COMPONENTS) - 1]

        def over_runs(fit_one, per_run, i):
            values = []
            for rows in per_run:
                try:
                    values.append(fit_one(self._with_average(rows)[i]))
                except Exception:
                    pass
            values = np.array([v for v in values if np.isfinite(v)])
            if values.size < 2:
                return None
            return float(np.std(values, ddof=1) / np.sqrt(values.size))

        results = {"Green-Kubo": [], "Helfand": []}
        for i in indices:
            try:
                fit = fit_green_kubo_integral(integral[i], ts, sigma=integral_err[i])
                gk = {
                    "value": fit.value,
                    "stderr": fit.stderr,
                    "stderr_runs": over_runs(
                        lambda y: fit_green_kubo_integral(y, ts).value,
                        self.integrals,
                        i,
                    ),
                    "tau": fit.tau,
                    "xs": fit.xs,
                    "ys": fit.ys,
                }
            except Exception as e:
                logger.warning(f"The Green-Kubo fit of {COMPONENTS[i]} failed: {e}")
                gk = None
            results["Green-Kubo"].append(gk)

            # The Helfand moment is linear once the correlations have decayed
            start = 1.0
            if gk is not None:
                start = max(start, 3 * max(gk["tau"]))
            try:
                fit = fit_helfand_moment(M[i], ts, sigma=M_err[i], start=start)
                helfand = {
                    "value": fit.slope,
                    "stderr": fit.stderr,
                    "stderr_runs": over_runs(
                        lambda y: fit_helfand_moment(y, ts, start=start).slope,
                        self.Ms,
                        i,
                    ),
                    "start": start,
                    "xs": fit.xs,
                    "ys": fit.ys,
                }
            except Exception as e:
                logger.warning(f"The Helfand fit of {COMPONENTS[i]} failed: {e}")
                helfand = None
            results["Helfand"].append(helfand)
        return results

    def analyze(self, indent="", P=None, style="full", configuration=None, **kwargs):
        """Fit, print, plot and store the viscosity.

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        P : dict
            The values of the parameters.
        configuration : molsystem._Configuration
            The configuration, for storing the results as properties.
        """
        results = self._fit()
        nruns = len(self.Ms)
        T = float(np.mean(self.T))
        V = float(np.mean(self.V))

        table = {"Method": [], "Component": [], "η": [], "±": [], "95%": []}
        for method in ("Green-Kubo", "Helfand"):
            for i, fit in enumerate(results[method]):
                table["Method"].append(method if i == 0 else "")
                table["Component"].append(COMPONENTS[i])
                if fit is None:
                    table["η"].append("failed")
                    table["±"].append("")
                    table["95%"].append("")
                else:
                    v, e = fmt_err(fit["value"], self._err95(fit))
                    table["η"].append(v)
                    table["±"].append("±" if e != "" else "")
                    table["95%"].append(e)

        text = tabulate(
            table,
            headers="keys",
            tablefmt="simple_outline",
            disable_numparse=True,
            colalign=("center", "center", "decimal", "center", "decimal"),
        )
        length = len(text.splitlines()[0])
        header = "Shear viscosity (mPa*s)".center(length)
        printer.normal(
            __("\n" + header + "\n" + text + "\n", indent=8 * " ", wrap=False)
        )

        notes = []
        if nruns < 2:
            notes.append(
                "With a single run there is no error estimate. Use several runs, or "
                "divide the run into blocks."
            )
        if results["Helfand"][-1] is None:
            notes.append(
                "The Helfand moments could not be fitted: the fit starts once the "
                "correlations have decayed, but the moments extend only to a tenth "
                "of the length of each run (or block). Use longer runs, or fewer "
                "blocks."
            )
        dt = Q_(self.timestep[0], "fs")
        if dt > MAX_INTERVAL:
            notes.append(
                f"The pressure tensor was sampled every {dt:~P}. It should be "
                f"sampled at least every {MAX_INTERVAL:~P} to follow its fastest "
                "fluctuations."
            )
        notes.append(f"The temperature was {T:.2f} K and the volume {V:.1f} Å^3.")
        for note in notes:
            printer.normal(__(note, indent=8 * " "))
        printer.normal("")

        self._plot(results)

        # The results, for variables, tables and properties
        data = {"T": T, "V": V, "nruns": nruns}
        for key, method in (
            ("viscosity", "Green-Kubo"),
            ("viscosity (Helfand)", "Helfand"),
        ):
            fit = results[method][-1]
            if fit is not None:
                data[key] = fit["value"]
                if fit["stderr_runs"] is not None:
                    data[key + ",stderr"] = fit["stderr_runs"]

        if "viscosity" in data:
            stderr = data.get("viscosity,stderr", None)
            publish_viscosity(
                self,
                Q_(data["viscosity"], "mPa*s").m_as("Pa*s"),
                stderr=None if stderr is None else Q_(stderr, "mPa*s").m_as("Pa*s"),
                T=T,
                method="Green-Kubo",
            )

        # The model fills in the names of the properties
        self.model = self._find_model()
        self.store_results(configuration=configuration, data=data)

    def _find_model(self):
        """The model (forcefield or model chemistry) of the MD, for properties."""
        try:
            if self.variable_exists("_model_chemistry"):
                return self.get_variable("_model_chemistry")["level"]
            ff = self.get_variable("_forcefield")
            if ff == "OpenKIM":
                return "OpenKIM/" + self.get_variable("_OpenKIM_Potential")
            if ff == "PyTorch":
                return "PyTorch/" + self.get_variable("_pytorch_model")
            return ff.current_forcefield
        except Exception:
            return "unknown"

    def _plot(self, results):
        """Plot the correlation functions, Green-Kubo integrals and Helfand moments."""
        dt = self.timestep[0] * 1.0e-3  # ps
        integral, integral_err = self._curves(self.integrals)
        M, M_err = self._curves(self.Ms)
        acf, acf_err = self._curves(self.acfs)
        ts = (np.arange(integral.shape[1]) * dt).tolist()
        end = len(ts) // 4

        for filename, name, rows, err, fits, args in (
            (
                "PP_correlation.graph",
                "ACF",
                acf[:, :end],
                acf_err[:, :end],
                None,
                ("P0Pt", "<P(0)P(t)> (MPa^2)", "MPa^2"),
            ),
            (
                "GK_integral.graph",
                "GKI",
                integral,
                integral_err,
                results["Green-Kubo"],
                ("η", "η (mPa*s)", "mPa*s"),
            ),
            (
                "HelfandMoments.graph",
                "M",
                M,
                M_err,
                results["Helfand"],
                ("M", "M (mPa*s*ps)", "mPa*s*ps"),
            ),
        ):
            figure = self.create_figure(
                module_path=("seamm",),
                template="line.graph_template",
                title=filename.split(".")[0].replace("_", " "),
            )
            fit = None
            if fits is not None:
                fit = [
                    f if f is not None else {"xs": [], "ys": [], "value": 0.0}
                    for f in fits
                ]
            ts_plot = ts[: rows.shape[1]]
            plot_components(
                figure,
                name,
                rows,
                ts_plot,
                *args,
                err=err,
                fit=fit,
                symbol="η",
                fit_units="mPa*s",
                fit_yunits=args[2],
            )
            figure.grid_plots(name)
            path = Path(self.directory) / filename
            figure.dump(path)
            if "html" in self.options and self.options["html"]:
                figure.template = "line.html_template"
                figure.dump(path.with_suffix(".html"))

    def md_requirements(self, P=None):
        """What the analysis needs from the MD: the pressure tensor every 5 fs."""
        return Requirements(
            needs=[
                Need(
                    "pressure_tensor",
                    max_time=MAX_INTERVAL,
                    why="the viscosity, to follow its fastest fluctuations",
                )
            ],
            periodic=True,
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed.

        Parameters
        ----------
        preceding : [seamm.Node]
            The steps that run before this one.

        Returns
        -------
        [str]
            A description of each problem; empty if none.
        """
        source = self._static("md source")
        if source is None:
            return []
        return check_md_flowchart(self, preceding, self.md_requirements(), source)

    def check_configuration(self, configuration):
        """Called by an MD step before it runs MD that this step will analyze.

        Raises
        ------
        RuntimeError
            If the system is not periodic.
        """
        if self._static("md source") == run_own_md:
            return
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: the pressure tensor is divided by the volume of the cell."
            )

    def _read_run(self, run, run_dir):
        """The pressure tensor of one run of the sub-flowchart."""
        paths = []
        for path in sorted(run_dir.glob("**/md.h5md")):
            with open_md(path) as md:
                if md.has_observable("pressure_tensor"):
                    paths.append(path)
        if len(paths) == 0:
            raise RuntimeError(
                f"There is no pressure tensor for run {run} in {run_dir}."
            )
        # The last MD run in the sub-flowchart, which sampled the pressure tensor
        return self._read_md_file(paths[-1])

    def _read_md_file(self, path):
        """The pressure tensor, interval, temperature and volume of one md.h5md.

        Returns
        -------
        (numpy.ndarray(n, 3, 3), pint.Quantity, float, float)
            The pressure tensor in Pa, the interval between samples, the mean
            temperature in K and the mean volume in Å^3.
        """
        with open_md(path) as md:
            if not md.periodic:
                raise RuntimeError(
                    f"The MD in {path} is not periodic, so it has no viscosity."
                )
            if not md.has_observable("pressure_tensor"):
                raise RuntimeError(f"There is no pressure tensor in {path}.")
            tensor = md.observable("pressure_tensor", units="Pa", time_units="fs")
            T = float(md.observable("temperature", units="K").values.mean())
            if md.has_observable("volume"):
                V = float(md.observable("volume", units="Å^3").values.mean())
            else:
                V = float(abs(np.linalg.det(md.cell(units="Å"))))
        return tensor.values, tensor.dt, T, V

    @staticmethod
    def _split_into_blocks(data, n_blocks):
        """Divide one run into consecutive blocks, analyzed as separate runs."""
        tensor, dt, T, V = data
        n_blocks = max(1, int(n_blocks))
        n = tensor.shape[0] // n_blocks
        if n < MIN_BLOCK_SAMPLES:
            raise RuntimeError(
                f"The run has only {tensor.shape[0]} samples of the pressure tensor, "
                f"too few for {n_blocks} blocks: each block needs at least "
                f"{MIN_BLOCK_SAMPLES} for the correlation functions to be fitted. "
                "Use fewer blocks or a longer run."
            )
        return [(tensor[i * n : (i + 1) * n], dt, T, V) for i in range(n_blocks)]

    def set_id(self, node_id=()):
        """Sequentially number the subnodes"""
        self.logger.debug("Setting ids for subflowchart {}".format(self))
        if self.visited:
            return None
        else:
            self.visited = True
            self._id = node_id
            self.set_subids(self._id)
            return self.next()

    def set_subids(self, node_id=()):
        """Set the ids of the nodes in the subflowchart"""
        node = self.subflowchart.get_node("1").next()
        n = 1
        while node is not None:
            node = node.set_id((*node_id, str(n)))
            n += 1
