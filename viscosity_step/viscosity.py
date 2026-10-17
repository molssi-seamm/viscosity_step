# -*- coding: utf-8 -*-

"""Non-graphical part of the Viscosity step in a SEAMM flowchart

The shear viscosity from the fluctuations of the pressure tensor in equilibrium
MD, by Green-Kubo and from the Helfand moments. The five independent shear
components of the pressure tensor (Pxy, Pxz, Pyz, (Pxx-Pyy)/2 and (Pyy-Pzz)/2)
are analyzed separately and averaged within each run.

The MD is run before this step -- in a Loop for independent runs -- or read from
files; the step analyzes it.
"""

import importlib.resources
import logging
from pathlib import Path

import numpy as np
from seamm_md_analysis import (
    SHEAR_COMPONENTS,
    Estimate,
    autocorrelation_functions,
    component_check,
    component_helfand_moments,
    fit_green_kubo_integral,
    fit_helfand_moment,
    isotropy,
    open_md,
    shear_stress_components,
)
from seamm_md_analysis.plots import plot_components
from seamm_md_analysis.requirements import Need, Requirements
from seamm_md_analysis.stats import fmt_err
from seamm_md_analysis.step_support import (
    check_md_flowchart,
    md_inputs,
    publish_viscosity,
    read_files,
    report_state,
    size_warnings,
    state_data,
)
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

#: The methods, in the order they are reported
METHODS = ("Green-Kubo", "Helfand")

#: The units of the viscosity in the output
UNITS = "mPa·s"


class Viscosity(seamm.Node):
    """
    The non-graphical part of a Viscosity step in a flowchart.

    Attributes
    ----------
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
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Viscosity {self}")
        super().__init__(
            flowchart=flowchart,
            title="Viscosity",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._metadata = viscosity_step.metadata
        self.parameters = viscosity_step.ViscosityParameters()
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
        self.paths = []  # The md.h5md files analyzed
        self.T = []  # The mean temperature of each run, K
        self.V = []  # The mean volume of each run, Å^3
        self.timestep = []  # The interval between samples of each run, fs
        self.n_samples = []  # The number of samples of each run or block
        self.acfs = []  # The autocorrelation functions per run, 5 x m, MPa^2
        self.integrals = []  # The Green-Kubo running integrals per run, 5 x m
        self.Ms = []  # The Helfand moments per run, 5 x m, mPa*s*ps

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = "viscosity-step"
        parser = seamm_util.getParser()

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
        return self.next()

    def description_text(self, P=None, short=False):
        """Create the text description of what this step will do.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        short : bool = False
            Unused: the description is always short.

        Returns
        -------
        str
            A description of the current step.
        """
        if P is None:
            P = self.parameters.values_to_dict()

        source = P["md source"]
        text = (
            "Calculate the shear viscosity by Green-Kubo and from the Helfand "
            "moments of "
        )
        if source == read_files:
            text += f"the MD in {P['md files']}."
        else:
            text += f"{source} earlier in the flowchart."
        return self.header + "\n" + __(text, indent=4 * " ").__str__()

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
        printer.important(__(self.description_text(P), indent=self.indent))

        _, configuration = self.get_system_configuration()

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

        self.paths = md_inputs(self, P, name="pressure_tensor")
        if len(self.paths) == 0:
            raise RuntimeError("No MD run finished, so there is nothing to analyze.")
        runs = [self._read_md_file(path) for path in self.paths]

        text = f"Analyzing {len(runs)} MD run{'s' if len(runs) > 1 else ''}"
        if len(runs) == 1:
            source = "blocks"
            runs = self._split_into_blocks(runs[0], P["number of blocks"])
            text += f", divided into {len(runs)} blocks to estimate the errors"
        else:
            source = "runs"
        text += ": " + ", ".join(str(p) for p in self.paths)
        printer.job(__(text, indent=self.indent + 4 * " "))
        printer.job("")

        for data in runs:
            self.process_run(data)

        self.analyze(P=P, configuration=configuration, source=source)

        return next_node

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
        self.n_samples.append(n)
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

    def _curves(self, per_run):
        """The mean over the runs of each component and their average, and the
        standard deviation over the runs (the error bands of the plots)."""
        rows = np.stack([self._with_average(r) for r in per_run])
        return rows.mean(axis=0), rows.std(axis=0)

    def _fit(self):
        """Fit the Green-Kubo integrals and Helfand moments.

        The value of each component and of their average is fitted to the mean
        over the runs. The average is that of the components within each run, so
        its error is that of the per-run averages. The error comes from the values
        fitted to each run (or block) separately: the errors of the fit itself
        treat strongly correlated points as independent, so they are far too small.

        Returns
        -------
        {str: [dict]}
            For "Green-Kubo" and "Helfand", a dict for each component and their
            average with the value, the value fitted to each run ("samples", NaN
            where the fit failed), their standard error ("stderr_runs") and the
            fitted curve. A component that could not be fitted has None.
        """
        dt = self.timestep[0] * 1.0e-3  # ps
        integral, integral_err = self._curves(self.integrals)
        M, M_err = self._curves(self.Ms)
        ts = (np.arange(integral.shape[1]) * dt).tolist()

        def over_runs(fit_one, per_run, i):
            values = []
            for rows in per_run:
                try:
                    value = fit_one(self._with_average(rows)[i])
                except Exception:
                    value = float("nan")
                values.append(value if np.isfinite(value) else float("nan"))
            values = np.array(values, dtype=float)
            finite = values[np.isfinite(values)]
            if finite.size < 2:
                return values, None
            return values, float(np.std(finite, ddof=1) / np.sqrt(finite.size))

        results = {"Green-Kubo": [], "Helfand": []}
        for i in range(len(COMPONENTS)):
            try:
                fit = fit_green_kubo_integral(integral[i], ts, sigma=integral_err[i])
                samples, stderr = over_runs(
                    lambda y: fit_green_kubo_integral(y, ts).value,
                    self.integrals,
                    i,
                )
                gk = {
                    "value": fit.value,
                    "stderr": fit.stderr,
                    "samples": samples,
                    "stderr_runs": stderr,
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
                samples, stderr = over_runs(
                    lambda y: fit_helfand_moment(y, ts, start=start).slope,
                    self.Ms,
                    i,
                )
                helfand = {
                    "value": fit.slope,
                    "stderr": fit.stderr,
                    "samples": samples,
                    "stderr_runs": stderr,
                    "start": start,
                    "xs": fit.xs,
                    "ys": fit.ys,
                }
            except Exception as e:
                logger.warning(f"The Helfand fit of {COMPONENTS[i]} failed: {e}")
                helfand = None
            results["Helfand"].append(helfand)
        return results

    @staticmethod
    def _formatted(estimate):
        """The value and the half-width of its 95% interval as text; a value
        without an error to three decimals, like the values of each run."""
        if estimate.stderr is None:
            return f"{estimate.value:.3f}", ""
        return estimate.formatted()

    @staticmethod
    def _estimate(fit, source):
        """The Estimate of one fitted value: the fit to the mean over the runs,
        with the standard error of the values fitted to each run."""
        if fit is None:
            return None
        samples = np.asarray(fit["samples"], dtype=float)
        finite = samples[np.isfinite(samples)]
        if finite.size < 2:
            return Estimate(float(fit["value"]), None, 1, "mPa*s", source)
        estimate = Estimate.from_samples(finite, units="mPa*s", source=source)
        estimate.value = float(fit["value"])
        return estimate

    def analyze(self, indent="", P=None, configuration=None, source="runs", **kwargs):
        """Fit, print, plot and store the viscosity.

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        P : dict
            The values of the parameters.
        configuration : molsystem._Configuration
            The configuration, for storing the results as properties.
        source : str
            What the samples are: "runs" or "blocks" of a single run.
        """
        results = self._fit()
        nruns = len(self.Ms)
        estimates = {
            method: [self._estimate(fit, source) for fit in results[method]]
            for method in METHODS
        }
        eta = estimates["Green-Kubo"][-1]
        isotropic, why = self._isotropy()

        # The main table: the average, and the components if not isotropic
        table = {"Property": [], "Value": [], "± (95%)": [], "Units": []}

        def row(label, estimate):
            table["Property"].append(label)
            if estimate is None:
                table["Value"].append("failed")
                table["± (95%)"].append("")
            else:
                v, e = self._formatted(estimate)
                table["Value"].append(v)
                table["± (95%)"].append(e)
            table["Units"].append(UNITS)

        for method in METHODS:
            row(f"η ({method})", estimates[method][-1])
            if not isotropic:
                for name, estimate in zip(SHEAR_COMPONENTS, estimates[method][:-1]):
                    row(f"η_{name} ({method})", estimate)
        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left", "decimal", "decimal", "left"),
        )
        length = len(text.splitlines()[0])
        header = "Shear viscosity".center(length)
        printer.normal(
            __("\n" + header + "\n" + text + "\n", indent=8 * " ", wrap=False)
        )

        # What the table is from, the errors, and the check of the components
        dt = Q_(self.timestep[0], "fs")
        if source == "runs":
            what = f"From {nruns} runs"
        elif nruns > 1:
            what = f"From 1 run divided into {nruns} blocks"
        else:
            what = "From 1 run"
        n_samples = int(np.mean(self.n_samples))
        headline = eta if eta is not None else estimates["Helfand"][-1]
        lines = [
            f"{what} of {n_samples} samples of the pressure tensor every "
            f"{dt.m_as('fs'):.3g} fs."
            + ("" if headline is None else " " + headline.describe_error())
        ]
        if isotropic:
            components = self._component_samples(results["Green-Kubo"])
            if components is None:
                lines.append("The components cannot be compared: a fit failed.")
            else:
                check = component_check(
                    components, names=[f"η_{c}" for c in SHEAR_COMPONENTS]
                )
                lines.append(check.text.replace("Isotropy:", "Isotropy (Green-Kubo):"))
        else:
            lines.append(
                f"The run is not treated as isotropic ({why}), so the five shear "
                "components are given; η is their average, which is only meaningful "
                "for an isotropic fluid."
            )
        if eta is not None:
            if eta.stderr is None:
                passed = f"η = {eta.value:.4g} {UNITS}, without an error"
            else:
                v, e = fmt_err(eta.value, eta.stderr)
                passed = f"η = {v} ± {e} {UNITS} (± the standard error)"
            lines.append(
                "The Green-Kubo value is passed to later steps, such as the "
                f"Yeh-Hummer correction in Diffusivity: {passed}."
            )

        notes = []
        if results["Helfand"][-1] is None:
            notes.append(
                "The Helfand moments could not be fitted: the fit starts once the "
                "correlations have decayed, but the moments extend only to a tenth "
                "of the length of each run (or block). Use longer runs, or fewer "
                "blocks."
            )
        if dt > MAX_INTERVAL:
            notes.append(
                f"The pressure tensor was sampled every {dt:~P}. It should be "
                f"sampled at least every {MAX_INTERVAL:~P} to follow its fastest "
                "fluctuations."
            )
        for line in lines + notes:
            printer.normal(__(line, indent=8 * " "))
        printer.normal("")

        # The details: each component, and each run or block
        self._print_components(estimates)
        self._print_samples(results, source)

        state = report_state(self, printer, self.paths)

        self._plot(results)

        # The results, for variables, tables and properties
        data = {"nruns": nruns}
        for key, method in (
            ("viscosity", "Green-Kubo"),
            ("viscosity (Helfand)", "Helfand"),
        ):
            estimate = estimates[method][-1]
            if estimate is not None:
                data[key] = estimate.value
                if estimate.stderr is not None:
                    data[key + ",stderr"] = estimate.stderr
        if state is not None:
            data.update(state_data(state))

        if eta is not None:
            if state is not None and state.temperature is not None:
                T = state.temperature.value
            else:
                T = float(np.mean(self.T))
            publish_viscosity(
                self,
                Q_(eta.value, "mPa*s").m_as("Pa*s"),
                stderr=(
                    None if eta.stderr is None else Q_(eta.stderr, "mPa*s").m_as("Pa*s")
                ),
                T=T,
                method="Green-Kubo",
            )

        # The model fills in the names of the properties
        self.model = self._find_model()
        self.store_results(configuration=configuration, data=data)

    @staticmethod
    def _component_samples(fits):
        """The five components fitted to each run or block, (n, 5), leaving out
        runs where a fit failed; None if fewer than two runs remain."""
        if any(fit is None for fit in fits[:-1]):
            return None
        rows = np.array([fit["samples"] for fit in fits[:-1]]).T
        rows = rows[np.all(np.isfinite(rows), axis=1)]
        return rows if len(rows) >= 2 else None

    def _print_components(self, estimates):
        """The value of each component by each method, in the detailed output."""
        table = {"Method": [], "Component": [], "η": [], "± (95%)": []}
        for method in METHODS:
            for i, estimate in enumerate(estimates[method]):
                table["Method"].append(method if i == 0 else "")
                table["Component"].append(COMPONENTS[i])
                if estimate is None:
                    table["η"].append("failed")
                    table["± (95%)"].append("")
                else:
                    v, e = self._formatted(estimate)
                    table["η"].append(v)
                    table["± (95%)"].append(e)
        text = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            disable_numparse=True,
            colalign=("left", "left", "decimal", "decimal"),
        )
        length = len(text.splitlines()[0])
        title = f"Shear viscosity by component ({UNITS})".center(length)
        printer.verbose(
            __("\n" + title + "\n" + text + "\n", indent=8 * " ", wrap=False)
        )

    def _print_samples(self, results, source):
        """The values fitted to each run or block, in the detailed output."""
        what = "Run" if source == "runs" else "Block"
        for method in METHODS:
            fits = results[method]
            table = {what: list(range(1, len(self.Ms) + 1))}
            for name, fit in zip(COMPONENTS, fits):
                if fit is None:
                    table[name] = ["failed"] * len(self.Ms)
                else:
                    table[name] = [
                        "" if not np.isfinite(v) else f"{v:.3f}" for v in fit["samples"]
                    ]
            text = tabulate(
                table,
                headers="keys",
                tablefmt="rounded_outline",
                disable_numparse=True,
                colalign=("center", *(["decimal"] * len(COMPONENTS))),
            )
            length = len(text.splitlines()[0])
            title = f"η by {what.lower()}, {method} ({UNITS})".center(length)
            printer.verbose(
                __("\n" + title + "\n" + text + "\n", indent=8 * " ", wrap=False)
            )

    def _isotropy(self):
        """Whether all the runs are isotropic, and why not if they are not."""
        reason = ""
        for path in self.paths:
            with open_md(path) as md:
                isotropic, reason = isotropy(md)
            if not isotropic:
                return False, reason
        return True, reason

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
        """What the analysis needs from the MD: the pressure tensor every 5 fs, in
        a periodic cell, with dynamics that are not perturbed by the thermostat."""
        return Requirements(
            needs=[
                Need(
                    "pressure_tensor",
                    max_time=MAX_INTERVAL,
                    why="the viscosity, to follow its fastest fluctuations",
                )
            ],
            periodic=True,
            dynamics=True,
        )

    def _static(self, key):
        """A parameter's value without evaluating variables, or None."""
        parameter = self.parameters[key]
        if parameter.is_expr:
            return None
        return parameter.get()

    def check_flowchart(self, preceding):
        """Check, before the job runs, that the MD will sample what is needed and
        is set up for the analysis.

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

        Returns
        -------
        [str]
            Warnings about the structure, e.g. that it is small.

        Raises
        ------
        RuntimeError
            If the system is not periodic.
        """
        if configuration.periodicity != 3:
            raise RuntimeError(
                f"'{self.title}' will analyze the MD, but the system is not "
                "periodic: the pressure tensor is divided by the volume of the cell."
            )
        return [f"'{self.title}': {w}" for w in size_warnings(configuration)]

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
