.. _user-guide:

**********
User Guide
**********

The Viscosity step calculates the shear viscosity of a liquid from the fluctuations
of the pressure tensor in equilibrium molecular dynamics. The pressure tensor must be
sampled often, at least every 5 fs, to follow its fastest fluctuations; the NVE, NVT
and NPT sub-steps of LAMMPS do this with **Sample the pressure tensor**.

For an isotropic fluid, the five independent shear components of the pressure tensor,
Pxy, Pxz, Pyz, (Pxx - Pyy)/2 and (Pyy - Pzz)/2, each give the viscosity, so the step
analyzes all five and averages them, which gives better statistics than the three
off-diagonal components alone (Daivis and Evans, 1994). The step reports the viscosity
of each component and their average in two ways, which should agree within their
errors:

Green-Kubo
    The plateau of the running integral of the autocorrelation function of the shear
    stress, multiplied by V / k\ :sub:`B`\ T, fitted by one or, if justified by the
    data, two rising exponentials.
Helfand
    The slope at long times of the Helfand moments, the mean square of the integral
    of the shear stress (Helfand, 1960), starting once the correlations have decayed.

The mean of the shear stress of a fluid is zero, so it is not subtracted from the
series: subtracting the mean of a finite run lowers the correlation function at every
time and biases the viscosity low. For the same reason a solid, or a system under
shear stress, does not have a viscosity that this step can calculate.

The step needs a periodic system, since the pressure is divided by the volume of the
cell; it stops before the MD runs if the structure is not periodic. It reads the
pressure tensor, the temperature and the volume from the ``md.h5md`` file that the
LAMMPS step writes for each MD run.


Where the MD comes from
=======================

The step analyzes MD that ran before it in the flowchart, or the ``md.h5md`` files of
earlier runs. For independent runs, put the MD in a Loop -- in parallel, with each
run starting from the same structure ("Discard the structure" in the MD sub-steps)
and new velocities -- and analyze "the last MD run of each loop iteration", the
default. A typical flowchart:

.. code-block:: text

    Model Chemistry (FF@...), Packmol, minimization
    NPT equilibration; Set Cell to its mean volume; NVT equilibration
    Loop over the replicas (parallel)
        LAMMPS: velocities, NVT decorrelation, NVT production
    Viscosity
    Diffusivity

**MD to analyze** chooses the runs:

the last MD run of each loop iteration
    The last MD run of each iteration of the loop around the MD, analyzed as
    independent runs; without a loop, the one MD run before the step. This is the
    default.
the last MD run
    The last MD run earlier in the flowchart that sampled the pressure tensor. One
    production run can feed several analysis steps, such as Viscosity and
    Diffusivity.
files
    The ``md.h5md`` files given, each an independent run; ``job://<n>/<path>`` reads a
    file of job ``n``.

The MD must sample the pressure tensor ("Sample the pressure tensor" in the LAMMPS
NVE, NVT and NPT sub-steps) at least every 5 fs. Checking the flowchart before it runs
reports what is missing, with a fix that turns the sampling on, and also warns when
the analysis is inside the loop holding its MD, when the loop changes the temperature
or pressure (the runs are then not replicas of one state), or when the thermostat
changes the dynamics (Langevin, Andersen) or there is none (NVE). Before the MD runs,
the step warns if the system is small for the analysis (fewer than about 200
molecules, or a cell narrow for the molecules).

The viscosity is a collective property: it does not average over the molecules, so a
larger box does not make it more precise. Use a box that is just large enough, and
spend the time on independent runs, each a few hundred picoseconds or more; for water,
runs adding up to a few nanoseconds give the viscosity to within a few percent.

Convergence and errors
======================

Each value is fitted to the average over the runs. Its error comes from the spread of
the values fitted to each run separately, so at least two runs are needed for errors,
and more give better ones. A single run is divided into **Blocks for errors from one
run** (5 by default), which are analyzed separately instead; each block must be much
longer than the correlation time, and have at least 2000 samples.

The table gives the 95% confidence interval: Student's t with n - 1 degrees of
freedom (2.78 for 5 runs) times the standard error over the runs or blocks, and says
so under the table. The value stored and passed to later steps carries the standard
error, which the output also gives.

For an isotropic run -- a cubic cell with nothing directional in it -- the viscosity
is one value, the average of the five shear components within each run, and a line
checks that the components agree within their errors (χ²). Otherwise the components
are given in the table too. The value of each component, and of each run or block,
are in the detailed output.

Results
=======

* the shear viscosity by Green-Kubo, with its standard error, stored as the property
  ``viscosity#Green-Kubo#<model>``; it can be put in variables and tables, and is
  passed to later steps: a Diffusivity step after it uses it for the Yeh-Hummer
  correction of the diffusion coefficients;
* the viscosity from the Helfand moments, stored as ``viscosity#Helfand#<model>``;
* the state of the runs: temperature, pressure, volume and density with their errors,
  and for runs at constant volume the equilibrium density at the target pressure when
  an NPT run before them gives the bulk modulus;
* the number of runs or blocks.

The step also writes plots of the correlation functions, the Green-Kubo integrals and
the Helfand moments.

Index
=====

* :ref:`genindex`
