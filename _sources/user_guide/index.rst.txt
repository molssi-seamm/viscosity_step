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

Each value is fitted to the average over the runs. Its error, given as a 95%
interval (twice the standard error), comes from the spread of the values fitted to
each run separately, so at least two runs (or blocks) are needed for errors, and more
give better ones.

The mean of the shear stress of a fluid is zero, so it is not subtracted from the
series: subtracting the mean of a finite run lowers the correlation function at every
time and biases the viscosity low. For the same reason a solid, or a system under
shear stress, does not have a viscosity that this step can calculate.

The Green-Kubo average is the result: it is stored as the property
``viscosity#Green-Kubo#<model>`` (with its standard error), can be put in variables
and tables, and is passed to later steps. A Diffusivity step after it uses it for the
Yeh-Hummer correction of the diffusion coefficients. The value from the Helfand
moments is stored as ``viscosity#Helfand#<model>``. The step also writes plots of the
correlation functions, the Green-Kubo integrals and the Helfand moments.

The step needs a periodic system, since the pressure is divided by the volume of the
cell; it stops before running any MD if the structure is not periodic. It reads the
pressure tensor, the temperature and the volume from the ``md.h5md`` file that the
LAMMPS step writes for each MD run.


Where the MD comes from
=======================

**MD to analyze** chooses where the MD comes from:

run the MD in the sub-flowchart
    The step runs its sub-flowchart, usually a LAMMPS step, the given number of times,
    and the errors come from the spread of the runs. If the last MD sub-step does not
    sample the pressure tensor often enough, the step turns it on, every 5 fs, with a
    note in the output.
the last MD run
    The last MD run earlier in the flowchart that sampled the pressure tensor. One
    production run can feed several analysis steps, such as Viscosity and Diffusivity,
    and the sub-flowchart is not used.
the last MD run of each loop iteration
    For a LAMMPS step inside a Loop: the run of each iteration, analyzed as
    independent runs.
files
    The ``md.h5md`` files given, each an independent run; ``job://<n>/<path>`` reads a
    file of job ``n``.

A single run is divided into **Blocks for errors from one run** (5 by default),
analyzed separately to estimate the errors. Each block must be much longer than the
correlation time, and have at least 2000 samples. Before the job runs, the step checks
that the MD before it samples the pressure tensor often enough, and the job does not
start if it does not.

A typical flowchart for a liquid equilibrates it with NPT, then runs a long NVT
production run that samples the pressure tensor every 5 fs and the centres of mass of
the molecules, followed by Viscosity and then Diffusivity, both analyzing that run.
For water, runs of a few nanoseconds are needed for the viscosity to within a few
percent.

Index
=====

* :ref:`genindex`
