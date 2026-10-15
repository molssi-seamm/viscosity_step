.. _regression:

Regression: SPC/E water
=======================

The tests in ``tests/`` check the analysis on synthetic data with a known viscosity.
They cannot catch problems in the MD itself, such as the wrong constraints, or in the
chain from the LAMMPS step through ``md.h5md`` to the Viscosity and Diffusivity steps.
This regression does, by computing the viscosity and diffusion coefficient of rigid
SPC/E water and comparing them with published values. It takes too long for CI, so it
is run by hand: after changes to this step, to Diffusivity, to ``seamm_md_analysis``,
or to how the LAMMPS step samples the pressure tensor or applies constraints.

The files are in ``regression/spce/``:

``spce.yaml``
    The flowchart, as a specification for ``seamm-flowchart build``.
``spce.flow``
    The flowchart built from it, ready to run.
``check.py``
    Reads the results from ``job.out`` and checks them against the expected ranges.

The calculation
---------------

* 512 SPC/E waters (``nacl_water.frc``), rigid, held by SHAKE with a tolerance of
  1e-6, in a cubic cell at 0.997 g/mL (24.86 Å), packed by Packmol.
* No minimization: the hydrogens of SPC/E have no Lennard-Jones repulsion, and
  without constraints the minimizer can collapse them onto oxygens.
* LAMMPS: PPPM (1e-5), a 10 Å cutoff with the long-range dispersion correction, Nosé–
  Hoover NVT at 298.15 K, 2 fs timestep. 100 ps of equilibration, then 3 ns of
  production sampling the pressure tensor every 4 fs and the centres of mass every
  1 ps.
* Viscosity and then Diffusivity (mean square displacement of the molecules, with the
  Yeh–Hummer correction from the Viscosity step) analyze the production run in 5
  blocks of 600 ps.

It takes about 35 minutes on 4 cores of an Apple M5 laptop, and 50 on an M3 Pro.

Running it
----------

From an empty directory::

    run_flowchart /path/to/viscosity_step/regression/spce/spce.flow
    python /path/to/viscosity_step/regression/spce/check.py .

``check.py`` takes one or more job directories (or ``job.out`` files); with several
it also checks their average, which reduces the noise by the square root of the
number of runs.

Expected values
---------------

The ranges are the references widened by the 95% error of a single 3 ns run, about
12% for the viscosity:

.. list-table::
   :header-rows: 1
   :widths: 40 20 40

   * - Quantity
     - Accepted
     - References
   * - Viscosity, Green–Kubo average (mPa s)
     - 0.62–0.82
     - 0.72 ± 0.02 (our LAMMPS baseline, 3 × 1.7 ns, 1 fs); 0.68 ± 0.02 (Tazi et
       al., J. Phys.: Condens. Matter 24, 284117 (2012))
   * - D\ :sub:`0`, corrected (10\ :sup:`-9` m\ :sup:`2`/s)
     - 2.75–3.05
     - 2.86 (baseline; 2.51 ± 0.05 before the correction); 2.97 ± 0.05 (Tazi et al.)

The baseline was computed outside SEAMM, directly with LAMMPS, at the NPT density of
the model (0.9995 g/mL) with a 1 fs timestep. The small differences in density and
timestep do not matter at this precision.

Results
-------

.. list-table::
   :header-rows: 1
   :widths: 22 30 16 16 16

   * - Date
     - Versions, machine
     - η, GK (mPa s)
     - η, Helfand
     - D\ :sub:`0`
   * - 2026-10-15
     - lammps-step 2026.10.15.4, viscosity-step 2026.10.15; Apple M5
     - 0.716 ± 0.085
     - 0.72 ± 0.12
     - 2.91 ± 0.08
   * - 2026-10-15
     - the same; Apple M3 Pro, a second run
     - 0.725 ± 0.036
     - 0.698 ± 0.087
     - 2.80 ± 0.08
   * - 2026-10-15
     - the average of the two runs
     - 0.720
     -
     - 2.86

The five shear components gave 0.74, 0.715, 0.76, 0.635 and 0.73 mPa s in the first
run and 0.79, 0.77, 0.709, 0.789 and 0.61 in the second; a spread of this size is
normal for runs of this length. The average of the two runs reproduces the baseline,
0.72 mPa s and 2.86 × 10\ :sup:`-9` m\ :sup:`2`/s.

History
-------

The first run of this regression, with lammps-step 2026.10.15.3, failed: η came out
1.42 mPa s, and D in y and z was up to 60 times too large. Rigid water was then held
by LAMMPS's RATTLE, which for this system gave a pressure about 850 atm too high, Pxx
about 150 atm above Pyy and Pzz for the whole run, and a center of mass that drifted
ever faster. Because the mean of the shear stress is not subtracted (it is zero for a
fluid), the persistent (Pxx - Pyy)/2 integrated into a viscosity of 6.3 for that
component. lammps-step 2026.10.15.4 uses SHAKE instead, and the regression passes.
