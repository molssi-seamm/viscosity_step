=======
History
=======
2026.10.17 -- Analyze MD run before the step; labelled errors and a component check
    * The step analyzes MD run before it in the flowchart -- the runs of a Loop around
      the MD, as independent replicas ("the last MD run of each loop iteration", the
      default) -- or md.h5md files; it no longer runs its own MD in a sub-flowchart.
      Older flowcharts that stored its settings for that can still be read; the check
      before the job runs says what replaces them.
    * The check before the job runs also warns when the analysis is inside the loop
      holding its MD, when the loop changes the temperature or pressure, and when the
      thermostat does not suit the analysis, with fixes; the step warns before the MD
      when the system is small.
    * Errors are the 95% confidence interval, Student's t for the number of runs or
      blocks, and the output says so; stored errors are standard errors.
    * The state of the runs -- temperature, pressure, volume and density with their
      errors, and the equilibrium density where an NPT run before them gives the bulk
      modulus -- is reported and stored. The values of each run or block are in
      step.out.
    * The 95% interval was 2 × the standard error, too small for few runs; the
      viscosity passed to later steps, e.g. for Diffusivity's box-size correction,
      carries the standard error, labelled as such.
    * One line checks that the five shear components agree; the components are given
      in full for a run that is not isotropic.

2026.10.15.1 -- Internal: a regression test against SPC/E water
    * ``regression/spce`` holds a flowchart computing the viscosity and diffusion
      coefficient of rigid SPC/E water, and a script checking them against published
      values. It is run by hand after changes that could affect the results; the new
      "Regression" page in the developer guide describes it and records its results.

2026.10.15 -- First release: the shear viscosity from MD
    * Calculates the shear viscosity of a liquid from equilibrium MD by Green-Kubo
      and from the Helfand moments, averaging the five independent shear components
      of the pressure tensor.
    * Runs its own MD several times, or analyzes MD run earlier in the flowchart, in
      each iteration of a loop, or in other jobs, from the md.h5md files the LAMMPS
      step writes. A single run is divided into blocks to estimate the errors.
    * Checks before the job runs that the MD samples the pressure tensor at least
      every 5 fs, and turns it on when it runs its own MD.
    * The result is stored as a property, can be saved in variables and tables, and
      is used by a following Diffusivity step for its finite-size correction.

