=======
History
=======
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

