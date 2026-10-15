=======
History
=======
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

