2026-10-15 Initial development: phase 4 of the MD analysis campaign
===================================================================

The Viscosity step was written in phase 4 of the campaign that moved SEAMM's MD
analysis onto the shared library ``seamm_md_analysis``. The design, the decisions
and what was done in each phase are recorded with that library, in its developer
guide (`MD analysis campaign
<https://molssi-seamm.github.io/seamm_md_analysis/developer_guide/campaigns/2026-10-14/index.html>`_).

In brief, the step:

* analyzes the five independent shear components of the pressure tensor, Pxy, Pxz,
  Pyz, (Pxx-Pyy)/2 and (Pyy-Pzz)/2, by Green-Kubo and from the Helfand moments, and
  fits their average;
* does not subtract the mean of the shear stress, which is zero for a fluid, since
  subtracting the mean of a finite run biases the viscosity low;
* reads the MD from ``md.h5md`` files, running its own MD or analyzing earlier MD,
  and checks before the job runs that the pressure tensor is sampled at least every
  5 fs;
* publishes the Green-Kubo average in the variable ``_viscosity`` for a following
  Diffusivity step.

Phase 8 (2026-10-17) retired the step's own MD: the MD now runs before the step, in a
Loop for independent runs, and the step only analyzes it.
