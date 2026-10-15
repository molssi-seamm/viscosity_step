***************
Getting Started
***************

Installation
============
The Viscosity step is installed with the rest of SEAMM by ``seamm-manager``. To add it
to an existing installation::

  seamm-manager install viscosity-step

See the `SEAMM installation`_ documentation for more.

.. _SEAMM installation: https://molssi-seamm.github.io/installation/index.html

An example
==========
The viscosity of water at 298 K, analyzing one production run:

#. **Forcefield**, choosing a water model, then **Packmol** to fill a cubic cell with
   a few hundred water molecules at the experimental density.
#. **LAMMPS**, with these sub-steps:

   * Initialization, Minimization and Velocities;
   * NPT dynamics for about 100 ps to equilibrate the density;
   * NVT dynamics for the production run, a few nanoseconds, with **Sample the
     pressure tensor** set to every 5 fs on its trajectory tab, and the centres of
     mass of the molecules sampled every 1 ps if the diffusivity is wanted too.

#. **Viscosity**, with **MD to analyze** set to "the last MD run". The run is divided
   into 5 blocks to estimate the errors.
#. Optionally **Diffusivity**, also analyzing the last MD run, which uses this
   viscosity for its finite-size correction.

If the pressure tensor is not sampled often enough, the flowchart is refused when it
is checked or submitted, with a message saying what to change.

That should be enough to get started. For more detail about the functionality in this
plug-in, see the :ref:`User Guide <user-guide>`.
