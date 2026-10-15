=======================
SEAMM Viscosity Plug-in
=======================

.. image:: https://img.shields.io/github/issues-pr-raw/molssi-seamm/viscosity_step
   :target: https://github.com/molssi-seamm/viscosity_step/pulls
   :alt: GitHub pull requests

.. image:: https://github.com/molssi-seamm/viscosity_step/workflows/CI/badge.svg
   :target: https://github.com/molssi-seamm/viscosity_step/actions
   :alt: Build Status

.. image:: https://codecov.io/gh/molssi-seamm/viscosity_step/branch/main/graph/badge.svg
   :target: https://codecov.io/gh/molssi-seamm/viscosity_step
   :alt: Code Coverage

.. image:: https://github.com/molssi-seamm/viscosity_step/workflows/CodeQL/badge.svg
   :target: https://github.com/molssi-seamm/viscosity_step/security/code-scanning
   :alt: Code Quality

.. image:: https://github.com/molssi-seamm/viscosity_step/workflows/Release/badge.svg
   :target: https://molssi-seamm.github.io/viscosity_step/index.html
   :alt: Documentation Status

.. image:: https://img.shields.io/pypi/v/viscosity_step.svg
   :target: https://pypi.python.org/pypi/viscosity_step
   :alt: PyPi VERSION

A SEAMM plug-in for the shear viscosity of liquids from equilibrium molecular
dynamics.

* Free software: BSD-3-Clause
* Documentation: https://molssi-seamm.github.io/viscosity_step/index.html
* Code: https://github.com/molssi-seamm/viscosity_step

Features
--------

* The shear viscosity by Green-Kubo and from the Helfand moments, averaging the
  five independent shear components of the pressure tensor.
* Runs its own MD several times, or analyzes MD run earlier in the flowchart or in
  other jobs, from the ``md.h5md`` files the LAMMPS step writes.
* Checks before the job runs that the MD samples the pressure tensor often enough.
* Passes its result to a later Diffusivity step for the Yeh-Hummer correction.

Acknowledgements
----------------

This package was created with the `molssi-seamm/cookiecutter-seamm-plugin`_ tool, which
is based on the excellent Cookiecutter_.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`molssi-seamm/cookiecutter-seamm-plugin`: https://github.com/molssi-seamm/cookiecutter-seamm-plugin

Developed by the Molecular Sciences Software Institute (MolSSI_),
which receives funding from the `National Science Foundation`_ under
award CHE-2136142.

.. _MolSSI: https://molssi.org
.. _`National Science Foundation`: https://www.nsf.gov
