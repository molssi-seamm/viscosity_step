# -*- coding: utf-8 -*-
"""
Control parameters for the Viscosity step in a SEAMM flowchart
"""

import logging
import seamm
from seamm_md_analysis.step_support import source_parameters

logger = logging.getLogger(__name__)


class ViscosityParameters(seamm.Parameters):
    """
    The control parameters for Viscosity.

    The MD to analyze is run earlier in the flowchart -- in a Loop for independent
    runs -- or in other jobs (see seamm_md_analysis.step_support.source_parameters).
    """

    # Dropped with the own-MD route (phase 8): ignored in older flowcharts
    retired = ("nruns", "errors")

    parameters = {
        **source_parameters("the viscosity"),
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("ViscosityParameters.__init__")

        super().__init__(
            defaults={**ViscosityParameters.parameters, **defaults}, data=data
        )
