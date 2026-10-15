# -*- coding: utf-8 -*-
"""
Control parameters for the Viscosity step in a SEAMM flowchart
"""

import logging
import seamm
from seamm_md_analysis.step_support import run_own_md, source_parameters

logger = logging.getLogger(__name__)


class ViscosityParameters(seamm.Parameters):
    """
    The control parameters for Viscosity.

    The MD to analyze comes either from this step's sub-flowchart, run several
    times as independent runs, or from MD run earlier in the flowchart or in other
    jobs (see seamm_md_analysis.step_support.source_parameters).
    """

    parameters = {
        **source_parameters("the viscosity"),
        "nruns": {
            "default": "5",
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Number of runs to average:",
            "help_text": "The number of independent MD runs to average.",
            "applies_when": {"md source": run_own_md},
        },
        "errors": {
            "default": "continue to next run",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "continue to next run",
                "exit the loop",
                "stop the job",
            ),
            "format_string": "",
            "description": "On errors",
            "help_text": "What to do if a run of the MD fails.",
            "applies_when": {"md source": run_own_md},
        },
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
