# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialog: create it and re-lay it out for every choice of
the MD to analyze. Skipped when no display is available."""

import pytest

from seamm_md_analysis.step_support import read_files, run_own_md, source_choices


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("Viscosity")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("Viscosity")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def test_dialog_layouts(tk_node):
    tk_node.create_dialog()
    tk_node.reset_dialog()
    notebook = tk_node["notebook"]
    tab = str(tk_node["flowchart frame"])
    for source in source_choices:
        tk_node["md source"].set(source)
        tk_node.reset_viscosity_frame()
        own = source == run_own_md
        assert (tk_node["nruns"].grid_info() != {}) == own
        assert (tk_node["errors"].grid_info() != {}) == own
        assert (tk_node["md files"].grid_info() != {}) == (source == read_files)
        assert (tk_node["number of blocks"].grid_info() != {}) == (not own)
        assert (notebook.tab(tab, "state") == "hidden") == (not own)
