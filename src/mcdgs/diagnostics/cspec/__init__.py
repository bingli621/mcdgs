"""The diagnostics of CSPEC, one module per page of its documentation.

In the order of the beam, which is also the order ``pixi run docs-diagnostics`` runs them in:

1. :mod:`.overview`: the printout, the instrument file and the drawings of CSPEC.
2. :mod:`.source`: the source, at the entrance of the guide system.
3. :mod:`.guides_diagnostics`: the beam along the guides, with McStasScript's diagnostics.
4. :mod:`.beam_monitor`: the beam at the beam monitor, after the BBG.
5. :mod:`.time_distance`: the time-distance diagram.
"""
