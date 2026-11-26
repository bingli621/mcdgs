"""The diagnostics of CSPEC, one module per page of its documentation.

In the order of the beam, which is also the order ``pixi run docs-diagnostics`` runs them in:

1. :mod:`.overview`: the printout, the instrument file and the drawings of CSPEC.
2. :mod:`.source`: the source, at the entrance of the guide system.
3. :mod:`.guides_diagnostics`: the beam along the guides, with McStasScript's diagnostics.
4. :mod:`.beam_monitor`: the beam at the beam monitor, after the BBG.
5. :mod:`.choppers`: each disk, and the acceptance diagram of the cascade.
6. :mod:`.time_distance`: the time-distance diagram.
7. :mod:`.neutrons_on_sample`: the neutrons at the sample monitor.
8. :mod:`.masked_source`: the plain and the masked source, at the sample monitor.
9. :mod:`.incident_energy`: the incident energy of each sub-pulse.
10. :mod:`.elastic_resolution_incoherent`: the elastic resolution, from vanadium.
11. :mod:`.magnon_powder`: the S(Q, E) of a powder magnon in every sub-pulse.
12. :mod:`.diffraction`: the horizontal scattering plane of a La2Ni7 crystal.

:mod:`.sample_runs` reduces the runs split at the sample of the diagnostics with a sample.
"""
