"""The diagnostics of T-Rex, one module per page of its documentation.

In the order of the beam, which is also the order ``pixi run docs-diagnostics`` runs them in:

1. :mod:`.overview`: the printout, the instrument file and the drawings of T-Rex.
2. :mod:`.source`: the source, at the entrance of guide 1.
3. :mod:`.extraction`: the extraction guides up to the beam monitor, with the bender out.
4. :mod:`.bender_scan`: the same with the bender in, and a scan of its rotation.
5. :mod:`.beam_monitor`: the beam monitor after guide 4.
6. :mod:`.choppers`: the disks and the chopper acceptance diagram.
7. :mod:`.time_distance`: the time-distance diagram.
8. :mod:`.neutrons_on_sample`: the neutrons at the sample monitor.
9. :mod:`.masked_source`: the source masked by the choppers, against the plain one.
10. :mod:`.incident_energy`: the incident energy of each sub-pulse.
11. :mod:`.elastic_plane`: the horizontal scattering plane of a La2Ni7 crystal. It writes the
    beam file that the next one reads.
12. :mod:`.elastic_scan`: a rotation scan of that crystal.
"""
