"""Ready-made samples, for any instrument.

:func:`aluminum` and :func:`la2ni7` return a :class:`.Single_crystal` with its own lattice,
orientation, size and mosaic, and its CIF file (if any) next to it: they scatter
elastically. :func:`lead` returns a :class:`.Phonon_simple`, which scatters from one phonon
branch. :func:`vanadium` returns an :class:`.Incoherent`, which scatters elastically in every
direction. Placement is up to the caller: ``at``, and ``relative`` to a component of the
instrument, e.g. the goniometer's ``"gonio_phi"``. Keyword arguments set the sample's other
parameters or change the defaults, e.g.
``la2ni7(name="sample", mosaic=60, relative="gonio_phi")``.

``reflections`` names a CIF file relative to this folder (``CIF_DIR`` in
core/components/crystal.py), so a saved parameters file is portable.

The samples are sorted by how they scatter and what they are: ``elastic/single_crystal.py``,
``elastic/powder.py`` and ``inelastic/single_crystal.py``.

To add a sample:

1. Put its CIF file next to its Python file (e.g. ``NaCl.cif``), or use a McStas built-in
   reflection list (e.g. ``"Al.lau"``), which needs no file.
2. Add a file in the fitting folder, with a function like :func:`la2ni7`, with
   ``lattice_params``, ``along_beam``, ``in_plane``, ``reflections`` (the CIF path relative to
   this folder) and any other ``Single_crystal`` value it needs (e.g. ``mosaic``), then
   ``**params`` last so a caller can change any of them.
3. Import it in that folder's ``__init__.py`` and here, and add its name to ``__all__``.
4. Add it to the loop of ``test_the_samples_are_crystals_placed_like_any_component``
   (tests/samples/test_samples.py).
"""

from mcdgs.samples.elastic.powder import vanadium
from mcdgs.samples.elastic.single_crystal import aluminum, la2ni7
from mcdgs.samples.inelastic.single_crystal import lead

__all__ = ["aluminum", "la2ni7", "lead", "vanadium"]
