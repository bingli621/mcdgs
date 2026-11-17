"""Setups that study one part of an instrument, such as a section of it with extra monitors.

They build on the instruments, mcstas/, tof/ and reduction/, and nothing else in mcdgs imports
them. One package per instrument, as in instruments/. The modules share :mod:`.compare`, which
puts two time and wavelength maps of a monitor side by side and draws them, :mod:`.drawing`,
which describes and draws any instrument, and :mod:`.cli`, their command line and files.
"""
