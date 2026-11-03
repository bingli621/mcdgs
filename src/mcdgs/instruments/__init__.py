"""The instruments, one package each.

Each is an Instrument subclass (core/) with its own components and parameters. Importing this
package adds every instrument to ``INSTRUMENTS``; io/ imports it, so a file can name any of them.

To add an instrument, add a package here like trex/, and import and export its class below.
"""

from .trex.model import TRex

__all__ = ["TRex"]
