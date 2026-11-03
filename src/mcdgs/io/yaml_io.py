"""YAML parameters files; the only module that knows YAML.

A parameters file holds the nested dict of :mod:`.dict_io` (whose docstring shows the layout),
for any :class:`.Instrument`.

- :func:`load` lays a file over the default instrument of the class it names, so a file need
  only hold what differs. Every file must name its ``instrument:``.
- :func:`save` writes a complete file.

Each instrument keeps its default saved as a reference of how a parameters file looks; for
T-Rex, instruments/trex/params.yaml (``PARAMS_FILE``).

Files are UTF-8 (a unit such as ``Å`` is written as it is), in the order of the dict, with the
innermost lists and dicts on one line (``at: [0, 0, 1.0]``,
``lambda0: {value: 2.0, unit: Å, comment: central wavelength}``).

A number in exponent form reads as a number whether or not it has a dot: ``1e-3`` as well as
``1.0e-3``. Plain YAML 1.1, which PyYAML reads, takes ``1e-3`` for a string, which a component
would take for a McStas expression.
"""

import re
from pathlib import Path
from typing import Any

import yaml

from mcdgs.core.base import Instrument

from .dict_io import from_overrides, to_dict


class _Loader(yaml.SafeLoader):
    """PyYAML's safe loader, reading ``1e-3`` as a float, as YAML 1.2 does."""


_Loader.add_implicit_resolver(  # copies SafeLoader's resolvers first: SafeLoader is not changed
    "tag:yaml.org,2002:float",
    re.compile(r"^[-+]?[0-9][0-9_]*(\.[0-9_]*)?[eE][-+]?[0-9]+$"),
    list("-+0123456789"),
)


def read_yaml(path: str | Path) -> dict[str, Any]:
    """Read a nested dict from a YAML file.

    Parameters
    ----------
    path:
        The YAML file.

    Returns
    -------
    :
        The dict, empty if the file is empty.
    """
    with open(path, encoding="utf-8") as f:
        return yaml.load(f, Loader=_Loader) or {}  # _Loader is a SafeLoader: no code is run


def write_yaml(data: dict[str, Any], path: str | Path) -> None:
    """Write a nested dict to a YAML file.

    Keys keep their order, and the innermost lists and dicts are written on one line, however
    long.

    Parameters
    ----------
    data:
        The dict.
    path:
        The YAML file, overwritten if it exists.
    """
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(
            data,
            f,
            sort_keys=False,
            default_flow_style=None,
            allow_unicode=True,
            width=float("inf"),  # never wrapped
        )


def load(path: str | Path, **kwargs: Any) -> Instrument:
    """Load an instrument from a parameters file.

    The file is laid over the default instrument of the class it names (see
    :func:`.from_overrides`), so it need only hold what differs.

    Parameters
    ----------
    path:
        The parameters file.
    **kwargs:
        Passed to the class to create the default instrument; for T-Rex, ``sample`` and
        ``params``.

    Returns
    -------
    :
        The instrument.

    Raises
    ------
    ValueError
        As :func:`.from_dict`.
    """
    return from_overrides(read_yaml(path), **kwargs)


def save(instrument: Instrument, path: str | Path) -> None:
    """Save the complete parameters file of an instrument, as YAML.

    The file names the class and holds every component, so :func:`load` gives the instrument
    back.

    Parameters
    ----------
    instrument:
        The instrument.
    path:
        The file, overwritten if it exists.
    """
    write_yaml(to_dict(instrument), path)
