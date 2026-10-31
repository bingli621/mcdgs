"""Griffe extension of the documentation site: Sphinx roles in docstrings as inline code.

The docstrings refer to each other with roles such as ``:func:`.params.param```. The site renders
them as Markdown, which has no roles and would show them as they are written.
"""

import re
from typing import Any

import griffe

ROLE = re.compile(r":(?:attr|class|func|meth):`\.?([^`]+)`")


class SphinxRoles(griffe.Extension):
    def on_instance(self, *, obj: griffe.Object, **kwargs: Any) -> None:
        if obj.docstring:
            obj.docstring.value = ROLE.sub(r"`\1`", obj.docstring.value)
