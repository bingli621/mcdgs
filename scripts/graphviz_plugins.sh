# Activation script of the pixi feature mcstas: register the plugins of graphviz, for
# sc.show_graph. conda-forge's graphviz does it in a post-link script, which pixi does not run, so
# without this `dot` has no layout engine. Only the first activation after an install runs it.
if ! ls "$CONDA_PREFIX"/lib/graphviz/config* > /dev/null 2>&1 && [ -x "$CONDA_PREFIX/bin/dot" ]; then
    "$CONDA_PREFIX/bin/dot" -c || true
fi
