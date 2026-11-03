"""The T-Rex instrument.

``PARAMS_FILE`` is the default T-Rex saved, ``save(TRex(), PARAMS_FILE)``: a reference of how
a parameters file looks. A test checks that it stays the instrument the code makes; after a
change to the model, write it again the same way.
"""

from pathlib import Path

PARAMS_FILE = Path(__file__).with_name("params.yaml")
