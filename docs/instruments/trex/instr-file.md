# McStas instrument file

This is the McStas instrument file the model generates, at the default
[control parameters](details.md#control-parameters). The chopper timings, `nu_<name>` and
`delay_<name>`, are the defaults of its instrument parameters.
[`diagnosis/trex/overview.py`](https://github.com/bingli621/mcdgs/blob/main/diagnosis/trex/overview.py)
writes it to `mcstas_input/TRex/` and prints it. [Download it](files/TRex.instr), or expand it
here:

??? abstract "TRex.instr"

    ```c
    --8<-- "docs/instruments/trex/files/TRex.instr"
    ```
