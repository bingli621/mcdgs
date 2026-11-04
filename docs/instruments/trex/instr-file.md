# McStas instrument file

This is the McStas instrument file the model generates, at the default
[control parameters](details.md#control-parameters). The chopper timings, `nu_<name>` and
`delay_<name>`, are the defaults of its instrument parameters.
[`notebooks/trex/overview.ipynb`](https://github.com/bingli621/mcdgs/blob/main/notebooks/trex/overview.ipynb)
writes it to `mcstas_input/TRex/` and prints it; `pixi run docs-mcstas` runs it again. The tests
check that this copy matches the current code. [Download it](files/TRex.instr), or expand it here:

??? abstract "TRex.instr"

    ```c
    --8<-- "docs/instruments/trex/files/TRex.instr"
    ```
