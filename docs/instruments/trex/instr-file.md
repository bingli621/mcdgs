# McStas instrument file

This is the McStas instrument file the model generates, at the default
[control parameters](details.md#control-parameters). The chopper timings, `nu_<name>` and
`delay_<name>`, are the defaults of its instrument parameters.
[`mcdgs.diagnostics.trex.overview`](../../reference/diagnostics.md) writes it to
`mcstas_input/TRex/` and copies it here; `pixi run docs-mcstas` runs it. The tests
check that this copy matches the current code. [Download it](figures/TRex.instr), or expand it here:

??? abstract "TRex.instr"

    ```c
    --8<-- "docs/instruments/trex/figures/TRex.instr"
    ```
