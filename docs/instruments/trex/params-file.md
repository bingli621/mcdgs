# Parameters file

This is the default T-Rex saved as a parameters file, `save(TRex(), PARAMS_FILE)`, the
reference of how a T-Rex file looks. `PARAMS_FILE` is `src/mcdgs/instruments/trex/params.yaml`,
and a test checks that it stays the instrument the code makes.

It holds the control parameters, then the McStas instrument parameters they give, then one
section per component in beam order: its McStas component, its parameters, where it sits (`at`,
`rotated`, `relative`) and its position in the lab frame. The instrument parameters and the
positions are for reading only, and loading the file does not read them. A file you load need
only hold what differs from this one; see [Parameters files](../../reference/io.md).
[View it on GitHub](https://github.com/bingli621/mcdgs/blob/main/src/mcdgs/instruments/trex/params.yaml),
or expand it here:

??? abstract "params.yaml"

    ```yaml
    --8<-- "src/mcdgs/instruments/trex/params.yaml"
    ```
