# Code layout

The code is the Python package `mcdgs`, in `src/mcdgs/`, one folder per layer.
[Core concepts](concepts.md) explains what the parts are.

## Layers

Dependencies point inwards, so an instrument never depends on how it is simulated:

| Package | Role |
| ------- | ---- |
| `core` | `Instrument`, the components, parameters and geometry. It imports nothing else from `mcdgs`. |
| `instruments`, `samples` | The instruments and ready-made crystals, built on `core`. |
| `mcstas` | An adapter that works with any instrument: builds it, runs it, reads the output. |
| `tof` | An adapter that works with any instrument: builds its chopper cascade as a tof model, gives its disks as scippneutron disk choppers, and draws its chopper acceptance diagram. |
| `io` | An adapter that saves any instrument to a YAML parameters file and loads it back. |
| `reduction` | Reserved for data reduction. It holds no code yet. |
| `diagnostics` | Setups that study one part of an instrument, such as T-Rex up to guide 4, and what reads their McStas and tof results. One module per setup, numbered per instrument, such as `trex_01_source`, built on `instruments`; nothing else imports it. |

Each box is a package, with the main things it offers, and each dashed arrow points to a package it
imports from:

```mermaid
classDiagram
    direction BT
    class core {
        <<package>>
        Instrument
        Component
        InstrParam
    }
    class instruments {
        <<package>>
        TRex
    }
    class samples {
        <<package>>
        aluminum()
    }
    class mcstas {
        <<package>>
        build_mcstas()
        Simulation
        read_pixels()
    }
    class tof {
        <<package>>
        build_tof()
        acceptance_diagram()
        time_distance_diagram()
    }
    class io {
        <<package>>
        save()
        load()
    }
    class reduction {
        <<package>>
    }
    class diagnostics {
        <<package>>
        trex_00_overview to trex_06_time_distance_diagram
    }
    instruments ..> core
    samples ..> core
    mcstas ..> core
    tof ..> core
    io ..> core
    diagnostics ..> instruments
    diagnostics ..> mcstas
    diagnostics ..> tof
```
