# QuantExQASM.jl

Circuit composition and OpenQASM translation package for QuantEx project.
Creates a backend-agnostic algorithm and circuit composition framework,
which can subsequently be used to generate OpenQASM, or tie directly to 
other quantum assembly languages.

Integrates with PicoQuant.jl to perform full end-to-end demonstrations
of quantum circuit simulation as part of the QuantEx project.

## Installation
QuantExQASM.jl is activated in the same manner as PicoQuant.jl, where the 
following instructions will mirror the details provided there.

The prototype comes in the form of a Julia package which is targeted to versions
of Julia from v1 on. Binaries and source for this can be downloaded from
[https://julialang.org/](https://julialang.org/).

Once installed, from the Julia REPL prompt navigate to the QuantExQASM folder
and activate the environment, instantiate it and then build QuantExQASM.
This should install dependencies specified in the `Project.toml` and
`Manifest.toml` files as well as carry out any package specific build tasks
detailed in `deps/build.jl`. To use a custom python environment see the section
below on using different python environments.

```
]activate .
]instantiate
]build QuantExQASM
```

## Running the unittests

Unittests can be run from the QuantExQASM root folder with

```
julia --project=. tests/runtests.jl
```

## Building the documentation

This package uses Documenter.jl to generate html documentation from the sources.
To build the documentation, run the make.jl script from the docs folder.

```
cd docs && julia make.jl
```

The documentation will be placed in the build folder and can be hosted locally
by starting a local http server with

```
cd build && python3 -m http.server
```

