# benchmarks

Benchmarks against [astropy](https://github.com/astropy/astropy)

## setup

1. (optional) instantiate python virtual environment
```
$ python -m virtualenv venv
$ source venv/bin/activate
```
2. install python dependencies
```
$ pip install -r bench/requirements.txt
```
3. instantiate julia environment and develop local project
```
$ julia --project=bench -e 'using Pkg; Pkg.instantiate(); Pkg.develop(path=pwd())'
```
or, from within the Julia REPL
```julia
pkg> activate bench
pkg> instantiate
pkg> dev .
```
4. ensure PyCall.jl is linked to appropriate python
```
$ PYTHON=$(which python) julia --project=bench -e 'using Pkg; Pkg.build("PyCall")`
```
or, from within the Julia REPL
```julia
shell> which python
julia> ENV["PYTHON"] = "/path/to/python"
pkg> build PyCall
```

## usage

To run the benchmark, first make sure you have appropriately linked the version of python you are using with PyCall, as shown in the setup above. Once PyCall is linked, simply activating the Julia environment and executing `benchmark.jl` will run the benchmark.

```
$ julia --project=bench bench/benchmark.jl
```
or, from within the Julia REPL
```julia
pkg> activate bench
julia> include("bench/benchmark.jl")
```
