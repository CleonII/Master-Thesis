# About

This file describes how to run the parameter estimation benchmark for the supported models (list of supported models below). 

## Setting up correct Julia environment

Assuming that you have cloned this directory the correct Julia environment for Linux can be set up by can be set up via a local installation;

1. Download the folder with [Julia 1.8.2](https://julialang.org/downloads/). 
2. Extract the [tar-archive](https://xkcd.com/1168/).
3. In the terminal go to the root directory of the Julia PeTab importer (where the Manifest.toml and Project.toml are) and launch Julia via in the terminal by: `path_julia_folder/bin/julia`
4. In the Julia prompt press **]** which will activate the [package manager](https://pkgdocs.julialang.org/v1/environments/). In the pacakge manager : 
    1. Write `activate .` (activate current environment)
    2. Write `instantiate` (will download and precompile all the required packages)

Now all the correct Julia packages are installed. To use the correct environment when running a script from the root-directory of the repository simply use the project flag;

~~~
path_julia_folder/bin/julia --project=. file.jl
~~~

In VScode you can also select the environment.

## Setting up correct Python environment (for Fides)

I have created a `.yml` which installs Fides and its dependencies. Given a [Conda installation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html) the Python packages can from the following terminal command;

~~~
conda env create -f PeTab.yml -n your_env_name
~~~

## Running ODE solver benchmark

To run the benchmark Boehm, Borghans, Bruno, Crauste, Elowitz, Fiedler, and Fujita models in the project root-directory run:

~~~
path_julia_folder/bin/julia --project=. Benchmarks/ODE_solvers/Benchmark_solvers.jl compare_fabian
~~~

The ODE solvers tested are KenCarp4, FBDF, QNDF, Rosenbrock23, TRBDF2, RadauIIA5, Rodas4, CVODE_BDF, "Rodas5" using two tolerance options; (abstol=1e-8, reltol=1e-8) and (abstol=1e-16, reltol=1e-8). For each solver and tolerance option we measure the run-time (wall-time) three times (three repetitions).

The parameter values used for the ODE benchmark are the same as the nominal values in the PEtab files from [here](https://github.com/Benchmarking-Initiative/Benchmark-Models-PEtab), and the ODE:s are solved for all experimental conditions for which there are measurement data. For each model the benchmark pipeline computes a high accuracy solution of the ODE using BigFloatm abstol=1e-16, and reltol=1e-16 which can take some time.

The result from the parameter estimation is automatically stored in *Intermediate/Benchmarks/ODE_solvers/file_name.csv* in a [Tidy-format](https://www.jstatsoft.org/article/view/v059i10).

## Running parameter estimation benchmark

Currently the benchmark can run for the **Boehm**, **Bachmann**, **Brännmark** and **Fujita** models. To run the benchmark in terminal from the root project directory run:

~~~
path_julia_folder/bin/julia --project=. Benchmarks/Parameter_estimation/Run_benchmark.jl model_name
~~~

So to run the benchmark for the Boehm model simply run;

~~~
path_julia_folder/bin/julia --project=. Benchmarks/Parameter_estimation/Run_benchmark.jl Boehm
~~~

If you want to use Fides optimizers you must also make a small change in the *Benchmarks/Parameter_estimation/Run_benchmark.jl* file. Specificailly, you must add the following line before the call to the benchmark function:

~~~
loadFidesFromPython("path_to_python_for_conda_env")
~~~

Where `path_to_python_for_conda_env` is the path to python executable for the Conda environment with the Fides installation. I know this is a bit suboptimal, but using Python within Julia can be tricky easy :).

The result from the parameter estimation is automatically stored in *Intermediate/Benchmarks/Parameter_estimation/Model_name/Benchmark_result.csv* in a [Tidy-format](https://www.jstatsoft.org/article/view/v059i10).

### Changing parameter estimation options

When running the parameter estimation you can choose ODE-solver, from a list of algorithms, and ODE-solver tolerances.

Consider the call to the **Boehm** model;

```
algsTest = [:IpoptAutoHess, :IpoptBlockAutoDiff, :IpoptLBFGS, :OptimIPNewtonAutoHess, :OptimIPNewtonBlockAutoDiff, :OptimLBFGS, :FidesAutoHess, :FidesBlockAutoHess]
    
# This is just for setting up the model
dirModel = pwd() * "/Intermediate/PeTab_models/model_Boehm_JProteomeRes2014/"
peTabModel = setUpPeTabModel("model_Boehm_JProteomeRes2014", dirModel)

tol = 1e-9 # Rel- and Abs-tol for ODE solver
benchmarkParameterEstimation(peTabModel, QNDF(), "QNDF", tol, 1000, algList=algsTest)

```

Here `QNDF()` is the ODE solver. Basically any Julia ODE solver can be used. So instead writing `Rodas5()` uses Rodas5 (if you change ODE-solver also change the string argument in the call to get everything correct in the result file). 

algsTest is a list with optimization algorithms. Supported algorithms are;

| Julia Syntax      | Description |
| ----------- | ----------- |
| :IpoptAutoHess      | Ipopt using full hessian (via autodiff)|
| :IpoptBlockAutoDiff   | Ipopt using block approximated hessian (via autodiff)        |
| :IpoptLBFGS   | Ipopt using L-BFGS hessian approximation        |
| :OptimIPNewtonAutoHess   | Optim.jl interior point Newton full hessian (via autodiff)   |
| :OptimIPNewtonBlockAutoDiff   | Optim.jl interior point Newton via bloack approximated hessian (via autodiff)        |
| :OptimLBFGS   | Optim.jl L-BFGS method|
| :FidesAutoHess      | Fides using full hessian (via autodiff)|
| :FidesBlockAutoHess   | Fides using block approximated hessian (via autodiff)        |