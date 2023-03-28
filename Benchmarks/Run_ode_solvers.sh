# Script must be run from the Root directory to correctly load Julia environment
currentDir=${PWD##*/}
if [ ! $currentDir == "Master-Thesis" ]; then
    >&2 echo "Error : Script must be run from directory Master-thesis to set it paths correctly"
    exit 1
fi

# Neede to load Conda environment with Fides 
eval "$(conda shell.bash hook)"
conda activate PeTab

runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"

if [ $1 == "Test_all" ];then
    echo "Testing all models"
    ${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Test_all
fi


if [ $1 == "Test_random_parameters" ];then
    echo "Testing random parameter for a subset of models"
    ${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Test_random_parameter
fi


if [ $1 == "Test_random_parameters_big_models" ];then
    echo "Testing random parameter for big models"
    ${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Large_models_random_p model_Bachmann_MSB2011
    ${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Large_models_random_p model_Lucarelli_CellSystems2018
    ${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Large_models_random_p model_Isensee_JCB2018
    ${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Large_models_random_p Smith_BMCSystBiol2013
    ${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Large_models_random_p model_Chen_MSB2009
fi
