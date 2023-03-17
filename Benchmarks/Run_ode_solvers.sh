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
echo "Testing all models"
${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Test_all
echo "Testing random parameters"
${runJulia} Benchmarks/ODE_solvers/Benchmark_solvers.jl Test_random_parameter