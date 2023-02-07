# $1 - Model name 
# $2 - Number of multistarts (Integer)
# $3 - List of optmizers to try with the Julia parameter estimation 

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
pathBenchmarkScript="./Benchmarks/Parameter_estimation/Run_benchmark.jl"

if [ $1 == "Boehm_JProteomeRes2014" ];then
    ${runJulia} ${pathBenchmarkScript} Boehm_JProteomeRes2014 $2 $3
fi

exit 0