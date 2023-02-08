# Script must be run from the Root directory to correctly load Julia environment
currentDir=${PWD##*/}
if [ ! $currentDir == "Master-Thesis" ]; then
    >&2 echo "Error : Script must be run from directory Master-thesis to set it paths correctly"
    exit 1
fi

# Neede to load Conda environment with Fides 
eval "$(conda shell.bash hook)"
conda activate PeTab

if [ $1 == "Test_flags" ];then
    runJulia1="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. -O1 --threads=1"
    runJulia3="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. -O3 --threads=1"
    pathBenchmarkScript="./Benchmarks/Cost_grad_hess/Cost_grad_hess.jl"
    
    ${runJulia1} ${pathBenchmarkScript} Test_flags O1
    ${runJulia3} ${pathBenchmarkScript} Test_flags O3
fi


if [ $1 == "Fix_parameters" ];then
    runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"
    pathBenchmarkScript="./Benchmarks/Cost_grad_hess/Cost_grad_hess.jl"
    
    # For Bachman 
    for i in {1..26}
    do
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Bachman $i
    done 
    for i in {1..70}
    do
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Lucarelli $i
    done 
fi


exit 0