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
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Bachman $i 321
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Bachman $i 123
    done 
    for i in {1..70}
    do
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Lucarelli $i 123
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Lucarelli $i 321
    done 
fi


if [ $1 == "Test_chunks_shuffle" ];then
    runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"
    pathBenchmarkScript="./Benchmarks/Cost_grad_hess/Cost_grad_hess.jl"
    
    # For Bachman 
    ${runJulia} ${pathBenchmarkScript} Test_chunks_random_p Bachman Shuffle 123
    ${runJulia} ${pathBenchmarkScript} Test_chunks_random_p Bachman Shuffle 321
fi


exit 0