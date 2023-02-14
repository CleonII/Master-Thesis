# $1 - Model name 
# $2 - Number of multistarts (Integer)
# $3 - List of optmizers to try with the Julia parameter estimation 

# Script must be run from the Root directory to correctly load Julia environment
currentDir=${PWD##*/}
if [ ! $currentDir == "Master-Thesis" ]; then
    >&2 echo "Error : Script must be run from directory Master-thesis to set it paths correctly"
    exit 1
fi

# load same venv we installed packages before, no need for second anaconda environment, next two lines need to be manually adjusted with the environment used before for installing packages and the julia binary of version 1.8.5
source /home/sgrein/julia_benchmark/pypesto_benchmark/venv/bin/activate
runJulia="/home/sgrein/julia/julia --project=. --threads=1"
pathBenchmarkScript="./Benchmarks/Parameter_estimation/Run_benchmark.jl"

if [ $1 == "Boehm_JProteomeRes2014" ];then
    ${runJulia} ${pathBenchmarkScript} Boehm_JProteomeRes2014 $2 $3
fi

if [ $1 == "Fiedler_BMC2016" ];then
    ${runJulia} ${pathBenchmarkScript} Fiedler_BMC2016 $2 $3
fi

if [ $1 == "Fujita_SciSignal2010" ];then
    ${runJulia} ${pathBenchmarkScript} Fujita_SciSignal2010 $2 $3
fi

if [ $1 == "Bachmann_MSB2011" ];then
    ${runJulia} ${pathBenchmarkScript} Bachmann_MSB2011 $2 $3
fi

if [ $1 == "Beer_MolBioSystems2014" ];then
    ${runJulia} ${pathBenchmarkScript} Beer_MolBioSystems2014 $2 $3
fi

exit 0
