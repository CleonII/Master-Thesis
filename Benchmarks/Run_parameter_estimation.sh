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

if [ $1 == "Brannmark_JBC2010" ];then
    ${runJulia} ${pathBenchmarkScript} Brannmark_JBC2010 $2 $3
fi

if [ $1 == "Bruno_JExpBot2016" ];then
    ${runJulia} ${pathBenchmarkScript} Bruno_JExpBot2016 $2 $3
fi

if [ $1 == "Weber_BMC2015" ];then
    ${runJulia} ${pathBenchmarkScript} Weber_BMC2015 $2 $3
fi

if [ $1 == "Zheng_PNAS2012" ];then
    ${runJulia} ${pathBenchmarkScript} Zheng_PNAS2012 $2 $3
fi

if [ $1 == "Sneyd_PNAS2002" ];then
    ${runJulia} ${pathBenchmarkScript} Sneyd_PNAS2002 $2 $3
fi

if [ $1 == "Lucarelli_CellSystems2018" ];then
    ${runJulia} ${pathBenchmarkScript} Lucarelli_CellSystems2018 $2 $3
fi

if [ $1 == "Schwen_PONE2014" ];then
    ${runJulia} ${pathBenchmarkScript} Schwen_PONE2014 $2 $3
fi

if [ $1 == "Elowitz_Nature2000" ];then
    ${runJulia} ${pathBenchmarkScript} Elowitz_Nature2000 $2 $3
fi

if [ $1 == "Crauste_CellSystems2017" ];then
    ${runJulia} ${pathBenchmarkScript} Crauste_CellSystems2017 $2 $3
fi

if [ $1 == "Isensee_JCB2018" ];then
    ${runJulia} ${pathBenchmarkScript} Isensee_JCB2018 $2 $3
fi

if [ $1 == "Borghans_BiophysChem1997" ];then
    ${runJulia} ${pathBenchmarkScript} Borghans_BiophysChem1997 $2 $3
fi

if [ $1 == "Okuonghae_ChaosSolitonsFractals2020" ];then
    ${runJulia} ${pathBenchmarkScript} Okuonghae_ChaosSolitonsFractals2020 $2 $3
fi

if [ $1 == "Oliveira_NatCommun2021" ];then
    ${runJulia} ${pathBenchmarkScript} Oliveira_NatCommun2021 $2 $3
fi

if [ $1 == "Rahman_MBS2016" ];then
    ${runJulia} ${pathBenchmarkScript} Rahman_MBS2016 $2 $3
fi


exit 0