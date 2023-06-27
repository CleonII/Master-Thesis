# Script must be run from the Root directory to correctly load Julia environment
currentDir=${PWD##*/}
if [ ! $currentDir == "Master-Thesis" ]; then
    >&2 echo "Error : Script must be run from directory Master-thesis to set it paths correctly"
    exit 1
fi

# Neede to load Conda environment with Fides 
eval "$(conda shell.bash hook)"
conda activate PeTab


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
    for i in {1..32}
    do
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Isensee $i 123
        ${runJulia} ${pathBenchmarkScript} Fix_parameters Isensee $i 321
    done 
fi

if [ $1 == "Test_chunks" ];then
    runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"
    pathBenchmarkScript="./Benchmarks/Cost_grad_hess/Cost_grad_hess.jl"
    
    ${runJulia} ${pathBenchmarkScript} Test_chunks model_Boehm_JProteomeRes2014
    ${runJulia} ${pathBenchmarkScript} Test_chunks model_Bachmann_MSB2011
    ${runJulia} ${pathBenchmarkScript} Test_chunks model_Lucarelli_CellSystems2018
    ${runJulia} ${pathBenchmarkScript} Test_chunks model_Isensee_JCB2018
fi


if [ $1 == "Test_adjoint_random" ];then
    runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"
    pathBenchmarkScript="./Benchmarks/Cost_grad_hess/Cost_grad_hess.jl"
    
    ${runJulia} ${pathBenchmarkScript} Test_adjoint_random_p model_Boehm_JProteomeRes2014
    ${runJulia} ${pathBenchmarkScript} Test_adjoint_random_p model_Bachmann_MSB2011
    ${runJulia} ${pathBenchmarkScript} Test_adjoint_random_p model_Lucarelli_CellSystems2018
    ${runJulia} ${pathBenchmarkScript} Test_adjoint_random_p Smith_BMCSystBiol2013
    ${runJulia} ${pathBenchmarkScript} Test_adjoint_random_p model_Chen_MSB2009
fi


if [ $1 == "Gradient_cost_small_models" ];then
    runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"
    pathBenchmarkScript="./Benchmarks/Cost_grad_hess/Cost_grad_hess.jl"
    
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Beer_MolBioSystems2014
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Boehm_JProteomeRes2014
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Bachmann_MSB2011
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Bruno_JExpBot2016
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Crauste_CellSystems2017
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Elowitz_Nature2000
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Fiedler_BMC2016
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Fujita_SciSignal2010
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Lucarelli_CellSystems2018
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Sneyd_PNAS2002
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Weber_BMC2015
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Brannmark_JBC2010
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Isensee_JCB2018
    ${runJulia} ${pathBenchmarkScript} Gradient_cost_small_models model_Zheng_PNAS2012
fi


if [ $1 == "Hessian_cost_small_models" ];then
    runJulia="/home/sebpe/julia-1.8.5-linux-x86_64/julia-1.8.5/bin/julia --project=. --threads=1"
    pathBenchmarkScript="./Benchmarks/Cost_grad_hess/Cost_grad_hess.jl"
    
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Boehm_JProteomeRes2014
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Beer_MolBioSystems2014
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Bachmann_MSB2011
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Bruno_JExpBot2016
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Crauste_CellSystems2017
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Elowitz_Nature2000
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Fiedler_BMC2016
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Fujita_SciSignal2010
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Lucarelli_CellSystems2018
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Sneyd_PNAS2002
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Weber_BMC2015
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Brannmark_JBC2010
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Isensee_JCB2018
    ${runJulia} ${pathBenchmarkScript} Hessian_cost_small_models model_Zheng_PNAS2012
fi

exit 0