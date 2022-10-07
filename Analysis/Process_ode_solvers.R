library(tidyverse)
library(ggthemes)


# General plotting parameters (plot using theme-tufte)
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
my_theme <- theme_hc(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
my_minimal <- theme_minimal(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=18))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0


get_solver_info <- function(solver_name)
{
  solver_list_short <- c("Vern6", 
                         "Vern7", 
                         "Vern8", 
                         "Tsit5", 
                         "DP5", 
                         "DP8", 
                         "Feagin14", 
                         "VCABM", 
                         "Rosenbrock23", 
                         "TRBDF2", 
                         "Rodas4", 
                         "Rodas4P", 
                         "Rodas4P2", 
                         "Rodas5", 
                         "QNDF", 
                         "FBDF", 
                         "Trapezoid", 
                         "ARKODE", 
                         "Kvaerno5", 
                         "RadauIIA3", 
                         "RadauIIA5", 
                         "Tsit5_Rosenbrock23", 
                         "Vern7_Rodas5", 
                         "Vern9_Rodas4P", 
                         "Vern9_Rodas5", 
                         "lsoda", 
                         "CVODE_BDF_N", # Dense
                         "CVODE_BDF_L", 
                         "CVODE_BDF_G", 
                         "CVODE_Admas_D", 
                         "CVODE_Admas_L", 
                         "ARKODE_E4", 
                         "ARKODE_E8", 
                         "ARKODE_I3", 
                         "ARKODE_I5",
                         "dopri5", 
                         "dop853", 
                         "radau", 
                         "radau5", 
                         "rodas", 
                         "auto", 
                         "nonstiff", 
                         "stiff")
  
  solver_type <- c("nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", 
                   "stiff", "stiff", "stiff", "stiff", "stiff", "stiff", "stiff", "stiff","stiff", "stiff", "stiff", 
                   "stiff", "stiff", "composite", "composite", "composite", "composite", "composite", "stiff", "stiff", 
                   "stiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "stiff", "stiff", "nonstiff", "nonstiff", 
                   "stiff", "stiff", "stiff", "hint", "hint", "hint")
  solver_library <- c("OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "LSODA.jl", "Sundials", "Sundials", "Sundials", "Sundials", "Sundials", "Sundials", "Sundials", 
                      "Sundials", "Sundials", "ODEInterface", "ODEInterface", "ODEInterface", "ODEInterface", 
                      "ODEInterface", "DifferentialEquations", "DifferentialEquations", "DifferentialEquations")
  
  i_solver <- which(solver_name == solver_list_short)[1]
  return(c(solver_type[i_solver], solver_library[i_solver]))
}


read_data <- function(path_file)
{

  
  data_raw <- read_csv(path_file, col_types = cols()) 

  
  solver_list <- unique(data_raw$solver)
  solver_list_short <- c("Vern6", 
                         "Vern7", 
                         "Vern8", 
                         "Tsit5", 
                         "DP5", 
                         "DP8", 
                         "Feagin14", 
                         "VCABM", 
                         "Rosenbrock23", 
                         "TRBDF2", 
                         "Rodas4", 
                         "Rodas4P", 
                         "Rodas4P2", 
                         "Rodas5", 
                         "QNDF", 
                         "FBDF", 
                         "Trapezoid", 
                         "ARKODE", 
                         "Kvaerno5", 
                         "RadauIIA3", 
                         "RadauIIA5", 
                         "Tsit5_Rosenbrock23", 
                         "Vern7_Rodas5", 
                         "Vern9_Rodas4P", 
                         "Vern9_Rodas5", 
                         "lsoda", 
                         "CVODE_BDF_N", 
                         "CVODE_BDF_L", 
                         "CVODE_BDF_G", 
                         "CVODE_Admas_D", 
                         "CVODE_Admas_L", 
                         "ARKODE_E4", 
                         "ARKODE_E8", 
                         "ARKODE_I3", 
                         "ARKODE_I5",
                         "dopri5", 
                         "dop853", 
                         "radau", 
                         "radau5", 
                         "rodas", 
                         "auto", 
                         "nonstiff", 
                         "stiff")
  
  model_list <- c("model_Alkan_SciSignal2018.jl", "model_Bachmann_MSB2011.jl",
                "model_Beer_MolBioSystems2014.jl", "model_Bertozzi_PNAS2020.jl",
                "model_Blasi_CellSystems2016.jl", "model_Boehm_JProteomeRes2014.jl",
                "model_Borghans_BiophysChem1997.jl", "model_Brannmark_JBC2010.jl",
                "model_Bruno_JExpBot2016.jl", "model_Chen_MSB2009.jl", 
                "model_Crauste_CellSystems2017.jl", "model_Elowitz_Nature2000.jl", 
                "model_Fiedler_BMC2016.jl", "model_Fujita_SciSignal2010.jl",
                "model_Giordano_Nature2020.jl", "model_Isensee_JCB2018.jl",
                "model_Laske_PLOSComputBiol2019.jl", "model_Lucarelli_CellSystems2018.jl",
                "model_Okuonghae_ChaosSolitonsFractals2020.jl", "model_Oliveira_NatCommun2021.jl",
                "model_Perelson_Science1996.jl", "model_Rahman_MBS2016.jl",
                "model_Raimundez_PCB2020.jl", "model_SalazarCavazos_MBoC2020.jl",
                "model_Schwen_PONE2014.jl", "model_Sneyd_PNAS2002.jl",
                "model_Weber_BMC2015.jl", "model_Zhao_QuantBiol2020.jl",
                "model_Zheng_PNAS2012.jl")
  
  n_states <- c(36, 25, 4,  3, 16, 8,  3,  9,  7, 500,  5,  8,  6,  9, 13, 25, 41,  33,  9,  9, 4,  7,  22, 75, 11,  6,  7,  5, 15)
  n_params <- c(52, 37, 8, 21, 34, 9, 20, 20, 17, 196, 12, 18, 18, 20, 72, 69, 90, 106, 17, 33, 5, 24, 101, 50, 14, 19, 38, 37, 47) 
  model_list_sort <- model_list[order(n_states)]
  model_list_sort_param <- model_list[order(n_params)]
  
  solver_type <- c("nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", 
                   "stiff", "stiff", "stiff", "stiff", "stiff", "stiff", "stiff", "stiff","stiff", "stiff", "stiff", 
                   "stiff", "stiff", "composite", "composite", "composite", "composite", "composite", "stiff", "stiff", 
                   "stiff", "nonstiff", "nonstiff", "nonstiff", "nonstiff", "stiff", "stiff", "nonstiff", "nonstiff", 
                   "stiff", "stiff", "stiff", "hint", "hint", "hint")
  solver_library <- c("OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", "OrdinaryDiffEq", 
                      "LSODA.jl", "Sundials", "Sundials", "Sundials", "Sundials", "Sundials", "Sundials", "Sundials", 
                      "Sundials", "Sundials", "ODEInterface", "ODEInterface", "ODEInterface", "ODEInterface", 
                      "ODEInterface", "DifferentialEquations", "DifferentialEquations", "DifferentialEquations")
  
  
  n_rows <- dim(data_raw)[1]
  data_raw <- data_raw %>%
    mutate(solver_type = "c", solver_lib = "c", n_param=0, n_states=0, model_short="c", model_short_param="p")
  for(i in 1:n_rows){
    i_solver <- which(data_raw$solver[i] == solver_list)
    i_model <- which(data_raw$model[i] == model_list)
    i_mod_sorted <- which(data_raw$model[i] == model_list_sort)
    i_param_sorted <- which(data_raw$model[i] == model_list_sort_param)
    data_raw$solver[i] <- solver_list_short[i_solver] 
    data_raw$solver_type[i] <- solver_type[i_solver]
    data_raw$solver_lib[i] <- solver_library[i_solver]
    data_raw$n_param[i] <- n_params[i_model]
    data_raw$n_states[i] <- n_states[i_model]
    data_raw$model_short[i] <- str_c("M", as.character(i_mod_sorted))
    data_raw$model_short_param[i] <- str_c("P", as.character(i_param_sorted))
  }
  
  return(data_raw)
}


calc_winner <- function(data, tol=1e-6)
{

  data_winner <- tibble()
  data <- data %>% filter(abstol == tol)
  model_list <- unique(data$model)
  
  for(i in 1:length(model_list)){
    
    
    data_mod_i <- data %>%
      filter(model == model_list[i]) %>%
      mutate(solver = as.factor(solver)) %>%
      group_by(solver) %>%
      summarise(median_time = median(runTime, na.rm = T), 
                median_sq = median(sqDiff, na.rm = T), 
                solver_type = median(solver_type), 
                solver_lib = median(solver_lib), 
                n_param = median(n_param), 
                n_states = median(n_states), 
                model_short = median(model_short), 
                model_short_param = median(model_short_param))
    # Add order
    data_mod_i <- data_mod_i[sort(data_mod_i$median_time, na.last = T, index.return=T, method="radix")$ix, ]
    data_mod_i <- data_mod_i %>% mutate(rank_time = 1:length(median_time))
    data_mod_i <- data_mod_i[sort(data_mod_i$median_sq, na.last = T, index.return=T, method="radix")$ix, ]
    data_mod_i <- data_mod_i %>% mutate(rank_sq = 1:length(median_time))
    data_mod_i <- data_mod_i %>% mutate(model = model_list[i])
    
    data_winner <- data_winner %>%
      bind_rows(data_mod_i)
  }
  
  return(data_winner)
}


calc_winner_category <- function(data, tol=1e-6, sq_res=F)
{
  
  data_winner <- tibble()
  model_list <- unique(data$model)
  solver_type_list <- unique(data$solver_type)
  data <- data %>% filter(abstol == tol)
  
  for(i in 1:length(model_list)){
    for(j in 1:length(solver_type_list)){
      
      data_mod_i <- data %>%
        filter(model == model_list[i]) %>%
        filter(solver_type == solver_type_list[j]) %>%
        mutate(solver = as.factor(solver)) %>%
        group_by(solver) %>%
        summarise(median_time = median(runTime, na.rm = T), 
                  median_sq = median(sqDiff, na.rm = T), 
                  solver_type = median(solver_type), 
                  solver_lib = median(solver_lib), 
                  n_param = median(n_param), 
                  n_states = median(n_states), 
                  model_short = median(model_short))
      if(sq_res == FALSE){
        data_mod_i <- data_mod_i[which.min(data_mod_i$median_time)[1], ]
      }else{
        data_mod_i <- data_mod_i[which.min(data_mod_i$median_sq)[1], ]
      }
      
      data_winner <- data_winner %>%
        bind_rows(data_mod_i)
    }
  }
  
  return(data_winner)
}


calc_score <- function(data, crit="scale", exclude_fail=FALSE, tol=1e-6)
{
  model_list <- unique(data$model)
  data <- data %>%
    filter(abstol == tol)
  
  if(crit == "scale"){
    data_use <- data %>%
      group_by(model) %>%
      summarise(mean_time = mean(runTime, na.rm = T), 
                mean_sq = mean(sqDiff, na.rm = T), 
                sd_time = sd(runTime, na.rm = T),
                sd_sq = sd(sqDiff, na.rm = T),
                max_time = max(runTime, na.rm = T),
                max_sq = max(sqDiff, na.rm = T))
  }else if(crit == "mad"){
    data_use <- data %>%
      group_by(model) %>%
      summarise(mean_time = median(runTime, na.rm = T), 
                mean_sq = median(sqDiff, na.rm = T), 
                sd_time = mad(runTime, na.rm = T, constant = 1),
                sd_sq = mad(sqDiff, na.rm = T, constant = 1),
                max_time = max(runTime, na.rm = T),
                max_sq = max(sqDiff, na.rm = T))
  }
  
  data_rel_score <- data %>%
    inner_join(data_use, by = "model") %>%
    mutate(rel_time = (runTime - mean_time) / sd_time, rel_sq = (sqDiff - mean_sq) / (sd_sq+1e-9), 
           max_time = (max_time - mean_time) / sd_time, max_sq = (max_sq - mean_sq) / sd_sq)
  
  if(exclude_fail == FALSE){
    data_rel_score$rel_sq[is.nan(data_rel_score$rel_sq)] <- data_rel_score$max_sq[is.nan(data_rel_score$rel_sq)]
    data_rel_score$rel_time[is.nan(data_rel_score$rel_time)] <- data_rel_score$max_time[is.nan(data_rel_score$rel_time)]
  }else{
    data_rel_score <- data_rel_score %>%
      drop_na()
  }
  
  data_rel_score_tmp <- data_rel_score %>%
    group_by(solver, solver_type, solver_lib) %>%
    summarise(rel_time = mean(rel_time, na.rm = T), 
              rel_sq = mean(rel_sq, na.rm = T)) %>%
    mutate(crit = crit, 
           exclude_fail = as.character(exclude_fail))
  
  data_rel_score_tmp <- data_rel_score_tmp[sort(data_rel_score_tmp$rel_time, na.last = T, index.return=T, method="radix")$ix, ]
  data_rel_score_tmp$rank_time <- 1:length(data_rel_score_tmp$rel_time)
  
  return(data_rel_score_tmp)
}


calc_fail <- function(data, tol=1e-6)
{

  model_list <- unique(data$model)
  solver_list <- unique(data$solver)
  data_fail <- tibble()
  data_tol <- data %>% filter(abstol == tol)
  for(i in 1:length(model_list)){
    mod_short <- (data %>% filter(model == model_list[i]))$model_short[1]
    for(j in 1:length(solver_list)){
      data_mod <- (data_tol %>%
                     filter(model == model_list[i] & solver == solver_list[j]))[1, ]
      
      if(is.na(data_mod$success) || is.nan(data_mod$sqDiff[1]) | data_mod$success == "FALSE"){
        
        solver_info <- get_solver_info(solver_list[j])
        data_tmp <- tibble(model = model_list[i], 
                           solver = solver_list[j],
                           n_param = data_mod$n_param, 
                           n_states = data_mod$n_states,
                           tol = tol, 
                           solver_type = solver_info[1], 
                           solver_lib = solver_info[2])
        data_fail <- data_fail %>% bind_rows(data_tmp)
        
      }
    }
  }
  
  return(data_fail)
}


plot_winner <- function(data, tol, sq_res=F)
{
  
  data_winner <- calc_winner(data, tol=tol)
  data_winner <- data_winner %>%
    mutate(model_short = factor(model_short, levels = c("M1", "M2", "M3", "M4", "M5", 
                                                        "M6", "M7", "M8", "M9", "M10", 
                                                        "M11", "M12", "M13", "M14", "M15", 
                                                        "M16", "M17", "M18", "M19", "M20", 
                                                        "M21", "M22", "M23", "M24", "M25", 
                                                        "M26", "M27", "M28")), 
           solver_lib = as.factor(solver_lib))
  
  if(sq_res == T){
    data_winner <- data_winner %>%
      rename("crit_plot" = "median_sq")
    ylab = "Median squared error"
  }else{
    data_winner <- data_winner %>%
      rename("crit_plot" = "median_time")
    ylab = "Median time [ns]"
  }
  
  # Plot data on stiff solvers 
  data_stiff <- data_winner %>%
    filter(solver_type == "stiff")
  
  solvers_col_map <- c("Rodas5" = cbPalette[2], "Rodas4P" = cbPalette[3], "Rosenbrock23" = cbPalette[4], 
                       "RadauIIA3" = cbPalette[5], "QNDF" = cbPalette[6])
  p_stiff <- ggplot(data_stiff, aes(x = model_short, crit_plot, group=solver)) + 
    geom_point(stat="summary", fun=sum, color = cbPalette[1]) + 
    stat_summary(fun=sum, geom="line", color = cbPalette[1]) +
    geom_point(data=filter(data_stiff, solver == "Rodas5"), stat="summary", fun=sum, mapping=aes(color = "Rodas5"), size=3.0) + 
    stat_summary(data=filter(data_stiff, solver == "Rodas5"), fun=sum, geom="line", mapping=aes(color = "Rodas5"), size=2.0) +
    geom_point(data=filter(data_stiff, solver == "Rodas4P"), stat="summary", fun=sum, mapping=aes(color = "Rodas4P"), size=3.0) + 
    stat_summary(data=filter(data_stiff, solver == "Rodas4P"), fun=sum, geom="line", mapping=aes(color = "Rodas4P"), size=2.0) +
    geom_point(data=filter(data_stiff, solver == "Rosenbrock23"), stat="summary", fun=sum, mapping=aes(color = "Rosenbrock23"), size=3.0) + 
    stat_summary(data=filter(data_stiff, solver == "Rosenbrock23"), fun=sum, geom="line", mapping=aes(color = "Rosenbrock23"), size=2.0) +
    geom_point(data=filter(data_stiff, solver == "RadauIIA3"), stat="summary", fun=sum, mapping=aes(color = "RadauIIA3"), size=3.0) + 
    stat_summary(data=filter(data_stiff, solver == "RadauIIA3"), fun=sum, geom="line", mapping=aes(color = "RadauIIA3"), size=2.0) +
    geom_point(data=filter(data_stiff, solver == "QNDF"), stat="summary", fun=sum, mapping=aes(color = "QNDF"), size=3.0) + 
    stat_summary(data=filter(data_stiff, solver == "QNDF"), fun=sum, geom="line", mapping=aes(color = "QNDF"), size=2.0) +
    scale_y_log10() + 
    labs(x = "Model (short name)", y = ylab, title = "Stiff solvers") +
    scale_color_manual(values = solvers_col_map, name = "Solver") + 
    my_minimal
  
  p_stiff_lib <- ggplot(data_stiff, aes(x = model_short, crit_plot, group=solver, color = solver_lib)) + 
    geom_point(stat="summary", fun=sum, size = 3.0, alpha=0.7) +
    stat_summary(fun=sum, geom="line", size = 1.5, alpha=0.7) +
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1], name = "Solver library") + 
    labs(x = "Model (short name)", y = ylab, title = "Library : Stiff solvers") +
    my_minimal
  
  data_non_stiff <-  data_winner %>%
    filter(solver_type == "nonstiff")
  
  solvers_col_map <- c("DP5" = cbPalette[2], "Tsit5" = cbPalette[3], "Vern7" = cbPalette[4])
  p_non_stiff <- ggplot(data_non_stiff, aes(x = model_short, crit_plot, group=solver)) + 
    geom_point(stat="summary", fun=sum, color = cbPalette[1]) + 
    stat_summary(fun=sum, geom="line", color = cbPalette[1]) + 
    geom_point(data=filter(data_non_stiff, solver == "DP5"), stat="summary", fun=sum, mapping=aes(color = "DP5"), size=3.0) + 
    stat_summary(data=filter(data_non_stiff, solver == "DP5"), fun=sum, geom="line", mapping=aes(color = "DP5"), size=2.0) +
    geom_point(data=filter(data_non_stiff, solver == "Tsit5"), stat="summary", fun=sum, mapping=aes(color = "Tsit5"), size=3.0) + 
    stat_summary(data=filter(data_non_stiff, solver == "Tsit5"), fun=sum, geom="line", mapping=aes(color = "Tsit5"), size=2.0) +
    geom_point(data=filter(data_non_stiff, solver == "Vern7"), stat="summary", fun=sum, mapping=aes(color = "Vern7"), size=3.0) + 
    stat_summary(data=filter(data_non_stiff, solver == "Vern7"), fun=sum, geom="line", mapping=aes(color = "Vern7"), size=2.0) +
    scale_y_log10() + 
    labs(x = "Model (short name)", y = ylab, title = "Non stiff solvers") +
    scale_color_manual(values = solvers_col_map, name = "Solver") + 
    my_minimal
  
  p_non_stiff_lib <- ggplot(data_non_stiff, aes(x = model_short, crit_plot, group=solver, color = solver_lib)) + 
    geom_point(stat="summary", fun=sum, size = 3.0, alpha=0.7) +
    stat_summary(fun=sum, geom="line", size = 1.5, alpha=0.7) +
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1], name = "Solver library") + 
    labs(x = "Model (short name)", y = ylab, title = "Library : Non stiff solvers") +
    my_minimal
  
  # Plot for sundials 
  data_sundials_stiff <- data_winner %>%
    filter(solver_lib == "Sundials") %>%
    filter(solver_type == "stiff")
  data_sundials_non_stiff <- data_winner %>%
    filter(solver_lib == "Sundials") %>% 
    filter(solver_type == "nonstiff")
  new_col <- c("#a0e85b", "#6314af", "#79c6c1", "#1f3ca6", "#bbc3fe")
  p1 <- ggplot(data_sundials_stiff, aes(x = model_short, crit_plot, group=solver, color = solver)) + 
    geom_point(stat="summary", fun=sum, size = 3.5) + 
    stat_summary(fun=sum, geom="line", size = 1.5) + 
    scale_y_log10() + 
    scale_color_manual(values = new_col, name = "Solver") + 
    labs(x = "Model (short name)", y = ylab, title = "Sundials : Stiff solvers") +
    my_minimal + theme(legend.position = "bottom") 
  p2 <- ggplot(data_sundials_non_stiff, aes(x = model_short, crit_plot, group=solver, color = solver)) + 
    geom_point(stat="summary", fun=sum, size = 3.5) + 
    stat_summary(fun=sum, geom="line", size = 1.5) + 
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1], name = "Solver") + 
    labs(x = "Model (short name)", y = "", title = "Sundials : Non stiff solvers") +
    my_minimal + theme(legend.position = "bottom")
  p_sundial <- ggpubr::ggarrange(p1, p2, ncol = 2)
  
  winner_category <- calc_winner_category(data, tol=tol, sq_res = sq_res) %>%
    mutate(model_short = factor(model_short, levels = c("M1", "M2", "M3", "M4", "M5", 
                                                        "M6", "M7", "M8", "M9", "M10", 
                                                        "M11", "M12", "M13", "M14", "M15", 
                                                        "M16", "M17", "M18", "M19", "M20", 
                                                        "M21", "M22", "M23", "M24", "M25", 
                                                        "M26", "M27", "M28"))) 
  
  if(sq_res == T){
    winner_category <- winner_category %>%
      rename("crit_plot" = "median_sq")
    ylab = "Median squared error"
  }else{
    winner_category <- winner_category %>%
      rename("crit_plot" = "median_time")
    ylab = "Median time [ns]"
  }
  
  p_winner <- ggplot(winner_category, aes(x = model_short, crit_plot, group=solver_type, color=solver_type)) + 
    geom_point(stat="summary", fun=sum, size=4.0) +
    stat_summary(fun=sum, geom="line", size=2.0) +
    labs(x = "Model (short name)", y = ylab, title = "Best solver each type") +
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1], name = "Solver type") + 
    my_minimal
  
  if(sq_res == T){
    dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Sq_res/")
  }else{
    dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Time/")
  }
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  
  ggsave(str_c(dir_save, "Stiff_solvers.svg"), p_stiff, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Stiff_solvers_lib.svg"), p_stiff_lib, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Non_stiff_solvers.svg"), p_non_stiff, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Non_stiff_solvers_lib.svg"), p_non_stiff_lib, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Sundials.svg"), p_sundial, width = BASE_WIDTH*4, height = BASE_HEIGHT*4)
  ggsave(str_c(dir_save, "Winner.svg"), p_winner, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  
}


plot_fail <- function(data, tol)
{
  
  data_fail <- calc_fail(data, tol) 
  
  p_fail <- ggplot(data_fail, aes(model, fill = solver_type)) + 
    geom_bar(position = position_dodge2(preserve = "single")) + 
    scale_fill_manual(values = cbPalette[-1], name = "Solver type") + 
    coord_flip() + 
    labs(y = "Count", x = "Model (short)", title = "Integration failures") + 
    my_minimal
  
  dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, "Solvers_fail.svg"), p_fail, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  
}


plot_score <- function(data, tol)
{

  data_score_mean <- calc_score(data, crit = "scale", tol=tol) %>%
    filter(rel_time < 1.5 && rel_sq < 1.5)
  data_score_mean_no_fail <- calc_score(data, crit = "scale", tol=tol, exclude_fail = T) %>%
    filter(rel_time < 1.5 && rel_sq < 1.5)
  p1 <- ggplot(data_score_mean, aes(rel_time, rel_sq, color = solver_lib, shape = solver_type)) + 
    geom_point(size = 6.0) + 
    scale_color_manual(values = cbPalette[-c(1, 5)], name = "Solver library") + 
    labs(x = "Relative time", y = "Relative squared error", title = "Mean normalisation") +
    my_minimal
  p2 <- ggplot(data_score_mean_no_fail, aes(rel_time, rel_sq, color = solver_lib, shape = solver_type)) + 
    geom_point(size = 6.0) + 
    scale_color_manual(values = cbPalette[-c(1, 5)], name = "Solver library") + 
    labs(x = "Relative time", y = "Relative squared error", title = "Mean normalisation (exclude fail)") +
    my_minimal
  
  dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, "Score.svg"), p1, width = BASE_WIDTH*1, height = BASE_HEIGHT*1)
  ggsave(str_c(dir_save, "Score_exclude.svg"), p2, width = BASE_WIDTH*1, height = BASE_HEIGHT*1)
  
}


plot_rank <- function(data, tol=1e-6, solver_type="all", sq_err=F)
{

  if(solver_type == "stiff"){
    data <- data %>%
      filter(solver_type == "stiff")
    name_save = "Position_stiff"
  }else if(solver_type == "not_stiff"){
    data <- data %>%
      filter(solver_type != "stiff")
    name_save = "Position_not_stiff"
  }else{
    data <- data
    name_save = "Position"
  }
  
  data_low_tol <- calc_winner(data, tol=1e-6) %>%
    group_by(solver) %>%
    summarise(mean_time = mean(rank_time, na.rm = T))
  solver_order <- data_low_tol$solver[order(data_low_tol$mean_time)]
  
  data_winner <- calc_winner(data, tol=tol) %>%
    mutate(model_short = factor(model_short, levels = c("M1", "M2", "M3", "M4", "M5", 
                                                        "M6", "M7", "M8", "M9", "M10", 
                                                        "M11", "M12", "M13", "M14", "M15", 
                                                        "M16", "M17", "M18", "M19", "M20", 
                                                        "M21", "M22", "M23", "M24", "M25", 
                                                        "M26", "M27", "M28")), 
           solver = factor(solver, levels = rev(solver_order)))
  
  if(sq_err == T){
    data_winner <- data_winner %>%
      rename("rank_plot" = "rank_sq") 
    ylab = "Ranking squared error"
    name_save <- str_c(name_save, "_sq")
     dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Sq_res/")
  }else if(sq_err == "Comb"){
    data_winner <- data_winner %>%
      mutate(rank_plot = (rank_sq + rank_time)*0.5)
    ylab = "Ranking squared error"
    name_save <- str_c(name_save, "_comb")
    dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Res_comb/")
  }else{
    data_winner <- data_winner %>%
      rename("rank_plot" = "rank_time")
    ylab = "Ranking squared error"
    name_save <- str_c(name_save, "_time")
    dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Time/")
  }
  
  p1 <- ggplot(data_winner, aes(x=model_short, y=solver)) + 
    geom_raster(aes(fill = rank_plot)) + 
    labs(x = "Model (increasing size ->)", y = ylab) + 
    scale_fill_viridis_c(direction = -1) + 
    theme_minimal()
  
  print(dir_save)
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, name_save, ".svg"), p1, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
}


plot_all_tol <- function(data)
{
  
  data_tot_stiff <- data %>%
    filter(solver_type == "stiff") %>%
    mutate(reltol = as.factor(reltol)) %>%
    group_by(reltol, model_short) %>%
    summarise(time = median(runTime, na.rm = T), 
              res_sq = median(sqDiff, na.rm = T)) %>%
    mutate(model_number = as.numeric(str_extract(model_short, "\\d+")))
  data_tot_not_stiff <- data %>%
    filter(solver_type != "stiff") %>%
    mutate(reltol = as.factor(reltol)) %>%
    group_by(reltol, model_short) %>%
    summarise(time = median(runTime, na.rm = T), 
              res_sq = median(sqDiff, na.rm = T)) %>%
    mutate(model_number = as.numeric(str_extract(model_short, "\\d+")))
  
  data_ref_stiff <- data_tot_stiff %>% filter(reltol == 1e-12) %>%
    select(model_short, time, res_sq) %>%
    rename("time_ref" = "time", "res_sq_ref" = "res_sq", "reltol_tmp" = "reltol")
  data_tot_stiff <- data_tot_stiff %>%
    inner_join(data_ref_stiff) %>%
    mutate(time_rel = time / time_ref, 
           res_sq_rel = (res_sq + 1e-100) / (res_sq_ref + 1e-100))
  
  data_ref_not_stiff <- data_tot_not_stiff %>% filter(reltol == 1e-12) %>%
    select(model_short, time, res_sq) %>%
    rename("time_ref" = "time", "res_sq_ref" = "res_sq", "reltol_tmp" = "reltol")
  data_tot_not_stiff <- data_tot_not_stiff %>%
    inner_join(data_ref_not_stiff) %>%
    mutate(time_rel = time / time_ref, 
           res_sq_rel = (res_sq + 1e-100) / (res_sq_ref + 1e-100))
  
  
  p1 <- ggplot(data_tot_stiff, aes(model_number, time_rel, color = reltol)) + 
    geom_line(size = 2.0) + 
    geom_point(size = 5.0) + 
    scale_color_manual(values = cbPalette[-1]) + 
    my_theme
  p2 <- ggplot(data_tot_stiff, aes(model_number, res_sq_rel, color = reltol)) + 
    geom_line(size = 2.0) + 
    geom_point(size = 5.0) + 
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1]) + 
    my_theme
  
  p3 <- ggplot(data_tot_not_stiff, aes(model_number, time_rel, color = reltol)) + 
    geom_line(size = 2.0) + 
    geom_point(size = 5.0) + 
    scale_color_manual(values = cbPalette[-1]) + 
    my_theme
  p4 <- ggplot(data_tot_not_stiff, aes(model_number, res_sq_rel, color = reltol)) + 
    geom_line(size = 2.0) + 
    geom_point(size = 5.0) + 
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1]) + 
    my_theme
  
  p5 <- ggplot(data_tot_not_stiff, aes(model_number, res_sq, color = reltol)) + 
    geom_line(size = 2.0) + 
    geom_point(size = 5.0) + 
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1]) + 
    my_theme
  p6 <- ggplot(data_tot_stiff, aes(model_number, res_sq, color = reltol)) + 
    geom_line(size = 2.0) + 
    geom_point(size = 5.0) + 
    scale_y_log10() + 
    scale_color_manual(values = cbPalette[-1]) + 
    my_theme
  
  dir_save <- str_c("../Results/ODE_solvers/All_tol/")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, "Stiff_time_rel.svg"), p1, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Stiff_sq_rel.svg"), p2, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Not_stiff_time_rel.svg"), p3, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Not_stiff_sq_rel.svg"), p4, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Not_stiff_sq.svg"), p5, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Stiff_sq.svg"), p6, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
}




data <- read_data("../Pipeline_ModelSolver/IntermediaryResults/benchmark_5.csv")


tol_list <- unique(data$abstol)
for(tol in tol_list){
  plot_winner(data, tol=tol)
  plot_winner(data, tol=tol, sq_res = T)
  plot_fail(data, tol=tol)
  plot_score(data, tol=tol)  
  plot_rank(data, tol=tol, solver_type = "stiff")
  plot_rank(data, tol=tol, solver_type = "not_stiff")
  plot_rank(data, tol=tol, solver_type = "all")
  plot_rank(data, tol=tol, solver_type = "stiff", sq_err = T)
  plot_rank(data, tol=tol, solver_type = "not_stiff", sq_err = T)
  plot_rank(data, tol=tol, solver_type = "all", sq_err = T)
  plot_rank(data, tol=tol, solver_type = "stiff", sq_err = "Comb")
  plot_rank(data, tol=tol, solver_type = "not_stiff", sq_err = "Comb")
  plot_rank(data, tol=tol, solver_type = "all", sq_err = "Comb")
}


plot_all_tol(data)


tol <- 1e-9
data_tmp <- data %>%
  filter(solver_type == "stiff")
data_low_tol <- calc_winner(data_tmp, tol=1e-6) %>%
  group_by(solver) %>%
  summarise(mean_time = mean(rank_time, na.rm = T))
solver_order <- data_low_tol$solver[order(data_low_tol$mean_time)]

data_winner <- calc_winner(filter(data, solver_type == "stiff"), tol=tol) %>%
  mutate(model_short_param = factor(model_short_param, levels = c("P1", "P2", "P3", "P4", "P5", 
                                                                  "P6", "P7", "P8", "P9", "P10", 
                                                                  "P11", "P12", "P13", "P14", "P15", 
                                                                  "P16", "P17", "P18", "P19", "P20", 
                                                                  "P21", "P22", "P23", "P24", "P25", 
                                                                  "P26", "P27", "P28")), 
         solver = factor(solver, levels = rev(solver_order)))


ggplot(data_winner, aes(x=model_short_param, y=solver)) + 
  geom_raster(aes(fill = rank_time)) + 
  labs(x = "Model (increasing param ->)", y = "Solver") + 
  scale_fill_viridis_c(direction = -1) + 
  theme_minimal()


data_plot1 <- tibble(P_vec = seq(0.0, 250), 
                     gamma = 1.0, 
                     D = 1.0 / (1.0 + P_vec*gamma), 
                     push = P_vec * gamma)
data_plot2 <- tibble(P_vec = seq(0.0, 250), 
                     gamma = 0.25, 
                     D = 1.0 / (1.0 + P_vec*gamma), 
                     push = P_vec * gamma)
data_plot3 <- tibble(P_vec = seq(0.0, 250), 
                     gamma = 0.1, 
                     D = 1.0 / (1.0 + P_vec*gamma), 
                     push = P_vec * gamma)
data_plot4 <- tibble(P_vec = seq(0.0, 250), 
                     gamma = 0.01, 
                     D = 1.0 / (1.0 + P_vec*gamma), 
                     push = P_vec * gamma)
data_plot <- data_plot1 %>% 
  bind_rows(data_plot2, data_plot3, data_plot4) %>%
  mutate(gamma = factor(gamma))
ggplot(data_plot, aes(P_vec, D, color = gamma)) + 
  geom_line(size = 2.0) + 
  geom_hline(yintercept = 0.05) + 
  scale_y_log10() + 
  geom_vline(xintercept = 90) + 
  scale_color_manual(values = cbPalette[-1]) + 
  my_theme
    
ggplot(data_plot, aes(P_vec, push, color = gamma)) + 
  geom_line(size = 2.0) + 
  geom_hline(yintercept = 10) + 
  scale_color_manual(values = cbPalette[-1]) + 
  my_theme
  

