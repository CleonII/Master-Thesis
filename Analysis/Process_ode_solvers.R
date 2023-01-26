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


read_data <- function(path_file)
{

  
  data_raw <- read_csv(path_file, col_types = cols()) 
  model_short = rep("", nrow(data_raw))
  for(i in 1:length(model_short)){
    iUse = str_locate(data_raw$model[i], "_[:alpha:]+_")
    model_short[i] = str_c(str_sub(data_raw$model[i], start=iUse[1]+1, end=iUse[2]-1), str_sub(data_raw$model[i], start=-4))
  }
  data_raw$model_short = model_short

  return(data_raw)
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



plot_rank <- function(data, tol=1e-6, solver_type="all", sq_err=F)
{

  if(solver_type == "stiff"){
    data <- data %>%
      filter(solverType == "stiff")
    name_save = "Position_stiff"
  }else if(solver_type == "nonstiff"){
    data <- data %>%
      filter(solver_type != "stiff")
    name_save = "Position_not_stiff"
  }else{
    data <- data
    name_save = "Position"
  }
  
  data_winner = calc_winner(data, tol)
  
  if(sq_err == T){
    data_winner <- data_winner |> 
      rename("rank_plot" = "rank_sq") 
    ylab = "Ranking squared error"
    name_save <- str_c(name_save, "_sq")
    dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Sq_res/")
  }else if(sq_err == "Comb"){
    data_winner <- data_winner |> 
      mutate(rank_plot = (rank_sq + rank_time)*0.5)
    ylab = "Ranking averaged time and error"
    name_save <- str_c(name_save, "_comb")
    dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Res_comb/")
  }else{
    data_winner <- data_winner |> 
      mutate(rank_plot = rank_time)
    ylab = "Ranking squared error"
    name_save <- str_c(name_save, "_time")
    dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/Time/")
  }
  
  posx = unique(data_winner$model_short[order(data_winner$n_states)])
  posy = data_winner |> group_by(solver) |> summarise(rank = median(rank_time))
  posy = posy$solver[order(posy$rank)]
  p1 = ggplot(data_winner, aes(model_short, solver)) + 
    geom_raster(aes(fill = rank_plot)) + 
    scale_fill_viridis_c(direction = -1) +
    scale_x_discrete(limits = posx) + 
    scale_y_discrete(limits = rev(posy)) + 
    labs(x = "Model (increasing size ->)", y = ylab, title = "Rank time") + 
    my_minimal +
    theme(plot.background = element_rect(fill = "white"), 
          axis.text.x = element_text(angle=90))
  
  print(dir_save)
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, name_save, ".svg"), p1, width = BASE_WIDTH*2.5, height = BASE_HEIGHT*3.0)
}


calc_fail <- function(data, tol=1e-8)
{

  model_list <- unique(data$model)
  solver_list <- unique(data$solver)
  data_use <- data |> filter(abstol == tol)
  
  data_fail = data_use |> 
    group_by(solver, model_short, solverType, solverLib) |> 
    summarise(nFail = sum(is.nan(runTime) / 3)) |> 
    filter(nFail != 0)
    
  return(data_fail)
}


plot_fail <- function(data, tol)
{
  
  data_fail <- calc_fail(data, tol) 
  data1 = data_fail |> 
    group_by(model_short, solverType) |> 
    summarise(nFail = as.integer(sum(nFail)))
  data2 = data_fail |> 
    group_by(solver, solverType) |> 
    summarise(nFail = as.integer(sum(nFail)))
  
  data1$solverType = factor(data1$solverType, levels = rev(c("nonstiff", "stiff", "hint", "composite")))
  pos1 = data1 |> group_by(model_short) |> summarise(nFail=sum(nFail)) 
  pos1 = pos1$model_short[order(pos1$nFail)]
  p1 = ggplot(data1, aes(model_short, nFail, fill = solverType)) + 
    geom_col(position=position_dodge(preserve = 'single')) + 
    scale_fill_manual(values = cbPalette[-1], name = "Solver type") + 
    scale_x_discrete(limits=pos1) + 
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data1$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures per model") + 
    my_minimal +
    theme(legend.position = "bottom")
  
  pos2 = unique(data2$solver[order(data2$nFail)])
  p2 = ggplot(data2, aes(solver, nFail, fill = solverType)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values = cbPalette[-1], name = "Solver type") + 
    scale_x_discrete(limits = pos2) +
    scale_y_continuous(breaks = seq(from=1, by=2, to=max(data2$nFail))) +
    coord_flip() + 
    labs(y = "Number of integration failures", x = "", title = "Number of integration failures (total of 27 models)") + 
    my_minimal + 
    theme(legend.position = "bottom")
  
  dir_save <- str_c("../Results/ODE_solvers/", as.character(tol), "/")
  if(!dir.exists(dir_save)) dir.create(dir_save, recursive = T)
  ggsave(str_c(dir_save, "Solvers_fail_solver.svg"), p1, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
  ggsave(str_c(dir_save, "Solvers_fail_model.svg"), p2, width = BASE_WIDTH*2, height = BASE_HEIGHT*2)
}


calc_winner <- function(data, tol=1e-6)
{
  
  return_first = function(value) return(value[1])
  
  data_winner <- tibble()
  data_use <- data |>  
    filter(abstol == tol) 
  model_list <- unique(data$model_short)
  n_solvers = length(unique(data$solver))
  
  for(i in 1:length(model_list)){
    
    data_model_i <- data_use |> 
      filter(model_short == model_list[i]) |> 
      group_by(solver, model_short) |> 
      summarise(median_time = median(runTime, na.rm = T), 
                median_sq = median(sqDiff, na.rm = T), 
                solver_type = return_first(solverType), 
                solver_lib = return_first(solverLib), 
                n_param = median(nParam), 
                n_states = median(nStates))
    # Add ranked index on time or sq-error
    data_model_i <- data_model_i[order(data_model_i$median_time), ]
    data_model_i$rank_time <- 1:nrow(data_model_i)
    data_model_i <- data_model_i[order(data_model_i$median_sq), ]
    data_model_i$rank_sq <- 1:nrow(data_model_i)
    
    # Account for integration failures     
    data_model_i$rank_time[is.na(data_model_i$median_time)] = n_solvers
    data_model_i$rank_sq[is.na(data_model_i$median_sq)] = n_solvers
    
    data_winner <- bind_rows(data_winner, data_model_i)
  }
  
  return(data_winner)
}



data <- read_data("../Intermediate/Benchmarks/ODE_solvers/Sparse_not_linsolvers_new.csv") |> 
  mutate(solverType =  case_when(solverType == "nonStiff" ~ "nonstiff", 
                                 T ~ solverType)) |> 
  filter(solverLib != "ODEInterface") |> 
  filter(solver != "lsoda")
tol = 1e-8



plot_fail(data, tol)
plot_rank(data, tol, solver_type = "all", sq_err = F)
plot_rank(data, tol, solver_type = "stiff", sq_err = F)
plot_rank(data, tol, solver_type = "nonstiff", sq_err = F)
plot_rank(data, tol, solver_type = "all", sq_err = T)
plot_rank(data, tol, solver_type = "stiff", sq_err = T)
plot_rank(data, tol, solver_type = "nonstiff", sq_err = T)



ggplot(data_winner, aes(model_short, solver)) + 
  geom_raster(aes(fill = rank_sq)) + 
  scale_fill_viridis_c(direction = -1) +
  scale_x_discrete(limits = posx) + 
  scale_y_discrete(limits = rev(posy)) + 
  labs(x = "Model (increasing size ->)", y = "", title = "Rank squared error") + 
  my_minimal


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


# -----------------------------------------------------------------------------------------------------------------
# Random parameters 
# -----------------------------------------------------------------------------------------------------------------
dir_data = "../Intermediate/Benchmarks/ODE_solvers/"
data_random_parameters = read_data(str_c(dir_data, "Random_parameters.csv")) |> 
  mutate(solverType = case_when(solverType == "nonStiff" ~ "nonstiff", 
                                T ~ solverType))

data_plot = data_random_parameters |> 
  filter(abstol == 1e-8) |> 
  mutate(solver = case_when(solver == "CVODE_BDF_default" ~ "CVODE_BDF", 
                            T ~ solver))

p1 = ggplot(data_plot, aes(solver, runTime, fill = solverType)) + 
  geom_violin(draw_quantiles = 0.5, linewidth=1.0) +
  geom_jitter(width = 0.1, size=0.2) + 
  facet_wrap(~model_short, scales="free_y") + 
  scale_y_log10() +
  labs(x = "", y = "Run time [s]", title = "Run time for 100 random paramter vectors", 
       subtitle = "On average stiff solvers have less variabillity") + 
  scale_fill_manual(values = cbPalette[-1], name = "Solver type") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")

data_fail = data_plot |> 
  group_by(solver, model_short, solverType) |> 
  summarise(n_fail = sum(is.nan(runTime)))
p2 = ggplot(data_fail, aes(solver, n_fail, fill = solverType)) + 
  geom_bar(stat="identity", position="dodge") + 
  facet_wrap(~model_short) + 
  ylim(0, 100) + 
  labs(x = "", y = "Percentage integration failures", title = "Percentage integration failure for 100 random parameter vectors", 
       subtitle = "On average nonstiff solvers have more integration failures") + 
  scale_fill_manual(values = cbPalette[-1], name = "Solver type") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "bottom")

ggsave("Random_parameter_run_time.png", p1, width = BASE_WIDTH*3.5, height = BASE_HEIGHT*2.5, dpi=300)
ggsave("Random_parameter_fail.png", p2, width = BASE_WIDTH*3.5, height = BASE_HEIGHT*2.5, dpi=300)
