library(tidyverse)
library(ggthemes)
library(gt)


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
col_highlight = c("#D0C0B0", "#B6C2CC", "#BEAAB4", "#ECE9CD", "#0A3D6B", "#0D5C3D", "#812F02")


# -------------------------------------------------------------------------------------------------------------
# Cost, grad and hess 
# -------------------------------------------------------------------------------------------------------------
dir_result <- "../Intermediate/Benchmarks/Cost_grad_hess/"
data <- read_csv(str_c(dir_result, "Cost_grad_hess_tol_low_composite.csv"), col_types = cols()) |> 
  group_by(model, solver) |> 
  summarise(T_cost = median(T_cost), 
            T_grad = median(T_grad), 
            T_grad_zygote = median(T_grad_zygote), 
            T_hess = median(T_hess)) |> 
  mutate(model = str_sub(str_extract(model, "model_([:alpha:]+)_"), 7, -2)) |> 
  mutate(n_param = case_when(model == "Bachmann" ~ 113, 
                             model == "Beer" ~ 72, 
                             model == "Boehm" ~ 9, 
                             model == "Bruno" ~ 13, 
                             model == "Crauste" ~ 12, 
                             model == "Elowitz" ~ 21, 
                             model == "Fiedler" ~ 22, 
                             model == "Fujita" ~ 19, 
                             model == "Lucarelli" ~ 84, 
                             model == "Sneyd" ~ 15))
  
data_min <- data |> 
  group_by(model) |> 
  summarise(T_cost_min = min(T_cost), 
            T_grad_min = min(T_grad), 
            T_hess_min = min(T_hess))
data_plot <- inner_join(data, data_min, by=c("model"))
pos = unique(data_plot$model[order(data_plot$n_param)])
data_plot$model = factor(data_plot$model, levels = pos)

ggplot(data_plot, aes(model, T_cost / T_cost_min, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_x_discrete(breaks = pos) + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Normalised run-time by fastest solver", title = "Cost") + 
  ylim(0, 10) + 
  geom_hline(yintercept = 1) +
  my_theme
ggplot(data_plot, aes(model, T_grad / T_grad_min, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Normalised run-time by fastest solver", title = "Gradient") + 
  ylim(0, 10) + 
  geom_hline(yintercept = 1) +
  my_theme
ggplot(data_plot, aes(model, T_hess / T_hess_min, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Normalised run-time by fastest solver", title = "Hessian") + 
  ylim(0, 10) + 
  geom_hline(yintercept = 1) +
  my_theme

ggpubr::ggarrange(p1, p2, p3, ncol = 2, nrow=2, common.legend = T)

data_param <- data_plot |> 
  group_by(model) |> 
  summarise(n_param = median(n_param))

ggplot(data_plot, aes(model)) + 
  geom_bar(aes(y = T_grad / T_cost, fill = solver), stat="identity", position = "dodge") +
  geom_text(aes(model, y= n_param, label = as.character(n_param)), size=3.0) + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Ratio Gradient / Cost", title = "Gradient", 
       subtitle = "Number = number of parameters to estimate") + 
  my_theme
ggplot(data_plot, aes(model)) + 
  geom_bar(aes(y = T_hess / T_cost, fill = solver), stat="identity", position = "dodge") +
  geom_text(aes(model, y= n_param^2, label = as.character(n_param^2)), size=3.0) + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Ratio Hessian / Cost", title = "Hessian", 
       subtitle = "Number = number of parameters to estimate squared") + 
  scale_y_log10() +
  my_theme
p = ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = T, legend = "bottom")
ggsave("Grad_hess.png", p, width = BASE_WIDTH*2.8, height = BASE_HEIGHT*1.2)

my_tab <- data_plot |> 
  select(model, solver, T_cost, T_grad, T_hess) |> 
  gt() |> 
  tab_header(title = "Run time for cost, gradient, and hessian") |> 
  tab_spanner(label = "ODE Solver", columns = c(solver)) |> 
  tab_spanner(label = "Run time [s]", 
              columns = c(T_cost, T_grad, T_hess)) |> 
  cols_label(solver= "", T_cost="Cost", T_grad="Gradient", T_hess="Hessian") |> 
  fmt(columns = starts_with("T_"), fns = function(x) sprintf("%.2e", x)) 
gtsave(my_tab, "Table_run_time.html")


# -------------------------------------------------------------------------------------------------------------
# Against AMICI 
# -------------------------------------------------------------------------------------------------------------
dir_result <- "../Intermediate/Benchmarks/Cost_grad_hess/"
if(file.exists(str_c(dir_result, "computation_times.csv"))){

  data <- read_csv(str_c(dir_result, "Cost_grad_hess.csv"), col_types = cols()) |> 
    group_by(model, solver) |> 
    summarise(T_cost = median(T_cost), 
              T_grad = median(T_grad), 
              T_grad_zygote = median(T_grad_zygote), 
              T_hess = median(T_hess)) |> 
    mutate(model = str_sub(str_extract(model, "model_([:alpha:]+)_"), 7, -2)) |> 
    mutate(n_param = case_when(model == "Bachmann" ~ 113, 
                               model == "Beer" ~ 72, 
                               model == "Boehm" ~ 9, 
                               model == "Bruno" ~ 13, 
                               model == "Crauste" ~ 12, 
                               model == "Elowitz" ~ 21, 
                               model == "Fiedler" ~ 22, 
                               model == "Fujita" ~ 19, 
                               model == "Lucarelli" ~ 84, 
                               model == "Sneyd" ~ 15))
  
  data_amici = read_csv(str_c(dir_result, "computation_times.csv"), col_types = cols())
  colnames(data_amici) = c("model", "T_cost", "AMICI_Forward", "AMICI_Adjoint", "n_param")
  data_amici_ = data_amici |> 
    pivot_longer(c("AMICI_Forward", "AMICI_Adjoint"), names_to = "solver", values_to = "T_grad") |> 
    mutate(model = str_extract(model, "^[:alpha:]+")) |> 
    select(-n_param)
  
  amici_best = data_amici_ |> 
    group_by(model) |> 
    summarise("Cost_best" = min(T_cost), 
              "Grad_best" = min(T_grad))
  
  data_plot = bind_rows(data, data_amici_) |> 
    inner_join(amici_best, by="model") |> 
    filter(!model %in% c("Borghans", "Brannmark", "Isensee", "Schwen", "Weber", "Zheng")) |> 
    filter(solver != "AMICI_Adjoint")
  pos = unique(data_plot$model[order(data_plot$n_param)])
  data_plot$model = factor(data_plot$model, levels = pos)
  
  p1 = ggplot(data_plot, aes(model, T_cost / Cost_best, fill = solver)) + 
    geom_bar(stat="identity", position="dodge") + 
    scale_fill_manual(values = cbPalette[-1]) + 
    my_theme
  p2 = ggplot(data_plot, aes(model, T_grad / Grad_best, fill = solver)) + 
    geom_bar(stat="identity", position="dodge") + 
    annotate("text", x = "Beer", y = 2.5, label = "Our code is bad at these kind of models") +
    scale_fill_manual(values = cbPalette[-1]) + 
    scale_y_continuous(breaks = 0:7) +
    my_theme
  
  data_amici_adjoint = data_amici_ |> 
    filter(solver == "AMICI_Adjoint")
  
  p_save = ggpubr::ggarrange(p1, p2, ncol=2, common.legend = T, legend="bottom")
  ggsave("Compare_grad_cost_AMICI.png", width = BASE_WIDTH*3.0, height = BASE_HEIGHT)
}


# -------------------------------------------------------------------------------------------------------------
# Bachman gradient scaling 
# -------------------------------------------------------------------------------------------------------------
dir_result <- "../Intermediate/Benchmarks/Cost_grad_hess/"
data  <- read_csv(str_c(dir_result, "Bachman_fix_param.csv"), col_types = cols()) |> 
  filter(!is.infinite(Time)) |> 
  mutate(N_param_est = 28 - N_param_fixed) 

#ForEq_AutoDiff
data_plot <- data |> 
  filter(solver == "QNDF") |> 
  filter(Method_info == "ForwardDiff")

p = ggplot(data_plot, aes(N_param_est, Time, color = chunk_size, fill = chunk_size))  +
  geom_point() + 
  geom_smooth() + 
  scale_fill_manual(values = cbPalette[-1], name = "Chunk size") + 
  scale_color_manual(values = cbPalette[-1], name = "Chunk size") + 
  scale_x_continuous(breaks = seq(from = 2, by = 2, to  = 28)) +
  labs(x = "Number of parameter to take gradient on", y = "Time [s]", 
       title = "Bachman model", 
       subtitle = "Chunking is the secret sauce behind the performance of Julia gradients") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
ggsave("Bachman_gradient_time.png", p, width = BASE_WIDTH*1.5, height = BASE_HEIGHT*1.5)

data_median = data_plot |> 
  filter(N_param_fixed == 3) |> 
  filter(chunk_size == "Default") |> 
  group_by(Method_info) |> 
  summarise(median = median(Time)) |> 
  pull(median)

# See effect of different chunks 
data_chunks  <- read_csv(str_c(dir_result, "Bachman_test_chunks.csv"), col_types = cols()) |> 
  filter(!is.infinite(Time)) |> 
  mutate(N_param_est = as.factor(28 - N_param_fixed))
data_plot = data_chunks |> 
  filter(solver == "QNDF") |> 
  filter(Method_info == "ForwardDiff") |> 
  filter(N_param_fixed != 0)

p = ggplot(data_plot, aes(chunk_size, Time, color = N_param_est, fill = N_param_est)) + 
  geom_point() + 
  geom_hline(yintercept = data_median) +   
  geom_smooth(method = "loess", formula = y ~ x, span=0.3) + 
  scale_fill_manual(values = cbPalette[-1], name = "Number of parameter to estimate") + 
  scale_color_manual(values = cbPalette[-1], name = "Number of parameter to estimate") + 
  scale_x_continuous(breaks = seq(from = 1, by = 2, to  = 22)) +
  annotate("text", x = 17, y = 1.3, label = "Run time default chunking (1.38s)", size=5.0) +
  annotate("text", x = 6, y = 1.0, label = "Best chunk size (1.07s)", size=5.0) +
  labs(x = "Number of chunks", y = "Run time [s]", title = "Testing the effect of different chunks for Bachman model", 
       subtitle = "The default value is slower than a tuned value") +
  ylim(0, 3.5) + 
  theme_bw(base_size = 16) + 
  theme(legend.position = "bottom")
ggsave("Bachman_gradient_chunks.png", p, width = BASE_WIDTH*1.5, height = BASE_HEIGHT*1.5)


# Try effect of random parameters on chunking 
data  <- read_csv(str_c(dir_result, "Bachman_test_chunks_random_p.csv"), col_types = cols()) |> 
  filter(!is.infinite(Time)) |> 
  filter(Method_info == "ForwardDiff") 

data_ret = tibble()
param_i = unique(data$i_random_parameter)
for(i in 1:length(param_i)){
  data_tmp = data |> filter(i_random_parameter == param_i[i])
  data_tmp = data_tmp[order(data_tmp$Time), ]
  data_tmp$rank = 1:nrow(data_tmp)
  data_ret = bind_rows(data_ret, data_tmp)
}
data_plot = data_ret |> mutate(chunk_size = factor(chunk_size, levels = 1:26))

ggplot(data_plot, aes(i_random_parameter, chunk_size, fill = rank)) + 
  geom_tile() + 
  scale_fill_viridis_c(direction = -1, name = "Rank") + 
  scale_x_continuous(expand=c(0, 0), breaks=seq(from=0, to = 100, by = 5)) +
  labs(x = "Index random parameter", y = "Chunk size", title = "Bachman - Run time ranking different chunk-sizes for random parameter vectors", 
       subtitle = "A chunk size of 7 or 9 consistently are top-performing across parameters") +
  my_theme

data_plot1 = data_plot |> 
  filter(chunk_size == 9 | chunk_size == 7)
data_plot2 = data_plot |> 
  filter(!(chunk_size == 9 | chunk_size == 7))

ggplot(data_plot1, aes(i_random_parameter, Time)) + 
  geom_point(data=data_plot2, mapping=aes(i_random_parameter, Time, group=chunk_size), color=cbPalette[1]) +
  geom_line(data=data_plot2, mapping=aes(i_random_parameter, Time, group=chunk_size), color=cbPalette[1]) +
  geom_point(aes(color = chunk_size), size=3.0) + 
  geom_line(aes(color = chunk_size), linewidth=2.0) + 
  scale_color_manual(values = cbPalette[-1], name = "Chunk size") +
  scale_x_continuous(breaks=seq(from=0, to = 100, by = 5)) +
  labs(x = "Index random parameter", y = "Run time [s]", title = "Bachman - Run time for different chunk-sizes", 
       subtitle = "Chunk size 7 consistently performs best (grey lines chunk-sizes from 1 - 26)") +
  my_theme

ggsave("Bachman_chunk_random_p_heat.png", p1, width = BASE_WIDTH*1.5, height = BASE_HEIGHT*1.5)
ggsave("Bachman_chunk_random_p_line.png", p2, width = BASE_WIDTH*1.5, height = BASE_HEIGHT*1.5)


# -------------------------------------------------------------------------------------------------------------
# Lucarelli gradient scaling 
# -------------------------------------------------------------------------------------------------------------
dir_result <- "../Intermediate/Benchmarks/Cost_grad_hess/"
data  <- read_csv(str_c(dir_result, "Lucarelli_fix_param.csv"), col_types = cols()) |> 
  filter(!is.infinite(Time)) |> 
  mutate(N_param_est = max(N_param_fixed) + 1- N_param_fixed)  

#ForEq_AutoDiff
data_plot <- data |> 
  filter(solver == "QNDF") |> 
  filter(Method_info == "ForwardDiff")

p = ggplot(data_plot, aes(N_param_est, Time, color = chunk_size, fill = chunk_size))  +
  geom_point() + 
  geom_smooth() + 
  scale_fill_manual(values = cbPalette[-1], name = "Chunk size") + 
  scale_color_manual(values = cbPalette[-1], name = "Chunk size") + 
  scale_x_continuous(breaks = seq(from = 1, by = 2, to  = 71)) +
  labs(x = "Number of parameter to take gradient on", y = "Time [s]", 
       title = "Lucarelli model", 
       subtitle = "Chunking is the secret sauce behind the performance of Julia gradients") +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom")
ggsave("Lucarelli_gradient_time.png", p, width = BASE_WIDTH*1.5, height = BASE_HEIGHT*1.5)


# See effect of different chunks 
data_chunks  <- read_csv(str_c(dir_result, "Lucarelli_chunk_size.csv"), col_types = cols()) |> 
  filter(!is.infinite(Time)) 
data_plot_ = data_chunks |> 
  filter(solver == "QNDF") |> 
  filter(Method_info == "ForwardDiff") |> 
  filter(N_param_fixed == 0)
data_plot = data_plot_ |> filter(chunk_size != "Default") |> mutate(chunk_size = as.integer(chunk_size))
default_val = data_plot_ |> filter(chunk_size == "Default") |> pull(Time)

p = ggplot(data_plot, aes(chunk_size, Time)) + 
  geom_point() + 
  geom_hline(yintercept = default_val) +   
  geom_smooth(method = "loess", formula = y ~ x, span=0.3) + 
  scale_x_continuous(breaks = seq(from = 1, by = 2, to  = 71)) +
  #annotate("text", x = 17, y = 1.3, label = "Run time default chunking (1.38s)", size=5.0) +
  #annotate("text", x = 6, y = 1.0, label = "Best chunk size (1.07s)", size=5.0) +
  labs(x = "Number of chunks", y = "Run time [s]", title = "Testing the effect of different chunks for Lucarelli  model", 
       subtitle = "The default value is slower than a tuned value, but we are not looking at a linear relationship") +
  theme_bw(base_size = 16) + 
  theme(legend.position = "bottom")
ggsave("Lucarelli_gradient_chunks.png", p, width = BASE_WIDTH*1.5, height = BASE_HEIGHT*1.5)
