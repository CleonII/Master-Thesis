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
  
data_min <- data |> 
  group_by(model) |> 
  summarise(T_cost_min = min(T_cost), 
            T_grad_min = min(T_grad), 
            T_hess_min = min(T_hess))
data_plot <- inner_join(data, data_min, by=c("model"))
pos = unique(data_plot$model[order(data_plot$n_param)])
data_plot$model = factor(data_plot$model, levels = pos)

p1 <- ggplot(data_plot, aes(model, T_cost / T_cost_min, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_x_discrete(breaks = pos) + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Normalised run-time by fastest solver", title = "Cost") + 
  my_theme
p2 <- ggplot(data_plot, aes(model, T_grad / T_grad_min, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Normalised run-time by fastest solver", title = "Gradient") + 
  my_theme
p3 <- ggplot(data_plot, aes(model, T_hess / T_hess_min, fill = solver)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Normalised run-time by fastest solver", title = "Hessian") + 
  my_theme

ggpubr::ggarrange(p1, p2, p3, ncol = 2, nrow=2, common.legend = T)

data_param <- data_plot |> 
  group_by(model) |> 
  summarise(n_param = median(n_param))

p1 <- ggplot(data_plot, aes(model)) + 
  geom_bar(aes(y = T_grad / T_cost, fill = solver), stat="identity", position = "dodge") +
  geom_text(aes(model, y= n_param, label = as.character(n_param)), size=3.0) + 
  scale_fill_manual(values = cbPalette[-1]) + 
  labs(x = "", y = "Ratio Gradient / Cost", title = "Gradient", 
       subtitle = "Number = number of parameters to estimate") + 
  my_theme
p2 <- ggplot(data_plot, aes(model)) + 
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