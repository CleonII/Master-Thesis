library(tidyverse)
library(ggthemes)
library(colorspace)


# General plotting parameters (plot using theme-tufte)
cbPalette <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
col_highlight = c("#D0C0B0", "#B6C2CC", "#BEAAB4", "#ECE9CD", "#0A3D6B", "#0D5C3D", "#812F02")
my_theme <- theme_hc(base_size = 22) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=22))
my_minimal <- theme_minimal(base_size = 22) + theme(plot.title = element_text(hjust = 0.5, size = 14, face="bold"), 
                                             plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.title=element_text(size=22))
BASE_HEIGHT <- 5
BASE_WIDTH <- 7.0



# Large scale model 
data_large = read_csv("../Intermediate/Benchmarks/ODE_solvers/Large_models.csv", 
                      col_types = cols()) |> 
  mutate(is_sparse = case_when(str_sub(solver, -1, -1) == "S" ~ T,
                               T ~ F)) |> 
  group_by(solver, is_sparse) |> 
  summarise(runTime = median(runTime))

min_time = min(data_large$runTime, na.rm = T)

pos = data_large$solver[order(data_large$runTime)]
p1 = ggplot(data_large, aes(solver, runTime, fill = is_sparse)) + 
  geom_bar(stat="identity") + 
  geom_text(aes(label = sprintf("%.2e", runTime)), stat="identity") +
  coord_flip() + 
  scale_fill_manual(values = cbPalette[-1], name = "Sparse Jacobian") + 
  scale_x_discrete(limits=rev(pos)) + 
  labs(x = "", y = "Run time [s]") + 
  my_minimal + 
  theme(legend.position = "bottom", 
        plot.background = element_rect(fill = "white"))
  
data_sparse = data_large |> 
  filter(str_sub(solver, -1, -1) == "S") |> 
  select(solver, runTime) |> 
  rename("runTimeSparse" = "runTime")
data_not_sparse = data_large |> 
  filter(str_sub(solver, -1, -1) != "S") |> 
  select(solver, runTime) |> 
  rename("runTimeDense" = "runTime")
data_plot = bind_cols(data_sparse, data_not_sparse)

p2 = ggplot(data_plot, aes(runTimeSparse, runTimeDense)) + 
  geom_abline(slope=1.0, intercept = 0.0) + 
  geom_abline(slope=2.0, intercept = 0.0, linetype=2) + 
  geom_abline(slope=4.0, intercept = 0.0, linetype=2) + 
  geom_abline(slope=8.0, intercept = 0.0, linetype=2) + 
  geom_point(size=3.0) + 
  geom_text(aes(label = solver...1), nudge_x = 1.5) +
  labs(x = "Run time sparse Jacobian", y = "Run time dense Jacobian", 
       title = "Dashed lines = Sparse Jacobian 2, 4, and 8 times faster") +
  xlim(0.0, 15.0) + 
  ylim(0.0, 15.0)  + 
  theme_classic()

ggsave("Run_time_chen_solve.png", p1, dpi = 300, width = BASE_WIDTH, height = BASE_HEIGHT)
ggsave("Chen_sparse_vs_dense.png", p2, dpi = 300, width = BASE_WIDTH, height = BASE_HEIGHT)
