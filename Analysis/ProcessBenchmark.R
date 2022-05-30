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


# Processing forward-gradient 
dir_result <- "..//Pipeline_ModelParameterEstimation/IntermediaryResults/ForwardGradient/"
overview_file <- read_csv(str_c(dir_result, "benchmark_model_Bachmann_MSB2011_1.csv"), col_types = cols()) %>%
  rename("solver" = "option1", 
         "n_steps" = "option2", 
         "balance" = "option3", 
         "step_range" = "option4") %>%
  mutate(n_steps = as.factor(n_steps), 
         step_range = as.factor(step_range), 
         balance = as.factor(balance), 
         startParameterIndex = as.factor(startParameterIndex)) %>%
  filter(cost != 0.0)

ggplot(overview_file, aes(n_steps, cost, fill = step_range)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = cbPalette[-1]) + 
  my_theme
ggplot(overview_file, aes(n_steps, cost, fill = balance)) + 
  geom_boxplot() + 
  scale_y_log10() + 
  scale_fill_manual(values = cbPalette[-1]) + 
  my_theme

data_10000 <- overview_file %>%
  filter(n_steps == 10000)
ggplot(data_10000, aes(startParameterIndex, cost, fill = step_range)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette[-1]) + 
  scale_y_log10()+ 
  my_theme
ggplot(data_10000, aes(startParameterIndex, cost, fill = balance)) + 
  geom_violin() + 
  scale_fill_manual(values = cbPalette[-1]) + 
  scale_y_log10()+ 
  my_theme
  

dir_result <- "..//Pipeline_ModelParameterEstimation/IntermediaryResults/ForwardAutomatic/"
overview_file1 <- read_csv(str_c(dir_result, "benchmark_model_Bachmann_MSB2011_3.csv"), col_types = cols()) %>%
  rename("solver" = "option1",
         "alg" = "option2",
         "balance" = "option3", 
         "step_range" = "option4") %>%
  mutate(step_range = as.factor(step_range), 
         balance = as.factor(balance), 
         startParameterIndex = as.factor(startParameterIndex)) %>%
  filter(cost != 0.0)

ggplot(overview_file1, aes(startParameterIndex, cost)) + 
  geom_point() + 
  my_theme
    

