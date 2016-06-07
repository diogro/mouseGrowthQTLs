source("./read_mouse_data.R")

power_Data = read_csv("./data/area traits/power_analysis.csv")

(power_plot <- power_Data %>% group_by(.id) %>% summarize(power = sum(V1 > 20)/1000) %>%
  ggplot(aes(.id*100, power)) + geom_smooth() + geom_point() + geom_hline(yintercept = 0.95) + 
  labs(x = "Effect size (% V_p)", y = "Proportion of significan tests"))
save_plot("./data/area traits/power_plot.png", power_plot, base_height = 5, base_aspect_ratio = 1.5)
