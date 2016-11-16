library(tidyverse)
library(viridis)
library(cowplot)

sigmoid <- function(pst, mu0, k, t0) 2 * mu0 / (1 + exp(-k*(pst - t0)))

N <- 100

pst <- runif(N)


# Discrete case -----------------------------------------------------------

df_discrete <- data_frame(
  pst = pst,
  x1 = sigmoid(pst, 2, 10, 0.7),
  x2 = sigmoid(pst, -2, 10, 0.7)
)

gather(df_discrete, covariate, y, -pst) %>% 
  arrange(pst) %>% 
  ggplot(aes(x = pst, y = y, color = covariate)) +
  geom_line(size = 1.5) +
  scale_color_brewer(name = "Discrete\nCovariate", palette = "Set1") +
  xlab("Patient trajectory") + ylab("Dynamic observable") +
  theme(legend.text = element_blank(), legend.direction = "horizontal",
        legend.position = "top", axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size = 11))

discrete_plot <- last_plot()



# Continuous case ---------------------------------------------------------

mu0 <- seq(-2, 2, length.out = 50)

df_conts <- lapply(mu0, function(m0) {
  data_frame(
    pst = pst,
    y = sigmoid(pst, m0, 10, 0.7),
    mu0 = m0
  )
})

dfc <- bind_rows(df_conts)

dfc %>% arrange(mu0, pst) %>% 
  ggplot(aes(x = pst, y = y, group = mu0, color = mu0)) +
  geom_line(size = 1.5) + scale_color_viridis(name = "Continuous\nCovariate") +
  theme(legend.text = element_blank(), legend.direction = "horizontal",
        legend.title = element_text(hjust = 0.5),
        legend.position = "top", axis.text = element_blank(),
        axis.ticks = element_blank(), axis.title = element_text(size = 11)) +
  xlab("Patient trajectory") + ylab("Dynamic observable")

continuous_plot = last_plot()

plot_grid(discrete_plot, continuous_plot)

phenotime_diagram <- last_plot()

save(phenotime_diagram, file = "../figs/phenotime_diagram.Rdata")
ggsave("../figs/phenotime_diagram.png", width = 9, height = 3)

  
  