setwd("~/Arbeit/DFG/Paper_SEMIFAR/Data")
library(ggplot2)
library(smootslm)

X = read.csv("london-air-quality.csv", sep = ";", header = TRUE)
y = X$AQI[!is.na(X$AQI)]
muy = mean(y)
n = length(y)
year = (1:n)/n * 138 + 1880
y = log(y)

result = tsmoothlm(y, pmax = 3, qmax = 3, p = 1, InfR = "Var")
result$iterations

result.der1 = dsmoothlm(y, pmax = 3, qmax = 3, pp = 3, d = 1, mu = 2, InfR.p = "Var")
result.der1$iterations

result.der2 = dsmoothlm(y, pmax = 3, qmax = 3, pp = 3, d = 2, mu = 3, InfR.p = "Var")
result.der2$iterations

g0 = result$ye
g0.der1 = result.der1$ye
g0.der2 = result.der2$ye

res = y - g0# + muy

df = data.frame(cbind(y, year, g0, g0.der1, g0.der2, res))

plot.trend <- ggplot(df, aes(x = year, y = y)) + 
  geom_line(aes(color = "Air Pollution"), size = 0.25, linetype = "dotted") +
  geom_line(aes(y = g0, color = "Trend (local cubic)"), size = 0.25) +
  labs(title = "(a) Log of air pollution (AQI) & estimated trend", y = "Log-Air Pollution - AQI", x = "") + 
  scale_x_continuous(name = "", breaks = seq(1880, 2020, 20)) +
  scale_color_manual(name = "Lines:",
                     breaks = c("Air Pollution", "Trend (local cubic)"),
                     values = c("Air Pollution" = "black", "Trend (local cubic)" = "red")) +
  theme(legend.position = c(0.1275, 0.845), legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', size = 0.25),
        plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

plot.der1 <- ggplot(df, aes(x = year, y = g0.der1)) + 
  geom_line(size = 0.25) +
  geom_hline(yintercept = 0, size = 0.25) +
  labs(title = "(c) Estimated first derivative", y = "1st derivative", x = "") +
  scale_x_continuous(name = "", breaks = seq(1880, 2020, 20)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25))

plot.res <- ggplot(df, aes(x = year, y = res)) +
  geom_line(size = 0.25) + 
  labs(title = "(b) Trend-adjusted Residuals", y = "Resiudals", x = "Year") +
  scale_x_continuous(name = "Year", breaks = seq(1880, 2020, 20)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 9),
        axis.text = element_text(size = 7), axis.ticks = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25))

plot.der2 <- ggplot(df, aes(x = year, y = g0.der2)) +
  geom_line(size = 0.25) +
  geom_hline(yintercept = 0, size = 0.25) +
  labs(title = "(d) Estimated second derivative", y = "2nd derivative", x = "Year") +
  scale_x_continuous(name = "Year", breaks = seq(1880, 2020, 20)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 9),
        axis.text = element_text(size = 7), axis.ticks = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25))

ggpubr::ggarrange(plot.trend, plot.der1, plot.res, plot.der2, heights = c(0.9, 1), nrow = 2, ncol = 2)
setwd("~/Arbeit/DFG/Paper_SEMIFAR/Latex_aktuell/Abb")
ggsave("London_AirQuality_AQI.pdf", height = 5, width = 9.5, dpi = 600)

