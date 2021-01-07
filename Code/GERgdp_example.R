library(ggplot2)
library(smootslm)
setwd("~/Arbeit/DFG/Paper_SEMIFAR/Data")
data = read.csv("GER_GDP_ann_1850_2016.csv", sep = ";", header = TRUE)


gdp = data$RealGDP


y <- log(gdp)
n = length(y)
year = (1:n)/n * 167 + 1850

result = tsmoothlm(y, q.min = 0, p.min = 0, p.max = 3, q.max = 3, p = 1, InfR = "Opt")
result$iterations

result.der1 = dsmoothlm(y, p.max = 3, q.max = 3, pp = 1, d = 1, mu = 2, InfR.p = "Opt")
result.der1$iterations

result.der2 = dsmoothlm(y, p.max = 3, q.max = 3, pp = 1, d = 2, mu = 3, InfR.p = "Opt")


result.ker = tsmoothlm(y, q.min = 0, p.min = 0, p.max = 3, q.max = 3, p = 1, InfR = "Opt", method = "kr")

g0 = result$ye
g0.der1 = result.der1$ye
g0.der2 = result.der2$ye
g.ker = result.ker$ye

res = y - g0

df = data.frame(cbind(year, g0, g0.der1, g0.der2, g.ker, res))

plot.trend <- ggplot(df, aes(x = year, y = y)) + 
  geom_line(aes(color = "Log of GDP"), size = 0.25) +
  geom_line(aes(y = g0, color = "Trend (local linear)"), size = 0.25) +
  geom_line(aes(y = g.ker, color = "Trend (kernel)"), size = 0.25, linetype = "dashed") +
  labs(title = "(a) Log of GER-GDP & estimated trends", y = "Log-GDP and trends") + 
  scale_x_continuous(name = "", breaks = seq(1850, 2020, 20)) +
  scale_color_manual(name = "Lines:",
                     breaks = c("Log of GDP", "Trend (local linear)", "Trend (kernel)"),
                     values = c("Log of GDP" = "black", "Trend (local linear)" = "red", "Trend (kernel)" = "blue")) +
  theme(legend.position = c(0.1275, 0.81), legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 6), legend.title = element_text(size = 8),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', size = 0.25),
        plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

plot.der1 <- ggplot(df, aes(x = year, y = g0.der1)) + 
  geom_line(size = 0.25) +
  labs(title = "(c) Estimated first derivative", y = "1st derivative", x = "") +
  scale_x_continuous(name = "", breaks = seq(1850, 2020, 20)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25))

plot.res <- ggplot(df, aes(x = year, y = res)) +
  geom_line(size = 0.25) + 
  labs(title = "(b) Trend-adjusted Residuals (local linear)", y = "Resiudals", x = "Year") +
  scale_x_continuous(name = "Year", breaks = seq(1850, 2020, 20)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 9),
        axis.text = element_text(size = 7), axis.ticks = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25))

plot.der2 <- ggplot(df, aes(x = year, y = g0.der2)) +
  geom_line(size = 0.25) +
  geom_hline(yintercept = 0, size = 0.25) +
  labs(title = "(d) Estimated second derivative", y = "2nd derivative", x = "Year") +
  scale_x_continuous(name = "Year", breaks = seq(1850, 2020, 20)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 9),
        axis.text = element_text(size = 7), axis.ticks = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25))

ggpubr::ggarrange(plot.trend, plot.der1, plot.res, plot.der2, heights = c(0.9, 1), nrow = 2, ncol = 2)
setwd("~/Arbeit/DFG/Paper_SEMIFAR/Latex_aktuell/Abb")
ggsave("GERgdp.pdf", height = 5, width = 9.5, dpi = 600)
