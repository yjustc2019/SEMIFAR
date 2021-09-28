library(ggplot2)
library(smootslm)

### downloading data from yahoo.finance
ticker = c("GDAXI", "GSPC", "DJI", "HSI", "N225", "FTSE", "STOXX50E")#SP500 is ganz gut #N225 und FTSE fehlen noch
idx = 2
start_date = "1990-01-01"
end_date = "2020-12-31"#Sys.Date()

data = quantmod::getSymbols(paste0("^",ticker[idx]), from = start_date, to = end_date, warnings = FALSE, auto.assign = FALSE)
IDX = data[which(!is.na(data[, paste0(ticker[idx],".Close")]))] ### deleting missing values
Clo = as.vector(IDX[, paste0(ticker[idx],".Close")]) ### Closing Price


ret = diff(log(Clo))
retc = ret - mean(ret)
lret = log(retc^2)
n = length(lret)

result1 = tsmoothlm(lret, p = 1, pmin = 1, pmax = 1, qmin = 1, qmax = 1, InfR = "Opt")
result3 = tsmoothlm(lret, p = 3, pmin = 1, qmin = 1, pmax = 1, qmax = 1, InfR = "Opt")

ye1 = result1$ye
ye3 = result3$ye
res1 = result1$res
res3 = result3$res

year = (1:n) / n * 31 + 1990



FARIMA11 <- fracdiff::fracdiff(res3, nar = result3$p.BIC, nma = result3$q.BIC)
Zt = FARIMA11$fitted

Csig = sd(retc / exp(ye3 / 2))
ye3sd = exp(ye3 / 2) * Csig 
res3sd = retc / ye3sd
h_star = exp(Zt)
Ch = var(res3sd / sqrt(h_star))    
cond.vol3 = sqrt(Ch * h_star)      ### nach Chefs Paper
tot.vol3 = cond.vol3 * ye3sd


## nach Smoots-Paper
# mulz = log(1/var(retc / exp(ye1/2)))
# mueps = -log(mean(exp(FARIMA11$residuals))) ## nach smoots Paper
# cond.vol = sqrt(exp(Zt + mulz - mueps))## nach smoots Paper


df = data.frame(year, retc, lret, ye1, ye3, res1, res3, cond.vol3, tot.vol3)

plot.ret = ggplot(data = df, aes(x = year, y = retc)) + 
  geom_line() + 
  labs(title = "(a) Centralized S&P 500 return series", y = "Returns") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25))

plot.trend = ggplot(data = df, aes(x = year, y = lret)) +
  geom_line(aes(color = "Log-data")) +
  geom_line(aes(y = ye3, color = "Trend 1 (local cubic)")) +
  geom_line(aes(y = ye1, color = "Trend 2 (local linear)"), linetype = "dashed") +
  labs(title = "(b) Log-transformed returns & estimated trends", y = "Log-data § trends") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  scale_color_manual(name = "Lines:",
                     breaks = c("Log-data", "Trend 1 (local cubic)", "Trend 2 (local linear)"),
                     values = c("Log-data" = "darkgrey", "Trend 1 (local cubic)" = "red", "Trend 2 (local linear)" = "blue")) +
  theme(legend.position = c(0.5, 0.1), legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', size = 0.25), 
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

plot.condv = ggplot(data = df, aes(x = year, y = cond.vol3)) +
  geom_line() +
  labs(title = "(c) Conditional volatility (obtained with trend 1)", y = "Cond. volatility") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25))

plot.totv = ggplot(data = df, aes(x = year, y = tot.vol3)) + 
  geom_line() +
  labs(title = "(d) Total volatility (obtained with trend 1)",y = "Total volatility", x = "Year") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text = element_text(size = 7), axis.ticks = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25))

SP500 = ggpubr::ggarrange(plot.ret, plot.trend, plot.condv, plot.totv, ncol = 1)
setwd("~/Arbeit/DFG/Paper_SEMIFAR/Latex_aktuell/Abb")
ggsave("SP500.pdf", plot = SP500, height = 9, width = 11, dpi = 600)

