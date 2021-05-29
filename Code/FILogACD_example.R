library(ggplot2)
library(smootslm)

### downloading data from yahoo.finance
ticker = c("GDAXI", "GSPC", "DJI", "HSI", "N225", "FTSE", "STOXX50E", "VIX")#SP500 is gut
idx = 2
start_date = "2000-01-01"
end_date = "2020-12-31"#Sys.Date()

data = quantmod::getSymbols(paste0("^",ticker[idx]), from = start_date, to = end_date, warnings = FALSE, auto.assign = FALSE)
IDX = data[which(!is.na(data[, paste0(ticker[idx],".Volume")])
                 & data[, paste0(ticker[idx],".Volume")] != 0)] ### deleting missing and zero values
Clo = as.vector(IDX[, paste0(ticker[idx],".Volume")]) ### Closing Price

retc = Clo
lret = log(Clo)
n = length(lret)

result1 = tsmoothlm(lret, p = 1, pmax = 3, qmax = 3, InfR = "Opt")
result3 = tsmoothlm(lret, p = 3, pmax = 3, qmax = 3, InfR = "Nai")

ye1 = result1$ye
ye3 = result3$ye
res1 = result1$res
res3 = result3$res

year = (1:n) / n * 21 + 2000

FARIMA11 <- fracdiff::fracdiff(res1, nar = result1$p.BIC, nma = result1$q.BIC)
FARIMA13 <- fracdiff::fracdiff(res3, nar = result3$p.BIC, nma = result3$q.BIC)
Zt1 = FARIMA11$fitted
Zt3 = FARIMA13$fitted

Csig = sqrt(mean(retc / exp(ye1)))
ye1sd = exp(ye1) * Csig 
res1sd = retc / ye1sd
h_star = exp(Zt1)
Ch = mean(res1sd / sqrt(h_star))    
cond.mean1 = Ch * h_star      ### nach Chefs Paper
tot.mean1 = cond.mean1 * ye1sd

Csig = sqrt(mean(retc / exp(ye3)))
ye3sd = exp(ye3) * Csig 
res3sd = retc / ye3sd
h_star = exp(Zt3)
Ch = mean(res3sd / sqrt(h_star))    
cond.mean3 = Ch * h_star      ### nach Chefs Paper
tot.mean3 = cond.mean3 * ye3sd

df = data.frame(year, retc, lret, ye1, ye3, res1, res3, cond.mean1, tot.mean1)

plot.ret = ggplot(data = df, aes(x = year, y = retc)) + 
  geom_line() + 
  labs(title = "(a) SP500 trading volume", y = "Cum. Volume") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25))

plot.trend = ggplot(data = df, aes(x = year, y = lret)) +
  geom_line(aes(color = "Log-data")) +
  geom_line(aes(y = ye1, color = "Trend 1 (local linear)")) +
  geom_line(aes(y = ye3, color = "Trend 2 (local cubic)"), linetype = "dashed") +
  labs(title = "(b) Log-transformed data & estimated trends", y = "Log-data & trends") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  scale_color_manual(name = "Lines:",
                     breaks = c("Log-data", "Trend 1 (local linear)", "Trend 2 (local cubic)"),
                     values = c("Log-data" = "darkgrey", "Trend 1 (local linear)" = "red", "Trend 2 (local cubic)" = "blue")) +
  theme(legend.position = c(0.5, 0.11), legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 8), legend.title = element_text(size = 10),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid', size = 0.25), 
        legend.direction = "horizontal",
        plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25)) +
  guides(color = guide_legend(override.aes = list(size = 0.2)))

plot.condv = ggplot(data = df, aes(x = year, y = cond.mean1)) +
  geom_line() +
  labs(title = "(c) Conditional means (obtained with local linear trend)", y = "Conditional means") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 7), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.ticks = element_line(size = 0.25), panel.grid.major = element_line(size = 0.25))

plot.totv = ggplot(data = df, aes(x = year, y = tot.mean1)) + 
  geom_line() +
  labs(title = "(d) Total means (obtained with local linear trend)",y = "Total means", x = "Year") +
  scale_x_continuous(name = "", breaks = seq(1990, 2021, 1)) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.y = element_text(size = 9), axis.title.x = element_blank(), 
        axis.text = element_text(size = 7), axis.ticks = element_line(size = 0.25),
        panel.grid.major = element_line(size = 0.25))

VOL = egg::ggarrange(plot.ret, plot.trend, plot.condv, plot.totv, ncol = 1)
setwd("~/Arbeit/DFG/Paper_SEMIFAR/Latex_aktuell/Abb")
ggsave("SP500VOL.pdf", plot = VOL, height = 9, width = 11, dpi = 600)

