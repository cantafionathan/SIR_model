library(ggplot2)
library(car)
library(MASS)
library(dplyr)
data <- read.csv("pandemic_data.csv")
head(data)
summary(data)

anova <- aov(peak.inf ~ soc.iso * rate.vac * quar.dur * num.daily, data = data)

# automatically save plot to latex folder so 
# that if the plots change we easily change in latex later
pdf("latex/soc-iso:vac-rate_interaction_plot.pdf", width = 8, height = 6)
interaction.plot(
  x.factor = data$soc.iso,
  trace.factor = data$rate.vac,
  response = data$peak.inf,
  col = 1:4,
  lty = 1,
  main = "Interaction Plot: Social Isolation x Vaccination Rate",
  xlab = "Social Isolation",
  ylab = "Peak Infected",
  trace.label = "Vaccination Rate"
)
dev.off()