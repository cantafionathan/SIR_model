library(ggplot2)
library(car)
library(MASS)
library(dplyr)
data <- read.csv("pandemic_data.csv")
head(data)
summary(data)

# Anova
anova_model <- aov(peak.inf ~ soc.iso * rate.vac * quar.dur * num.daily, data = data)
summary(anova_model)

# Regression and Box-Cox
lm_model <- lm(peak.inf ~ soc.iso * rate.vac * quar.dur * num.daily, data = data)
summary(lm_model)
confint(lm_model, level = 0.95)

# Box-Cox Transformation
pdf("latex/boxcox.pdf", width = 8, height = 6)
boxcox_result <- boxcox(lm_model, lambda = seq(-2, 2, by = 0.1))
dev.off()
optimal_lambda <- boxcox_result$x[which.max(boxcox_result$y)]
optimal_lambda # no transformation

# Residual plots
residuals_plot <- ggplot(data.frame(fitted = lm_model$fitted.values, 
                                    residuals = lm_model$residuals), 
                         aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted Values", y = "Residuals")

pdf("latex/residuals.pdf", width = 8, height = 6)
print(residuals_plot)
dev.off()

# Q-Q Plot
pdf("latex/qqplot.pdf", width = 8, height = 6)
qqnorm(lm_model$residuals)
qqline(lm_model$residuals)
dev.off()


# Graphics
# Interaction plots
pdf("latex/interactionplot1.pdf", width = 8, height = 6)
interaction.plot(
  x.factor = data$soc.iso,
  trace.factor = data$rate.vac,
  response = data$peak.inf,
  col = 1:4,
  lty = 1,
  xlab = "Social Isolation",
  ylab = "Peak Infected",
  trace.label = "Vaccination Rate"
)
dev.off()
pdf("latex/interactionplot2.pdf", width = 8, height = 6)
interaction.plot(
  x.factor = data$soc.iso,
  trace.factor = data$quar.dur,
  response = data$peak.inf,
  col = 1:4,
  lty = 1,
  xlab = "Social Isolation",
  ylab = "Peak Infected",
  trace.label = "Quarantine Duration"
)
dev.off()
pdf("latex/interactionplot3.pdf", width = 8, height = 6)
interaction.plot(
  x.factor = data$soc.iso,
  trace.factor = data$num.daily,
  response = data$peak.inf,
  col = 1:4,
  lty = 1,
  xlab = "Social Isolation",
  ylab = "Peak Infected",
  trace.label = "Daily Interactions"
)
dev.off()
pdf("latex/interactionplot4.pdf", width = 8, height = 6)
interaction.plot(
  x.factor = data$rate.vac,
  trace.factor = data$quar.dur,
  response = data$peak.inf,
  col = 1:4,
  lty = 1,
  xlab = "Vaccination Rate",
  ylab = "Peak Infected",
  trace.label = "Quarantine Duration"
)
dev.off()
pdf("latex/interactionplot5.pdf", width = 8, height = 6)
interaction.plot(
  x.factor = data$rate.vac,
  trace.factor = data$num.daily,
  response = data$peak.inf,
  col = 1:4,
  lty = 1,
  xlab = "Vaccination Rate",
  ylab = "Peak Infected",
  trace.label = "Daily Interactions"
)
dev.off()
pdf("latex/interactionplot6.pdf", width = 8, height = 6)
interaction.plot(
  x.factor = data$quar.dur,
  trace.factor = data$num.daily,
  response = data$peak.inf,
  col = 1:4,
  lty = 1,
  xlab = "Quarantine Duration",
  ylab = "Peak Infected",
  trace.label = "Daily Interactions"
)
dev.off()

# Boxplots
pdf("latex/boxplot1.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(soc.iso), y = peak.inf)) +
  geom_boxplot() +
  labs(x = "Social Isolation Level", y = "Peak Infected")
dev.off()

pdf("latex/boxplot2.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(rate.vac), y = peak.inf)) +
  geom_boxplot() +
  labs(x = "Vaccination Rate", y = "Peak Infected")
dev.off()

pdf("latex/boxplot3.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(quar.dur), y = peak.inf)) +
  geom_boxplot() +
  labs(x = "Quarantine Duration (Days)", y = "Peak Infected")
dev.off()

pdf("latex/boxplot4.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(num.daily), y = peak.inf)) +
  geom_boxplot() +
  labs(x = "Daily Interactions", y = "Peak Infected")
dev.off()



#Optimal Combination
aggregated_results <- aggregate(peak.inf ~ soc.iso + rate.vac + quar.dur + num.daily, data = data, mean)
optimal_combination <- aggregated_results[which.min(aggregated_results$peak.inf), ]
optimal_combination

# Treatment Contrasts
data$soc.iso <- as.factor(data$soc.iso)
data$rate.vac <- as.factor(data$rate.vac)
data$quar.dur <- as.factor(data$quar.dur)
data$num.daily <- as.factor(data$num.daily)
head(data)

contrasts(data$soc.iso) <- contr.treatment(levels(data$soc.iso), base = 1) 
contrasts(data$rate.vac) <- contr.treatment(levels(data$rate.vac), base = 1) 
contrasts(data$quar.dur) <- contr.treatment(levels(data$quar.dur), base = 1) 
contrasts(data$num.daily) <- contr.treatment(levels(data$num.daily), base = 1) 

# Fit linear model with treatment contrasts
lm_treatment <- lm(peak.inf ~ soc.iso + rate.vac + quar.dur + num.daily, data = data)
summary(lm_treatment)

confint(lm_treatment, level = 0.95)
coef_lm <- coef(lm_treatment)
coef_lm
confint_lm <- confint(lm_treatment)
confint_lm