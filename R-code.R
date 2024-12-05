library(ggplot2)
library(car)
library(MASS)
library(dplyr)
data <- read.csv("pandemic_data.csv")
head(data)
summary(data)

# Anova
anova_model <- aov(load ~ soc.iso * rate.vac * quar.dur * num.daily, data = data)
summary(anova_model)

# Regression and Box-Cox
lm_model <- lm(load ~ soc.iso * rate.vac * quar.dur * num.daily, data = data)
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
# Boxplots
pdf("latex/boxplot_soc.iso.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(soc.iso), y = load)) +
  geom_boxplot() +
  labs(x = "Social Isolation Level", y = "Load")
dev.off()

pdf("latex/boxplot_vac.rate.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(rate.vac), y = load)) +
  geom_boxplot() +
  labs(x = "Vaccination Rate", y = "Load")
dev.off()

pdf("latex/boxplot_quar.dur.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(quar.dur), y = load)) +
  geom_boxplot() +
  labs(x = "Quarantine Duration (Days)", y = "Load")
dev.off()

pdf("latex/boxplot_num.daily.pdf", width = 8, height = 6)
ggplot(data, aes(x = as.factor(num.daily), y = load)) +
  geom_boxplot() +
  labs(x = "Daily Interactions", y = "Load")
dev.off()



# Optimal Combination
aggregated_results <- aggregate(load ~ soc.iso + rate.vac + quar.dur + num.daily, data = data, mean)
aggregated_results
# Optimal Levels for each level of daily interaction
optimal_levels <- aggregated_results %>%
  group_by(num.daily) %>%
  slice_min(load, n = 1, with_ties = FALSE)
print(optimal_levels)

# For the most optimal combination of all factors levels, we somehow see 7 day quarantine
# producing least load (less than 14 day quarantine)
# 1	0.02	14	15	0.01031777
# 1	0.02	7	15	0.01001565
optimal_combination <- aggregated_results[which.min(aggregated_results$load), ]
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
lm_treatment <- lm(load ~ soc.iso + rate.vac + quar.dur + num.daily, data = data)
summary(lm_treatment)

confint(lm_treatment, level = 0.95)
coef_lm <- coef(lm_treatment)
coef_lm
confint_lm <- confint(lm_treatment)
confint_lm


