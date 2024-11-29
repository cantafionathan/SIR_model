pandemic <- read.csv("fractional_pandemic.csv")
summary(aov(peak_infected ~ ., data = pandemic))
MASS::boxcox(peak_infected ~., data = pandemic)