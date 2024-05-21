setwd("C:/Users/mathi/Documents/ENSAE/2A/2A S2/Séries temporelles linéaires/Projet séries temporelles")
data <- read.csv("pommes_de_terre.csv", header = TRUE, sep = ";", skip = 4, stringsAsFactors = FALSE)

head(data)

dates <- data[, 1]
indices <- as.numeric(data[, 2])  
indices <- na.omit(indices)

print(head(dates))
print(head(indices))


# Partie 1 :

ts_data <- ts(indices, start = c(1990, 1), frequency = 12)

print(ts_data)

library(ggplot2)
library(forecast)

autoplot(ts_data) + ggtitle("Série temporelle des indices de production des pommes de terre")

# Différenciation simple
diff_ts_data <- diff(ts_data)
autoplot(diff_ts_data) + ggtitle("Série temporelle différenciée des indices de production des pommes de terre")

# Test de Dickey-Fuller sur la série différenciée

library(tseries)
adf_test <- adf.test(diff_ts_data)
print(adf_test)

library(astsa)

# Estimation des autocorrélations et autocorrélations partielles
acf_diff <- acf(diff_ts_data, lag.max = 20)
pacf_diff <- pacf(diff_ts_data, lag.max = 20)
