setwd("C:/Users/mathi/Documents/ENSAE/2A/2A S2/Séries temporelles linéaires/Projet séries temporelles")
data <- read.csv("pommes_de_terre.csv", header = FALSE, sep = ";", skip = 52, stringsAsFactors = FALSE)
#On skip les 51 premières lignes, qui correspondent à la fois au header et aux données postérieures à décembre 2019 pour éviter les données pendant le Covid

colnames(data) <- c("dates", "indices")
head(data)

# Partie 1 :

ts_data <- ts(indices, start = c(1997, 1), end = c(2019, 12), frequency = 12)
print(ts_data)

library(ggplot2)
library(forecast)

autoplot(ts_data) + ggtitle("Série temporelle des indices de production des pommes de terre")

# Différenciation à l'ordre 1
diff_ts_data <- diff(ts_data,1)
autoplot(diff_ts_data) + ggtitle("Série temporelle différenciée des indices de production des pommes de terre")

# Test de Durbin Watson
library(car)

model <- lm(indices ~ dates, data=data)
dwtest <- durbinWatsonTest(model)
dwtest

# Fonction Qtests pour tester l'autocorrélation

require(fUnitRoots) #tests de racine unitaire plus modulables

adfTest(ts_data, lag=0, type="ct")




# Test de l'autocorrélation des résidus jusqu'à l'ordre 24

# Ajuster votre modèle de régression (remplacez "modele" par votre modèle de régression)
modele <- lm(indices ~ dates, data=data)

residus <- residuals(modele)
print(residus)
print(modele)

autocorr <- acf(residus, plot=FALSE)
plot(autocorr, main="Autocorrelogramme des résidus")


# Test de Dickey-Fuller sur la série différenciée

library(tseries)
adf_test <- adf.test(diff_ts_data)
print(adf_test)

library(astsa)

# Estimation des autocorrélations et autocorrélations partielles
acf_diff <- acf(diff_ts_data, lag.max = 24, main = "Autocorrélation des données différenciées")
pacf_diff <- pacf(diff_ts_data, lag.max = 24, main = "Autocorrélation partielle des données différenciées")

# L'ACF est significative jusqu'à l'ordre 3 au maximum
# La PACF est significative jusqu'à l'ordre 12 au maximum

# On a choisi le modèle ARMA(3,12)
