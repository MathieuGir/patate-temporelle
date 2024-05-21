setwd("C:/Users/mathi/Documents/ENSAE/2A/2A S2/Séries temporelles linéaires/Projet séries temporelles")
serie <- read.csv("pommes_de_terre.csv", header = FALSE, sep = ";", skip = 52, stringsAsFactors = FALSE)
#On skip les 51 premières lignes, qui correspondent à la fois au header et aux données postérieures à décembre 2019 pour éviter les données pendant le Covid

library(zoo)

colnames(serie) <- c("date", "indices")
serie$date <- as.yearmon(serie$date, "%Y-%m")
str(serie)

# Partie 1 :

ts_serie <- ts(serie$indices, start = c(2005, 1), end = c(2019, 12), frequency = 12)
print(ts_serie)

library(ggplot2)
library(forecast)

autoplot(ts_serie) + ggtitle("Série temporelle des indices de production des pommes de terre")

# Différenciation à l'ordre 1
diff_ts_serie <- diff(ts_serie,1)
autoplot(diff_ts_serie) + ggtitle("Série temporelle différenciée des indices de production des pommes de terre")

# Fonction de test de l'autocorrélation des résidus

model <- lm(indices ~ date, data=serie)
residuals <- residuals(model)
test_result <- Box.test(residuals, lag = 1, type = "Ljung-Box")
print(test_result)
print(ts_serie)

require(fUnitRoots)
adf <- adfTest(ts_serie, lag=0, type="ct")
print(adf)

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}
residuals <- residuals(adf@test$lm)

Qtests(residuals, length(adf@test$lm$coefficients))

adfTest_valid <- function(series,kmax,type){ #tests ADF jusqu’`a des r´esidus non autocorr´el´es
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("nope \n")
    k <- k + 1
  }
  return(adf)
}

# Test sur la série non différenciée (lag de 12)

adf <- adfTest_valid(ts_serie,24,"ct")
adf

#Maintenant pour la série différenciée (lag de 6)
adf <- adfTest_valid(diff_ts_serie,24, type="nc")
adf

# Test de Dickey-Fuller sur la série différenciée
library(tseries)
adf_test <- adf.test(diff_ts_serie)
print(adf_test)

library(astsa)

##PARTIE II

#Question 4

# Estimation des autocorrélations et autocorrélations partielles
acf_diff <- acf(diff_ts_serie, lag.max = 24, main = "Autocorrélation des données différenciées")
pacf_diff <- pacf(diff_ts_serie, lag.max = 24, main = "Autocorrélation partielle des données différenciées")

# L'ACF est significative jusqu'à l'ordre 5 au maximum
# La PACF est significative jusqu'à l'ordre 7 au maximum

# On a choisi le modèle ARMA(5,7)

library(stats)

arma_model <- arima(diff_ts_serie, order=c(5,0,7))

summary(arma_model)

arima(x = diff_ts_serie, order = c(5,0,7), include.mean=T)


#Question 5

