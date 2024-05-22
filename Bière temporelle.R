##### Packages ####
install.packages("zoo")
install.packages("fUnitRoots")
install.packages("tseries")
install.packages("purrr")
install.packages("urca")

library(zoo)
library(fUnitRoots)
library(tseries)
library(purrr)
library(fUnitRoots)
library(urca)  
library(ggplot2)
library(forecast)


#### PART I: THE DATA ####
##### Import & Processing #####

### Q1: Que représente la série choisie ? ###

setwd("C:/Users/mathi/Documents/ENSAE/2A/2A S2/Séries temporelles linéaires/Projet séries temporelles")
serie <- read.csv("biere.csv", header = FALSE, sep = ";", skip = 4, stringsAsFactors = FALSE)
serie <- serie[,-3]
colnames(serie) <- c("date", "indices")

serie$indices <- as.numeric(serie$indices)
str(serie)

ts_serie <- ts(serie$indices, start = c(1990, 1), end = c(2018, 12), frequency = 12)
autoplot(ts_serie) + ggtitle("Série temporelle des indices de production de la bière")


### Q2 : Rendre la série stationnaire : ###

regression <- lm(ts_serie ~ index(ts_serie))
summary(regression)
confint(regression)


adf <- adfTest(ts_serie, lag=0, type="ct")
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

# Regression model
regression <- lm(ts_serie ~ index(ts_serie))
residuals <- residuals(regression)

Qtests(adf@test$lm$residuals[0:24], length(adf@test$lm$coefficients))
# Absence of residual autocorrelation: the test is valid
adf


## Philip-Perron test:
pp.test(ts_serie)

# The unit root hypothesis is rejected by both tests, so it seems stationary, which is inconsistent with the form of the series when plotted
acf(ts_serie)
# Too many autocorrelations, the series is very persistent, which can pose a problem

# KPSS Test for the original series xm
kpss_test <- ur.kpss(ts_serie, type = "mu")
cat("KPSS Test for the original series xm:\n")
print(summary(kpss_test))

##KPSS rejects the stationarity of xm at 1%!


# On différencie la série à l'ordre 1

ts_serie_diff <- diff(ts_serie,1) # First difference
adf_diff <- adfTest(ts_serie_diff, lag=0, type="ct")

# The series is stationary

Qtests(adf_diff@test$lm$residuals[0:24], length(adf_diff@test$lm$coefficients))
# Hypothesis of no correlation not rejected, test valid
adf_diff


pp.test(ts_serie_diff)
# The unit root hypotheses are rejected, so we will say the series is stationary
acf(ts_serie_diff)

kpss_test_diff <- ur.kpss(ts_serie_diff, type = "mu")

cat("\nKPSS Test for the differenced series:\n")
summary(kpss_test_diff)

### Q3 :
autoplot(ts_serie)+ ggtitle("Série temporelle de l'IPI de la bière")
autoplot(ts_serie_diff) + ggtitle("Série temporelle différenciée de l'IPI de la bière")

#### PART II: ARMA MODELS ####

#Question 4

acf(ts_serie_diff) # The first-order autocorrelation (total or partial, it’s the same thing) is about -0.35, which is small and far from being equal to 1. The series seems stationary
pacf(ts_serie_diff)

pmax <- 9
qmax <- 2

# Estimation des autocorrélations et autocorrélations partielles
acf_diff <- acf(ts_serie_diff, lag.max = 24, main = "Autocorrélation des données différenciées")
pacf_diff <- pacf(ts_serie_diff, lag.max = 24, main = "Autocorrélation partielle des données différenciées")

# L'ACF est significative jusqu'à l'ordre 9 au maximum 
# La PACF est significative jusqu'à l'ordre 2 au maximum

auto_model <- auto.arima(ts_serie_diff)
summary(auto_model)

pmax <- 9
qmax <- 2


## On évalue toutes les ARMA(p,q) dans cette plage et on les compare avec l'AIC et BIC
mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) 
rownames(mat) <- paste0("p=",0:pmax) 
colnames(mat) <- paste0("q=",0:qmax) 
AICs <- mat #matrice des AIC à remplir
BICs <- mat #matrice des BIC à remplir
pqs <- expand.grid(0:pmax,0:qmax) #l'ens des combinaisons de p et q
for (row in 1:dim(pqs)[1]){ 
  p <- pqs[row,1] 
  q <- pqs[row,2] 
  estim <- try(arima(ts_serie_diff,c(p,0,q),include.mean = F)) #tente d’estimer l’ARIMA
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigne l’AIC
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigne le BIC
}

AICs
BICs


#Test ARMA(1,1)
arma11 <- arima(ts_serie_diff, order=c(1,0,1), include.mean=F)
arma11

#Test ARMA(2,1)
arma21 <- arima(ts_serie_diff, order=c(2,1,1), include.mean = F)
arma21


#Test ARMA(2,2)
arma22 <- arima(ts_serie_diff, order=c(2,1,2))
arma22 

#En ayant comparé les Sigma carré, la vraisemblance, on retient le modèle 2,2
# L'auto arima nous donne un ARMA(2,1) ; à voir si on change pas pour avoir un modèle plus parcimonieux


#Question 5

#On choisit un ARIMA(2,1,2)
