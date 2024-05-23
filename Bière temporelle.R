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
library(ellipse)


#### PART I: THE DATA ####
##### Import & Processing #####

### Q1: Que représente la série choisie ? ###

setwd("C:/Users/mathi/Documents/ENSAE/2A/2A S2/Séries temporelles linéaires/Projet séries temporelles")
serie <- read.csv("biere.csv", header = FALSE, sep = ";", skip = 4, stringsAsFactors = FALSE)
serie <- serie[,-3]
colnames(serie) <- c("date", "indices")

serie$indices <- as.numeric(serie$indices)
str(serie)

ts_serie <- ts(serie$indices, start = c(1990, 1), end = c(2019, 12), frequency = 12)
autoplot(ts_serie) + ggtitle("Série temporelle des indices de production de la bière")


### Q2 : Rendre la série stationnaire : ###

regression <- lm(ts_serie ~ index(ts_serie))
summary(regression)
confint(regression)
# Le coefficient associé à la tendance linéaire n'est spas significatif. 
# On fait donc des tests de racine unitaire avec une constante 

adf <- adfTest(ts_serie, lag=0, type="ct")
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

# Régression
regression <- lm(ts_serie ~ index(ts_serie))
residuals <- residuals(regression)

Qtests(adf@test$lm$residuals[0:24], length(adf@test$lm$coefficients))
# Le test conclut à l'absence d'autocorrelation résiduelle 
adf

## Test de Philip-Perron :
pp.test(ts_serie)

acf(ts_serie, lag = 24)
# Les deux tests rejettent l'hypothèse de racine unitaire : la série devrait donc être stationnaire.
# Mais en regardant les autocorrélations, même à lag 24, les résidus sont auocorrélés

kpss_test <- ur.kpss(ts_serie, type = "mu")
summary(kpss_test)

## Le test de KPSS rejette la stationarité de la série à 1%

# On va donc différencier la série à l'ordre 1

ts_serie_diff <- diff(ts_serie,1) # First difference
autoplot(ts_serie_diff) + ggtitle("Série temporelle différenciée de l'IPI de la bière")

# La série semble être stationnaire 

Qtests(adf_diff@test$lm$residuals[0:24], length(adf_diff@test$lm$coefficients))

adf_diff <- adfTest(ts_serie_diff, lag=0, type="ct")
adf_diff

pp.test(ts_serie_diff)

# L'hypothèse de racine unitaire est rejetée : on considère que la série est bien stationnaire 
acf(ts_serie_diff)

kpss_test_diff <- ur.kpss(ts_serie_diff, type = "mu")
summary(kpss_test_diff)

# On ne rejette pas à 10% l'hypothèse nulle, qui est que la série est stationnaire.
# La série différenciée à l'ordre 1 est donc stationnaire

### Q3 :
autoplot(ts_serie)+ ggtitle("Série temporelle de l'IPI de la bière")
autoplot(ts_serie_diff) + ggtitle("Série temporelle différenciée de l'IPI de la bière")

#### PART II: ARMA MODELS ####

#Question 4

# Estimation des autocorrélations et autocorrélations partielles
acf_diff <- acf(ts_serie_diff, lag.max = 24, main = "Autocorrélation des données différenciées")
pacf_diff <- pacf(ts_serie_diff, lag.max = 24, main = "Autocorrélation partielle des données différenciées")

# L'ACF est significative jusqu'à l'ordre 9 au maximum 
# La PACF est significative jusqu'à l'ordre 2 au maximum


pmax <- 9
qmax <- 2


# On va évaluer tous les ARMA(p,q) dans cette plage et les comparer avec l'AIC et BIC
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
arma21 <- arima(ts_serie_diff, order=c(2,0,1), include.mean = F)
arma21

#Test ARMA(2,2)
arma22 <- arima(ts_serie_diff, order=c(2,0,2))
arma22 

#En ayant comparé les Sigma carré, la vraisemblance, on retient le modèle ARMA(2,2)

#On fait le test de portmaneau pour les résidus
residuals_arma22 <- residuals(arma22)
ljung_box_test <- Box.test(residuals_arma22, lag = 12, type = "Ljung-Box")
print(ljung_box_test)
acf(residuals_arma22)
autoplot(residuals_arma22)


#Question 5

#On choisit un ARIMA(2,1,2)


### Partie III

var_res <- arma22$sigma2
phi1 <- -arma22$coef[1] #coeff phi1
theta1 <- -arma22$coef[3] #coeff theta1
sigma_eps <- var_res
Sigma <- sigma_eps * matrix(c(1, phi1 + theta1, phi1 + theta1, (1 + phi1 + theta1)^2), nrow = 2)
Sigma #matrice de variance covariance

residuals <- residuals(arma22)
residuals
residuals[1]
residuals[2]
errors <- c(residuals[2], residuals[2] + (1 + phi1 + theta1) * residuals[1])
errors

# On calcule le quantile d'ordre 1-alpha d'un chi2
alpha <- 0.05 # Confidence level
q <- qchisq(1 - alpha, df = 2)


# Région de confiance : 

predictions <- predict(arma22, n.ahead = 2)
pred_values <- predictions$pred


# Générer l'ellipse de confiance

mean_vector <- as.numeric(pred_values)  # centre de l'ellipse = prédictions
conf_ellipse <- ellipse(Sigma, level = 0.95, centre = mean_vector, npoints = 10000)

# Tracé avec limites spécifiées pour les axes x et y
plot(conf_ellipse, type = "l", xlab = "X1", ylab = "X2", main = "Région de confiance au niveau 95%", xlim = c(-20, 20), ylim = c(-50, 50))
points(mean_vector[1], mean_vector[2], col = "red", pch = 19)  # Ajouter le centre
text(mean_vector[1], mean_vector[2], labels = expression(hat(X)), pos = 3)
