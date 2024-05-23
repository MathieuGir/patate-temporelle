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

residuals_arma22 <- residuals(arma22)

# Test de portemanteau pour les résidus
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
sigma_epsilon <- var_residus
Sigma <- sigma * matrix(c(1, phi1 + theta1, phi1 + theta1, (1 + phi1 + theta1)^2), nrow = 2) # Covariance Matrix
Sigma


residuals <- residuals(arma22)
residuals
residuals[1]
residuals[2]
errors <- c(residuals[2], residuals[2] + (1 + phi1 + theta1) * residuals[1])
errors

t <- t(errors) %*% solve(Sigma) %*% errors

# Calculate the quantile of the chi-squared distribution
alpha <- 0.05 # Confidence level
q <- qchisq(1 - alpha, df = 2)


# Confidence region
confidence_region <- qchisq(1 - alpha, df = 2)
confidence_region


predictions <- predict(arma22, n.ahead = 2)
pred_values <- predictions$pred
pred_se <- predictions$se
pred_interval_upper <- pred_values + qnorm(1 - alpha / 2) * pred_se
pred_interval_lower <- pred_values - qnorm(1 - alpha / 2) * pred_se


# Générer les points de l'ellipse
library(MASS)
theta <- seq(0, 2*pi, length.out=100)
circle <- cbind(cos(theta), sin(theta))


eigen_Sigma <- eigen(Sigma)
values <- eigen_Sigma$values
vectors <- eigen_Sigma$vectors
values
vectors

# Générer les points de l'ellipse
theta <- seq(0, 2*pi, length.out = 5000)
circle <- cbind(cos(theta), sin(theta))

center <- as.numeric(pred_values)
center
sqrt_values <- sqrt(values * q)
sqrt_values
transformation_matrix <- vectors %*% diag(sqrt_values)
transformation_matrix
ellipse_points <- t(transformation_matrix %*% t(circle)) + center

# Tracé
plot(ellipse_points, type = 'l', xlab = "Prédiction de X1", ylab = "Prédiction de X2", main = "Région de Confiance au Niveau Alpha",xlim = c(-20, 20), ylim = c(-50, 50))
points(center[1], center[2], col = "red", pch = 19)  # Ajouter le centre
text(center[1], center[2], labels = "valeurs prédites", pos = 3)


library(ellipse)

# Paramètres
alpha <- 0.05
mean_vector <- as.numeric(pred_values)  # Extraire le centre de l'ellipse des prédictions
Sigma

# Générer l'ellipse de confiance
conf_ellipse <- ellipse(Sigma, level = 0.95, centre = mean_vector, npoints = 10000)
plot(conf_ellipse)

# Tracé avec limites spécifiées pour les axes x et y
plot(conf_ellipse, type = "l", xlab = "X1", ylab = "X2", main = "Ellipse de Confiance au Niveau Alpha", xlim = c(-10, 10), ylim = c(-10, 10))
points(mean_vector[1], mean_vector[2], col = "red", pch = 19)  # Ajouter le centre
text(mean_vector[1], mean_vector[2], labels = expression(hat(X)), pos = 3)


## SUpplisson
# Construction de la matrice de covariance Sigma pour un modèle ARMA(2,2)
sigma_x <- sqrt(arma22$fit$sigma2)
sigma_y <- sqrt(arma22$fit$sigma2 * (1 + sum(arma22$fit$coef[1:2]^2)))
rho <- arma22$fit$sigma2 * sum(arma22$fit$coef[1:2]) / (sigma_x * sigma_y)
Sigma <- matrix(c(sigma_x^2, rho * sigma_x * sigma_y, rho * sigma_x * sigma_y, sigma_y^2), nrow = 2)

# Calcul des valeurs propres et vecteurs propres de la matrice de covariance
eigen_Sigma <- eigen(Sigma)
values <- eigen_Sigma$values
vectors <- eigen_Sigma$vectors

# Génération de l'ellipse de confiance
library(ellipse)
ell <- ellipse(rho, scale = c(sigma_x, sigma_y), centre = c(Prevision$pred[1], Prevision$pred[2]), level = 0.95, npoints = 10000, xlab = "Prévision à T+1", ylab = "Prévision à T+2")

# Tracé de l'ellipse
plot(ell)
points(x = Prevision$pred[1], y = Prevision$pred[2], type = "p", lwd = 7)
