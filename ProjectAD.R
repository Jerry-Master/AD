
library("lmtest")

#################Validation#################################
validation=function(model,dades){
  s=frequency(get(model$series))
  resid=model$residuals
  par(mfrow=c(2,2),mar=c(3,3,3,3))
  #Residuals plot
  plot(resid,main="Residuals")
  abline(h=0)
  abline(h=c(-3*sd(resid),3*sd(resid)),lty=3,col=4)
  #Square Root of absolute values of residuals (Homocedasticity)
  scatter.smooth(sqrt(abs(resid)),main="Square Root of Absolute residuals",
                 lpars=list(col=2))

  #ACF & PACF of square residuals 
  acf(resid^2,ylim=c(-1,1),lag.max=60,col=c(2,rep(1,s-1)),lwd=1)
  pacf(resid^2,ylim=c(-1,1),lag.max=60,col=c(rep(1,s-1),2),lwd=1)
  
  par(mfrow=c(1,2))
  #Normal plot of residuals
  qqnorm(resid)
  qqline(resid,col=2,lwd=2)
  
  ##Histogram of residuals with normal curve
  hist(resid,breaks=20,freq=FALSE)
  curve(dnorm(x,mean=mean(resid),sd=sd(resid)),col=2,add=T)
  
  
  #ACF & PACF of residuals
  par(mfrow=c(1,2))
  acf(resid,ylim=c(-1,1),lag.max=100,col=c(2,rep(1,s-1)),lwd=1)
  pacf(resid,ylim=c(-1,1),lag.max=100,col=c(rep(1,s-1),2),lwd=1)
  par(mfrow=c(1,1))
  
  
  
  #Ljung-Box p-values
  par(mar=c(2,2,1,1))
  tsdiag(model,gof.lag=7*s)
  cat("\n--------------------------------------------------------------------\n")
  print(model)
  
  #Stationary and Invertible
  cat("\nModul of AR Characteristic polynomial Roots: ", 
      Mod(polyroot(c(1,-model$model$phi))),"\n")
  cat("\nModul of MA Characteristic polynomial Roots: ",
      Mod(polyroot(c(1,model$model$theta))),"\n")
  
  #Model expressed as an MA infinity (psi-weights)
  psis=ARMAtoMA(ar=model$model$phi,ma=model$model$theta,lag.max=36)
  names(psis)=paste("psi",1:36)
  cat("\nPsi-weights (MA(inf))\n")
  cat("\n--------------------\n")
  print(psis[1:20])
  
  #Model expressed as an AR infinity (pi-weights)
  pis=-ARMAtoMA(ar=-model$model$theta,ma=-model$model$phi,lag.max=36)
  names(pis)=paste("pi",1:36)
  cat("\nPi-weights (AR(inf))\n")
  cat("\n--------------------\n")
  print(pis[1:20])
  
  cat("\nNormality Tests\n")
  cat("\n--------------------\n")
  
  ##Shapiro-Wilks Normality test
  print(shapiro.test(resid(model)))
  
  suppressMessages(require(nortest,quietly=TRUE,warn.conflicts=FALSE))
  ##Anderson-Darling test
  print(ad.test(resid(model)))
  
  suppressMessages(require(tseries,quietly=TRUE,warn.conflicts=FALSE))
  ##Jarque-Bera test
  print(jarque.bera.test(resid(model)))
  
  cat("\nHomoscedasticity Test\n")
  cat("\n--------------------\n")
  suppressMessages(require(lmtest,quietly=TRUE,warn.conflicts=FALSE))
  ##Breusch-Pagan test
  obs=get(model$series)
  print(bptest(resid(model)~I(obs-resid(model))))
  
  cat("\nIndependence Tests\n")
  cat("\n--------------------\n")
  
  ##Durbin-Watson test
  print(dwtest(resid(model)~I(1:length(resid(model)))))
  
  ##Ljung-Box test
  cat("\nLjung-Box test\n")
  print(t(apply(matrix(c(1:4,(1:4)*s)),1,function(el) {
    te=Box.test(resid(model),type="Ljung-Box",lag=el)
    c(lag=(te$parameter),statistic=te$statistic[[1]],p.value=te$p.value)})))
  
  
  #Sample ACF vs. Teoric ACF
  par(mfrow=c(2,2),mar=c(3,3,3,3))
  acf(dades, ylim=c(-1,1) ,lag.max=36,main="Sample ACF")
  
  plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36),ylim=c(-1,1), 
       type="h",xlab="Lag",  ylab="", main="ACF Teoric")
  abline(h=0)
  
  #Sample PACF vs. Teoric PACF
  pacf(dades, ylim=c(-1,1) ,lag.max=36,main="Sample PACF")
  
  plot(ARMAacf(model$model$phi,model$model$theta,lag.max=36, pacf=T),ylim=c(-1,1),
       type="h", xlab="Lag", ylab="", main="PACF Teoric")
  abline(h=0)
  par(mfrow=c(1,1))
}
################# Fi Validaci? #################################

prevision <- function(model, data) {
  ultim <- c(2018,11)
  data_2 <-window(data, end=ultim)
  
  mod2 <- update(model, x = data_2)
  
  #prediction 12 next observations
  pred <- predict(mod2, n.ahead=12)
  point <- exp(pred$pred)
  
  up <- point*exp(pred$se * 1.96)
  lo <- point/exp(pred$se * 1.96)
  
  ts.plot(data,lo,up,point,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2013,2019), type="o",
          main="Predicciones de 2018 según el modelo 1") 
  abline(v=2013:2019,lty=3,col=4)
  
  #RMSPE and MAPE
  obs <- window(data, start=c(2018, 1), end=c(2018, 12))
  
  print(paste0("RMSPE: ", RMSPE<-sqrt(mean(((obs-point)/obs)^2))))
  print(paste0("MAPE: ", MAPE<-mean(abs(obs-point)/obs)))
  
  print(paste0("Mean conf. interval length: ", mean(up-lo)))
}


setwd("~/Dropbox/CFIS/q4/AnalisisDatos/Lab/AD/")

peajes <- ts(read.table("./peajes.dat"), start = 1990, frequency = 12)

plot(peajes, main="Intensidad Media Diaria en peajes en España", ylab="miles de vehiculos")
abline(v=1990:2020,lty=3,col=4)
var(peajes)
monthplot(peajes)

# varianza constante ?
boxplot(peajes~floor(time(peajes)), main="peajes")
plot.mean.var <- function(serie, title){
  m=apply(matrix(serie,nr=frequency(serie)),2,mean)
  v=apply(matrix(serie,nr=frequency(serie)),2,var)
  plot(m,v,xlab="Medias anuales",ylab="Varianzas anuales",main=title, asp = 1)
  abline(lm(v~m),col=2,lty=3,lwd=2)
}

plot.mean.var(peajes, "peajes")

boxplot(log(peajes)~floor(time(peajes)), main="peajes")
plot.mean.var(log(peajes), "peajes")

#hacemos transformacion logaritmica

peajes_log <- log(peajes)

#diferenciacion regular orden 12
plot(peajes_12 <- diff(peajes_log, 12))
abline(h = 0, lty = 2, col = 4)
(var_12 <- var(peajes_12))
monthplot(peajes_12)
abline(v=1990:2020,lty=2,col=4)


# dif reg 1
plot(peajes_12_1 <- diff(peajes_12, 1))
abline(h = 0, lty = 2, col = 4)
(var_12_1 <- var(peajes_12_1))

#dif reg 2
plot(peajes_12_2 <- diff(peajes_12_1, 1))
abline(h = 0, lty = 2, col = 4)
(var_12_2 <- var(peajes_12_2))


# vemos si es mejor peajes_12_1 o peajes_12 mirando intercept

(arima(peajes_12, order=c(4,0,0), seasonal=list(order=c(2,0,0), period=12)))

#serie definitiva
serie <- peajes_12


#ACF & PACF of residuals
s=frequency(serie)
par(mfrow=c(1,2))
acf(serie,ylim=c(-1,1),lag.max=100,col=c(2,rep(1,s-1)),lwd=1, main = "")
pacf(serie,ylim=c(-1,1),lag.max=100,col=c(rep(1,s-1),2),lwd=1, main = "")
par(mfrow=c(1,1))

# modelos propuestos
(mod1 <- arima(peajes_log, order=c(4,0,0), seasonal=list(order=c(2,1,1), period=12)))
(mod2 <- arima(peajes_log, order=c(1,0,1), seasonal=list(order=c(2,1,1), period=12)))


# validacion (SALEN MAL)
validation(mod1, peajes_log)
validation(mod2, peajes_log)

AIC(mod1)
BIC(mod1)
AIC(mod2)
BIC(mod2)


# PREDICCIONES CON MOD 1
ultim <- c(2018,11)
data_2 <-window(peajes_log, end=ultim) 

(mod1_2 <- update(mod1, x = data_2))

#prediction 12 next observations
pred <- predict(mod1_2, n.ahead=12)
point <- pred$pred

up <- point + (pred$se * 1.96)
lo <- point - (pred$se * 1.96)

ts.plot(peajes_log,lo,up,point,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2013,2020), type="o",
        main="") 
abline(v=2013:2019,lty=3,col=4)

#RMSPE and MAPE
obs <- window(peajes_log, start=c(2018, 12), end=c(2019, 11))

print(paste0("RMSPE: ", RMSPE<-sqrt(mean(((obs-point)/obs)^2))))
print(paste0("MAPE: ", MAPE<-mean(abs(obs-point)/obs)))

print(paste0("Mean conf. interval length: ", mean(up-lo)))


# PREDICCIONES CON MOD 2
(mod2_2 <- update(mod2, x = data_2))

#prediction 12 next observations
pred2 <- predict(mod2_2, n.ahead=12)
point2 <- pred2$pred

up2 <- point2 + (pred2$se * 1.96)
lo2 <- point2 - (pred2$se * 1.96)

ts.plot(peajes_log,lo2,up2,point2,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2013,2020), type="o",
        main="") 
abline(v=2013:2019,lty=3,col=4)

#RMSPE and MAPE
obs <- window(peajes_log, start=c(2018, 12), end=c(2019, 11))

print(paste0("RMSPE: ", RMSPE<-sqrt(mean(((obs-point2)/obs)^2))))
print(paste0("MAPE: ", MAPE<-mean(abs(obs-point2)/obs)))

print(paste0("Mean conf. interval length: ", mean(up2-lo2)))


# Escogemos modelo 1

# Tratamiento outliers

source("atipics2.R")

(atip <- outdetec(mod1,dif=c(0,12), crit=2.8, LS=T))

lnserie.lin=lineal(peajes_log, atip$atip)

ts.plot(peajes,exp(lnserie.lin), lty=c(1,2),col=c(1,3),xlim=c(1990,2019), type="o",
        main="")
legend("bottomright", legend = c("original", "linearizada"), fill = c(1,4))
abline(v=1990:2019,lty=3,col=4)


ts.plot(peajes- exp(lnserie.lin), lty=c(1,2),col=c(1,4),xlim=c(1990,2019), type="l",
        main="")
legend("bottomright", legend = c("original", "linearizada"), fill = c(1,4))
abline(v=1990:2019,lty=3,col=4)



monthplot(lnserie.lin)


#diferenciacion regular orden 12
plot(lnserie.lin_12 <- diff(lnserie.lin, 12), ylab = "lnserie.lin_12")
abline(h = 0, lty = 2, col = 4)
(var_12 <- var(lnserie.lin_12))
monthplot(lnserie.lin_12)
abline(v=1990:2020,lty=2,col=4)

# dif reg 1
plot(lnserie.lin_12_1 <- diff(lnserie.lin_12, 1))
abline(h = 0, lty = 2, col = 4)
(var_12_1 <- var(lnserie.lin_12_1))

#dif red 2
plot(lnserie.lin_12_2 <- diff(lnserie.lin_12_1, 1))
abline(h = 0, lty = 2, col = 4)
(var_12_2 <- var(lnserie.lin_12_2))


#serie definitiva
serie <- lnserie.lin_12_1


#ACF & PACF of residuals
s=frequency(serie)
par(mfrow=c(1,2))
acf(serie,ylim=c(-1,1),lag.max=100,col=c(2,rep(1,s-1)),lwd=1, main = "")
pacf(serie,ylim=c(-1,1),lag.max=100,col=c(rep(1,s-1),2),lwd=1, main = "")
par(mfrow=c(1,1))

# miramos intercept
(arima(serie, order=c(3,0,0), seasonal=list(order=c(1,0,0), period=12)))

(modl1 <- arima(lnserie.lin, order=c(3,1,0), seasonal=list(order=c(2,1,1), period=12)))
(modl2 <- arima(lnserie.lin, order=c(0,1,2), seasonal=list(order=c(2,1,1), period=12)))

validation(modl1, lnserie.lin)
validation(modl2, lnserie.lin)


AIC(modl1)
BIC(modl1)
AIC(modl2)
BIC(modl2)

# PREDICCION modl1
ultim <- c(2018,11)
data_2 <-window(lnserie.lin, end=ultim) 

(modl1_2 <- update(modl1, x = data_2))

#prediction 12 next observations
pred <- predict(modl1_2, n.ahead=12)
point <- pred$pred

up <- point + (pred$se * 1.96)
lo <- point - (pred$se * 1.96)

ts.plot(lnserie.lin,lo,up,point,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2013,2020), type="o",
        main="") 
abline(v=2013:2020,lty=3,col=4)

#RMSPE and MAPE
obs <- window(lnserie.lin, start=c(2018, 12), end=c(2019, 11))

print(paste0("RMSPE: ", RMSPE<-sqrt(mean(((obs-point)/obs)^2))))
print(paste0("MAPE: ", MAPE<-mean(abs(obs-point)/obs)))

print(paste0("Mean conf. interval length: ", mean(up-lo)))


# PREDICCION modl2
ultim <- c(2018,11)
data_2 <-window(lnserie.lin, end=ultim) 

(modl2_2 <- update(modl2, x = data_2))

#prediction 12 next observations
pred <- predict(modl2_2, n.ahead=12)
point <- pred$pred

up <- point + (pred$se * 1.96)
lo <- point - (pred$se * 1.96)

ts.plot(lnserie.lin,lo,up,point,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2013,2020), type="o",
        main="") 
abline(v=2013:2020,lty=3,col=4)

#RMSPE and MAPE
obs <- window(lnserie.lin, start=c(2018, 12), end=c(2019, 11))

print(paste0("RMSPE: ", RMSPE<-sqrt(mean(((obs-point)/obs)^2))))
print(paste0("MAPE: ", MAPE<-mean(abs(obs-point)/obs)))

print(paste0("Mean conf. interval length: ", mean(up-lo)))


#####################################################################
# Predicciones para la serie linealizada

ultim <- c(2019,11)
data_2 <-window(lnserie.lin, end=ultim) 

(modl2_1 <- update(modl1, x = data_2))

#prediction 12 next observations
pred <- predict(modl2_1, n.ahead=12)

## Efecto de los cambios de nivel detectados
wLS=sum(atip$atip[atip$atip$type_detected=="LS" & atip$atip$Obs<=length(serie)-12,3])

point <- exp(pred$pred + wLS)


up <- point*exp((pred$se * 1.96))
lo <- point/exp((pred$se * 1.96))

ts.plot(peajes,lo,up,point,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2013,2021), type="o",
        main="") 
abline(v=2013:2021,lty=3,col=4)



