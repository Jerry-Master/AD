
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
  
  #Normal plot of residuals
  qqnorm(resid)
  qqline(resid,col=2,lwd=2)
  
  ##Histogram of residuals with normal curve
  hist(resid,breaks=20,freq=FALSE)
  curve(dnorm(x,mean=mean(resid),sd=sd(resid)),col=2,add=T)
  
  
  #ACF & PACF of residuals
  par(mfrow=c(1,2))
  acf(resid,ylim=c(-1,1),lag.max=60,col=c(2,rep(1,s-1)),lwd=1)
  pacf(resid,ylim=c(-1,1),lag.max=60,col=c(rep(1,s-1),2),lwd=1)
  par(mfrow=c(1,1))
  
  #ACF & PACF of square residuals 
  par(mfrow=c(1,2))
  acf(resid^2,ylim=c(-1,1),lag.max=60,col=c(2,rep(1,s-1)),lwd=1)
  pacf(resid^2,ylim=c(-1,1),lag.max=60,col=c(rep(1,s-1),2),lwd=1)
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


setwd("~/Desktop")

peajes <- ts(read.table("./peajes.dat"), start = 1990, frequency = 12)

plot(peajes, main="Intensidad Media Diaria en peajes en España", ylab="miles de vehiculos")
abline(v=1990:2020,lty=3,col=4)
var(peajes)


monthplot(peajes)

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

peajes_log <- log(peajes)

#diferenciacion regular orden 12
plot(peajes_12 <- diff(peajes_log, 12))
abline(v=1990:2020,lty=2,col=4)
(var_12 <- var(peajes_12))
monthplot(peajes_12)

# dif reg 1
plot(peajes_12_1 <- diff(peajes_12, 1))
abline(h = 0, lty = 2, col = 4)
(var_12_1 <- var(peajes_12_1))

#dif red 2
plot(peajes_12_2 <- diff(peajes_12_1, 1))
abline(h = 0, lty = 2, col = 4)
(var_12_2 <- var(peajes_12_2))

#serie definitiva
serie <- peajes_12_1


#ACF & PACF of residuals
s=frequency(serie)
par(mfrow=c(1,2))
acf(serie,ylim=c(-1,1),lag.max=1000,col=c(2,rep(1,s-1)),lwd=1, main = "")
pacf(serie,ylim=c(-1,1),lag.max=1000,col=c(rep(1,s-1),2),lwd=1, main = "")
par(mfrow=c(1,1))


mod1 <- arima(peajes_log, order=c(0,1,2), seasonal=list(order=c(2,1,0), period=12))
mod2 <- arima(peajes_log, order=c(4,1,0), seasonal=list(order=c(2,1,0), period=12))
mod3 <- arima(peajes_log, order=c(1,1,1), seasonal=list(order=c(2,1,0), period=12))
mod4 <- arima(peajes_log, order=c(23,1,0), seasonal=list(order=c(0,1,0), period=12))
mod5 <- arima(peajes_log, order=c(23,1,0), fixed = c(rep(0,10), NA, rep(0,11), NA), seasonal=list(order=c(0,1,0), period=12))




res1 <- resid(mod1)
res2 <- resid(mod2)

validation(mod1, peajes_log)
validation(mod2, peajes_log)
validation(mod3, peajes_log)
validation(mod4, peajes_log)
validation(mod5, peajes_log)

AIC(mod1)
BIC(mod1)
AIC(mod2)
BIC(mod2)


ultim <- c(2018,11)
data_2 <-window(peajes_log, end=ultim)

mod2 <- update(mod1, x = data_2)

#prediction 12 next observations
pred <- predict(mod2, n.ahead=12)
point <- exp(pred$pred)

up <- point*exp(pred$se * 1.96)
lo <- point/exp(pred$se * 1.96)

ts.plot(peajes,lo,up,point,lty=c(1,2,2,1),col=c(1,4,4,2),xlim=c(2013,2020), type="o",
        main="Predicciones de 2018 según el modelo 1") 
abline(v=2013:2019,lty=3,col=4)

#RMSPE and MAPE
obs <- window(data, start=c(2018, 1), end=c(2018, 12))

print(paste0("RMSPE: ", RMSPE<-sqrt(mean(((obs-point)/obs)^2))))
print(paste0("MAPE: ", MAPE<-mean(abs(obs-point)/obs)))

print(paste0("Mean conf. interval length: ", mean(up-lo)))


