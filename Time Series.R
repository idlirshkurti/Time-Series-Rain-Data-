packages <- c('ggplot2', 'nlme', 'mgcv', 'TSA', 'forecast', 'devtools', 'ggfortify',
 'reshape2', 'stats', 'dlm', 'grid', 'gridExtra','scales')
lapply(packages, require, character.only = TRUE)
install_github('sinhrks/ggfortify')


ow = read.table("http://www.stats.ox.ac.uk/~reinert/dtc/ow.txt",header=T)
rain = ts(ow$rain, start=ow$yyyy[1], frequency=12)
tma = ts(ow$tmax, start=ow$yyyy[1], frequency=12)
# ------------------ Exploratory Analysis ------------------------------ #

# -- Histograms -- #
hist1 <- qplot(ow$rain, geom="histogram", main = "Histogram for rainfall in Oxford",
xlab = "Rainfall", ylab="Frequency", fill=I("blue"), col=I("red"), alpha=I(.8))
hist1 <- hist1 + labs(title = "Histogram for rainfall in Oxford") +
theme(legend.position = "none", title = element_text(size=15))
authist1 <- autoplot(hist1)
#Same for ow$tmax

#Call the autoplotted histogram authist2
grid.arrange(authist1, authist2, ncol=2)

### Monthly histograms
par(mfrow=c(3,4))
for (i in 1:12) { hist(rain[seq(i,1853,12)],
main=month.name[i],xlab="Rainfall")}
par(mfrow=c(3,4))

for (i in 1:12) { hist(tma[seq(i,1853,12)], main=month.name[i], xlab="Maximum temperature")}

##Lagplots
gglagplot(rain, lags = 12) + labs(title = "Monthly lagplots for rainfall data") + 
theme(legend.position = "none", title = element_text(size=15))

gglagplot(tma, lags = 12) + labs(title = "Monthly lagplots for maximum temperature data") + 
theme(legend.position = "none", title = element_text(size=15))

####### Observed rainfall(zooming in seasonality detection) #########
auto1 <- autoplot(rain, ylab= "Rainfall", xlab="Year")
auto1 <- auto1 + labs(title = "Rainfall in Oxford from 1853 to 2014")
auto1 <- auto1 + theme(legend.position = "none", title = element_text(size=15))
auto1
rain1 <- ts(ow$rain,start=ow$yyyy[1],end=ow$yyyy[565], frequency= 12)

auto2 <- autoplot(rain1, ylab= "Rainfall", xlab="Year")
auto2 <- auto2 + labs(title = "1853 to 1900")
auto2 <- auto2 + theme(legend.position = "none", title = element_text(size=15))
auto2
rain2<-ts(ow$rain,start=ow$yyyy[1],end=ow$yyyy[240], frequency= 12)

auto3 <- autoplot(rain2, ylab= "Rainfall", xlab="Year")
auto3 <- auto3 + labs(title = "1853 to 1873")
auto3 <- auto3 + theme(legend.position = "none", title = element_text(size=15))
auto3
# same for ow$tmax
# name the plots aut1,aut2,aut3
grid.arrange(auto1,aut1,auto2,aut2, auto3,aut3, nrow = 3)

##### Timeplots for separate months
freq1 <- ggfreqplot(rain, ncol = 12) + scale_x_date(labels=date_format("%Y"),
breaks = as.Date(c("1900-1-1","1980-12-31"))) +
labs(title = "Rainfall separate for each month. 1-12 denoting months Jan-Dec") + theme(legend.position = "none",
title = element_text(size=15))

freq2 <- ggfreqplot(tma, ncol = 12) + scale_x_date(labels=date_format("%Y"),
breaks=as.Date(c("1900-1-1","1980-12-31"))) + labs(title =
"Maximum temperature separate for each month. 1-12 denoting months Jan-Dec")+ 
theme(legend.position = "none", title = element_text(size=15))
grid.arrange(freq1,freq2,nrow=2)

###### Seasonal decomposition and Autocorrelation function#######
acfr <- autoplot(acf(ow$rain,lag=50)) + labs(title = "ACF for the rain data") +
theme(legend.position = "none", title = element_text(size=15))

pacfr <- autoplot(pacf(ow$rain,lag=50)) + labs(title = "PACF for the rain data") +
labs(y="PACF") + theme(legend.position = "none", title = element_text(size=15))
acft <- autoplot(acf(ow$tma,lag=50)) + labs(title = "ACF for the maximum temperature data") + 
theme(legend.position = "none", title = element_text(size=15))

pacft <- autoplot(pacf(ow$tma,lag=50)) + labs(title = "PACF for the maximum temperature data") +
labs(y="PACF")+ theme(legend.position = "none", title = element_text(size=15))

grid.arrange(acfr,acft,pacfr,pacft, ncol=2)
seas <- autoplot(stl(rain,s.window='periodic'), ts.colour='deepskyblue4')
seas <- seas + labs(title = "Rainfall data")
seas <- seas + theme(legend.position = "none",
title = element_text(size=15))
seas
# same for tma
# name the seasonal decomposition plot seast
grid.arrange(seas,seast, ncol=2)

########## Raw spectogram ###########
raw.spec <- spec.pgram(ow$rain,plot=F)
spec.df <- data.frame(freq = raw.spec$freq, spec = raw.spec$spec)

per1 <- ggplot(data = subset(spec.df)) + geom_line(aes(x = freq, y = spec)) + 
labs(title = "Periodogram for the rain data",y="Periodogram",x="Frequency") + 
theme(legend.position = "none", title = element_text(size=15))

raw.spec1 <- spec.pgram(ow$tmax,plot=F)
spec.df1 <- data.frame(freq = raw.spec1$freq, spec = raw.spec1$spec)
per2 <- ggplot(data = subset(spec.df1)) + geom_line(aes(x = freq, y = spec)) + 
labs(title = "Periodogram for the maximum temperature data", y="Periodogram",x="Frequency") + 
theme(legend.position = "none", title = element_text(size=15))
grid.arrange(per1,per2,ncol=1)

######### Smoothing ##########
k = kernel("daniell", c(9, 9, 9))
smooth.spec <- spec.pgram(ow$rain, kernel = k, taper = 0,log="no")
spec.df <- data.frame(freq = smooth.spec$freq, `c(9,9,9)` = smooth.spec$spec)
names(spec.df) <- c("freq", "c(9,9,9)")
k <- kernel("daniell", c(9, 9))
spec.df[, "c(9,9)"] <- spec.pgram(ow$rain, kernel = k, taper = 0, plot = FALSE)$spec
k <- kernel("daniell", c(9))
spec.df[, "c(9)"] <- spec.pgram(ow$rain, kernel = k, taper = 0, plot = FALSE)$spec
spec.df <- melt(spec.df, id.vars = "freq", value.name = "spec", variable.name = "kernel")
smgg1 <- ggplot(data = subset(spec.df)) + geom_path(aes(x = freq, y = spec, color = kernel)) + 
labs(title = "Smoothed periodogram for the rainfall data", y="Periodogram",x="Frequency") + 
theme(title = element_text(size=15))
smgg1

# repeat code for maximum temperature and call the plot smgg2

k = kernel("daniell", c(9, 9, 9))
smooth.spec <- spec.pgram(ow$tmax, kernel = k, taper = 0,log="no")
spec.df1 <- data.frame(freq = smooth.spec1$freq, 'c(9,9,9)' = smooth.spec1$spec)
names(spec.df1) <- c("freq", "c(9,9,9)")
k <- kernel("daniell", c(9, 9))
spec.df1[, "c(9,9)"] <- spec.pgram(ow$tmax, kernel = k, taper = 0, plot = FALSE)$spec
k <- kernel("daniell", c(9))
spec.df1[, "c(9)"] <- spec.pgram(ow$tmax, kernel = k, taper = 0, plot = FALSE)$spec
spec.df1 <- melt(spec.df1, id.vars = "freq", value.name = "spec", variable.name = "kernel")
smgg2 <- ggplot(data = subset(spec.df1)) + geom_path(aes(x = freq, y = spec, color = kernel)) + 
labs(title = "Smoothed periodogram for the rainfall data", y="Periodogram",x="Frequency") + 
theme(title = element_text(size=15))
smgg2
grid.arrange(smgg1,smgg2,ncol=1)

########################## AR(2) model #####################################
a.spec <- spec.ar(ow$rain,log="no", plot=F)
#Assumes data comes from an AR model
aspe1 <- autoplot(a.spec)+ labs(title = "Rain AR(2) spectrum", y="Spectrum",x="Frequency") + 
theme(title = element_text(size=15))
a.spec1 <- spec.ar(ow$tmax,log="no"), plot=F)#Assumes data comes from an AR model
aspe2 <- autoplot(a.spec1)+ labs(title = "tmax AR(25) spectrum", y="Spectrum",x="Frequency") + 
theme(title = element_text(size=15))
grid.arrange(aspe1,aspe2,ncol=1)

###Fit rain into AR(2)
a1 <- Arima(rain,order=c(2,0,0))
# fitted vs observed plots
plot(a1$x,col="lightblue",ylab="rainfall",xlab="Year", main = 
"Time plot for AR(2) rainfall fitted values(red) vs observed rainfall(blue)") + 
lines(fitted(a1),col="red")

# Residuals plot
autoplot(residuals(a1), colour="darkred")
#Residuals test
ggtsdiag(a1)
#auto.arima model fitting
a2t <- auto.arima(tma)
a2t$coef
# fitted vs observed plots
plot(a2t$x, col="lightblue", main="Observed temperature values (blue) vs fitted values (red)", 
ylab = "Maximum temperature")
lines(fitted(a2t),col="red")

# Residuals plot
autoplot(residuals(a2t),colour="darkred")
#Residuals test
ggtsdiag(a2t)
Box.test(residuals(a2t),lag=10,fitdf=2,type="Ljung")

########### Seasonal Means Model ###################
#Rain Data
mon = month.abb[((1:length(rain))-1) %% 12+1]
model = lm(rain~mon)
autoplot(model)
#tmax Data
model1 = lm(tma~mon+ow$yyyy)
autoplot(model1)
############ Question 3 - Square rooted values #############
rainsqrt <- sqrt(rain)
acf(rainsqrt,lag=50)
pacf(rainsqrt,lag=50)
model = lm(rainsqrt~mon)
autoplot(model)
########### Seasonal State Space Model(rain) ##############
build = function(parm) {
r = parm[1]; q = parm[2]
mod_s = dlmModSeas(12, dV=r^2, dW=c(q^2, rep(0,10)))
mod_s
}
# Find the MLE
mle = dlmMLE(rain, c(1,1), build, lower=rep(0.001,2), hessian=T)
mod = build(mle$par)
filt = dlmFilter(rain, mod)
smth = dlmSmooth(rain, mod)
p <- autoplot(filt)+labs(x= "Year", y = "Rainfall") + 
ggtitle("Kalman filtered model (red) and rainfall observations (black)") +
theme(plot.title = element_text(size= 20, lineheight=8.8, face = "bold"))
m_t = smth$s[,1]
d_t = 1.96*sqrt(simplify2array(dlmSvd2var(smth$U.S, smth$D.S))[1,1,])
m_s = smth$s[,3]

d_s = 1.96*sqrt(simplify2array(dlmSvd2var(smth$U.S, smth$D.S))[3,3,])
autoplot(smth$s)
# Residual diagnostics
tsdiag(filt)
########### Seasonal State Space Model(tmax) ##############
build = function(parm) {
r = parm[1]; q = parm[2]
mod_s = dlmModSeas(12, dV=r^2, dW=c(q^2, rep(0,10)))
mod_s
}
# Find the MLE
mle = dlmMLE(tma, c(1,1), build, lower=rep(0.001,2), hessian=T)
mod = build(mle$par)
filt = dlmFilter(tma, mod)
smth = dlmSmooth(tma, mod)
p1 <- autoplot(filt)+labs(x= "Year", y = "Maximum Temperature") + 
ggtitle("Kalman filtered model (red) and maximum temperature observations (black)") +
theme(plot.title = element_text(size= 20, lineheight=8.8, face="bold"))
m_t = smth$s[,1]
d_t = 1.96*sqrt(simplify2array(dlmSvd2var(smth$U.S, smth$D.S))[1,1,])
m_s = smth$s[,3]
d_s = 1.96*sqrt(simplify2array(dlmSvd2var(smth$U.S, smth$D.S))[3,3,])
grid.arrange(p1,p,ncol=1)

######## Seasonal State Space Model With Linear Trend (rain) #################
# Build a state-space model with polynomial trend
# plus a seasonal component
build = function(parm) {
a = parm[1]; r = parm[2]; q = parm[3]; b = parm[4]
mod = dlmModPoly(2, dV=r^2)
W(mod) = matrix(c(a^2,0,0,b^2), nrow=2)
mod_s = dlmModSeas(12, dV=r^2/2, dW=c(q^2, rep(0,10)))
mod + mod_s
}

# Find the MLE
mle = dlmMLE(rain, c(1,1,1,1), build, lower=rep(0.001,4), hessian=T)
mod = build(mle$par)
filt = dlmFilter(rain, mod)
smth = dlmSmooth(rain, mod)
p<-autoplot(filt)+labs(x= "Year", y = "Rainfall") + 
ggtitle("Kalman filtered model (red) and rainfall observations (black)") +
theme(plot.title = element_text(size= 20, lineheight=8.8, face="bold"))
autoplot(p)

m_t = smth$s[,1]
d_t = 1.96*sqrt(simplify2array(dlmSvd2var(smth$U.S, smth$D.S))[1,1,])
m_s = smth$s[,3]
d_s = 1.96*sqrt(simplify2array(dlmSvd2var(smth$U.S, smth$D.S))[3,3,])

# Forecasting with the model
fc = dlmForecast(filt, nAhead=24)
m_t_fc = fc$a[,1]
d_t_fc = 1.96*sqrt(simplify2array(fc$R)[1,1,])
m_s_fc = fc$a[,3]
d_s_fc = 1.96*sqrt(simplify2array(fc$R)[3,3,])
x_fc = fc$f; d_fc = 1.96*sqrt(unlist(fc$Q))

# Plot the data + prediction, and states + predictions
par(mfrow=c(3,1))
ts.plot(rain, x_fc, x_fc+d_fc, x_fc-d_fc, col=c(1,2,2,2), lty=c(1,2,2,2))
ts.plot(m_t, m_t_fc, m_t_fc+d_t_fc, m_t_fc-d_t_fc, col=c(1,2,2,2),lty=c(1,2,2,2))
ts.plot(m_s, m_s_fc, m_s_fc+d_s_fc, m_s_fc-d_s_fc, col=c(1,2,2,2),lty=c(1,2,2,2))
autoplot(smth$s)
# Residual diagnostics
tsdiag(filt)