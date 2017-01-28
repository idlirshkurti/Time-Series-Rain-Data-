# --------------- Seasonal State Space Model(rain) ---------------- #

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
