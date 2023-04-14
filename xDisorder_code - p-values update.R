library(survival)
library(ggplot2)

dx1.time <- dx2.time <- comorbid.time <- seq(from=78, to=1144, by=52)

# ---- Estimating ----

fit.dx1 = coxph(Surv(time = dx1_follow_duration, event = dx1) ~ exposure + confounders,
                weights = weight, data = df)
fit.dx2 = coxph(Surv(time = dx2_follow_duration, event = dx2) ~ exposure + confounders,
                weights = weight, data = df)
fit.comorbid = coxph(Surv(time = comorbid_follow_duration, event = comorbid) ~ exposure + confounders,
                     weights = weight, data = df)

# Create high and low exposure dataframe:

predict.high <- predict.low <- df
predict.high$exposure = 1
predict.low$exposure = 0

# In each group the base cumulative hazard is determined and saved as a function of the time:

tmp.dx1 = basehaz(fit.dx1)
H0 = stepfun(tmp.dx1$time, c(0, tmp.dx1$hazard))
tmp.dx2 = basehaz(fit.dx2)
H1 = stepfun(tmp.dx2$time, c(0, tmp.dx2$hazard))
tmp.comorbid = basehaz(fit.comorbid)
H2 = stepfun(tmp.comorbid$time, c(0, tmp.comorbid$hazard))

# This is saved as in a matrix to prepare for calculating the survival function at each time point

H0.Mat = matrix(-H0(dx1.time))
H1.Mat = matrix(-H1(dx2.time))
H2.Mat = matrix(-H2(comorbid.time))

# In each group the hazard is calculated for the high and low exposure data
# Then the probability of having been diagnosed prior to each time point is calculated for both exposures

tmp.dx1.high = predict(fit.dx1, newdata = predict.high, type='risk')
tmp.dx1.low = predict(fit.dx1, newdata = predict.low, type='risk')
org.dx1.high = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx1.high)), na.rm = T)
org.dx1.low = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx1.low)), na.rm = T)
tmp.dx2.high = predict(fit.dx2, newdata = predict.high, type='risk')
tmp.dx2.low = predict(fit.dx2, newdata = predict.low, type='risk')
org.dx2.high = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx2.high)), na.rm = T)
org.dx2.low = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx2.low)), na.rm = T)
tmp.comorbid.high = predict(fit.comorbid, newdata = predict.high, type='risk')
tmp.comorbid.low = predict(fit.comorbid, newdata = predict.low, type='risk')
org.comorbid.high = 1 - rowMeans(exp(H0.Mat %*% t(tmp.comorbid.high)), na.rm = T)
org.comorbid.low = 1 - rowMeans(exp(H0.Mat %*% t(tmp.comorbid.low)), na.rm = T)

# Calculating the Risk Ratios:

org.lnrr_dx1 = org.dx1.high/org.dx1.low
org.dx1.high.comorbid = org.dx1.high - org.comorbid.high
org.dx1.low.comorbid = org.dx1.low - org.comorbid.low
org.lnrr_dx1_comorbid = org.dx1.high.comorbid/org.dx1.low.comorbid

org.lnrr_dx2 = org.dx2.high/org.dx2.low
org.dx2.high.comorbid = org.dx2.high - org.comorbid.high
org.dx2.low.comorbid = org.dx2.low - org.comorbid.low

org.lnrr_dx2_comorbid = org.dx2.high.comorbid/org.dx2.low.comorbid
org.lnrr_comorbid = org.comorbid.high/org.comorbid.low

# Making Plots

tmp.plot = as.data.frame(cbind(dx1.time, org.lnrr_dx1_comorbid, org.lnrr_dx2_comorbid, org.lnrr_comorbid))
colnames(tmp.plot) = c('time', 'dx1_Only', 'dx2_Only', 'Comorbid')

org.tmp <- ggplot() +
  geom_point(data=tmp.plot, aes(time, dx1_Only), color='blue') +
  geom_line(data=tmp.plot, aes(time, dx1_Only), color='blue') +
  geom_point(data=tmp.plot, aes(time, dx2_Only), color='green') +
  geom_line(data=tmp.plot, aes(time, dx2_Only), color='green') +
  geom_point(data=tmp.plot, aes(time, Comorbid), color='purple') +
  geom_line(data=tmp.plot, aes(time, Comorbid), color='purple') +
  geom_hline(yintercept = 1, color = 'red', size = 1.5) +
  labs(title = 'Risk Ratio for Original HR') +
  xlab('Follow Duration in Days') +
  ylab('Risk Ratio')


# ---- Bootstrap ----

B = 100

# matrices to contain results:

CI_matrix_dx1_high <- CI_matrix_dx1_low <- matrix(data = NA, nrow = B, ncol = length(dx1.time))
CI_matrix_dx2_high <- CI_matrix_dx2_low <- matrix(data = NA, nrow = B, ncol = length(dx2.time))
CI_matrix_comorbid_high <- CI_matrix_comorbid_low <- matrix(data = NA, nrow = B, ncol = length(comorbid.time))

for (b in 1:B){
  set.seed(1111+b)

  # Sampling the weights:

  df.temp = as.data.frame(df)
  df.temp$weightexp <- df.temp$weightbs <- NA

  index_control = which(df.temp$dx1 == 0 & df.temp$dx2 == 0)
  df.temp$weightexp[index_control] = rexp(length(index_control))
  df.temp$weightexp[index_control] = df.temp$weightexp[index_control]/mean(df.temp$weightexp[index_control])

  index_dx1 = which(df.temp$dx1 == 1 & df.temp$dx2 == 0)
  df.temp$weightexp[index_dx1] = rexp(length(index_dx1))
  df.temp$weightexp[index_dx1] = df.temp$weightexp[index_dx1]/mean(df.temp$weightexp[index_dx1])

  index_dx2 = which(df.temp$dx1 == 0 & df.temp$dx2 == 1)
  df.temp$weightexp[index_dx2] = rexp(length(index_dx2))
  df.temp$weightexp[index_dx2] = df.temp$weightexp[index_dx2]/mean(df.temp$weightexp[index_dx2])

  index_comorbid = which(df.temp$dx1 == 1 & df.temp$dx2 == 1)
  df.temp$weightexp[index_comorbid] = rexp(length(index_comorbid))
  df.temp$weightexp[index_comorbid] = df.temp$weightexp[index_comorbid]/mean(df.temp$weightexp[index_comorbid])

  df.temp$weightbs = df.temp$weight * df.temp$weightexp
  df.temp$weightbs = round(df.temp$weightbs*(10^9), 0) # rounding the weights like this makes little to no difference and makes the code much faster to run due to R not having to handle long/small floats

  # the remainder of the code in the for-loop is identical to other case:

  # Fitting models

  fit.dx1 = coxph(Surv(time = dx1_follow_duration, event = dx1) ~ exposure + confounders,
                  weights = weightbs, data = df.temp)
  fit.dx2 = coxph(Surv(time = dx2_follow_duration, event = dx2) ~ exposure + confounders,
                  weights = weightbs, data = df.temp)
  fit.comorbid = coxph(Surv(time = comorbid_follow_duration, event = comorbid) ~ exposure + confounders,
                       weights = weightbs, data = df.temp)

  # Create high and low exposure dataframe:

  predict.high <- predict.low <- df
  predict.high$exposure = 1
  predict.low$exposure = 0

  # In each group the base cumulative hazard is determined and saved as a function of the time:

  tmp.dx1 = basehaz(fit.dx1)
  H0 = stepfun(tmp.dx1$time, c(0, tmp.dx1$hazard))

  tmp.dx2 = basehaz(fit.dx2)
  H1 = stepfun(tmp.dx2$time, c(0, tmp.dx2$hazard))

  tmp.comorbid = basehaz(fit.comorbid)
  H2 = stepfun(tmp.comorbid$time, c(0, tmp.comorbid$hazard))

  # This is saved as in a matrix to prepare for calculating the survival function at each time point

  H0.Mat = matrix(-H0(dx1.time))
  H1.Mat = matrix(-H1(dx2.time))
  H2.Mat = matrix(-H2(comorbid.time))

  # In each group the hazard is calculated for the high and low exposure data
  # Then the probability of having been diagnosed prior to each time point is calculated for both exposures

  tmp.dx1.high = predict(fit.dx1, newdata = predict.high, type='risk')
  tmp.dx1.low = predict(fit.dx1, newdata = predict.low, type='risk')
  CI_matrix_dx1_high[b,] = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx1.high)), na.rm = T)
  CI_matrix_dx1_low[b,] = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx1.low)), na.rm = T)
  tmp.dx2.high = predict(fit.dx2, newdata = predict.high, type='risk')
  tmp.dx2.low = predict(fit.dx2, newdata = predict.low, type='risk')
  CI_matrix_dx2_high[b,] = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx2.high)), na.rm = T)
  CI_matrix_dx2_low[b,] = 1 - rowMeans(exp(H0.Mat %*% t(tmp.dx2.low)), na.rm = T)
  tmp.comorbid.high = predict(fit.comorbid, newdata = predict.high, type='risk')
  tmp.comorbid.low = predict(fit.comorbid, newdata = predict.low, type='risk')
  CI_matrix_comorbid_high[b,] = 1 - rowMeans(exp(H0.Mat %*% t(tmp.comorbid.high)), na.rm = T)
  CI_matrix_comorbid_low[b,] = 1 - rowMeans(exp(H0.Mat %*% t(tmp.comorbid.low)), na.rm = T)

}

# Calculating the Risk Ratios:

lnrr_dx1_comorbid = (CI_matrix_dx1_high - CI_matrix_comorbid_high)/(CI_matrix_dx1_low - CI_matrix_comorbid_low)
sd_dx1_comorbid = apply(log(lnrr_dx1_comorbid), 2, sd, na.rm=T)
CI_U_dx1_comorbid = log(lnrr_dx1_comorbid) + qnorm(.975) * sd_dx1_comorbid
CI_L_dx1_comorbid = log(lnrr_dx1_comorbid) + qnorm(.025) * sd_dx1_comorbid
plot_dx1_comorbid = as.data.frame(cbind(dx1.time, lnrr_dx1_comorbid, exp(CI_U_dx1_comorbid), exp(CI_L_dx1_comorbid)))
colnames(plot_dx1_comorbid) = c('dx1.time', 'RR', 'RR_Upper', 'RR_Lower')


lnrr_dx2_comorbid = (CI_matrix_dx2_high - CI_matrix_comorbid_high)/(CI_matrix_dx2_low - CI_matrix_comorbid_low)
sd_dx2_comorbid = apply(log(lnrr_dx2_comorbid), 2, sd, na.rm=T)
CI_U_dx2_comorbid = log(lnrr_dx2_comorbid) + qnorm(.975) * sd_dx2_comorbid
CI_L_dx2_comorbid = log(lnrr_dx2_comorbid) + qnorm(.025) * sd_dx2_comorbid
plot_dx2_comorbid = as.data.frame(cbind(dx2.time, lnrr_dx2_comorbid, exp(CI_U_dx2_comorbid), exp(CI_L_dx2_comorbid)))
colnames(plot_dx2_comorbid) = c('dx2.time', 'RR', 'RR_Upper', 'RR_Lower')

lnrr_comorbid = (CI_matrix_comorbid_high)/(CI_matrix_comorbid_low)
sd_comorbid = apply(log(lnrr_comorbid), 2, sd, na.rm=T)
CI_U_comorbid = log(lnrr_comorbid) + qnorm(.975) * sd_comorbid
CI_L_comorbid = log(lnrr_comorbid) + qnorm(.025) * sd_comorbid
plot_comorbid = as.data.frame(cbind(comorbid.time, lnrr_comorbid, exp(CI_U_comorbid), exp(CI_L_comorbid)))
colnames(plot_dx1_comorbid) = c('comorbid.time', 'RR', 'RR_Upper', 'RR_Lower')

# Making Plots

plot_all = as.data.frame(cbind(plot_dx1_comorbid, plot_dx2_comorbid, plot_comorbid))
plot_all = as.data.frame(plot_all[,-c(5,9)]) # remove duplicates of time
colnames(plot_all) = c('time',
                       'dx1.comorbid.RR', 'dx1.comorbid.RR.U', 'dx1.comorbid.RR.L',
                       'dx2.comorbid.RR', 'dx2.comorbid.RR.U', 'dx2.comorbid.RR.L',
                       'comorbid.RR', 'comorbid.RR.U', 'comorbid.RR.L')

p.all <- ggplot() +
  geom_line(data=plot_all, aes(x = time, y=dx1.comorbid.RR, color='dx1')) +
  geom_line(data=plot_all, aes(x = time, y=dx2.comorbid.RR, color='dx2')) +
  geom_line(data=plot_all, aes(x = time, y=comorbid.RR, color='comorbid')) +
  scale_colour_manual(name = 'Condition', values=c('dx1=blue', 'dx2'='green', 'comorbid=purple'),
                                                   labels = c('dx1', 'dx2', 'dx1+dx2')) +
  ylab('Risk Ratio and 95% Confidence Interval') +
  coord_trans(y='log10') +
  scale_x_continuous(name = 'Follow-up Time (Years)') +
  theme(axis.text.x = element_text(size=12, face='bold', angle=75)) +
  ggtitle('Risk Ratio and 95% Confidence Interval for Exposure') +
  geom_hline(yintercept = 1, color = 'red', size=1.5)

p.all <- p.all + geom_point(data=plot_all, aes(x=time, y=dx1.comorbid.RR), color='blue') +
  geom_point(data=plot_all, aes(x=time, y=dx2.comorbid.RR), color='green') +
  geom_point(data=plot_all, aes(x=time, y=comorbid.RR), color='purple') +
  geom_ribbon(data=plot_all, aes(x=time, ymin=dx1.comorbid.RR.L, ymax=dx1.comorbid.RR.U), alpha=0.3, fill='blue') +
  geom_ribbon(data=plot_all, aes(x=time, ymin=dx2.comorbid.RR.L, ymax=dx2.comorbid.RR.U), alpha=0.3, fill='green') +
  geom_ribbon(data=plot_all, aes(x=time, ymin=comorbid.RR.L, ymax=comorbid.RR.U), alpha=0.3, fill='purple')

# ---- p-values: average and whole curve ----

avgrr_dx1 = rowMeans(log(lnrr_dx1_comorbid))
avgrr_dx2 = rowMeans(log(lnrr_dx2_comorbid))
avgrr_comorbid = rowMeans(log(CI_matrix_comorbid_high/CI_matrix_comorbid_low))

org.avgrr_dx1 = mean(log(org.lnrr_dx1_comorbid))
org.avgrr_dx2 = mean(log(org.lnrr_dx2_comorbid))
org.avgrr_comorbid = mean(log(org.lnrr_comorbid))

p_vals = data.frame(matrix(NA, 2, 4))
colnames(p_vals) = c('all', 'dx1 vs dx2', 'dx1 vs comorbid', 'dx2 vs comorbid')
rownames(p_vals) = c('avg eq', 'whole curve eq')

# average equal:

Sigmahat = var(cbind(avgrr_dx1, avgrr_dx2, avgrr_comorbid))

Lmat = rbind(c(1,-1,0),
             c(0,1,-1))
est = Lmat %*% c(org.avgrr_dx1, org.avgrr_dx2, org.avgrr_comorbid)
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% solve(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['avg eq', 'all'] = pval

Lmat = rbind(c(1,-1,0))
est = Lmat %*% c(org.avgrr_dx1, org.avgrr_dx2, org.avgrr_comorbid)
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% solve(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['avg eq', 'dx1 vs dx2'] = pval

Lmat = rbind(c(1,0,-1))
est = Lmat %*% c(org.avgrr_dx1, org.avgrr_dx2, org.avgrr_comorbid)
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% solve(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['avg eq', 'dx1 vs comorbid'] = pval

Lmat = rbind(c(0,1,-1))
est = Lmat %*% c(org.avgrr_dx1, org.avgrr_dx2, org.avgrr_comorbid)
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% solve(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['avg eq', 'dx2 vs comorbid'] = pval

# whole curve equal

Sigmahat = var(cbind(log(lnrr_dx1_comorbid)[,3:20],
                     log(lnrr_dx2_comorbid)[,3:20],
                     log(CI_matrix_comorbid_high/CI_matrix_comorbid_low)[,3:20]))

Lmat = cbind(diag(18), -diag(18), diag(rep(0, 18)))
est = Lmat %*% c(log(org.lnrr_dx1_comorbid)[3:20], log(org.lnrr_dx2_comorbid)[3:20], log(org.lnrr_comorbid)[3:20])
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% ginv(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['whole cure eq', 'dx1 vs dx2'] = pval

Lmat = cbind(diag(18), diag(rep(0, 18)),  -diag(18))
est = Lmat %*% c(log(org.lnrr_dx1_comorbid)[3:20], log(org.lnrr_dx2_comorbid)[3:20], log(org.lnrr_comorbid)[3:20])
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% ginv(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['whole cure eq', 'dx1 vs comorbid'] = pval

Lmat = cbind(diag(rep(0, 18)), diag(18),  -diag(18))
est = Lmat %*% c(log(org.lnrr_dx1_comorbid)[3:20], log(org.lnrr_dx2_comorbid)[3:20], log(org.lnrr_comorbid)[3:20])
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% ginv(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['whole cure eq', 'dx2 vs comorbid'] = pval

Lmat = rbind(cbind(diag(rep(0, 18)), diag(18),  -diag(18)),
             cbind( diag(18), diag(rep(0, 18)), -diag(18)))
est = Lmat %*% c(log(org.lnrr_dx1_comorbid)[3:20], log(org.lnrr_dx2_comorbid)[3:20], log(org.lnrr_comorbid)[3:20])
estV = est %*% Sigmahat %*% est
chisq = t(est) %*% ginv(estV) %*% est
pval = 1 - pchisq(chisq, df = Matrix::rankMatrix(estV))
p_vals_urban['whole cure eq', 'all'] = pval
