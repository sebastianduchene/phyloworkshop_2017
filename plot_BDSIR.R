library(NELSI)

data <- read.table('NorthAm.BDSIR.log', head = T)
most_recent_sample <- 2009
burnin <- 3000
data <- data[burnin:nrow(data), ]

susceptibles_cols <-   grep('dSEs.+', colnames(data), value = F)
recovereds_cols <- grep('dREs.+', colnames(data), value = F)
intervals <- length(susceptibles_cols)

I <- vector(length = intervals)
S <- vector(length = intervals)
R <- vector(length = intervals)
t <- vector(length = intervals)
incidence <- vector(length = intervals)
I[1] <- 1
R[1] <- 0
S[1] <- median(data$S0Es)
t[1] <- most_recent_sample - median(data$originEs)
incidence[1] <- 1
step <- median(data$originEs) / intervals
for(i in 2:100){
    S[i] <- S[i-1] - median(data[, susceptibles_cols[i]])
    R[i] <- R[i-1] + median(data[, recovereds_cols[i]])
    I[i] <- I[i-1] + median(data[, susceptibles_cols[i]]) - median(data[, recovereds_cols[i]])
    incidence[i] <- median(data[, susceptibles_cols[i]])
    t[i] <- t[i-1] + step
}

par(mfrow = c(2, 2))
plot(t, I, type = 'l', ylim = c(0, max(S)), col = 'orange', ylab = 'Number of individuals', main = 'SIR trajectories', xlab = 'Time')
lines(t, R, col = 'blue')
lines(t, S, col = 'green')
legend(2008, 100000, legend = c('S', 'I', 'R'), fill  = c('green', 'orange', 'blue'))
# Prevalence is number of infected at a time
# Incidence is the number of new cases
plot(t, I, type = 'l', main = 'Prevalence and incidence', ylab = 'Number of individuals', col = 'orange', xlab = 'Time')
legend(2008, 400, legend = c('Prevalence', 'Incidence'), fill = c('orange', 'red'))
lines(t, incidence, col = 'red')

# Plot R over time
birth <- median(data$birthEs)
birth_lower <- quantile(data$birthEs, 0.025)
birth_upper <- quantile(data$birthEs, 0.975)
become_uninfectious <- median(data$becomeUninfectiousRateEs)
become_uninfectious_lower <- quantile(data$becomeUninfectiousRateEs, 0.025)
become_uninfectious_upper <- quantile(data$becomeUninfectiousRateEs, 0.975)
Re <- vector(length = intervals)
Re_lower <- vector(length = intervals)
Re_upper <- vector(length = intervals)
for(i in 1:intervals){
    Re[i] <- birth * S[i] / become_uninfectious
    Re_lower[i] <- birth_lower * S[i] / become_uninfectious_lower
    Re_upper[i] <- birth_upper * S[i] / become_uninfectious_upper
}
plot(t, Re, type = 'l', main = 'Reproductive ratio over time', xlab = 'Time')
lines(t, rep(1, length(t)), lty = 2)
