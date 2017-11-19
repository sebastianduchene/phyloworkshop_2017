

```{r}
#options(repr.plot.width=4, repr.plot.height=4)
```

## Opening the MCMC blackbox: a Bayesian linear model

- In this example, we will simulate some data according to a linear model. The data will follow the regression line 

y = mx + b + e

Where *m* is the slope, *b* is the intercept, and *e* is the error term. In particular, we expect the points to be normally distributed along this line, so we will model the error term as a normal distribution with an unknown standard deviation.

Our simulation will take the form:

y = 3*x + 15 + N(0, 7)

The slope is 3, the *y* intercept is 15, and the error is described by a normal distribution with mean of 0 and standard deviation of 7.

We will first generate the data and plot it:


```{r , dev='pdf'}
set.seed(121432)
x <- runif(100, 0, 100)
y <- 3*x + 15 + rnorm(length(x), sd = 7)
plot(x, y)
```


Next we need to define some functions for the MCMC. If you are well acquainted with R, you can follow the code and ask any questions. 

**The Likelihood**

- Notice that the likelihood is taking the a set of given values for the slope, intercept, and error to to predict y. It then calculates the difference between the predicted value of y and the actual value of y and calculates their density according to a normal distribution with the given standard deviation:


```{r , dev='pdf'}
likelihood <- function(par){
    y_predicted <- par[1] * x + par[2]
    single_likelihoods <- dnorm(y, mean = y_predicted, sd = par[3], log = T)
    return(sum(single_likelihoods))
}
```

**The prior**

- For the prior, we simply calculate the density of the three parameters in our model given any distribution of our choice. Our prior will initially consist or normal distributions with means of 3, 15, and 10 for the slope, intercept, and error term, respectively. All with a standard deviation of 10. In a later stage we will modify the prior. Note that we do not use the data to calculate the prior.


```{r}
prior <- function(par){
    prior_slope <- dnorm(par[1], mean = 3, sd = 10, log = T)
    prior_intercept <- dnorm(par[2], mean = 15, sd = 10, log = T)
    prior_error <- dnorm(par[3], mean = 10, sd = 10, log = T)
    return(prior_slope+prior_intercept+prior_error)
}
```

**The proposal function**
    
- This is a function that will modify our parameter values at each MCMC step. It dictates how well the chain will move through parameter space. Concretely, we will add random numbers from a standard normal distribution for each parameter, with the condition that the error term cannot be negative because it corresponds to a standard devation, which always has positive values.


```{r}
proposal_function <- function(par){
    slope_proposal <- par[1] + rnorm(1, mean = 0, sd = 1)
    intercept_proposal <- par[2] + rnorm(1, mean = 0, sd = 1)
    error_proposal <- abs(par[3] + rnorm(1, mean = 0, sd = 1))
    return(c(slope_proposal, intercept_proposal, error_proposal))
}
```

Next we will define a function to run the MCMC. We have not defined the posterior, but we will do so internally by multiplying the prior and the posterior. The MCMC will require two inputs, a set of parameter starting values, and the number of iterations. We will modify these through the prac to understand how the algorithm works.


```{r}
start_values <- c(1, 1, 1)
iterations <- 100000

run_mcmc <- function(start_values, iterations, sample_from_prior = F){
    chain <- matrix(NA, iterations, 6)
    prior_start <- prior(start_values)
    likelihood_start <- likelihood(start_values)
    if(sample_from_prior) likelihood_start <- 0
    posterior_start <- prior_start+likelihood_start
    chain[1, ] <- c(prior_start, likelihood_start, posterior_start, start_values)
    colnames(chain) <- c('prior', 'likelihood', 'posterior', 'slope', 'intercept','error')
    for(i in 2:iterations){
        proposal <- proposal_function(chain[i-1, 4:6])
        prior_temp <- prior(proposal)
        likelihood_temp <- likelihood(proposal)
        if(sample_from_prior) likelihood_temp <- 0
        posterior_temp <- prior_temp + likelihood_temp
        mh_ratio <- exp(posterior_temp - chain[i-1, 'posterior'])
        if(mh_ratio > 1){ # Accept if the new proposal has a higher posterior than the previous step
            chain[i, ] <- c(prior_temp, likelihood_temp, posterior_temp, proposal)
        }else if(mh_ratio > runif(1)){ # Accept if the new proposal has a higher posterior than a random number from a uniform distribuiton
            chain[i, ] <- c(prior_temp, likelihood_temp, posterior_temp, proposal)
        }else{ # Reject proposal and return to previous step.
            chain[i, ] <- chain[i-1, ]
        }
    }
    return(chain)
}
```

Assessing sufficient sampling of an MCMC is very difficult, and it is often helpful to visualise the run. Below we will define a function to plot the output of the MCMC. It requires the user to specify an MCMC chain, a burnin, and whether it should use colours for plotting.


```{r}
plot_mcmc <- function(chain, burnin, use_heat_colors = T, col_line = rgb(0, 0, 0, 0.5)){
    ml_fit <- summary(lm(y ~ x))
    ml_estimates <- ml_fit$coefficients[, 'Estimate']

    samples <- (1:nrow(chain))[-c(1:burnin)]

    par(mfrow = c(3, 3))
    plot(chain[-c(1:burnin), 'slope'], type = 'l', xlab = 'step', ylab = 'Slope')
    lines( samples, rep(ml_estimates[2], length(samples)), col = 'red', lwd = 2)
    plot(chain[-c(1:burnin), 'intercept'], type = 'l', xlab = 'step', ylab = 'Intercept')
    lines( samples, rep(ml_estimates[1], length(samples)), col = 'red', lwd = 2)
    plot(chain[-c(1:burnin), 'error'], type = 'l', xlab = 'step', ylab = 'Error')

    hist(chain[-c(1:burnin), 'slope'], xlab = 'Slope', main = '')
    lines(x = rep(ml_estimates[2], 2), y = c(0, 10000), col = 'red', lwd = 3)
    hist(chain[-c(1:burnin), 'intercept'], xlab = 'Intercept', main = '')
    lines(x = rep(ml_estimates[1], 2), y = c(0, 10000), col = 'red', lwd = 3)

    hist(chain[-c(1:burnin), 'error'], xlab = 'Error', main = '')
    
    cols <-  rev( heat.colors(nrow(chain), alpha = 1.0))

    if(use_heat_colors){
       plot(chain[order(chain[,'posterior'],decreasing =F),  
                  c('slope', 'intercept')],pch = 20, type = 'b', col = cols)#0, 0, 0.2))
    }else{
        plot(chain[order(chain[,'posterior'],decreasing =F),  
               c('slope', 'intercept')],pch = 20, type = 'l', col = rgb(0, 0, 0, 0.5))
    }

    plot(x, y)
    for(i in 1:length(samples)){
    abline(b =  mean( chain[samples[i], 'slope']), 
           a = mean( chain[samples[i], 'intercept']),col = rgb(1, 0, 0, 0.2), lwd = 0.5)


    }    
}
```
**Sufficient sampling**

- In our first example, we will run the chain with 100 iterations. The starting values for the parameters will be 1, 1, 1, for the slope, intercept, and error. Because this is a stochastic algorithm, the starting values should not matter. Feel free to try differnet values to see this. We will then plot the results from the chain, with a burnin of 10 samples.


```{r}
#options(repr.plot.width=7,repr.plot.height=7)
```


```{r , dev='pdf'}
chain <- run_mcmc(start_values = c(1, 1, 1), iterations = 100, sample_from_prior = F)
plot_mcmc(chain, 10, use_heat_colors = F)

```


The first row of plots are the traces for each parameter, the second row is a histogram of the post-burnin samples, which correspond the the posterior distribution for each parameter. The red lines correspond to the maximum-likelihood esimates, for comparison.

The third row shows the trace of the intercept vs. the slope. If you used use_heat_colors =T,  the brighter colours corresponding to regions with higher posterior probability. The last plot is our original data, and regression lines sampled from the posterior, which is the posterior distribution of the regression. 

** Question 5.1**
Judging by the traces, do you think that we have a reliable sample from the posterior?

- We will now run the chain with 50,000 iterations. We will also increase the burnin accordingly.


```{r , dev='pdf'}
chain <- run_mcmc(start_values = c(1, 1, 1), iterations = 50000, sample_from_prior = F)
plot_mcmc(chain, 5000, use_heat_colors = F)
```


**Question 5.2** Has increasing the number of iterations improved our confidence in our estimates? How do we know this?

**The proposal function**

- We will now tweak the proposal function to see its effect on the MCMC. For example, we will allow it to make very large moves by increasing the standard deviation from 1 to 10. To do this, we need to redefine the proposal function:


```{r , dev='pdf'}
proposal_function <- function(par){
    slope_proposal <- par[1] + rnorm(1, mean = 0, sd = 10)
    intercept_proposal <- par[2] + rnorm(1, mean = 0, sd = 10)
    error_proposal <- abs(par[3] + rnorm(1, mean = 0, sd = 10))
    return(c(slope_proposal, intercept_proposal, error_proposal))
}

chain <- run_mcmc(start_values = c(1, 1, 1), iterations = 50000, sample_from_prior = F)
plot_mcmc(chain, 5000, use_heat_colors = F)
```

**Question 5.3** Why has increasing the size of the moves in the proposal function made the MCMC less efficient?


It should be clear that our previous proposal function was better. We will redefine it for the rest of the prac as we did before:


```{r}
proposal_function <- function(par){
    slope_proposal <- par[1] + rnorm(1, mean = 0, sd = 1)
    intercept_proposal <- par[2] + rnorm(1, mean = 0, sd = 1)
    error_proposal <- abs(par[3] + rnorm(1, mean = 0, sd = 1))
    return(c(slope_proposal, intercept_proposal, error_proposal))
}


```

**Prior sensitivity: a misleading prior**

- An other factor that is worth inspecting is the influence of the priors. We will set some clearly misleading and very informative priors for the slope and the intercept (we can also do this for the error term, but it is more tractable to do this for a few parameters at a time) and re run our MCMC.

    - slope: N(10, 0.01)
    - intercept: N(1, 0.01)


```{r}
prior <- function(par){
    prior_slope <- dnorm(par[1], mean = 10, sd = 0.01, log = T)
    prior_intercept <- dnorm(par[2], mean = 1, sd = 0.01, log = T)
    prior_error <- dnorm(par[3], mean = 10, sd = 10, log = T)
    return(prior_slope+prior_intercept+prior_error)
}
```


```{r , dev='pdf'}
chain <- run_mcmc(start_values = c(1, 1, 1), iterations = 50000, sample_from_prior = F)
plot_mcmc(chain, 5000, use_heat_colors = F, col_line = 0.1)
```

**Prior sensitivity: uninformative priors**

The MCMC has struggled to obtain sufficient samples from the posterior. The error term appears to increase indefinitely. Although we could possibly address this by using more iterations, note that the regression line deviates substantially from what we obtained using more reasonable priors. Our priors here are overly influential and might mislead our inferences.

- It is often the case, that the researcher has no grounds for choosing one prior over an other. In such cases, it might be appropriate to choose an 'uninformative' prior. Redefine the priors using uniform distributions with -100 and 100 and minimum and maximum bounds respectively.


```{r}
prior <- function(par){
    prior_slope <- dunif(par[1], min = -100, max = 100,log = T)
    prior_intercept <- dunif(par[2], min = -100, max = 100, log = T)
    prior_error <- dunif(par[3], min = -100, max = 100, log = T)
    return(prior_slope+prior_intercept+prior_error)
}

```

```{r , dev='pdf'}
chain <- run_mcmc(start_values = c(1, 1, 1), iterations = 100000, sample_from_prior = F)
plot_mcmc(chain, 5000, use_heat_colors = F, col_line = 0.1)
```

**Prior sensitivity: sampling from the uniform prior**

- It appears that our data are sufficiently informative to recover the correct estimates, even with uninformative priors. In practice, however, it is often difficult to determine whether the prior is overly influential, and the extent to which the data are driving the estimates. An easy way to do this is to sample from the prior. This consists in running the MCMC with no data, or by fixing the likelihood to an arbitrary value. To do this, we will set the option sample_from_prior to TRUE. This simply assigns a likelihood of 0, regardless of the parameter values. We will also use the prior from our second set of analyses. This analysis usually runs much faster becuase it involves fewer calculations:




```{r}
chain <- run_mcmc(start_values = c(1, 1, 1), iterations = 150000, sample_from_prior = T)
```


```{r , dev='pdf'}
plot_mcmc(chain, 5000, use_heat_colors = F, col_line = 0.1)
```

**Prior sensitivity: sampling from the normal prior**

- For comparisson, we will also sample from the prior that we used initially.

```{r}
prior <- function(par){
    prior_slope <- dnorm(par[1], mean = 3, sd = 10, log = T)
    prior_intercept <- dnorm(par[2], mean = 15, sd = 10, log = T)
    prior_error <- dnorm(par[3], mean = 10, sd = 10, log = T)
    return(prior_slope+prior_intercept+prior_error)
}
```

```{r , dev='pdf'}
chain <- run_mcmc(start_values = c(1, 1, 1), iterations = 150000, sample_from_prior = T)
plot_mcmc(chain, 5000, use_heat_colors = F, col_line = 0.1)
```
