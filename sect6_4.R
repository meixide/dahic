# This R script replicates Model 25 from Section 6.2 of the "Causal Dantzig" paper by D. Rothenhäusler et al.
# with an additional test environment e = 2. It also reproduces the results in Section 6.4 of the "Generative
# Invariance" paper by Meixide and Insua. The script serves to compare the performance of Generative Invariance (GI),
# Causal Dantzig, and Instrumental Variables (IV) in terms of their prediction accuracy in an unseen domain e = 2.
# The script runs num_runs iterations to calculate the (MSEs) for each method in e=2.

library(ivreg)

source('hd.R')

t=2 # Index of test environment (could be non-integer)

# Parameters
n <- 1000 # Number of samples
num_runs <- 100 # Number of iterations

# Initialize vectors to store MSEs for each run
mseiv_vector <- numeric(num_runs)
msegi_vector <- numeric(num_runs)
msedantz_vector <- numeric(num_runs)

for (run in 1:num_runs) {
  # Simulate exogenous random variables
  H <- rnorm(n, mean = 0, sd = 1)         # H ~ N(0, 1)
  epsilon1 <- rnorm(n, mean = 0, sd = 1)  # eps1 ~ N(0, 1)
  epsilon2 <- rnorm(n, mean = 0, sd = 1)  # eps2 ~ N(0, 1)
  epsilon3 <- rnorm(n, mean = 0, sd = 1)  # eps3 ~ N(0, 1)
  ENV <- sample(c(0, 1,t), n, replace = TRUE) # ENV  in 0, 1, 2
  
  # Compute X and Y
  X <- H + 2 * ENV * (0.25 + epsilon3) + epsilon1
  # or X <- H + 2 *( ENV * 0.25 + 1)*epsilon3 + epsilon1 as in page 1701
  Y <- 2 * X + H + 2 * epsilon2
  
  # Store the results in a data frame
  datasim <- data.frame(H, epsilon1, epsilon2, epsilon3, ENV, X, Y)
  datatrain <- datasim[datasim$ENV != t, ]
  
  mean(datatrain$X[datatrain$ENV==0])
  mean(datatrain$X[datatrain$ENV==1])
  
  
  # Linear and instrumental variable models
  mlm <- lm(Y ~ X, data = datatrain)
  miv <- ivreg(Y ~ X | ENV, data = datatrain)
  
  # Generate predictions for ENV == 2
  xnew <- as.numeric(X[ENV == t])
  gi <- hd(datatrain$Y, data.frame(X = datatrain$X, ENV = datatrain$ENV))
  ygi <- generator(gi$betahat, gi$khat, data.frame(xnew), 0, distributional = FALSE)
  
  # Prepare data for Causal Dantzig estimation
  data_drig <- list()
  unique_indicators <- c(0, 1) # Training environments
  
  for (value in unique_indicators) {
    rows <- which(datasim$ENV == value)
    combined_matrix <- as.matrix(cbind(Y[rows], cbind(1, X[rows])))
    data_drig[[as.character(value)]] <- combined_matrix
  }
  
  # Causal Dantzig 
  dantz <- est_drig(data_drig, gamma = Inf, y_idx = 1, del_idx = NULL, unif_weight = FALSE)

  # Predicted Y for ENV == 2
  ydantz <- dantz[1] + dantz[2] * xnew
  yiv <- coef(miv)[1] + coef(miv)[2] * xnew
  ytrue <- Y[ENV == t]
  
  # Compute Mean Squared Errors
  mseiv_vector[run] <- mean((yiv - ytrue)^2)
  msegi_vector[run] <- mean((ygi - ytrue)^2)
  msedantz_vector[run] <- mean((ydantz - ytrue)^2)
  
  if(!(run %% 10)) print(run)
}

# Display the results
q=0.75
median(mseiv_vector)
c(quantile(mseiv_vector,1-q),quantile(mseiv_vector,q))
median(msegi_vector)
c(quantile(msegi_vector,1-q),quantile(msegi_vector,q))
median(msedantz_vector)
c(quantile(msedantz_vector,1-q),quantile(msedantz_vector,q))

# Mean and standard deviation
mean(mseiv_vector)
sd(mseiv_vector)
mean(msegi_vector)
sd(msegi_vector)
mean(msedantz_vector)
sd(msedantz_vector)











