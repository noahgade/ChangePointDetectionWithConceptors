library(parallel)
library(hdbinseg)
library(ecp)
library(tidyverse)
library(mclust)
library(xtable)

# Load conceptorCP package locally
library(conceptorCP)

##################################################
# Functions for data generation
##################################################

# Gaussian Data
genGaussianData <- function(Length, change_point, mean_change = c(0, 0), var_change = c(1, 1), cov_change = 0) {
  if(!is.null(change_point)) {
    Section1 <- dplyr::as_tibble(matrix(MASS::mvrnorm(change_point, c(0, 0), matrix(c(1, 0, 0, 1), nrow = 2)), nrow = change_point))
    Section2 <- dplyr::as_tibble(matrix(MASS::mvrnorm(Length - change_point, mean_change, matrix(c(var_change[1], cov_change, cov_change, var_change[2]), nrow = 2)), nrow = Length - change_point))
    
    dataset <- list(data = dplyr::bind_rows(Section1, Section2), 
                    parameters = list(Length = Length, change_point = change_point, mean_change = mean_change, var_change = var_change, cov_change = cov_change))
  } else {
    Section1 <- dplyr::as_tibble(matrix(MASS::mvrnorm(Length, c(0, 0), matrix(c(1, 0, 0, 1), nrow = 2)), nrow = Length))
    dataset <- list(data = Section1, 
                    parameters = list(Length = Length, change_point = change_point, mean_change = mean_change, var_change = var_change, cov_change = cov_change))
  }
  return(dataset)
}

# Vector Autoregressive Data
genVARData <- function(Length, change_point, lag, spec_rad1, spec_rad2, noise = 1) {
  repeat{
    A1 <- matrix(runif(4 * lag, -1, 1), nrow = 2)
    A2 <- matrix(0, nrow = 2 * lag - 2, ncol = 2 * lag)
    diag(A2) <- 1
    A <- rbind(A1, A2)
    eA <- max(abs(eigen(A)$values))
    if(eA < 1 && abs(eA - spec_rad1) < 0.02){
      break
    }
  }
  
  repeat{
    B1 <- matrix(runif(4 * lag, -1, 1), nrow = 2)
    B2 <- matrix(0, nrow = 2 * lag - 2, ncol = 2 * lag)
    diag(B2) <- 1
    B <- rbind(B1, B2)
    eB <- max(abs(eigen(B)$values))
    if(eB < 1 && abs(eB - spec_rad2) < 0.02){
      break
    }
  }
  
  washoutperiod <- 100
  Series <- dplyr::as_tibble(rbind(matrix(rnorm(2 * lag, 0, 1), nrow = 1), matrix(0, ncol = 2 * lag, nrow = Length + washoutperiod)))
  Noise <- dplyr::as_tibble(cbind(matrix(rnorm(2 * (Length + washoutperiod), 0, noise), ncol = 2), matrix(0, ncol = 2 * lag - 2, nrow = Length + washoutperiod)))
  
  if(!is.null(change_point)) {
    for(t in 1:(Length + washoutperiod)) {
      if(t <= (change_point + washoutperiod)) {
        Series[t + 1,] <- as.matrix(Series[t,]) %*% t(A) + Noise[t,]
      } else {
        Series[t + 1,] <- as.matrix(Series[t,]) %*% t(B) + Noise[t,]
      }
    }
  } else {
    for(t in 1:(Length + washoutperiod)) {
      Series[t + 1,] <- as.matrix(Series[t,]) %*% t(A) + Noise[t,]
    }
  }
  
  dataset <- list(data = dplyr::select(Series[-c(1:101),], c(1, 2)), 
                  parameters = list(Length = Length, change_point = change_point, lag = lag, spec_rad1 = eA, spec_rad2 = eB, Matrix1 = A, Matrix2 = B, noise = noise))
  return(dataset)
}

# Periodic Data
genPeriodicData <- function(Length, change_point, mean_change = c(0, 0), amplitude = c(1, 1), frequency = c(2 * pi, 2 * pi), noise = 1) {
  if(!is.null(change_point)) {
    Section1 <- dplyr::as_tibble(cbind(amplitude[1] * cos(frequency[1] * 1:change_point) + rnorm(change_point, 0, noise), amplitude[1] * sin(frequency[1] * 1:change_point) + rnorm(change_point, 0, noise)))
    Section2 <- dplyr::as_tibble(cbind(amplitude[2] * cos(frequency[2] * (change_point + 1):Length) + rnorm(Length - change_point, 0, noise) + mean_change[1],
                                       amplitude[2] * sin(frequency[2] * (change_point + 1):Length) + rnorm(Length - change_point, 0, noise) + mean_change[2]))
    dataset <- list(data = dplyr::bind_rows(Section1, Section2), 
                    parameters = list(Length = Length, change_point = change_point, mean_change = mean_change, amplitude = amplitude, frequency = frequency, noise = noise))
  } else {
    Section1 <- dplyr::as_tibble(cbind(amplitude[1] * cos(frequency[1] * 1:Length) + rnorm(Length, 0, noise), amplitude[1] * sin(frequency[1] * 1:Length) + rnorm(Length, 0, noise)))
    dataset <- list(data = Section1, 
                    parameters = list(Length = Length, change_point = change_point, mean_change = mean_change, amplitude = amplitude, frequency = frequency, noise = noise))
  }
  return(dataset)
}

# Gaussian Process Data
genGPData <- function(Length, change_point, mean = c(0, 0), reverting = c(1, 1), volatility = c(1, 1)) {
  Series1 <- vector('numeric', Length)
  Series2 <- vector('numeric', Length)
  Series1[1] <- rnorm(1, mean[1], volatility[1])
  Series2[1] <- rnorm(1, mean[1], volatility[1])
  
  if(!is.null(change_point)) {
    for(i in 1:(change_point - 1)) {
      Series1[i + 1] <- Series1[i] + reverting[1] * (mean[1] - Series1[i]) + rnorm(1, 0, volatility[1])
      Series2[i + 1] <- Series2[i] + reverting[1] * (mean[1] - Series2[i]) + rnorm(1, 0, volatility[1])
    }
    for(i in change_point:(Length - 1)) {
      Series1[i + 1] <- Series1[i] + reverting[2] * (mean[2] - Series1[i]) + rnorm(1, 0, volatility[2])
      Series2[i + 1] <- Series2[i] + reverting[2] * (mean[2] - Series2[i]) + rnorm(1, 0, volatility[2])
    }
    dataset <- list(data = dplyr::bind_cols(V1 = Series1, V2 = Series2), 
                    parameters = list(Length = Length, change_point = change_point, mean = mean, reverting = reverting, volatility = volatility))
  } else {
    for(i in 1:(Length - 1)) {
      Series1[i + 1] <- Series1[i] + reverting[1] * (mean[1] - Series1[i]) + rnorm(1, 0, volatility[1])
      Series2[i + 1] <- Series2[i] + reverting[1] * (mean[1] - Series2[i]) + rnorm(1, 0, volatility[1])
    }
    dataset <- list(data = dplyr::bind_cols(V1 = Series1, V2 = Series2), 
                    parameters = list(Length = Length, change_point = change_point, mean = mean, reverting = reverting, volatility = volatility))
  }
  return(dataset)
}

##################################################
# Simulation Data Generation
##################################################

set.seed(10242022)
num_sims <- 300
Length <- 1000
lower_limit <- 181
upper_limit <- 999

ID1 <- replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 1, 0.5, 0.5, noise = 0.5), simplify = F) %>% 
  append(replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 1, 0.5, 0.8, noise = 0.5), simplify = F)) %>% 
  append(replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 1, 0.8, 0.5, noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 1, 0.8, 0.8, noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genVARData(Length, change_point = NULL, lag = 1, 0.5, 0.5, noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genVARData(Length, change_point = NULL, lag = 1, 0.8, 0.8, noise = 0.5), simplify = F))

ID2 <- replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 2, 0.5, 0.5, noise = 0.5), simplify = F) %>% 
  append(replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 2, 0.5, 0.8, noise = 0.5), simplify = F)) %>% 
  append(replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 2, 0.8, 0.5, noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genVARData(Length, sample(lower_limit:upper_limit, 1), lag = 2, 0.8, 0.8, noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genVARData(Length, change_point = NULL, lag = 2, 0.5, 0.5, noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genVARData(Length, change_point = NULL, lag = 2, 0.8, 0.8, noise = 0.5), simplify = F))

ID3 <- replicate(num_sims, genPeriodicData(Length, sample(lower_limit:upper_limit, 1), frequency = c(1, 0.5), noise = 0.5), simplify = F) %>%
  append(replicate(num_sims, genPeriodicData(Length, sample(lower_limit:upper_limit, 1), frequency = c(1, 0.8), noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genPeriodicData(Length, sample(lower_limit:upper_limit, 1), frequency = c(1, 1.2), noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genPeriodicData(Length, sample(lower_limit:upper_limit, 1), frequency = c(1, 1.5), noise = 0.5), simplify = F)) %>%
  append(replicate(num_sims, genPeriodicData(Length, change_point = NULL, frequency = c(1, 1), noise = 0.5), simplify = F))

ID4 <- replicate(num_sims, genGPData(Length, sample(lower_limit:upper_limit, 1), reverting = c(0.5, 0), volatility = c(0.5, 0.5)), simplify = F) %>%
  append(replicate(num_sims, genGPData(Length, sample(lower_limit:upper_limit, 1), reverting = c(0.5, 1), volatility = c(0.5, 0.5)), simplify = F)) %>%
  append(replicate(num_sims, genGPData(Length, sample(lower_limit:upper_limit, 1), reverting = c(1, 0), volatility = c(0.5, 0.5)), simplify = F)) %>%
  append(replicate(num_sims, genGPData(Length, sample(lower_limit:upper_limit, 1), reverting = c(1, 0.5), volatility = c(0.5, 0.5)), simplify = F)) %>%
  append(replicate(num_sims, genGPData(Length, sample(lower_limit:upper_limit, 1), reverting = c(0.5, 0.5), volatility = c(0.5, 0.2)), simplify = F)) %>%
  append(replicate(num_sims, genGPData(Length, sample(lower_limit:upper_limit, 1), reverting = c(0.5, 0.5), volatility = c(0.5, 0.8)), simplify = F)) %>%
  append(replicate(num_sims, genGPData(Length, sample(lower_limit:upper_limit, 1), reverting = c(0.5, 0.5), volatility = c(0.5, 1)), simplify = F)) %>%
  append(replicate(num_sims, genGPData(Length, change_point = NULL, reverting = c(0.5, 0.5), volatility = c(0.5, 0.5)), simplify = F)) %>%
  append(replicate(num_sims, genGPData(Length, change_point = NULL, reverting = c(1, 1), volatility = c(0.5, 0.5)), simplify = F))

ID5 <- replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), mean_change = c(0.5, 0.5)), simplify = F) %>%
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), mean_change = c(0.8, 0.8)), simplify = F)) %>%
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), mean_change = c(1.0, 1.0)), simplify = F)) %>%
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), var_change = c((0.5)^2, (0.5)^2)), simplify = F)) %>%
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), var_change = c((0.8)^2, (0.8)^2)), simplify = F)) %>%
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), var_change = c((1.2)^2, (1.2)^2)), simplify = F)) %>%
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), var_change = c((1.5)^2, (1.5)^2)), simplify = F)) %>%
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), cov_change = -0.8), simplify = F)) %>% 
  append(replicate(num_sims, genGaussianData(Length, sample(lower_limit:upper_limit, 1), cov_change = 0.8), simplify = F)) %>% 
  append(replicate(num_sims, genGaussianData(Length, change_point = NULL), simplify = F))

save(ID1, file = "ID1.RData")
save(ID2, file = "ID2.RData")
save(ID3, file = "ID3.RData")
save(ID4, file = "ID4.RData")
save(ID5, file = "ID5.RData")

# Simulation data is available upon request.

##################################################
# conceptorCP method on simulated data
##################################################

CCP02 <- function(data) {
  result <- ccp(data$data, washoutL = 60, trainL = 120, tol = 0.02, plot.it = F)
  out <- list(estimate = result$estimate, MBBquant = result$MBBquant, MBBnull = result$MBBnull, MBBblockL = result$MBBblockL, angles = result$angles, statSeries = result$statSeries,
              netParams = list(washoutL = result$netParams$washoutL, trainL = result$netParams$trainL))
  return(out)
}

CCP04 <- function(data) {
  result <- ccp(data$data, washoutL = 60, trainL = 120, tol = 0.04, plot.it = F)
  out <- list(estimate = result$estimate, MBBquant = result$MBBquant, MBBnull = result$MBBnull, MBBblockL = result$MBBblockL, angles = result$angles, statSeries = result$statSeries,
              netParams = list(washoutL = result$netParams$washoutL, trainL = result$netParams$trainL))
  return(out)
}

CCP08 <- function(data) {
  result <- ccp(data$data, washoutL = 60, trainL = 120, tol = 0.08, plot.it = F)
  out <- list(estimate = result$estimate, MBBquant = result$MBBquant, MBBnull = result$MBBnull, MBBblockL = result$MBBblockL, angles = result$angles, statSeries = result$statSeries,
              netParams = list(washoutL = result$netParams$washoutL, trainL = result$netParams$trainL))
  return(out)
}

CCP16 <- function(data) {
  result <- ccp(data$data, washoutL = 60, trainL = 120, tol = 0.16, plot.it = F)
  out <- list(estimate = result$estimate, MBBquant = result$MBBquant, MBBnull = result$MBBnull, MBBblockL = result$MBBblockL, angles = result$angles, statSeries = result$statSeries,
              netParams = list(washoutL = result$netParams$washoutL, trainL = result$netParams$trainL))
  return(out)
}

load("ID1.RData")

ccp02_1 <- mclapply(ID1, CCP02, mc.cores = numCores)
save(ccp02_1, file = "ccp02_1.RData")

ccp04_1 <- mclapply(ID1, CCP04, mc.cores = numCores)
save(ccp04_1, file = "ccp04_1.RData")

ccp08_1 <- mclapply(ID1, CCP08, mc.cores = numCores)
save(ccp08_1, file = "ccp08_1.RData")

ccp16_1 <- mclapply(ID1, CCP16, mc.cores = numCores)
save(ccp16_1, file = "ccp16_1.RData")

load("ID2.RData")

ccp02_2 <- mclapply(ID2, CCP02, mc.cores = numCores)
save(ccp02_2, file = "ccp02_2.RData")

ccp04_2 <- mclapply(ID2, CCP04, mc.cores = numCores)
save(ccp04_2, file = "ccp04_2.RData")


ccp08_2 <- mclapply(ID2, CCP08, mc.cores = numCores)
save(ccp08_2, file = "ccp08_2.RData")

ccp16_2 <- mclapply(ID2, CCP16, mc.cores = numCores)
save(ccp16_2, file = "ccp16_2.RData")

load("ID3.RData")

ccp02_3 <- mclapply(ID3, CCP02, mc.cores = numCores)
save(ccp02_3, file = "ccp02_3.RData")

ccp04_3 <- mclapply(ID3, CCP04, mc.cores = numCores)
save(ccp04_3, file = "ccp04_3.RData")

ccp08_3 <- mclapply(ID3, CCP08, mc.cores = numCores)
save(ccp08_3, file = "ccp08_3.RData")

ccp16_3 <- mclapply(ID3, CCP16, mc.cores = numCores)
save(ccp16_3, file = "ccp16_3.RData")

load("ID4.RData")

ccp02_4 <- mclapply(ID4, CCP02, mc.cores = numCores)
save(ccp02_4, file = "ccp02_4.RData")

ccp04_4 <- mclapply(ID4, CCP04, mc.cores = numCores)
save(ccp04_4, file = "ccp04_4.RData")

ccp08_4 <- mclapply(ID4, CCP08, mc.cores = numCores)
save(ccp08_4, file = "ccp08_4.RData")

ccp16_4 <- mclapply(ID4, CCP16, mc.cores = numCores)
save(ccp16_4, file = "ccp16_4.RData")

load("ID5.RData")

ccp02_5 <- mclapply(ID5, CCP02, mc.cores = numCores)
save(ccp02_5, file = "ccp02_5.RData")

ccp04_5 <- mclapply(ID5, CCP04, mc.cores = numCores)
save(ccp04_5, file = "ccp04_5.RData")

ccp08_5 <- mclapply(ID5, CCP08, mc.cores = numCores)
save(ccp08_5, file = "ccp08_5.RData")

ccp16_5 <- mclapply(ID5, CCP16, mc.cores = numCores)
save(ccp16_5, file = "ccp16_5.RData")

##################################################
# comparison methods on simulated data
##################################################

load("ID1.RData")
load("ID2.RData")
load("ID3.RData")
load("ID4.RData")
load("ID5.RData")

EDivisive <- function(data) {
  L <- nrow(data$data)
  result <- e.divisive(data$data[181:L,], sig.lvl = 0.05)
  return(result)
}

ed1 <- mclapply(ID1, EDivisive, mc.cores = numCores)
save(ed1, file = "ed1.RData")

ed2 <- mclapply(ID2, EDivisive, mc.cores = numCores)
save(ed2, file = "ed2.RData")

ed3 <- mclapply(ID3, EDivisive, mc.cores = numCores)
save(ed3, file = "ed3.RData")

ed4 <- mclapply(ID4, EDivisive, mc.cores = numCores)
save(ed4, file = "ed4.RData")

ed5 <- mclapply(ID5, EDivisive, mc.cores = numCores)
save(ed5, file = "ed5.RData")

KCP <- function(data) {
  Len <- nrow(data$data)
  result <- kcpa(as.matrix(data$data[181:Len,]), L = 1, C = 2) - 1
  return(result)
}

kcp1 <- mclapply(ID1, KCP, mc.cores = numCores)
save(kcp1, file = "kcp1.RData")

kcp2 <- mclapply(ID2, KCP, mc.cores = numCores)
save(kcp2, file = "kcp2.RData")

kcp3 <- mclapply(ID3, KCP, mc.cores = numCores)
save(kcp3, file = "kcp3.RData")

kcp4 <- mclapply(ID4, KCP, mc.cores = numCores)
save(kcp4, file = "kcp4.RData")

kcp5 <- mclapply(ID5, KCP, mc.cores = numCores)
save(kcp5, file = "kcp5.RData")

SBS1 <- function(data) {
  L <- nrow(data$data)
  result01 <- sbs.alg(t(data$data[181:L,]), cp.type = 1, q = 0.01, do.parallel = 0)
  result05 <- sbs.alg(t(data$data[181:L,]), cp.type = 1, q = 0.05, do.parallel = 0)
  result10 <- sbs.alg(t(data$data[181:L,]), cp.type = 1, q = 0.10, do.parallel = 0)
  result20 <- sbs.alg(t(data$data[181:L,]), cp.type = 1, q = 0.20, do.parallel = 0)
  result <- list(result01 = result01, result05 = result05, result10 = result10, result20 = result20)
  return(result)
}

sbsmean1 <- mclapply(ID1, SBS1, mc.cores = numCores)
save(sbsmean1, file = "sbsmean1.RData")

sbsmean2 <- mclapply(ID2, SBS1, mc.cores = numCores)
save(sbsmean2, file = "sbsmean2.RData")

sbsmean3 <- mclapply(ID3, SBS1, mc.cores = numCores)
save(sbsmean3, file = "sbsmean3.RData")

sbsmean4 <- mclapply(ID4, SBS1, mc.cores = numCores)
save(sbsmean4, file = "sbsmean4.RData")

sbsmean5 <- mclapply(ID5, SBS1, mc.cores = numCores)
save(sbsmean5, file = "sbsmean5.RData")

SBS2 <- function(data) {
  L <- nrow(data$data)
  result01 <- sbs.alg(t(data$data[181:L,]), cp.type = 2, q = 0.01, do.parallel = 0)
  result05 <- sbs.alg(t(data$data[181:L,]), cp.type = 2, q = 0.05, do.parallel = 0)
  result10 <- sbs.alg(t(data$data[181:L,]), cp.type = 2, q = 0.10, do.parallel = 0)
  result20 <- sbs.alg(t(data$data[181:L,]), cp.type = 2, q = 0.20, do.parallel = 0)
  result <- list(result01 = result01, result05 = result05, result10 = result10, result20 = result20)
  return(result)
}

sbs1 <- mclapply(ID1, SBS2, mc.cores = numCores)
save(sbs1, file = "sbs1.RData")

sbs2 <- mclapply(ID2, SBS2, mc.cores = numCores)
save(sbs2, file = "sbs2.RData")

sbs3 <- mclapply(ID3, SBS2, mc.cores = numCores)
save(sbs3, file = "sbs3.RData")

sbs4 <- mclapply(ID4, SBS2, mc.cores = numCores)
save(sbs4, file = "sbs4.RData")

sbs5 <- mclapply(ID5, SBS2, mc.cores = numCores)
save(sbs5, file = "sbs5.RData")

# naive conceptor method included for study

NaiveConceptor <- function(data, trainL = "", washoutL = "", nboots = 240, plot.it = TRUE) {
  input <- as.matrix(scalein(data))
  L <- nrow(input)
  dim <- ncol(input)
  aperture <- 1
  C <- conceptorCalc(array(t(as.matrix(input[(washoutL + 1):(washoutL + trainL),])), dim = c(dim, 1, trainL)), aperture)
  projectedStates <- C[,,1] %*% t(as.matrix(input))
  angles <- angleCalc(array(t(as.matrix(input)), dim = c(dim, 1, L)), 
                      array(as.matrix(projectedStates), dim = c(dim, 1, L)))
  KSseries <- KSstatCalc(angles[(washoutL + trainL + 1):L,], 0.01)
  statistic <- max(KSseries)
  estimate <- which.max(KSseries) + washoutL + trainL
  
  Lstar <- washoutL
  M <- 40
  B <- 5
  L <- nrow(data)
  Lb <- ceiling(seq(L^(1/5), L^(1/2), length.out = B))
  split_input <- array(dim = c(L - M + 1, ncol(data), M * B))
  for (b in 1:B) {
    for (m in 1:M) {
      split_input[,,(b - 1) * M + m] <- bootdata(rbind(data[1:(washoutL + trainL),], data[(washoutL + trainL + m):(L - M + m), ]), washoutL, trainL, Lb[b]) 
    }
  }
  base_input <- replicate(M, bootdata(data, washoutL, trainL, Lstar))
  
  MBBests <- rep(NA, M*B)
  for(iter in 1:(M*B)) {
    tempC <- conceptorCalc(array(t(as.matrix(split_input[(washoutL + 1):(washoutL + trainL),,iter])), dim = c(dim, 1, trainL)), aperture)
    tempPS <- tempC[,,1] %*% t(as.matrix(split_input[,, iter]))
    tempA <- angleCalc(array(t(as.matrix(split_input[,, iter])), dim = c(dim, 1, L - M + 1)), 
                       array(as.matrix(tempPS), dim = c(dim, 1, L - M + 1)))
    MBBests[iter] <- max(KSstatCalc(tempA[(washoutL + trainL + 1):(L - M + 1),], 0.01))
  }
  
  MBBbase <- rep(NA, M)
  for(iter in 1:M) {
    tempC <- conceptorCalc(array(t(as.matrix(base_input[(washoutL + 1):(washoutL + trainL),, iter])), dim = c(dim, 1, trainL)), aperture)
    tempPS <- tempC[,,1] %*% t(as.matrix(base_input[,, iter]))
    tempA <- angleCalc(array(t(as.matrix(base_input[,, iter])), dim = c(dim, 1, L)), 
                       array(as.matrix(tempPS), dim = c(dim, 1, L)))
    MBBbase[iter] <- max(KSstatCalc(tempA[(washoutL + trainL + 1):L,], 0.01))
  }
  
  MBBblockL <- ceiling((L/(L - M + 1))^(1/3) * Lb[which.min(colSums(matrix((MBBests - mean(MBBbase))^2, nrow = M)))])
  
  binput <- replicate(nboots, bootdata(data, washoutL, trainL, MBBblockL))
  MBBnull <- rep(NA, nboots)
  for(iter in 1:nboots) {
    tempC <- conceptorCalc(array(t(as.matrix(binput[(washoutL + 1):(washoutL + trainL),, iter])), dim = c(dim, 1, trainL)), aperture)
    tempPS <- tempC[,,1] %*% t(as.matrix(binput[,, iter]))
    tempA <- angleCalc(array(t(as.matrix(binput[,, iter])), dim = c(dim, 1, L)), 
                       array(as.matrix(tempPS), dim = c(dim, 1, L)))
    MBBnull[iter] <- max(KSstatCalc(tempA[(washoutL + trainL + 1):L,], 0.01))
  }
  MBBquant <- sum(statistic <= MBBnull) / nboots
  output <- list("estimate" = estimate,
                 "statistic" = statistic,
                 "MBBquant" = MBBquant,
                 "MBBnull" = MBBnull,
                 "MBBblockL" = MBBblockL,
                 "statSeries" = KSseries,
                 "angles" = angles,
                 "netParams" = list(washoutL = washoutL, trainL = trainL))
  
  if(plot.it == TRUE) {
    CPplot <- plotCP(output)
    print(CPplot)
  }
  return(output)
}

NC <- function(data) {
  result <- NaiveConceptor(data$data, washoutL = 60, trainL = 120, plot.it = F)
  return(result)
}

nc1 <- mclapply(ID1, NC, mc.cores = numCores)
save(nc1, file = "nc1.RData")

nc2 <- mclapply(ID2, NC, mc.cores = numCores)
save(nc2, file = "nc2.RData")

nc3 <- mclapply(ID3, NC, mc.cores = numCores)
save(nc3, file = "nc3.RData")

nc4 <- mclapply(ID4, NC, mc.cores = numCores)
save(nc4, file = "nc4.RData")

nc5 <- mclapply(ID5, NC, mc.cores = numCores)
save(nc5, file = "nc5.RData")

##################################################
# simulation analysis, plot generation, etc.
##################################################

L_MAX = 1000
PAD = 180
neach <- 5

simdata <- c()
L <- c()
ID <- c()
Key <- c()

for(simID in 1:neach) {
  load(paste0("ID", simID,".RData"))
  simdata <- c(simdata, eval(parse(text = paste0("ID", simID))))
  L <- c(L, length(eval(parse(text = paste0("ID", simID)))))
  ID <- c(ID, rep(simID, L[simID]))
  Key <- c(Key, seq(L[simID]))
}
simresults <- dplyr::tibble(ID = ID, Key = Key, 
                            Class = c(rep("VAR(1)", 1800), rep("VAR(2)", 1800), rep("Periodic", 1500), rep("Ornstein-Uhlenbeck", 2700), rep("Gaussian", 3000)),
                            Details = c(rep("\u03c1 = 0.5 \u2192 0.5", 300), rep("\u03c1 = 0.5 \u2192 0.8", 300), rep("\u03c1 = 0.8 \u2192 0.5", 300), rep("\u03c1 = 0.8 \u2192 0.8", 300), rep("\u03c1 = 0.5 \u2192 NC", 300), rep("\u03c1 = 0.8 \u2192 NC", 300),
                                        rep("\u03c1 = 0.5 \u2192 0.5", 300), rep("\u03c1 = 0.5 \u2192 0.8", 300), rep("\u03c1 = 0.8 \u2192 0.5", 300), rep("\u03c1 = 0.8 \u2192 0.8", 300), rep("\u03c1 = 0.5 \u2192 NC", 300), rep("\u03c1 = 0.8 \u2192 NC", 300),
                                        rep("\u03c9 = 1 \u2192 0.5", 300), rep("\u03c9 = 1 \u2192 0.8", 300), rep("\u03c9 = 1 \u2192 1.2", 300), rep("\u03c9 = 1 \u2192 1.5", 300), rep("\u03c9 = 1 \u2192 NC", 300),
                                        rep("\u03b8 = 0.5 \u2192 0; \u03bb = 0.5", 300), rep("\u03b8 = 0.5 \u2192 1; \u03bb = 0.5", 300), rep("\u03b8 = 1 \u2192 0; \u03bb = 0.5", 300), rep("\u03b8 = 1 \u2192 0.5; \u03bb = 0.5", 300),
                                        rep("\u03b8 = 0.5; \u03bb = 0.5 \u2192 0.2", 300), rep("\u03b8 = 0.5; \u03bb = 0.5 \u2192 0.8", 300), rep("\u03b8 = 0.5; \u03bb = 0.5 \u2192 1", 300),
                                        rep("\u03b8 = 0.5, \u03bb = 0.5 \u2192 NC", 300), rep("\u03b8 = 1, \u03bb = 0.5 \u2192 NC", 300),
                                        rep("\u03bc = 0 \u2192 0.5", 300), rep("\u03bc = 0 \u2192 0.8", 300), rep("\u03bc = 0 \u2192 1.0", 300),
                                        rep("\u03c3 = 1 \u2192 0.5", 300), rep("\u03c3 = 1 \u2192 0.8", 300), rep("\u03c3 = 1 \u2192 1.2", 300), rep("\u03c3 = 1 \u2192 1.5", 300),
                                        rep("\u03c1 = 0 \u2192 -0.8", 300), rep("\u03c1 = 0 \u2192 0.8", 300), rep("\u03bc, \u03c1 = 0, \u03c3 = 1 \u2192 NC", 300)),
                            change_point = unlist(lapply(simdata, function(D) ifelse(!is.null(D$parameters$change_point), D$parameters$change_point, L_MAX))))

ggplot(simresults %>% filter(change_point < L_MAX)) + 
  geom_histogram(aes(x = change_point), binwidth = 20, center = 11, 
                 color = "gray50", fill = "lightblue", alpha = 0.7) +
  geom_hline(aes(yintercept = nrow(simresults[simresults$change_point < 1000,]) * 20 / (1000 - 181)), linetype = "dashed", color = "black") + 
  scale_x_continuous("Random Change Point", expand = c(0.02, 0.02)) +
  scale_y_continuous(paste0("# Simulations (out of ", nrow(simresults[simresults$change_point < 1000,]), ")"), expand = c(0.02, 0.02)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size  = 12))

eddata <- c()

for(simID in 1:neach) {
  load(paste("ed", simID, ".RData", sep = ""))
  eddata <- c(eddata, eval(parse(text = paste0("ed", simID))))
}
simresults <- simresults %>% mutate(ed_est = unlist(lapply(eddata, function(D) ifelse(!is.na(D$order.found[3]), D$order.found[3] + PAD, L_MAX))),
                                    ed_q = unlist(lapply(eddata, function(D) D$p.values[1])))

kcpdata <- c()

for(simID in 1:neach) {
  load(paste("kcp", simID, ".RData", sep = ""))
  kcpdata <- c(kcpdata, eval(parse(text = paste0("kcp", simID))))
}
simresults <- simresults %>% mutate(kcp_est = unlist(lapply(kcpdata, function(D) D[2] + PAD)))

sbsmeandata <- c()

for(simID in 1:neach) {
  load(paste("sbsmean", simID, ".RData", sep = ""))
  sbsmeandata <- c(sbsmeandata, eval(parse(text = paste0("sbsmean", simID))))
}
simresults <- simresults %>% mutate(sbs1_est01 = unlist(lapply(sbsmeandata, function(D) ifelse(D$result01$tree[[1]][3,1] != 0, D$result01$tree[[1]][3,1] + PAD, L_MAX))),
                                    sbs1_est05 = unlist(lapply(sbsmeandata, function(D) ifelse(D$result05$tree[[1]][3,1] != 0, D$result05$tree[[1]][3,1] + PAD, L_MAX))),
                                    sbs1_est10 = unlist(lapply(sbsmeandata, function(D) ifelse(D$result10$tree[[1]][3,1] != 0, D$result10$tree[[1]][3,1] + PAD, L_MAX))),
                                    sbs1_est20 = unlist(lapply(sbsmeandata, function(D) ifelse(D$result20$tree[[1]][3,1] != 0, D$result20$tree[[1]][3,1] + PAD, L_MAX))))

sbsdata <- c()

for(simID in 1:neach) {
  load(paste("sbs", simID, ".RData", sep = ""))
  sbsdata <- c(sbsdata, eval(parse(text = paste0("sbs", simID))))
}
simresults <- simresults %>% mutate(sbs2_est01 = unlist(lapply(sbsdata, function(D) ifelse(D$result01$tree[[1]][3,1] != 0, D$result01$tree[[1]][3,1] + PAD, L_MAX))),
                                    sbs2_est05 = unlist(lapply(sbsdata, function(D) ifelse(D$result05$tree[[1]][3,1] != 0, D$result05$tree[[1]][3,1] + PAD, L_MAX))),
                                    sbs2_est10 = unlist(lapply(sbsdata, function(D) ifelse(D$result10$tree[[1]][3,1] != 0, D$result10$tree[[1]][3,1] + PAD, L_MAX))),
                                    sbs2_est20 = unlist(lapply(sbsdata, function(D) ifelse(D$result20$tree[[1]][3,1] != 0, D$result20$tree[[1]][3,1] + PAD, L_MAX))))

nccpdata <- c()

for(simID in 1:neach) {
  load(paste("nc", simID, ".RData", sep = ""))
  nccpdata <- c(nccpdata, eval(parse(text = paste0("nc", simID))))
}
simresults <- simresults %>% mutate(nccp_est = unlist(lapply(nccpdata, function(D) D$estimate)),
                                    nccp_q = unlist(lapply(nccpdata, function(D) D$MBBquant)))

ccp02data <- c()

for(simID in 1:neach) {
  load(paste("ccp02_", simID, ".RData", sep = ""))
  ccp02data <- c(ccp02data, eval(parse(text = paste0("ccp02_", simID))))
}
simresults <- simresults %>% mutate(ccp02_est = unlist(lapply(ccp02data, function(D) D$estimate)),
                                    ccp02_q = unlist(lapply(ccp02data, function(D) D$MBBquant)))

ccp04data <- c()

for(simID in 1:neach) {
  load(paste("ccp04_", simID, ".RData", sep = ""))
  ccp04data <- c(ccp04data, eval(parse(text = paste0("ccp04_", simID))))
}
simresults <- simresults %>% mutate(ccp04_est = unlist(lapply(ccp04data, function(D) D$estimate)),
                                    ccp04_q = unlist(lapply(ccp04data, function(D) D$MBBquant)))

ccp08data <- c()
for(simID in 1:neach) {
  load(paste("ccp08_", simID, ".RData", sep = ""))
  ccp08data <- c(ccp08data, eval(parse(text = paste0("ccp08_", simID))))
}
simresults <- simresults %>% mutate(ccp08_est = unlist(lapply(ccp08data, function(D) D$estimate)),
                                    ccp08_q = unlist(lapply(ccp08data, function(D) D$MBBquant)))

ccp16data <- c()

for(simID in 1:neach) {
  load(paste("ccp16_", simID, ".RData", sep = ""))
  ccp16data <- c(ccp16data, eval(parse(text = paste0("ccp16_", simID))))
}
simresults <- simresults %>% mutate(ccp16_est = unlist(lapply(ccp16data, function(D) D$estimate)),
                                    ccp16_q = unlist(lapply(ccp16data, function(D) D$MBBquant)))

simresults_adj <- simresults %>% mutate(ed_est = ifelse(ed_q <= 0.05, ed_est, L_MAX),
                                        ccp02_est = ifelse(ccp02_q <= 0.05, ccp02_est, L_MAX),
                                        ccp04_est = ifelse(ccp04_q <= 0.05, ccp04_est, L_MAX),
                                        ccp08_est = ifelse(ccp08_q <= 0.05, ccp08_est, L_MAX),
                                        ccp16_est = ifelse(ccp16_q <= 0.05, ccp16_est, L_MAX)) %>% 
  select(ID, Key, Class, Details, change_point, ed_est, kcp_est, sbs1_est05, sbs2_est05, ccp02_est, ccp04_est, ccp08_est, ccp16_est, nccp_est)

ARI <- function(change_point, estimate) {
  tab <- matrix(c(min(change_point, estimate) - PAD, max(0, estimate - change_point), max(0, change_point - estimate), L_MAX - max(change_point, estimate)), nrow = 2)
  if(all(dim(tab) == c(1, 1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  out <- (a - (a + b) * (a + c) / (a + b + c + d)) / ((a + b + a + c) / 2 - (a + b) * (a + c) / (a + b + c + d))
  return(out)
}

ARItable <- simresults_adj %>% 
  filter(change_point < 1000) %>%
  rowwise() %>%
  mutate(across(6:14, function(x) ARI(change_point, x))) %>% 
  select(Class, Details, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>% #, nccp_est) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  group_by(Class, Details) %>%
  summarise(across(everything(), function(x) format(round(mean(x), 3), nsmall = 3)))

create_table <- function(class) {
  tab <- ARItable %>% filter(Class == class) %>% ungroup() %>% select(3:10)
  tabltx <- xtable(tab, type = "latex", align = "ccccccccc")
  return(tabltx)
}

xtable(ncpTable)

ncpTable <- simresults_adj %>% 
  filter(change_point == 1000) %>%
  select(Class, Details, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>% #, nccp_est) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  group_by(Class, Details) %>%
  rowwise() %>%
  mutate(across(everything(), function(x) x != 1000)) %>%
  group_by(Class, Details) %>%
  summarise(across(everything(), mean))
ncpTable

CDFs <- c()
for(diff in 0:105) {
  CDFs <- simresults_adj %>% 
    filter(change_point < 1000) %>%
    rowwise() %>%
    mutate(across(6:14, function(x) abs(x - change_point))) %>%
    select(Class, Details, ed_est, kcp_est, sbs1_est05, sbs2_est05, ccp02_est, ccp04_est, ccp08_est, ccp16_est, nccp_est) %>%
    group_by(Class, Details) %>%
    summarise(across(everything(), function(x) sum(x <= diff) / n())) %>%
    mutate(Diff = diff) %>%
    bind_rows(CDFs)
}

VARdata <- CDFs %>% filter(Class == "VAR(1)" | Class == "VAR(2)") %>% 
  select(Class, Details, Diff, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  pivot_longer(4:11, names_to = "Method", values_to = "CDF")

VARplot <- ggplot(VARdata) +
  facet_grid(vars(Class), vars(Details)) +
  geom_line(aes(x = Diff, y = CDF, color = Method, linetype = Method)) +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dashed", "solid")) +
  scale_x_continuous(name = "|Estimate - Actual|", limits = c(0, 105), expand = c(0.008, 0.008), breaks = seq(0, 100, 20)) +
  scale_y_continuous(expand = c(0.008, 0.008), breaks = seq(0.2, 1.0, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))

quartz(file = "VARplot.pdf", type = "pdf", width = 15, height = 7)
VARplot
dev.off()

PERdata <- CDFs %>% filter(Class == "Periodic") %>% 
  select(Class, Details, Diff, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  pivot_longer(4:11, names_to = "Method", values_to = "CDF")

PERplot <- ggplot(PERdata) +
  facet_wrap(vars(Details), nrow = 1) +
  geom_line(aes(x = Diff, y = CDF, color = Method, linetype = Method)) +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dashed", "solid")) +
  scale_x_continuous(name = "|Estimate - Actual|", limits = c(0, 105), expand = c(0.008, 0.008), breaks = seq(0, 100, 20)) +
  scale_y_continuous(expand = c(0.008, 0.008), breaks = seq(0, 1.0, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))

quartz(file = "PERplot.pdf", type = "pdf", width = 15, height = 4.3)
PERplot
dev.off()

OUdetails <- unique(filter(CDFs, Class == "Ornstein-Uhlenbeck")$Details)

OUdata1 <- CDFs %>% filter(Class == "Ornstein-Uhlenbeck") %>% 
  filter(Details %in% OUdetails[c(1, 2, 6, 7)]) %>%
  select(Class, Details, Diff, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  pivot_longer(4:11, names_to = "Method", values_to = "CDF")

OUplot1 <- ggplot(OUdata1) +
  facet_wrap(vars(Details), nrow = 1) +
  geom_line(aes(x = Diff, y = CDF, color = Method, linetype = Method)) +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dashed", "solid")) +
  scale_x_continuous(name = "|Estimate - Actual|", limits = c(0, 105), expand = c(0.008, 0.008), breaks = seq(0, 100, 20)) +
  scale_y_continuous(expand = c(0.008, 0.008), breaks = seq(0, 1.0, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))

quartz(file = "OUplot1.pdf", type = "pdf", width = 15, height = 4.3)
OUplot1
dev.off()

OUdata2 <- CDFs %>% filter(Class == "Ornstein-Uhlenbeck") %>% 
  filter(Details %in% OUdetails[c(3, 4, 5)]) %>%
  select(Class, Details, Diff, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  pivot_longer(4:11, names_to = "Method", values_to = "CDF")

OUplot2 <- ggplot(OUdata2) +
  facet_wrap(vars(Details), nrow = 1) +
  geom_line(aes(x = Diff, y = CDF, color = Method, linetype = Method)) +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dashed", "solid")) +
  scale_x_continuous(name = "|Estimate - Actual|", limits = c(0, 105), expand = c(0.008, 0.008), breaks = seq(0, 100, 20)) +
  scale_y_continuous(expand = c(0.008, 0.008), breaks = seq(0, 1.0, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))

quartz(file = "OUplot2.pdf", type = "pdf", width = 13, height = 4.3)
OUplot2
dev.off()

Gdetails <- unique(filter(CDFs, Class == "Gaussian")$Details)
Gdata1 <- CDFs %>% filter(Class == "Gaussian") %>% 
  filter(Details %in% Gdetails[1:3]) %>%
  select(Class, Details, Diff, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  pivot_longer(4:11, names_to = "Method", values_to = "CDF")

Gplot1 <- ggplot(Gdata1) +
  facet_wrap(vars(Details), nrow = 1) +
  geom_line(aes(x = Diff, y = CDF, color = Method, linetype = Method)) +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dashed", "solid")) +
  scale_x_continuous(name = "|Estimate - Actual|", limits = c(0, 105), expand = c(0.008, 0.008), breaks = seq(0, 100, 20)) +
  scale_y_continuous(expand = c(0.008, 0.008), breaks = seq(0, 1.0, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))

quartz(file = "Gplot1.pdf", type = "pdf", width = 13, height = 4.3)
Gplot1
dev.off()

Gdata2 <- CDFs %>% filter(Class == "Gaussian") %>% 
  filter(Details %in% Gdetails[4:5]) %>%
  select(Class, Details, Diff, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  pivot_longer(4:11, names_to = "Method", values_to = "CDF")

Gplot2 <- ggplot(Gdata2) +
  facet_wrap(vars(Details), nrow = 1) +
  geom_line(aes(x = Diff, y = CDF, color = Method, linetype = Method)) +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dashed", "solid")) +
  scale_x_continuous(name = "|Estimate - Actual|", limits = c(0, 105), expand = c(0.008, 0.008), breaks = seq(0, 100, 20)) +
  scale_y_continuous(expand = c(0.008, 0.008), breaks = seq(0, 1.0, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))

quartz(file = "Gplot2.pdf", type = "pdf", width = 9, height = 4.3)
Gplot2
dev.off()

Gdata3 <- CDFs %>% filter(Class == "Gaussian") %>% 
  filter(Details %in% Gdetails[6:9]) %>%
  select(Class, Details, Diff, ccp02_est, ccp04_est, ccp08_est, ccp16_est, ed_est, kcp_est, sbs1_est05, sbs2_est05) %>%
  rename(CCP02 = ccp02_est, CCP04 = ccp04_est, CCP08 = ccp08_est, CCP16 = ccp16_est, EDiv = ed_est, KCP = kcp_est, SBS1 = sbs1_est05, SBS2 = sbs2_est05) %>%
  pivot_longer(4:11, names_to = "Method", values_to = "CDF")

Gplot3 <- ggplot(Gdata3) +
  facet_wrap(vars(Details), nrow = 1) +
  geom_line(aes(x = Diff, y = CDF, color = Method, linetype = Method)) +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dashed", "solid")) +
  scale_x_continuous(name = "|Estimate - Actual|", limits = c(0, 105), expand = c(0.008, 0.008), breaks = seq(0, 100, 20)) +
  scale_y_continuous(expand = c(0.008, 0.008), breaks = seq(0, 1.0, 0.2), limits = c(0, 1)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))

quartz(file = "Gplot3.pdf", type = "pdf", width = 15, height = 4.3)
Gplot3
dev.off()

qdata1 <- simresults %>% 
  filter(change_point == 1000) %>%
  select(Class, Details, ccp02_q, ccp04_q, ccp08_q, ccp16_q, ed_q) %>%
  group_by(Class, Details) %>%
  mutate(across(everything(), sort)) %>%
  mutate(Reference = rep(seq(1 / n(), 1, 1 / n()), length(unique(Details)))) %>%
  rename(CCP02 = ccp02_q, CCP04 = ccp04_q, CCP08 = ccp08_q, CCP16 = ccp16_q, EDiv = ed_q) %>%
  pivot_longer(3:7, names_to = "Method", values_to = "T1E")

qdata2 <- simresults %>% 
  filter(change_point == 1000) %>%
  select(Class, Details, kcp_est, sbs1_est01, sbs1_est05, sbs1_est10, sbs1_est20, sbs2_est01, sbs2_est05, sbs2_est10, sbs2_est20) %>%
  group_by(Class, Details) %>%
  rowwise() %>%
  mutate(across(everything(), function(x) x != 1000)) %>%
  group_by(Class, Details) %>%
  summarise(across(everything(), mean)) %>%
  mutate(kcp2 = kcp_est) %>%
  pivot_longer(3:12, names_to = "Method", values_to = "Reference") %>%
  mutate(T1E = rep(c(0, 0.01, 0.05, 0.10, 0.20, 0.01, 0.05, 0.10, 0.20, 1), length(unique(Details)))) %>%
  mutate(Method = ifelse(str_detect(Method, "kcp"), "KCP", ifelse(str_detect(Method, "sbs1"), "SBS1", "SBS2")))

qdata <- bind_rows(qdata1, qdata2)

NCplot <- ggplot(qdata) +
  facet_wrap(vars(Class, Details), nrow = 2) +
  geom_line(aes(x = T1E, y = Reference, color = Method, linetype = Method)) +
  geom_point(aes(x = T1E, y = Reference, shape = Method, color = Method), size = 6) + 
  geom_abline(slope = 1, intercept = 0, color = "gray") +
  scale_color_manual("", values = c("blue", "blue", "blue", "blue", "firebrick2", "green3", "darkgoldenrod1", "darkgoldenrod1")) +
  scale_linetype_manual("", values = c("solid", "longdash", "dashed", "dotted", "solid", "solid", "dotted", "solid")) +
  scale_shape_manual("", values = c("", "", "", "", "", "", "*", "*")) + 
  scale_x_continuous(name = "Target Type 1 Error", expand = c(0.008, 0.008), breaks = seq(0.1, 0.9, 0.2)) +
  scale_y_continuous(name = "Observed Type 1 Error", expand = c(0.008, 0.008), breaks =seq(0, 1, 0.2)) + 
  theme_bw() +
  theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 20), 
        legend.text = element_text(size = 20), legend.position = "bottom",
        axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        panel.grid = element_blank(), strip.background = element_blank(), 
        strip.text.x = element_text(size = 20), 
        strip.text.y = element_text(size = 20, angle = -90))


quartz(file = "NCplot.pdf", type = "pdf", width = 15, height = 9)
NCplot
dev.off()
