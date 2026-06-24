

# function to load required packages 

load_libraries <- function(pack) {
  for (p in pack) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
    }
    library(p, character.only = TRUE)
  }
}


#initialise individuals in patches 
create_n_per_patch <- function(patches, carrying_capacity) {
  
  if (patches < 1) {
    stop("Number of patches must be at least 1.")
  }
  if (carrying_capacity < 0) {
    stop("Carrying capacity must not be < 0")
  }
  n_per_patch <- rep(0, patches)
  n_per_patch[1] <- carrying_capacity
  return(n_per_patch)
}





# # function to make autosome
# 
# make_autosome <- function(individual, prefix, n_loci) {
#  
#    matrix(
#   c(paste0(rep(prefix, individual * (n_loci)), rep(1:(n_loci), each = individual))),
#   nrow = individual, ncol = n_loci
# )
# 
# }
# 
# #make_autosome(5, "A", 5)



# function to make allosome (sex chromosome)

# make_allosome <- function(sex_alleles, individual) {
#   
#   matrix(
#     c(sample(sex_alleles, size = individual, replace = TRUE)),
#   )
#   
# }


make_allosome <- function(sex_alleles, individual) {
  
  cbind(
    sample(sex_alleles, size = individual, replace = TRUE),
    rep(0L, individual)
  )
}
  

# Loci selection matrix: function to place loci at random on the genome (of size = 1)
# also takes exponential decay and variance to produce variance-covariance matrix

place_loci_mat <- function(loci, genome.size = 1, var = 1, decay){
  loci_positions <- sort((runif(loci, max = genome.size)))
  loci_dist_matrix <- as.matrix(dist(loci_positions))^2 
  loci_cov_matrix <- var*exp(-decay*loci_dist_matrix)
  return(loci_cov_matrix)
}


# growth degree day estimation (Abbasi et al., Environmental Entomology, 2023, Vol. 52, No. 6)

# function to calculate growth degree-days accumulation 

cal_dd <- function(daily_max_temp, daily_min_temp, T_base) {
  max(0, (daily_max_temp+daily_min_temp)/2 - T_base)
}


# function to estimate the probability of transition based on mean growth degree 
# days at different percentile using inverse cumulative distribution function* 

erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)

sigma_etimate <- function(x, mu, p) {
  (x - mu)/(sqrt(2) * erfinv(2*p-1))
}

prob_trans <- function(dd, mu, sigma) {
  pnorm(dd, mu, sigma)
}


# estimated survival based on temperature and population density (aquatic stage),  
# and temperature and humidity (adult stage) using daily mortality hazard, and 
# developmental rate (see Golding et al., unpublished data)

ensure_positive <- function(x) {
  x * as.numeric(x > 0)
}

# reload lifehistory functions from saved objects (RDS file) to used for survival. 
# Adapted from Golding et al., unpublished)

rehydrate_lifehistory_function <- function(path_to_object) {
  object <- readRDS(path_to_object)
  do.call(`function`,
          list(object$arguments,
               body(object$dummy_function)))
}


aquatic_stage <- "R/das_temp_dens_As.RDS"
adult_stage <- "R/ds_temp_humid.RDS"

das_temp_dens_As <- rehydrate_lifehistory_function(aquatic_stage)
ds_temp_humid_As <- rehydrate_lifehistory_function(adult_stage)


# function to simulate oviposition frequency and batch sizes. mean eggs per female
# per day (EFD) and mean temperature were estimated from Villena et al., https://doi.org/10.1002/ecy.3685

# return the parameters of lognormal with specified mean and variance
lognormal_params <- function(mean, sd) {
  var <- sd ^ 2
  list(
    meanlog = log((mean ^ 2) / sqrt(var + mean ^ 2)),
    sdlog = sqrt(log(1 + var / (mean ^ 2)))
  )
}

# simulate from a lognormal, given the mean and sd of the distribution (not the
# meanlog and sdlog parameters)
rlnorm_mean_var <- function(n, mean, sd) {
  params <- lognormal_params(mean, sd)
  rlnorm(n, params$meanlog, params$sdlog)
}

# simulate delays between egg batches, in days


# We simulate the expected batch sizes based on Suleman, 1990 https://doi.org/10.1093/jmedent/27.5.819
# to match the mean and SD but modelled as negative binomial
sim_batch_sizes <- function(n) {
  rnbinom(n, mu = 96.8, size = 1 / 0.16)
}

# calculate the frequency of oviposition or expected delay between batches,
# given the expected batch size
# 
egg_laying_rate <- function(temp) {
peak_temp <- 28
peak_val <- 26.2
temp_sd <- 6

unscaled_value <- dnorm(temp,
                        mean = peak_temp,
                        sd = temp_sd)
normalisation <- dnorm(peak_temp,
                       mean = peak_temp,
                       sd = temp_sd)
peak_val * unscaled_value / normalisation

}

expected_egg_laying_delay <- function(temp, expected_batch_size = 96.8) {
  expected_batch_size / egg_laying_rate(temp)
}

sim_delays <- function(n, temp) {
  expected_delay <- expected_egg_laying_delay(temp)
  delays_continuous <- rlnorm_mean_var(n,
                                       expected_delay,
                                       sd = 0.5)
  delays <- pmax(1, round(delays_continuous))
  delays
}



#### negative exponential dispersal kernel  

metapop <- function(coords, lambda, disp_prob) {
  # dispersal matrix 
  dist_matrix <- as.matrix(dist(coords, method = "euclidean"))
  
  #exponential dispersal kernel
  dispersal_kernel <- exp(-lambda * dist_matrix)
  
  # set the diagonal elements to 0 to prevent self-dispersal
  diag(dispersal_kernel) <- 0
  
  
  # make these rows sum to 1 to get probability of moving to other patch
  # *if* they left. This dispersal matrix gives the probability of the vector
  # vector moving between patches
  rel_dispersal_matrix <- sweep(dispersal_kernel, 1,
                                rowSums(dispersal_kernel), FUN = "/")
  
  # normalise these to have the overall probability of dispersing to that patch,
  # and add back the probability of remaining
  dispersal_matrix <- disp_prob * rel_dispersal_matrix +
    (1 - disp_prob) * diag(nrow(dispersal_kernel))
  
  return(dispersal_matrix)
}


# adjacency matrix 

step_stone <- function(n_patches, disp_prob) {
  
  matrix_landscape <- matrix(0, n_patches, n_patches)
  adjacency <- abs(row(matrix_landscape) - col(matrix_landscape)) == 1
  adjacency[] <- as.numeric(adjacency)

  # make these rows sum to 1 to get probability of moving to other patch
  # *if* they left. This dispersal matrix gives the probability of the vector
  # vector moving between patches
  rel_dispersal_matrix <- sweep(adjacency, 1,
                                rowSums(adjacency), FUN = "/")

  # normalise these to have the overall probability of dispersing to that patch,
  # and add back the probability of remaining
  dispersal_matrix <- disp_prob * rel_dispersal_matrix +
    (1 - disp_prob) * diag(nrow(adjacency))

  return(dispersal_matrix)
}






#Homing gene drive function (conversion mechanism)

home_drive_conv <- function(parent, prob1, prob2) {
  # browser()
  
  loci1 <- parent$autosome1 
  loci2 <- parent$autosome2
  
  # if (any(is.na(loci1)) | any(loci2)) {
  #     warning("NA detected in allele input!")
  # }
  
  drive_wt <- ((loci1 == 0) & (loci2 == 1)) | ((loci1 == 1) & (loci2 == 0))
  
  #cleavage
  cleavage  <- matrix(rbinom(nrow(loci1), 1, prob1), ncol(loci1), # drive cleavage at each locus
                      nrow = nrow(loci1), ncol = ncol(loci1))
  # homing
  homing  <- matrix(rbinom(nrow(loci1), 1, prob2), ncol(loci1), # drive conversion at each locus
                    nrow = nrow(loci1), ncol = ncol(loci1)) 
  
  conv_event <- homing*cleavage # conversion event?
  conv_heterozygous <- conv_event*drive_wt  #This is where the trick is....
  
  #successful homing (0 to 1 )
  loci1[loci1 == 0 & conv_event == 1 & conv_heterozygous == 1] <- 1 # successful conversions
  loci2[loci2 == 0 & conv_event == 1 & conv_heterozygous == 1] <- 1
  
  # Resistance development if homing fails (0 to 2)
  # failed_conv <- drive_wt & conv_event == 0
  # resistance_event <- rbinom(length(parent$chromosome1), 1, prob2)
  loci1[loci1 == 0 & cleavage == 1 & conv_event == 0 & conv_heterozygous == 0] <- 2  #Thoughts/To do: individuals that did not develop resistance, yet heterozygous can be can be designated as those with functional resistance and resistant to future Cas9 cutting
  loci2[loci2 == 0 & cleavage == 1 & conv_event == 0 & conv_heterozygous == 0] <- 2
  
  parent$autosome1  <- loci1
  parent$autosome2  <- loci2
  
  return(parent)
}


