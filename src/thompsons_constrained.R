
# Loading packages
library(dplyr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(tidyverse)
library(geodist)
library(countrycode)
library(geodist)
library(dplyr)
library(maps)
library(ggmap)
library(igraph)
library(geosphere)
library(GNAR)
library(MCMCpack)
library(glue)
library(here)
library(gridExtra)


# Read function definitions
source(here("R", "thompsons_functions.R"))

# Load data and networks
unemployment_ts <- read.csv(here("data", "unemployment_ts.csv"))
unemployment_vts <- ts(subset(unemployment_ts, select = -Month_Year))
load(here('data', 'MST_Distance_Networks.RData'))


# Finding the number of unique edges in the proximity inferred matrix
MST_37_edges <- adj_matrix_to_edge_list(MST_37)
distance_matrix_600_edges <- adj_matrix_to_edge_list(distance_matrix_600)

# Find unique edges in distance_matrix_600 that are not in MST_37
unique_edges <- find_unique_edges(MST_37_edges, distance_matrix_600_edges)

# Count the number of unique edges
num_unique_edges <- length(unique_edges)
print(glue('Number of additional edges versus MST; {num_unique_edges}'))

# Paarameters
n_iterations <- 1000
n_nodes <- 37
n_arms <- n_nodes*(n_nodes -1)*0.5 # Number of arms in upper triangle
upper_tri_indices <- get_upper_triangular_indices(n_nodes) # Finds the indices of a num_nodes x num_nodes sized matrix
mu_0 <- 0   # prior mean
sigma <- 1 # prior variance - Should this be uninformative? Where do we actually use this?
kappa_0 <- 1
alpha_0 <- 1
beta_0 <- 1
init_indices <- find_upper_triangular_indices(MST_37) # Initial MST_37 Indices
init_indices_flat <- which(MST_37[upper.tri(MST_37)] == 1) # Flattening the initial MST_37 indices

# Define file paths
models_dir <- here("models")

mu_file <- file.path(models_dir, "la_1_3_mu_history.RData")
sigma_file <- file.path(models_dir, "la_1_3_sigma_history.RData")
counts_file <- file.path(models_dir, "la_1_3_counts.RData")

# Thompsons Sampling Function - Negative Inverse Gamma Prior
thompson_sampling_constrained <- function(n_nodes, mu_0, kappa_0, alpha_0, beta_0, init_indices, n_arms, n_iterations, num_unique_edges) {
  # Initialise parameters for each arm
  mu <- rep(mu_0, n_arms)
  kappa <- rep(kappa_0, n_arms)
  alpha <- rep(alpha_0, n_arms)
  beta <- rep(beta_0, n_arms)
  n <- rep(0, n_arms)
  x_bar <- rep(0, n_arms)
  S <- rep(0, n_arms)
  chosen_arms <- numeric(n_iterations * num_unique_edges)  # Vector to store chosen arms

  # Containers for plotting later
  mu_history <- list()
  sigma_history <- list()
  chosen_arms_count <- matrix(0, ncol = n_arms, nrow = ceiling(n_iterations))

  # Initialise the network as MST_37
  last_net_indices <- init_indices

  start_time <- proc.time()
  timing <- system.time({

    for (t in 1:n_iterations) {
      chosen_this_iter <- numeric(0)

      # Initialize the network with the initial MST
      last_net_indices <- init_indices

      for (edge in 1:num_unique_edges) {
        sampled_means <- numeric(n_arms)

        for (i in 1:n_arms) {
          sampled_sigma2 <- rinvgamma(1, alpha[i], beta[i])  # Sample a variance
          sampled_means[i] <- rnorm(1, mu[i], sqrt(sampled_sigma2 / kappa[i]))
        }

        available_arms <- setdiff(1:n_arms, chosen_this_iter)
        sampled_means <- sampled_means[available_arms]
        chosen_arm_idx <- which.min(abs(sampled_means))
        chosen_arm <- available_arms[chosen_arm_idx]
        chosen_arms[(t - 1) * num_unique_edges + edge] <- chosen_arm
        chosen_this_iter <- c(chosen_this_iter, chosen_arm)

        # Select the arm and update the last generated network with the selected indices
        last_net_indices <- rbind(upper_tri_indices[[chosen_arm]], last_net_indices)

        # Reconstruct the network from the updated indices
        net <- as.GNARnet(reconstruct_matrix(n_nodes, last_net_indices))
        fit <- GNARfit(vts = unemployment_vts[1:275,], net = net, alphaOrder = 1, betaOrder = c(3), globalalpha = FALSE)

        # Find error (GNAR prediction accuracy)
        error <- sum(unemployment_vts[276,] - predict(fit))^2

        # Update the posterior parameters for the chosen arm
        n_chosen <- n[chosen_arm]
        n[chosen_arm] <- n_chosen + 1
        x_bar[chosen_arm] <- (x_bar[chosen_arm] * n_chosen + error) / n[chosen_arm]
        S[chosen_arm] <- S[chosen_arm] + (error - x_bar[chosen_arm])^2 * n_chosen / (n_chosen + 1)
        kappa[chosen_arm] <- kappa[chosen_arm] + 1
        mu[chosen_arm] <- (kappa[chosen_arm] * mu[chosen_arm] + error) / (kappa[chosen_arm] + 1)
        alpha[chosen_arm] <- alpha[chosen_arm] + 0.5
        beta[chosen_arm] <- beta[chosen_arm] + 0.5 * ((error - mu[chosen_arm])^2 / (kappa[chosen_arm] + 1) + S[chosen_arm])
      }
      
      # Store mu, sigma, and counts for every iteration
      mu_history[[t]] <- mu 
      sigma_history[[t]] <- sqrt(beta / (alpha - 1))
      chosen_arms_count[t, ] <- table(factor(chosen_arms[1:(t * num_unique_edges)], levels = 1:n_arms))
      

      if (t %% 5 == 0) {

        elapsed_time <- proc.time() - start_time  # Calculate elapsed time
        print(paste("Iteration:", t, "- Unique arms sampled:", length(unique(chosen_arms[1:(t * num_unique_edges)])),
                    "Time Elapsed:", elapsed_time[3], "seconds"))
      }
    }})

  # Print the elapsed time
  print(paste("Elapsed time:", timing["elapsed"], "seconds"))
  return(list(mu = mu, mu_history = mu_history, sigma_history = sigma_history, chosen_arms_count = chosen_arms_count))
}


# Check if results are already saved
results <- load_results()

if (is.null(results$mu_history) || is.null(results$sigma_history) || is.null(results$counts)) {
  # Run the Thompson Sampling Algorithm if results are not found
  thompson_constrained_results <- thompson_sampling_constrained(
    n_nodes = n_nodes,
    mu_0 = mu_0,
    kappa_0 = kappa_0,
    alpha_0 = alpha_0,
    beta_0 = beta_0,
    init_indices = init_indices,
    n_arms = n_arms,
    n_iterations = n_iterations,
    num_unique_edges = num_unique_edges
  )
  
  # Save the results
  save_results(
    mu_history = thompson_constrained_results$mu_history,
    sigma_history = thompson_constrained_results$sigma_history,
    counts = thompson_constrained_results$chosen_arms_count
  )
} else {
  # Use the loaded results
  thompson_constrained_results <- list(
    mu_history = results$mu_history,
    sigma_history = results$sigma_history,
    chosen_arms_count = results$counts
  )
}


