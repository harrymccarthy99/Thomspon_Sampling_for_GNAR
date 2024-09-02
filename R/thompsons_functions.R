
# Helper function for plotting

plot_from_matrix <- function(mat, title){
  g <- graph.adjacency(mat, mode = "undirected", weighted = TRUE)
  
  # Setting node labels
  V(g)$label <- rownames(mat)
  
  # Plot network, adjust parameters as needed
  plot(g, layout = layout_with_fr,
       vertex.size = 10,          
       vertex.label.cex = 0.8,    
       vertex.label.dist = 0.8,    
       vertex.color = "lightblue",
       vertex.frame.color = "gray",
       edge.width = 2,      
       edge.color = "black",  
       main = title)
}


# Function to get upper triangular indices of an n x n matrix
get_upper_triangular_indices <- function(n) {
  # Empty n x n matrix
  mat <- matrix(0, n, n)
  
  # Get the upper triangular indices
  indices <- which(upper.tri(mat, diag = FALSE), arr.ind = TRUE)
  
  # Convert to a list of pairs
  indices_list <- split(indices, row(indices))
  
  return(indices_list)
}

# Finds Row/Col where indices equal 1
find_upper_triangular_indices <- function(adj_matrix) {
  # Get indices where adjacency matrix equals 1 i.e. an edge exists
  indices <- which(adj_matrix == 1, arr.ind = TRUE)
  
  # Keep only upper triangular elements indices 
  upper_triangular_indices <- indices[indices[,1] < indices[,2], ]
  
  return(upper_triangular_indices)
}

# Function to reconstruct a matrix from upper triangular indices
reconstruct_matrix <- function(n, indices_mat) {
  # Placeholder n x n martrix
  mat <- matrix(0, n, n)
  
  # Use matrix indexing to fill in the upper triangular part
  mat[cbind(indices_mat[,1], indices_mat[,2])] <- rep(1, nrow(indices_mat))
  
  # Mirror the upper triangular part to the lower triangular part to keep the matrix symmetric
  mat <- mat + t(mat)
  
  return(mat)
}

# Function to plot means std devs and arm counts for each selected arm
plot_histograms <- function(mu_history, sigma_history, chosen_arms_count) {
  # Create a sequence of iteration indices for every hundredth iteration
  iteration_indices <- seq(10, length(mu_history), by = 10)
  
  for (i in iteration_indices) {
    mu_df <- data.frame(arm = 1:length(mu_history[[i]]), mu = mu_history[[i]])
    sigma_df <- data.frame(arm = 1:length(sigma_history[[i]]), sigma = sigma_history[[i]])
    count_df <- data.frame(arm = 1:ncol(chosen_arms_count), count = chosen_arms_count[i, ])
    
    p1 <- ggplot(mu_df, aes(x = arm, y = mu)) +
      geom_bar(stat = "identity", fill = "blue") +
      ggtitle(paste("Posterior Means at Iteration", i*5))
    
    p2 <- ggplot(sigma_df, aes(x = arm, y = sigma)) +
      geom_bar(stat = "identity", fill = "red") +
      ggtitle(paste("Posterior Std Devs at Iteration", i*5))
    
    p3 <- ggplot(count_df, aes(x = arm, y = count)) +
      geom_bar(stat = "identity", fill = "green") +
      ggtitle(paste("Chosen Arms Count at Iteration", i*5))
    
    grid.arrange(p1, p2, p3, nrow = 3)
  }
}


# Function to convert adjacency matrix to edge list
adj_matrix_to_edge_list <- function(adj_matrix) {
  edges <- which(adj_matrix != 0, arr.ind = TRUE)
  # Ensure each edge is represented as (min, max) to handle undirected edges correctly
  edges <- t(apply(edges, 1, function(x) sort(x)))
  # Remove duplicate edges (since undirected)
  edges <- unique(edges)

  return(edges)
}

# Function to find unique edges
find_unique_edges <- function(edges1, edges2) {
  edges1_str <- apply(edges1, 1, paste, collapse = "-")
  edges2_str <- apply(edges2, 1, paste, collapse = "-")
  unique_edges <- setdiff(edges2_str, edges1_str)
  return(unique_edges)
}


# Function to save results
save_results <- function(mu_history, sigma_history, counts) {
  if (!dir.exists(models_dir)) {
    dir.create(models_dir)
  }
  save(mu_history, file = mu_file)
  save(sigma_history, file = sigma_file)
  save(counts, file = counts_file)
}

load_results <- function() {
  # Check if files exist and load them if they do
  if (file.exists(mu_file) && file.exists(sigma_file) && file.exists(counts_file)) {
    load(mu_file)
    load(sigma_file)
    load(counts_file)
    return(list(mu_history = mu_history, sigma_history = sigma_history, counts = chosen_arms_count))
  } else {
    return(list(mu_history = NULL, sigma_history = NULL, counts = NULL))
  }
}

