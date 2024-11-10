hcmv <- read.table("hcmv.txt", header = TRUE)

# Question 1 - Random Scatter ---------------------------------------------

# Parameters
sequence_length <- 229354
num_palindromes <- 296
interval_size <- 1000

# Histogram of Real Data Palindromic Site Location
hist(hcmv$location, breaks = 25, probability = TRUE, col = rgb(1, 0, 0, 0.25),
     main = "Palindromic Sites of Real Data",
     xlab = "Location in DNA Sequence")
lines(density(hcmv$location, adjust = 1), col = 2)

# Simulation of Uniform Distribution of Palindromic Sites
set.seed(123)
simulated_palindromes <- runif(num_palindromes, min = 1, max = sequence_length)

# Histogram of Simulated Data Palindromic Site Location
hist(simulated_palindromes, breaks = 25, probability = TRUE, col = rgb(0,0,1,0.25), main = "Palindromic Sites of Simulated Data",
     xlab = "Location in DNA Sequence")
lines(density(simulated_palindromes, adjust = 1), col = 4)

simulated_palindromes_repeated <- replicate(4, sort(sample(1:sequence_length, num_palindromes, replace = TRUE)))

# Plots a few random simulations for qualitative comparison
hist(simulated_palindromes_repeated[,1], breaks = 100, probability = TRUE,
     main = "Simulated Random Scatter of Palindromic Sites",
     xlab = "Location in DNA Sequence")
lines(density(simulated_palindromes_repeated[,1], adjust = 1))
hist(simulated_palindromes_repeated[,2], breaks = 100, probability = TRUE,
     main = "Simulated Random Scatter of Palindromic Sites",
     xlab = "Location in DNA Sequence")
lines(density(simulated_palindromes_repeated[,2], adjust = 1))
hist(simulated_palindromes_repeated[,3], breaks = 100, probability = TRUE,
     main = "Simulated Random Scatter of Palindromic Sites",
     xlab = "Location in DNA Sequence")
lines(density(simulated_palindromes_repeated[,3], adjust = 1))
hist(simulated_palindromes_repeated[,4], breaks = 100, probability = TRUE,
     main = "Simulated Random Scatter of Palindromic Sites",
     xlab = "Location in DNA Sequence")
lines(density(simulated_palindromes_repeated[,4], adjust = 1))

# Layered Histogram for Graphical Comparison
hist(hcmv$location, breaks = 25, probability = TRUE, col = rgb(1, 0, 0, 0.25),
     main = "Palindromic Sites of Real & Simulated Data",
     xlab = "Location in DNA Sequence")
lines(density(hcmv$location, adjust = 1), col = 2)
hist(simulated_palindromes, breaks = 25, probability = TRUE, col = rgb(0,0,1,0.25), add = TRUE)
lines(density(simulated_palindromes, adjust = 1), col = 4)
legend("topright", legend = c("Real", "Simulated"), fill = c(rgb(1, 0, 0, 0.25), rgb(0,0,1,0.25)), border = NA)

# Grouping palindrome location in bins of 1000
tab <- table(cut(hcmv$location, 
                     breaks = seq(from = 1, to = 230000, by = 1000), 
                     include.lowest = TRUE))
counts <- as.vector(tab)
print(table(counts))

sim_tab <- table(cut(simulated_palindromes, 
                 breaks = seq(from = 1, to = 230000, by = 1000), 
                 include.lowest = TRUE))
sim_counts <- as.vector(sim_tab)
print(table(sim_counts))

# Histogram of Counts in Each Interval
hist(sim_counts, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), probability = TRUE, col = rgb(0, 0, 1, 0.25), main = "Real & Simulated Interval Counts", xlab = "number of points inside an interval", ylim = c(0,0.5))
lines(density(sim_counts, adjust = 2), col = rgb(0,0,1,0.5))

hist(counts, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), probability = TRUE, col = rgb(1, 0, 0, 0.25), add = TRUE)
lines(density(counts, adjust = 2), col = rgb(1,0,0,0.5))

legend("topright", legend = c("Real", "Simulated"), fill = c(rgb(1, 0, 0, 0.25), rgb(0,0,1,0.25)), border = NA)

# Calculate spacing for each simulation

# Real data
real_spacing <- diff(hcmv$location)
real_counts <- table(cut(hcmv$location, breaks = seq(0, sequence_length, by = interval_size), include.lowest = TRUE))

# Simulated data
simulated_spacings <- diff(sort(simulated_palindromes))
simulated_counts <- table(cut(simulated_palindromes, breaks = seq(0, sequence_length, by = interval_size), include.lowest = TRUE))

# Spacing histograms
hist(real_spacing, breaks = 30, col = rgb(1, 0,0,0.25), probability = TRUE, main = "Real & Simulated Data Spacing", xlab = "Spacing")
hist(unlist(simulated_spacings), breaks = 20, probability = TRUE, col = rgb(0, 0, 1, 0.25), add = TRUE)
lines(density(real_spacing, adjust = 1), col = rgb(1,0,0,0.5))
lines(density(unlist(simulated_spacings), adjust = 1), col = rgb(0,0,1,0.5))
legend("topright", legend = c("Real Data", "Simulated Data"),
       fill = c(rgb(1, 0, 0, 0.25), rgb(0, 0, 1, 0.25)), border = NA)

# Quantitative comparison
mean_real_spacing <- mean(real_spacing)
mean_simulated_spacing <- mean(unlist(simulated_spacings))
mean(diff(simulated_palindromes_repeated[,1]))
mean(diff(simulated_palindromes_repeated[,2]))
mean(diff(simulated_palindromes_repeated[,3]))
mean(diff(simulated_palindromes_repeated[,4]))

var_real_spacing <- var(real_spacing)
var_simulated_spacing <- var(unlist(simulated_spacings))

std_real_spacing <- sqrt(var_real_spacing)
std_simulated_spacing <- sqrt(var_simulated_spacing)

mean_real_counts <- mean(real_counts)
mean_simulated_counts <- mean(unlist(simulated_counts))

var_real_counts <- var(real_counts)
var_simulated_counts <- var(unlist(simulated_counts))

std_real_counts <- sqrt(var_real_counts)
std_simulated_counts <- sqrt(var_simulated_counts)

# Question 2 - Locations and Spacings -------------------------------------

hist(real_spacing, breaks = 30, probability = TRUE, main = "Real Data Spacing", xlab = "Spacing")


# intervals of 1000
tab <- table(cut(hcmv$location, 
                 breaks = seq(from = 1, to = 230000, by = 1000), 
                 include.lowest = TRUE))
counts <- as.vector(tab)
table(counts)

h <- hist(counts, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), ylim = c(0, 90))
text(h$mids,h$counts,labels=h$counts, adj=c(0.5, -0.5))

# Function for analyzing palindromes
analyze_palindromes <- function(region_length) {
  # Total number of non-overlapping regions
  num_regions <- ceiling(sequence_length / region_length)
  
  # Define region breaks and count palindromes in each interval
  region_breaks <- seq(1, sequence_length, by = region_length)
  
  # Count occurrences in each interval
  region_indices <- cut(hcmv$location, breaks = region_breaks, include.lowest = TRUE, labels = FALSE)
  counts <- table(region_indices)
  
  # Initialize counts vector with zeros for all regions
  counts_vector <- rep(0, num_regions)
  
  # Fill counts_vector based on the numeric indices
  counts_vector[as.numeric(names(counts))] <- as.vector(counts)
  
  # Expected count per region under uniform distribution
  expected_count <- length(hcmv$location) / num_regions
  
  hist(counts_vector, 
       breaks = seq(-0.5, max(counts_vector) + 0.5, by = 1),
       col = "blue", 
       border = "black", 
       main = paste("Palindrome Count Distribution for Region Length =", region_length),
       xlab = "Number of Palindromes per Region", 
       ylab = "Frequency")
  abline(v = expected_count, col = "red", lty = "dashed")
  
  # Classify regions by number of palindromes
  classified_counts <- table(cut(counts_vector, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, Inf), include.lowest = TRUE))
  print("Classified Counts by Number of Palindromes:")
  print(classified_counts)
  
  # Chi-square test comparing observed to expected uniform distribution
  chi_square_test <- chisq.test(counts_vector, p = rep(1 / num_regions, num_regions))
  print(paste("Chi-square test p-value for region length", region_length, ":", chi_square_test$p.value))
  
  return(counts_vector)
}

# Analyze for different region lengths
region_lengths <- c(500, 1000, 5000, 10000, 20000, 50000)
results <- lapply(region_lengths, analyze_palindromes)


# Question 4 - The Biggest Cluster ----------------------------------------

# Function to identify high palindrome density regions at different interval lengths
analyze_palindromic_clusters <- function(region_length) {
  # Define the region breaks
  region_breaks <- seq(1, sequence_length, by = region_length)
  
  # Initialize counts for each interval to zero
  num_intervals <- length(region_breaks) - 1
  counts_vector <- rep(0, num_intervals)
  
  # Count palindromes in each interval, using labels = FALSE to get numeric indices
  region_indices <- cut(hcmv$location, breaks = region_breaks, include.lowest = TRUE, labels = FALSE)
  counts <- table(region_indices)
  
  # Populate counts_vector with actual counts based on region indices
  counts_vector[as.numeric(names(counts))] <- as.vector(counts)
  
  # Find the interval with the highest palindrome count
  max_count <- max(counts_vector)
  max_interval <- which(counts_vector == max_count)
  
  # Expected count per interval under uniform distribution
  expected_count <- length(hcmv$location) / num_intervals
  
  # Print the interval and count information
  print(paste("Region length:", region_length))
  print(paste("Max palindrome count:", max_count))
  print(paste("Expected count:", expected_count))
  print(paste("Interval with max count:", max_interval))
  
  # Statistical test to check if max count significantly deviates from expectation
  chi_square_test <- chisq.test(counts_vector, p = rep(1 / num_intervals, num_intervals), simulate.p.value = TRUE)
  print(paste("Chi-square p-value for region length", region_length, ":", chi_square_test$p.value))
  
  return(list(max_count = max_count, max_interval = max_interval, p_value = chi_square_test$p.value))
}

# Analyze at multiple region lengths
region_lengths <- c(500, 1000, 2000, 5000)
results <- lapply(region_lengths, analyze_palindromic_clusters)

# Advanced Analysis -------------------------------------------------------

interval_sizes <- c(500, 1000, 5000, 10000)

# Chi-square and Poisson test functions
perform_chisq_test <- function(observed_counts, expected_counts) {
  chisq_test <- chisq.test(observed_counts, p = expected_counts / sum(expected_counts), simulate.p.value = TRUE)
  return(chisq_test$p.value)
}

p_values_chi <- numeric(length(interval_sizes))
p_values_poisson <- numeric(length(interval_sizes))
results_summary <- data.frame()

# Loop over different interval sizes
for (interval_size in interval_sizes) {
  
  # Real data: group counts by interval
  real_counts <- table(cut(hcmv$location, 
                           breaks = seq(from = 1, to = sequence_length, by = interval_size), 
                           include.lowest = TRUE))

  # Simulated data: group counts by interval
  simulated_counts <- table(cut(simulated_palindromes, 
                                breaks = seq(1, sequence_length, by = interval_size), 
                                include.lowest = TRUE))
  
  # output results
  print(interval_size)
  print(table(real_counts))
  print(table(simulated_counts))
  
  # Chi-square test on real vs expected under uniform distribution
  expected_count <- rep(num_palindromes * interval_size / sequence_length, length(real_counts))
  chi_p_value <- perform_chisq_test(as.vector(real_counts), expected_count)
  
  # Poisson test (mean assumed as expected count - Poisson)
  observed_counts <- as.numeric(real_counts)
  lambda <- mean(observed_counts)
  poisson_p_value <- poisson.test(sum(observed_counts), T = length(real_counts), r = lambda)$p.value
  
  # results
  results_summary <- rbind(results_summary, 
                           data.frame(interval_size = interval_size,
                                      mean_real_counts = mean(observed_counts),
                                      mean_sim_counts = mean(simulated_counts),
                                      var_real_counts = var(observed_counts),
                                      var_sim_counts = var(simulated_counts),
                                      chi_p_value = chi_p_value,
                                      expected_count = expected_count[1],
                                      poisson_p_value = poisson_p_value))
  
  # Histogram plot for counts
  sim_density <- density(simulated_counts)
  obs_density <- density(observed_counts)
  y_max <- max(sim_density$y, obs_density$y) * 2
  x_max <- max(sim_density$x, obs_density$x) * 1.1
  
  if (interval_size == 500) {
    hist(simulated_counts, 
         breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), 
         probability = TRUE, 
         col = rgb(0, 0, 1, 0.25),
         ylim = c(0, 0.6),
         main = paste("Counts Distribution for Interval Size:", interval_size),
         xlab = "Counts per Interval",
         xlim = c(-1,x_max))
    hist(observed_counts, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), col = rgb(1, 0, 0, 0.25), probability = TRUE, add = TRUE)
  } else if (interval_size == 1000) {
  # Plot the histograms with the same y-axis limit
    hist(simulated_counts, ylim = c(0,0.4), breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), probability = TRUE, col = rgb(0, 0, 1, 0.25), main = paste("Counts Distribution for Interval Size:", interval_size),
         xlab = "Counts per Interval")
    hist(counts, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), probability = TRUE, col = rgb(1, 0, 0, 0.25), add = TRUE)
  } else if (interval_size == 5000) {
    hist(simulated_counts, 
         probability = TRUE, 
         col = rgb(0, 0, 1, 0.25),
         breaks = 5,
         main = paste("Counts Distribution for Interval Size:", interval_size),
         xlab = "Counts per Interval",
         xlim = c(-1,20))
    hist(observed_counts, breaks = 10, col = rgb(1, 0, 0, 0.25), probability = TRUE, add = TRUE)
  } else {
    hist(simulated_counts, 
         probability = TRUE, 
         col = rgb(0, 0, 1, 0.25),
         ylim = c(0, 0.2),
         breaks = 5,
         main = paste("Counts Distribution for Interval Size:", interval_size),
         xlab = "Counts per Interval",
         xlim = c(0,25))
    hist(observed_counts, breaks = 8, col = rgb(1, 0, 0, 0.25), probability = TRUE, add = TRUE)
    
  }
  abline(v = expected_count[1], col = "red", lwd = 1, lty = 2)
  legend("topright", legend = c("Real Data", "Simulated Data"),
         fill = c(rgb(1, 0, 0, 0.25), rgb(0, 0, 1, 0.25)), border = NA)
}

# Display results summary
print(results_summary)
