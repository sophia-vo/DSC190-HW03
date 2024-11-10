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
lines(density(hcmv$location, adjust = 2), col = 2)

# Simulation of Uniform Distribution of Palindromic Sites
set.seed(123)
simulated_palindromes <- runif(num_palindromes, min = 1, max = sequence_length)

# Histogram of Simulated Data Palindromic Site Location
hist(simulated_palindromes, breaks = 25, probability = TRUE, col = rgb(0,0,1,0.25), main = "Palindromic Sites of Simulated Data",
     xlab = "Location in DNA Sequence")
lines(density(simulated_palindromes, adjust = 2), col = 4)

# Layered Histogram for Graphical Comparison
hist(hcmv$location, breaks = 25, probability = TRUE, col = rgb(1, 0, 0, 0.25),
     main = "Palindromic Sites of Real & Simulated Data",
     xlab = "Location in DNA Sequence")
lines(density(hcmv$location, adjust = 2), col = 2)
hist(simulated_palindromes, breaks = 25, probability = TRUE, col = rgb(0,0,1,0.25), add = TRUE)
lines(density(simulated_palindromes, adjust = 2), col = 4)
legend("topright", legend = c("Real", "Simulated"), fill = c(rgb(1, 0, 0, 0.25), rgb(0,0,1,0.25)), border = NA)

tab <- table(cut(hcmv$location, 
                     breaks = seq(from = 1, to = 230000, by = 1000), 
                     include.lowest = TRUE))
counts <- as.vector(tab)

sim_tab <- table(cut(simulated_palindromes, 
                 breaks = seq(from = 1, to = 230000, by = 1000), 
                 include.lowest = TRUE))
sim_counts <- as.vector(sim_tab)

hist(sim_counts, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), probability = TRUE, col = rgb(0, 0, 1, 0.25), main = "Real & Simulated Interval Counts", xlab = "number of points inside an interval", ylim = c(0,0.5))
lines(density(sim_counts, adjust = 2), col = rgb(0,0,1,0.5))

hist(counts, breaks = c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8), probability = TRUE, col = rgb(1, 0, 0, 0.25), add = TRUE)
lines(density(counts, adjust = 2), col = rgb(1,0,0,0.5))

legend("topright", legend = c("Real", "Simulated"), fill = c(rgb(1, 0, 0, 0.25), rgb(0,0,1,0.25)), border = NA)

# Calculate spacings and interval counts for each simulation

# Real data analysis
real_spacing <- diff(hcmv$location)
real_counts <- table(cut(hcmv$location, breaks = seq(0, sequence_length, by = interval_size), include.lowest = TRUE))

# Simulated data analysis
simulated_spacings <- diff(sort(simulated_palindromes))
simulated_counts <- table(cut(simulated_palindromes, breaks = seq(0, sequence_length, by = interval_size), include.lowest = TRUE))

# Quantitative comparison
mean_real_spacing <- mean(real_spacing)
mean_simulated_spacing <- mean(unlist(simulated_spacings))

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

# Visualization
# Spacing histograms
hist(real_spacing, breaks = 30, col = rgb(1, 0,0,0.25), probability = TRUE, main = "Real & Simulated Data Spacing", xlab = "Spacing")
hist(unlist(simulated_spacings), breaks = 20, probability = TRUE, col = rgb(0, 0, 1, 0.25), add = TRUE)
lines(density(real_spacing, adjust = 2), col = rgb(1,0,0,0.5))
lines(density(unlist(simulated_spacings), adjust = 2), col = rgb(0,0,1,0.5))
legend("topright", legend = c("Real Data", "Simulated Data"),
       fill = c(rgb(1, 0, 0, 0.25), rgb(0, 0, 1, 0.25)), border = NA)

# Question 2 - Locations and Spacings -------------------------------------

hist(real_spacing, breaks = 30, col = rgb(1, 0,0,0.25), probability = TRUE, main = "Real Data Spacing", xlab = "Spacing")

# Question 3 - Counts -----------------------------------------------------

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
