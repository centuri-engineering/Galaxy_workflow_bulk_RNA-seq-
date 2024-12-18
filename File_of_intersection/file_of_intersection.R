# Avoid crashing Galaxy with a UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Take in trailing command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Function to get the value of an option from command line arguments
getOptionValue <- function(option) {
  idx <- match(option, args)
  if (!is.na(idx) && (idx + 1) <= length(args)) {
    return(args[idx + 1])
  }
  return(NULL)
}

# Input
input_file <- getOptionValue("--galaxy_input")

# Output
output_file <- getOptionValue("--galaxy_output")

data <- read.table(input_file, sep = "\t", header = FALSE)

conditions <- data[-1, 2]

# Base list of conditions
conditions_unique <- unique(conditions)

# Step 1: Generate base elements with "up" and "down"
generate_elements <- function(conditions_unique) {
  comb <- expand.grid(cond1 = conditions_unique, cond2 = conditions_unique)  # All possible pairs
  comb <- comb[comb$cond1 != comb$cond2, ]  # Remove pairs where cond1 == cond2
  
  # Add "up" and "down" to each pair
  elements <- c(
    paste(comb$cond1, "vs", comb$cond2, "up"),
    paste(comb$cond1, "vs", comb$cond2, "down")
  )
  
  return(elements)
}

# Step 2: Generate all combinations without inverted duplicates
generate_combinations <- function(elements) {
  comb <- combn(elements, 2)  # Generate all unique pairs (no inverted duplicates)
  results <- apply(comb, 2, function(x) paste(x[1], "AND", x[2]))  # Concatenate with "AND"
  
  return(results)
}

# Call functions to get the final result
elements <- generate_elements(conditions_unique)
final_result <- generate_combinations(elements)

# Save the final result to a TSV file
write.table(final_result, file = output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
