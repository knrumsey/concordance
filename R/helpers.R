# This function is inspired by the 'abind' function from the 'abind' package
# Source: https://cran.r-project.org/web/packages/abind/index.html
abind_custom <- function(..., along){
  arrays <- list(...)
  # Check that all arrays have the same dimensions except for the 'along' dimension
  # Implement the binding logic here
  # This is a simplified example and may need more logic to handle all edge cases
  if (length(arrays) == 0) stop("No arrays provided")

  dims <- dim(arrays[[1]])
  for (array in arrays) {
    if (!all(dim(array)[-along] == dims[-along])) {
      stop("All arrays must have the same dimensions except for the 'along' dimension")
    }
  }

  result_dims <- dims
  result_dims[along] <- sum(sapply(arrays, function(x) dim(x)[along]))

  result <- array(NA, dim = result_dims)
  current_pos <- 1
  for (i in seq_along(arrays)){
    array <- arrays[[i]]
    slice <- seq(current_pos, length.out = dim(array)[along])
    index <- rep(list(TRUE), length(result_dims))
    index[[along]] <- slice
    result <- do.call(`[<-`, c(list(result), index, list(array)))
    current_pos <- current_pos + length(slice)
  }
  return(result)
}


find_duplicate_groups <- function(A) {
  # Convert rows to strings for easy comparison
  row_strings <- apply(A, 1, paste, collapse = ",")

  # Find unique strings and their indices
  unique_strings <- unique(row_strings)
  duplicate_groups <- list()

  # Loop through each unique string and find all matching row indices
  for (unique_string in unique_strings) {
    matching_indices <- which(row_strings == unique_string)
    if (length(matching_indices) > 1) {
      duplicate_groups <- c(duplicate_groups, list(matching_indices))
    }
  }
  return(duplicate_groups)
}
