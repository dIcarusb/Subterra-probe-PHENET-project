# Useful functions script
# Here are the functions used as source in other scripts, thus, this script
# MUST be in the same directory as the other R scripts, in order to get use of
# the functions written in it!!



# Extract data from a plot into 2D matrix
extract_data <- function(file, plot_name) {
  # Initialize an empty list to store all data points
  all_data <- list()
  
  # Check for the existence of the date attribute at the plot level
  date <- if (file[[plot_name]]$attr_exists("date")) {
    file[[plot_name]]$attr_open("date")$read()
  } else {
    NA
  }
  
  # Identify insertions in the specified plot
  insertions <- list.groups(file[[plot_name]], recursive = FALSE)
  
  # Loop through each insertion point
  for (insertion in insertions) {
    # Get the list of datasets for the current insertion point
    insertion_path <- paste0(plot_name, "/", insertion)
    
    depths <- file[[paste0(insertion_path, "/depths")]][]
    forces <- file[[paste0(insertion_path, "/forces")]][]
    VISNIR_spectra <- file[[paste0(insertion_path, "/VISNIR_spectra")]][,]
    NIRSWIR_spectra2 <- file[[paste0(insertion_path, "/NIRSWIR_spectra2")]][,]
    
    # Get SOC
    soc <- if (file[[insertion_path]]$exists("SOC")) {
      file[[paste0(insertion_path, "/SOC")]][]
    } else {
      NA
    }
    
    # Get water content
    water_content <- if (file[[insertion_path]]$exists("water_content")) {
      file[[paste0(insertion_path, "/water_content")]][]
    } else {
      NA
    }
    
    # Get the GPS attribute for the current insertion
    gps <- if (file[[insertion_path]]$attr_exists("GPS")) {
      file[[insertion_path]]$attr_open("GPS")$read()
    } else {
      NA
    }
    
    # Get the wavelength names
    wavelengths_names <- h5attr(file[[plot_name]], "wavelengths_names")
    wavelengths_names2 <- h5attr(file[[plot_name]], "wavelengths_names2")
    
    # Get the number of data points (rows)
    num_data_points <- length(depths)
    
    # Loop through each data point and construct the data rows
    for (i in 1:num_data_points) {
      # Create a data frame for each data point
      data_point <- data.frame(
        PlotName = plot_name,
        Insertion = insertion,
        GPS = gps,
        Date = date,  # Convert date to Date type if not NA
        Depth = depths[i],
        Force = forces[i],
        SOC = soc[i],
        Water = water_content[i]
        
      )
      
      # Add VISNIR spectra columns
      for (j in 1:length(wavelengths_names)) {
        col_name <- paste0("VISNIR_", wavelengths_names[j])
        data_point[[col_name]] <- VISNIR_spectra[i, j]
      }
      
      # Add NIRSWIR spectra columns
      for (j in 1:length(wavelengths_names2)) {
        col_name <- paste0("NIRSWIR_", wavelengths_names2[j])
        data_point[[col_name]] <- NIRSWIR_spectra2[i, j]
      }
      
      # Append the data point to the all_data list
      all_data <- append(all_data, list(data_point))
    }
  }
  
  # Combine all data points into a single data frame
  final_data <- do.call(rbind, all_data)
  return(final_data)
}



# Outliers filtering, per insertion
outlier_filter <- function(Dataset, thresh1=0.8, thresh2=2.2) {
  ins = Dataset[, 9:2756]
  # Spectral Information Divergence
  sidD <- sid(Xr = as.matrix(ins),
              Xu = as.matrix(ins),
              mode = "density",
              center = FALSE, scale = TRUE)
  distances <- apply(sidD$sid, 1, function(row) mean(row[row != 0]))
  # Once the distances are calculated we perform Robust Location and Scatter Estimation via MCD
  # check ?CovMcd for more details
  # With this we obtain a center of the distances, which later can be used to determine the oultiers
  #cov_model1 <- CovMcd(distances) # Robust Covariance Estimation (Elliptic Envelope)
  
  # Local Outlier Factor (LOF)
  P <- prcomp(ins, scale. = TRUE)
  lof_scores <- lof(P$x[, 1:10], minPts = 30)
  #cov_model2 <- CovMcd(lof_scores) # Robust Covariance Estimation
  
  par(mfrow=c(1,2))
  plot(1:length(distances), distances, xlab="Data points", ylim=c(0,1.6))
  title("Spectral Information Divergence spectra mean")
  abline(h=thresh1, col=2)
  plot(1:dim(P$x)[1], lof_scores, xlab="Data points", ylim=c(0,12))
  title("Local Outlier Factor (LOF)")
  abline(h=thresh2, col=2)
  
  #outliers <- sort(unique(c(which( distances > quantile(distances, quantil)), 
                            #which( lof_scores > quantile(lof_scores, quantil)))))
  # outliers are detected as the points which are out of the center multiplied by the magnitude
  outliers <- sort(unique(c(which(distances > thresh1), 
                            which(lof_scores > thresh2))))
  
  #new_ins <- ins[-outliers,]
  if(length(outliers) == 0) {
    new_dataset <- Dataset
  }
  else {
    new_dataset <- Dataset[-outliers,]
  }
  
  return(new_dataset)
  
}


# Filter duplicate depths
filter_duplicates_by_insertion <- function(data) {
  # Split data by the "Insertion" column
  split_data <- split(data, data$Insertion)
  
  # Define a helper function to filter each subset
  filter_subset <- function(subset) {
    subset[!duplicated(round(subset$Depth)), ]
  }
  
  # Apply the filtering function to each subset
  filtered_list <- lapply(split_data, filter_subset)
  
  # Combine the filtered subsets back into a single data frame
  filtered_data <- do.call(rbind, filtered_list)
  
  return(filtered_data)
}


# filter the first n rows by insertion
filter_nrows_by_insertion <- function(data, n) {
  # Split data by the "Insertion" column
  split_data <- split(data, data$Insertion)
  
  # Define a helper function to filter each subset
  filter_subset <- function(subset) {
    if (nrow(subset) > n) {
      subset <- subset[-(1:n), ]
    } else {
      subset <- subset[0, ]
    }
    return(subset)
  }
  
  # Apply the filtering function to each subset
  filtered_list <- lapply(split_data, filter_subset)
  
  # Combine the filtered subsets back into a single data frame
  filtered_data <- do.call(rbind, filtered_list)
  
  return(filtered_data)
}



# smoothing
smooth_data <- function(data, method = "sgolay", ...) {
  if (method == "sgolay") {
    t(apply(data, 1, function(row) {
      sgolayfilt(row, ...)
    }))
  } else if (method == "whittaker") {
    t(apply(data, 1, function(row) {
      whittaker(row, ...)
    }))
  } else {
    stop("Invalid smoothing method. Choose 'sgolay' or 'whittaker'.")
  }
}



# normalize one spectrum
normalize<-function(reflection, white, dark,  epsilon = 1e-8){ #white/dark normalization
  denominator <- white - dark
  result <- (reflection - dark) / (denominator + (denominator == 0) * epsilon)
  return(result)
}

# function to normalize the entire dataset
normalize_data <- function(data, white, dark) {
  # Ensure white and dark references are vectors of the same length as the number of columns in data
  if (length(white) != ncol(data) || length(dark) != ncol(data)) {
    stop("White and dark references must have the same length as the number of columns in data.")
  }
  
  # Apply the normalize function to each row
  normalized_data <- t(apply(data, 1, function(row) normalize(row, white, dark)))
  return(normalized_data)
}


# Interpolation function for spectra resampling
resample_spectra <- function(data, oldWavs, By) {
  require(prospectr)
  newWavs <- seq(floor(min(oldWavs)), ceiling(max(oldWavs)), by = By)
  resampled_spectral_data <- prospectr::resample(X = data,
                                                wav = oldWavs,
                                                new.wav = newWavs,
                                                interpol = "spline")
  return(resampled_spectral_data)
}


# Remove noise from your spectra using the signal-to-noise ratio (SNR)
# Compute the SNR as the ratio of the mean signal to the standard deviation
calculate_snr <- function(data) {
  # Calculate mean and standard deviation for each wavelength
  mean_signal <- apply(data, 2, mean)
  std_noise <- apply(data, 2, sd)
  
  # Calculate SNR
  snr <- mean_signal / std_noise
  return(snr)
}

# Exclude the wavelengths (or regions) where the SNR is below this threshold
# Use the identified non-noisy wavelengths to filter the spectral data
filter_by_snr <- function(data, snr, threshold, waves) {
  # Identify wavelengths with SNR above the threshold
  valid_wavelengths <- which(snr > threshold) #< or > depends
  
  # Filter the data by valid wavelengths
  filtered_data <- data[, valid_wavelengths, drop = FALSE]
  wav_filtered <- waves[valid_wavelengths, drop = FALSE]
  return(list(dat=filtered_data, wav=wav_filtered))
}


# filter by pca
filter_by_pca <- function(data, threshold, wavelength) {
  pca_result <- prcomp(t(as.matrix(data)), center = TRUE, scale. = TRUE)
  # Assuming the first few principal components capture the signal
  noise_components <- pca_result$x[, -(1:2)]
  # Reconstruct noise-only signal
  noise_signal <- noise_components %*% t(pca_result$rotation[, -(1:2)])
  good_indices <- which(rowSums(abs(noise_signal)) < threshold)
  filtered_data <- data[, good_indices, drop = FALSE]
  wavelength_filtered <- wavelength[ good_indices, drop = FALSE]
  
  return(list(dat=filtered_data, wav=wavelength_filtered))
}

# get insertion depths in one session
get_insertion_depths <- function(sess01, f) {
  depths_list <- list()
  
  # Loop through each insertion
  for (ins in seq_along(sess01)) {
    # Get the insertion path
    insertion_path <- sess01[ins]
    
    # Extract the depth (number of rows in the 'spectra' dataset)
    depth <- dim(f[[insertion_path]][["spectra"]][,])[2]
    
    # Create a data frame for the current insertion and its depth
    depth_data <- data.frame(
      Insertion = insertion_path,
      Depth = depth
    )
    
    # Append the data frame to the list
    depths_list <- append(depths_list, list(depth_data))
  }
  
  # Combine all data frames into a single data frame
  depths_table <- do.call(rbind, depths_list)
  
  return(depths_table)
}


# loop to add insertions into a new file based on the old file
# specify : new file name, insertion list seq(1,15,1), plot name, spectra1 names, spectra2 names, wavelength1 and 2 names, finally the order either starting from 1 to 15 or from 15 to 1
process_insertions <- function(file, insertions, plot_name, spectrometer1_groups, spectrometer2_groups, wavelengths_name, wavelengths_name2, order=TRUE,lab=FALSE) {
  # Reverse the insertions if order is FALSE
  if (!order) { insertions <- rev(insertions) }
  
  for (i in seq_along(insertions)) {
    # Create the insertion path with sequential names starting from ins001
    insertion_name <- sprintf("ins%03d", i)
    path <- paste0(plot_name, insertion_name)
    file$create_group(path)

    # Extract depths
    if (!lab) { 
      depths <- f[[spectrometer1_groups[insertions[i]]]][["spectrum_depths"]][]
    }
    else {
      depths <- rep(0, depths_table[1,2])
    }
    
    # Extract forces
    if (!lab) {
      forces <- f[[spectrometer1_groups[insertions[i]]]][["spectrum_forces"]][]
    }
    else {
      forces <- rep(0, depths_table[1,2])
    }
    
    # Extract VISNIR spectra
    VISNIR_spectra <- t(f[[spectrometer1_groups[insertions[i]]]][["spectra"]][,])
    colnames(VISNIR_spectra) <- wavelengths_name
    
    # Extract NIRSWIR spectra
    NIRSWIR_spectra2 <- t(f[[spectrometer2_groups[insertions[i]]]][["spectra"]][,])
    colnames(NIRSWIR_spectra2) <- wavelengths_name2
    
    # Add datasets to HDF5 file
    file[[path]]$create_dataset(name="depths", robj=depths)
    file[[path]]$create_dataset(name="forces", robj=forces)
    file[[path]]$create_dataset(name="VISNIR_spectra", robj=VISNIR_spectra)
    file[[path]]$create_dataset(name="NIRSWIR_spectra2", robj=NIRSWIR_spectra2)
  }
}




    
    
# clean insertions before loop 
clean_insertions <- function(spectrometer1_groups, spectrometer2_groups, insertions_to_remove) {
  # Ensure the insertions_to_remove is a vector of indices
  insertions_to_remove <- as.integer(insertions_to_remove)
  
  # Remove the specified insertions from spectrometer1_groups and spectrometer2_groups
  cleaned_spectrometer1_groups <- spectrometer1_groups[-insertions_to_remove]
  cleaned_spectrometer2_groups <- spectrometer2_groups[-insertions_to_remove]
  
  return(list(
    spectrometer1_groups = cleaned_spectrometer1_groups,
    spectrometer2_groups = cleaned_spectrometer2_groups
  ))
}



# data split into train, validation and testing sets
# seed is a random seed for sampling, 
# train_prop, valid_prop, test_prop are proportions from 0 to 1 (100%)
split_data <- function(data, seed, train_prop, valid_prop, test_prop) {
  # Check if the proportions sum to 1
  if (train_prop + valid_prop + test_prop != 1) {
    stop("The sum of train, validation, and test proportions must equal 1.")
  }
  
  # Set seed for reproducibility
  set.seed(seed)
  
  # Determine the number of observations in the data
  n <- nrow(data)
  
  # Create indices for the splits
  train_indices <- sample(seq_len(n), size = train_prop * n)
  remaining_indices <- setdiff(seq_len(n), train_indices)
  validation_indices <- sample(remaining_indices, size = valid_prop * n)
  test_indices <- setdiff(remaining_indices, validation_indices)
  
  # Create the datasets
  train_data <- data[train_indices, ]
  validation_data <- data[validation_indices, ]
  test_data <- data[test_indices, ]
  
  return(list(train = train_data, validation = validation_data, test = test_data))
}


# data split of multiple datasets, it requieres split_data
# similar arguments as split_data are requiered
process_multiple_datasets <- function(datasets, seed, train_prop, valid_prop, test_prop) {
  # Initialize empty lists to store combined datasets
  combined_train <- list()
  combined_validation <- list()
  combined_test <- list()
  
  # Loop through each dataset and apply the split_data function
  for (dataset in datasets) {
    split_result <- split_data(dataset, seed, train_prop, valid_prop, test_prop)
    combined_train <- append(combined_train, list(split_result$train))
    combined_validation <- append(combined_validation, list(split_result$validation))
    combined_test <- append(combined_test, list(split_result$test))
  }
  
  # Combine all datasets in each list into one data frame
  final_train <- do.call(rbind, combined_train)
  final_validation <- do.call(rbind, combined_validation)
  final_test <- do.call(rbind, combined_test)
  
  # Shuffle the datasets
  final_train <- final_train[sample(nrow(final_train)), ]
  final_validation <- final_validation[sample(nrow(final_validation)), ]
  final_test <- final_test[sample(nrow(final_test)), ]
  
  return(list(train = final_train, validation = final_validation, test = final_test))
}

