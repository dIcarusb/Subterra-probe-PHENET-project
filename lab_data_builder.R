# Check the structure of the h5 file
library(hdf5r)
setwd("D:/Soil_spectra_h5files") # your directory
filename <- "lab_ages.h5" #AGES_exp
file.h5 <- H5File$new(filename, mode = "r")

# print the file, it shows the attributes and groups that it contains
file.h5


getMeanSpectrum <- function(spectra) {
  # Ensure the input is a matrix
  if (!is.matrix(spectra)) {
    spectra <- as.matrix(spectra)
  }
  
  
  # Compute the mean of each column (i.e., mean for each wavelength)
  mean_spectrum <- rowMeans(spectra)
  
  # Return the mean spectrum (a vector of length 613)
  return(mean_spectrum)
} # get mean spectrum out of 3

# the spectrometer number must be changed in order to get one or another
getInsertionSpectraMatrix <- function(session, insertions, call_number) {
  # Check that the insertions vector is not empty
  if (length(insertions) == 0) {
    stop("The insertions vector is empty.")
  }
  
  # Use the first insertion to determine the number of wavebands
  first_ins_str <- sprintf("%03d", insertions[1])
  # Build the path: "cal001/insXXX/spectrometer1/spectra"
  first_path <- paste0(call_number, "/ins", first_ins_str, "/spectrometer1/spectra")
  
  # Retrieve the spectra matrix for the first insertion (assumed to be 8x613)
  spectra_first <- session[[first_path]][, ]
  
  # Compute the mean spectrum using your helper function
  mean_spec_first <- getMeanSpectrum(spectra_first)
  
  # Determine the number of wavebands (i.e., length of the mean spectrum)
  n_wavebands <- length(mean_spec_first)
  
  # Initialize a result matrix with one row per insertion and columns equal to wavebands
  result_matrix <- matrix(NA, nrow = length(insertions), ncol = n_wavebands)
  
  # Loop over each insertion in the provided vector
  for (i in seq_along(insertions)) {
    # Format the insertion number to a three-digit string (e.g., "007", "012")
    ins_str <- sprintf("%03d", insertions[i])
    
    # Build the path for the current insertion using the provided call number
    path <- paste0(call_number, "/ins", ins_str, "/spectrometer1/spectra")
    
    # Retrieve the spectra matrix for this insertion
    spectra <- session[[path]][, ]
    
    # Compute the mean spectrum for this insertion
    mean_spec <- getMeanSpectrum(spectra)
    
    # Add the mean spectrum as a row in the result matrix
    result_matrix[i, ] <- mean_spec
  }
  
  return(result_matrix)
}



#####
## DRY
RT <- getInsertionSpectraMatrix(file.h5[["session002"]], 1:3, "cal001")
RT30_40 <- getInsertionSpectraMatrix(file.h5[["session007"]], 1, "cal001")
RT <- rbind(RT,RT30_40)
rt <- getInsertionSpectraMatrix(file.h5[["session002"]], 4:5, "cal001")
RT <- rbind(RT, rt)

CT <- getInsertionSpectraMatrix(file.h5[["session007"]], 2, "cal001")
ct <- getInsertionSpectraMatrix(file.h5[["session003"]], 1:5, "cal001")
CT <- rbind(CT, ct)

NT <- getInsertionSpectraMatrix(file.h5[["session004"]], 1:6, "cal002")

GL <- getInsertionSpectraMatrix(file.h5[["session005"]], c(7,12), "cal001")

DRY <- rbind(RT, CT, NT, GL)

#####
## WET
filename <- "AT_lab.h5" #AGES_exp
file2.h5 <- H5File$new(filename, mode = "r")

# print the file, it shows the attributes and groups that it contains
file2.h5

RT <- getInsertionSpectraMatrix(file2.h5[["session012"]], 1:6, "cal002")

CT <- getInsertionSpectraMatrix(file2.h5[["session012"]], 7:12, "cal002")

NT <- getInsertionSpectraMatrix(file2.h5[["session012"]], 13:18, "cal002")

GL <- getInsertionSpectraMatrix(file2.h5[["session012"]], 19:20, "cal002")

WET <- rbind(RT, CT, NT, GL)

####
# Normalization by dark and white

# dry
white <- h5attr(file.h5[["session002"]][["cal001"]],"white_spectrum")
dark <- h5attr(file.h5[["session002"]][["cal001"]],"dark_spectrum")

# wet
white <- h5attr(file2.h5[["session012"]][["cal001"]],"white_spectrum")
dark <- h5attr(file2.h5[["session012"]][["cal001"]],"dark_spectrum")

# Normalization
source("Useful_functions.R")
norm_data <- normalize_data(DRY, white, dark)
plot(1:dim(norm_data)[2], norm_data[10,], "l") # example

norm_data2 <- normalize_data(WET, white, dark)
plot(1:dim(norm_data2)[2], norm_data2[10,], "l") # example

# Trimmer < ask thresholds Akos
#dry
oldWavs <- h5attr(file.h5[["session001"]][["cal001"]],"spec1_wavelengths_vector")
# wet
oldWavs <- h5attr(file2.h5[["session012"]][["cal001"]],"spec1_wavelengths_vector")

new_wav <- oldWavs[200:1700] # exclude noisy extremes
filtered <- norm_data2[,200:1700]

plot(1:dim(filtered)[2], filtered[10,], "l")

# Interpolation
resampled_spectral_data <- resample_spectra(filtered, new_wav, By=1)
dim(resampled_spectral_data)
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")

# Smoothing
library(signal)
smoothy <- smooth_data(resampled_spectral_data, method = "sgolay")
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")
lines(1:dim(resampled_spectral_data)[2], smoothy[10,], col=3)


#####
## EPO correction

epo0_ref<- smoothy # dry
epo1_ref <- smoothy # wet

epo <- function(D, npc){
  # npc is the number of components to use
  # return: P: the projection matrix
  D <- as.matrix(D)
  n <- nrow(D)
  p <- ncol(D)
  dtd <- t(D) %*% (D)
  # singular value decomposition of the D (n x n) matrix
  s <- svd(dtd)
  # extract the no. factors
  ld <- s$v[, 1:npc]
  # projection matrix
  P <- diag(p) - ld %*% t(ld)
  return(P)
}

D = as.matrix(epo0_ref - epo1_ref)
SC <- (epo0_ref + epo1_ref)/2
B <- cov(SC)
TrB <- sum(diag(B))
T <- cov(epo1_ref)
nc <- 20
xc <- seq(0:nc)
wilks <- matrix(0, nrow = nc + 1, ncol = 1) # create
# no transformation
wilks[1] <- TrB/sum(diag(T))
# make a for loop to estimate the Wilks" lambda
for (i in 1:nc) {
  npc <- i
  D = as.matrix(epo0_ref - epo1_ref)
  P <- epo(D, npc)
  Z0 <- as.matrix(epo0_ref) %*% P
  Z1 <- as.matrix(epo1_ref) %*% P
  T <- cov(Z1)
  SC <- (Z0+Z1)/2
  B <- cov(SC)
  TrB <-sum(diag(B))
  wilks[i+1] <- TrB/sum(diag(T))
}

wilks

plot(1:length(wilks[,1]), wilks[,1], "l")

which.max(wilks)




epo_correct <- function(epo0_ref, epo1_ref, npc=3) {
  D = as.matrix(epo0_ref - epo1_ref)
  P <- epo(D, npc)
  return(P)
}

EPO_matrix_sp1 <- epo_correct(epo0_ref, epo1_ref, npc=5) #5
save(EPO_matrix_sp1, file = "EPO_matrix_sp1.RData")


# example
# visnir
EPO_corrected2 <- as.matrix(visnir_dat) %*% EPO_matrix_sp1
colnames(EPO_corrected) <- colnames(resampled_spectral_data)


#####
# SWIR
# the spectrometer number must be changed in order to get one or another
getInsertionSpectraMatrix <- function(session, insertions, call_number) {
  # Check that the insertions vector is not empty
  if (length(insertions) == 0) {
    stop("The insertions vector is empty.")
  }
  
  # Use the first insertion to determine the number of wavebands
  first_ins_str <- sprintf("%03d", insertions[1])
  # Build the path: "cal001/insXXX/spectrometer1/spectra"
  first_path <- paste0(call_number, "/ins", first_ins_str, "/spectrometer2/spectra")
  
  # Retrieve the spectra matrix for the first insertion (assumed to be 8x613)
  spectra_first <- session[[first_path]][, ]
  
  # Compute the mean spectrum using your helper function
  mean_spec_first <- getMeanSpectrum(spectra_first)
  
  # Determine the number of wavebands (i.e., length of the mean spectrum)
  n_wavebands <- length(mean_spec_first)
  
  # Initialize a result matrix with one row per insertion and columns equal to wavebands
  result_matrix <- matrix(NA, nrow = length(insertions), ncol = n_wavebands)
  
  # Loop over each insertion in the provided vector
  for (i in seq_along(insertions)) {
    # Format the insertion number to a three-digit string (e.g., "007", "012")
    ins_str <- sprintf("%03d", insertions[i])
    
    # Build the path for the current insertion using the provided call number
    path <- paste0(call_number, "/ins", ins_str, "/spectrometer2/spectra")
    
    # Retrieve the spectra matrix for this insertion
    spectra <- session[[path]][, ]
    
    # Compute the mean spectrum for this insertion
    mean_spec <- getMeanSpectrum(spectra)
    
    # Add the mean spectrum as a row in the result matrix
    result_matrix[i, ] <- mean_spec
  }
  
  return(result_matrix)
}




## DRY
RT <- getInsertionSpectraMatrix(file.h5[["session002"]], 1:3, "cal001")
RT30_40 <- getInsertionSpectraMatrix(file.h5[["session007"]], 1, "cal001")
RT <- rbind(RT,RT30_40)
rt <- getInsertionSpectraMatrix(file.h5[["session002"]], 4:5, "cal001")
RT <- rbind(RT, rt)

CT <- getInsertionSpectraMatrix(file.h5[["session007"]], 2, "cal001")
ct <- getInsertionSpectraMatrix(file.h5[["session003"]], 1:5, "cal001")
CT <- rbind(CT, ct)

NT <- getInsertionSpectraMatrix(file.h5[["session004"]], 1:6, "cal002")

GL <- getInsertionSpectraMatrix(file.h5[["session005"]], c(7,12), "cal001")

DRY <- rbind(RT, CT, NT, GL)


## WET
filename <- "AT_lab.h5" #AGES_exp
file2.h5 <- H5File$new(filename, mode = "r")

# print the file, it shows the attributes and groups that it contains
file2.h5

RT <- getInsertionSpectraMatrix(file2.h5[["session012"]], 1:6, "cal002")

CT <- getInsertionSpectraMatrix(file2.h5[["session012"]], 7:12, "cal002")

NT <- getInsertionSpectraMatrix(file2.h5[["session012"]], 13:18, "cal002")

GL <- getInsertionSpectraMatrix(file2.h5[["session012"]], 19:20, "cal002")

WET <- rbind(RT, CT, NT, GL)

####
# Normalization by dark and white

# dry
white <- h5attr(file.h5[["session002"]][["cal001"]],"white_spectrum2")
dark <- h5attr(file.h5[["session002"]][["cal001"]],"dark_spectrum2")

# wet
white <- h5attr(file2.h5[["session012"]][["cal001"]],"white_spectrum2")
dark <- h5attr(file2.h5[["session012"]][["cal001"]],"dark_spectrum2")

# Normalization
source("Useful_functions.R")
norm_data <- normalize_data(DRY, white, dark)
plot(1:dim(norm_data)[2], norm_data[10,], "l") # example

norm_data2 <- normalize_data(WET, white, dark)
plot(1:dim(norm_data2)[2], norm_data2[10,], "l") # example

# Trimmer < ask thresholds Akos
#dry
oldWavs <- h5attr(file.h5[["session001"]][["cal001"]],"spec2_wavelengths_vector")
new_wav <- oldWavs[80:670] # exclude noisy extremes
filtered <- norm_data[,80:670]
# wet
oldWavs <- h5attr(file2.h5[["session012"]][["cal001"]],"spec2_wavelengths_vector")
new_wav <- oldWavs[80:670] # exclude noisy extremes
filtered <- norm_data2[,80:670]

plot(1:dim(filtered)[2], filtered[10,], "l")

# Interpolation
resampled_spectral_data <- resample_spectra(filtered, new_wav, By=1)
dim(resampled_spectral_data)
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")

# Smoothing
library(signal)
smoothy <- smooth_data(resampled_spectral_data, method = "sgolay")
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")
lines(1:dim(resampled_spectral_data)[2], smoothy[10,], col=3)


#####
## EPO correction

epo0_ref<- smoothy # dry
epo1_ref <- smoothy # wet

epo <- function(D, npc){
  # npc is the number of components to use
  # return: P: the projection matrix
  D <- as.matrix(D)
  n <- nrow(D)
  p <- ncol(D)
  dtd <- t(D) %*% (D)
  # singular value decomposition of the D (n x n) matrix
  s <- svd(dtd)
  # extract the no. factors
  ld <- s$v[, 1:npc]
  # projection matrix
  P <- diag(p) - ld %*% t(ld)
  return(P)
}

D = as.matrix(epo0_ref - epo1_ref)
SC <- (epo0_ref + epo1_ref)/2
B <- cov(SC)
TrB <- sum(diag(B))
T <- cov(epo1_ref)
nc <- 30
xc <- seq(0:nc)
wilks <- matrix(0, nrow = nc + 1, ncol = 1) # create
# no transformation
wilks[1] <- TrB/sum(diag(T))
# make a for loop to estimate the Wilks" lambda
for (i in 1:nc) {
  npc <- i
  D = as.matrix(epo0_ref - epo1_ref)
  P <- epo(D, npc)
  Z0 <- as.matrix(epo0_ref) %*% P
  Z1 <- as.matrix(epo1_ref) %*% P
  T <- cov(Z1)
  SC <- (Z0+Z1)/2
  B <- cov(SC)
  TrB <-sum(diag(B))
  wilks[i+1] <- TrB/sum(diag(T))
}

wilks

plot(1:length(wilks[,1]), wilks[,1], "l")

which.max(wilks)




epo_correct <- function(epo0_ref, epo1_ref, npc=3) {
  D = as.matrix(epo0_ref - epo1_ref)
  P <- epo(D, npc)
  return(P)
}

EPO_matrix_sp2 <- epo_correct(epo0_ref, epo1_ref, npc=21) #5
save(EPO_matrix_sp2, file = "EPO_matrix_sp2.RData")

