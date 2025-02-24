# Extract information from the h5 file
library(signal)
library(hdf5r)
library(pracma)
library(lifecycle)
library(prospectr)
setwd("D:/Soil_spectra_h5files/")
filename <- "AGES_exp.h5"
file.h5 <- H5File$new(filename, mode = "r")
source("Useful_functions.R")

plot_name <- "NT_I"
grassland_data <- extract_data(file.h5, plot_name)

#write.csv(grassland_data, file = "CT_I_data.csv", row.names = FALSE)


# Filter the unnecessary rows based on depth
dat <- filter_duplicates_by_insertion(grassland_data)
dim(dat)


# Remove first depths (optional) because they are too noisy
# dat <- filter_nrows_by_insertion(dat, 2) # first 2 depths
# Remove outliers instead of first nrows
#split_data <- split(dat, dat$Insertion)

depths <- dat$Depth

# Create bins for every 10 cm interval
bins <- cut(depths, breaks = seq(0, max(depths) + 10, by = 10), right = FALSE)
# View how the depths are grouped into bins
table(bins)
# Split the original data into a list based on these bins
split_data <- split(dat, bins)

#For dry lab insertions, submit depth and ins
#dat$Depth[5] <- 65

#dat$Insertion[5] <- 'ins015'


require(resemble)
require(dbscan)
require(rrcov)
## methods used are Local Outlier Factor (LOF) on PCA result and Spectral Information Divergence
## the methods are combined, 
## examples
#new_ins1 <- outlier_filter(split_data$ins001[,9:2756], 2, 0.95) # 2 is the magnitude by which is the elliptic center is multiplied, for selection outliers in spectral information divergence
#new_ins9 <- outlier_filter(split_data$ins009[,9:2756],2, 0.95) # 0.95 is the quantile selected for outliers determination in LOF
#new_ins11 <- outlier_filter(split_data$ins011[,9:2756], 2, 0.95)

### the magnitude can be changed (0.1 to inf), quantile as well (0 to 0.9999) so the selection of the outliers can be adjusted
## Loop through the whole data, for every insertions, after adjusting quantile
filtered_list <- list()

for (insertion in names(split_data)) {
  
  filtered_data <- outlier_filter(split_data[[insertion]], 0.8, 2.2) # specify the magnitude and quantile
  mtext(insertion, side = 4, line = -1.8, outer = TRUE, cex = 1.5)
  # Store the filtered data in the list
  filtered_list[[insertion]] <- filtered_data
}


dat <- do.call(rbind, filtered_list)
dim(dat)


# Filter the spectra
# exact column numbers are requiered!!
#####
## VIS-NIR , 9:2056
visnir_dat <- dat[, 9:2056]
white <- h5attr(file.h5[[plot_name]], "white_ref")
dark <- h5attr(file.h5[[plot_name]], "dark_ref")
oldWavs <- h5attr(file.h5[[plot_name]], "wavelengths_names")
plot(1:dim(visnir_dat)[2],visnir_dat[10,], "l") # example of raw spectrum
  
### 1-Normalization
norm_data <- normalize_data(visnir_dat, white, dark)
plot(1:dim(norm_data)[2], norm_data[10,], "l") # example

### 2-SNR (signal to noise ratio) filtering --- THIS should be done once in one plot to determine
# the wavehttp://127.0.0.1:47705/graphics/plot_zoom_png?width=1188&height=756lgths range, then when building whole data, just use the wavenumbers
snr <- calculate_snr(norm_data)
plot(1:length(snr),snr, "l")
abline(h=4, col=3)
#### This is one way to do it
# res <- filter_by_snr(norm_data, snr, 4, oldWavs)
# filtered <- res$dat
# dim(filtered)
# new_wav <- res$wav
# plot(1:dim(filtered)[2], filtered[20,], "l") #check the filtering
#### Another way is to trimmer vertically
plot(1:length(snr),snr, "l")
abline(v=200, col=2) # for example 
abline(v=1700, col=2) # for example 
new_wav <- oldWavs[200:1700] # exclude noisy extremes
filtered <- norm_data[,200:1700]

### 3-Resample, interpolation
resampled_spectral_data <- resample_spectra(filtered, new_wav, By=1)
dim(resampled_spectral_data)
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")

### 4-Smoothing, methods: "sgolay" , "whittaker"
library(signal)
smoothy <- smooth_data(resampled_spectral_data, method = "sgolay")
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")
lines(1:dim(resampled_spectral_data)[2], smoothy[10,], col=3)

# Ensure visnir_dat has the same column names as resampled_spectral_data
colnames(smoothy) <- colnames(resampled_spectral_data)

### Data substition
visnir_dat <- smoothy

### EPO correction
load("EPO_matrix_sp1.RData")
EPO_corrected <- as.matrix(visnir_dat) %*% EPO_matrix_sp1
colnames(EPO_corrected) <- colnames(resampled_spectral_data)


#####
## SWIR , 2057:2756
swir_dat <- dat[,2057:2756]
white <- h5attr(file.h5[[plot_name]], "white_ref2")
dark <- h5attr(file.h5[[plot_name]], "dark_ref2")
oldWavs <- h5attr(file.h5[[plot_name]], "wavelengths_names2")
plot(1:dim(swir_dat)[2],swir_dat[10,], "l") # example of raw spectrum

### 1-Normalization
norm_data <- normalize_data(swir_dat, white, dark)
plot(1:dim(norm_data)[2], norm_data[10,], "l") # example

### 2-PCA filtering --- THIS should be done once in one plot to determine
# the wavelgths range, then when building whole data, just use the wavenumbers
pca_result <- prcomp(t(as.matrix(norm_data)), center = TRUE, scale. = TRUE)
noise_components <- pca_result$x[, -(1:2)] 
noise_signal <- noise_components %*% t(pca_result$rotation[, -(1:2)])
noiz<-rowSums(abs(noise_signal))
plot(1:700, noiz, "l")
abline(h=150,col=3)
#### THis is one way to do it, it will cut by horizontal line
# res <- filter_by_pca(norm_data, 150, oldWavs)
# filtered <- res$dat
# dim(filtered)
# new_wav <- res$wav
# plot(1:dim(filtered)[2], filtered[20,], "l") #check the filtering
# #### However, another way is to Trimmer the spectrum based on what is visible
# in the plot(1:700, noiz, "l")
plot(1:700, noiz, "l")
abline(v=80, col=2) # for example
abline(v=670, col=2) # for example
new_wav <- oldWavs[80:670] # exclude noisy extremes
filtered <- norm_data[,80:670]


### 3-Resample, interpolation
resampled_spectral_data <- resample_spectra(filtered, new_wav, By=1)
dim(resampled_spectral_data)
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")

### 4-Smoothing, methods: "sgolay" , "whittaker"
smoothy <- smooth_data(resampled_spectral_data, method = "sgolay")
plot(1:dim(resampled_spectral_data)[2], resampled_spectral_data[10,], "l")
lines(1:dim(resampled_spectral_data)[2], smoothy[10,], col=3)

# Ensure visnir_dat has the same column names as resampled_spectral_data
colnames(smoothy) <- colnames(resampled_spectral_data)

### Data substition
swir_dat <- smoothy

### EPO correction
load("EPO_matrix_sp2.RData")
EPO_corrected2 <- as.matrix(swir_dat) %*% EPO_matrix_sp2
colnames(EPO_corrected2) <- colnames(resampled_spectral_data)
plot(1:dim(EPO_corrected2)[2], EPO_corrected2[10,], "l")
lines(1:dim(EPO_corrected2)[2], EPO_corrected2[100,], col=3)
lines(1:dim(EPO_corrected2)[2], EPO_corrected2[1000,], col=2)
lines(1:dim(EPO_corrected2)[2], EPO_corrected2[1100,], col=4)

pca_result <- prcomp(t(as.matrix(EPO_corrected2)), center = TRUE, scale. = TRUE)
noise_components <- pca_result$x[, -(1:2)] 
noise_signal <- noise_components %*% t(pca_result$rotation[, -(1:2)])
noiz<-rowSums(abs(noise_signal))
plot(1:length(noiz), noiz, "l")




#####
# Finally filtered spectra build a new dataset, of the plot
filtered_dat <- cbind(dat[,1:8], visnir_dat, swir_dat)
# Save it !!
write.csv(filtered_dat, file = "Dry_RT_III.csv", row.names = FALSE)

