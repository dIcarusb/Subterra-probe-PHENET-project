## Image feature extraction

library(imager)
library(moments)
library(entropy)
library(moments)


# Loading
setwd("D:/Soil_spectra_h5files")
# image
original_image <- load.image("soil_foto/A3_T001_L005_2019.12.15_100633_001_ABC.jpg")
original_image2 <- load.image("soil_foto/GN2020M05D23_T010_L002_2020.05.26_165656_002_PS.jpg") # control
original_image2 <-imrotate(original_image2, 90)
light_image <- load.image("soil_foto/5.60-80 Grassland b.jpg")
plot(original_image2)
# mask
masking <- function(im, t=0.5) {
  g<- grayscale(im)
  th <- g < t
  plot(th)
  return(th)
}

mask <- masking(original_image2)
mask <- masking(light_image, 0.6)


library(abind)
mask_3d <- abind(mask, mask, mask, along = 4) 
# semgented image
segmented_image <- original_image2 * mask_3d
segmented_image <- light_image * mask_3d
plot(segmented_image)

# size reduction, image compression
# Resize using bicubic interpolation
compressed_bicubic <- imresize(segmented_image, scale = 1/8, interpolation = "bicubic")

# Resize using Lanczos interpolation
compressed_lanczos <- imresize(original_image, scale = 1/8, interpolation = "lanczos")

segmented_image1 <- compressed_lanczos
segmented_image2 <- compressed_lanczos # light image
control_im <- compressed_lanczos


# LAB
example <- segmented_image1 #original_image2
w <- dim(example)[1]
h <- dim(example)[2]

rgb_matrix <- matrix(segmented_image1[,,,-4, drop = FALSE], ncol = 3)
rgb_matrix <- matrix(original_image2, ncol = 3)

lab_matrix <- convertColor(rgb_matrix, from = "sRGB", to = "Lab", scale.in = 1)

L_img <- as.cimg(matrix(lab_matrix[,1], nrow = w, ncol = h))
a_img <- as.cimg(matrix(lab_matrix[,2], nrow = w, ncol = h))
b_img <- as.cimg(matrix(lab_matrix[,3], nrow = w, ncol = h))

# Optionally, you can display the three channels side by side:
par(mfrow = c(1, 3))
plot(L_img, main = "L Channel")
plot(a_img, main = "a Channel")
plot(b_img, main = "b Channel")

# HSV
rgb_matrix_t <- t(rgb_matrix)

#   Row 3: Value
hsv_matrix <- rgb2hsv(rgb_matrix_t)

# Extract the individual HSV components as vectors:
H_values <- hsv_matrix[1, ]
S_values <- hsv_matrix[2, ]
V_values <- hsv_matrix[3, ]

### 3. Reshape and plot each HSV component

# For plotting, reshape each component into a matrix with dimensions w x h.
# (Be sure to use the same ordering as when flattening.)
H_img <- as.cimg(matrix(H_values, nrow = w, ncol = h, byrow = FALSE))
S_img <- as.cimg(matrix(S_values, nrow = w, ncol = h, byrow = FALSE))
V_img <- as.cimg(matrix(V_values, nrow = w, ncol = h, byrow = FALSE))

par(mfrow = c(1, 3))
plot(H_img, main = "Hue")
plot(S_img, main = "Saturation")
plot(V_img, main = "Value")



#####
# Munsell colors
library(munsell)
library(colorspace)


image_array <- as.array(segmented_image1)
rgb_matrix <- col2rgb(as.raster(segmented_image1))#cbind(as.vector(image_array[,,,1]), as.vector(image_array[,,,2]), as.vector(image_array[,,,3]))
set.seed(42)

# Apply K-means clustering with 2 clusters (2 dominant colors)
kmeans_result <- kmeans(t(rgb_matrix), centers = 2)

# Get the cluster assignment for each pixel
clustered_pixels <- kmeans_result$cluster
# Replace each pixel in the image with the centroid color for its cluster
clustered_image <- t(rgb_matrix)
clustered_image[clustered_pixels == 1, ] <- kmeans_result$centers[1, ]  # First color
clustered_image[clustered_pixels == 2, ] <- kmeans_result$centers[2, ]  # Second color

# Reshape the clustered image matrix back into an image array (original dimensions)
height <- dim(image_array)[1]
width <- dim(image_array)[2]
clustered_image_array <- array(clustered_image, dim = c(height, width, 3))

# Convert the array back into a cimg object
clustered_image_cimg <- as.cimg(clustered_image_array)

# Plot the segmented image
plot(clustered_image_cimg)



compressed_lanczos <- imresize(segmented_image2, scale = 1/80, interpolation = "lanczos")
plot(compressed_lanczos)



rgb_values <- col2rgb(as.raster(compressed_lanczos))
munsell_col <- rgb2mnsl(rgb_values[1,], rgb_values[2,], rgb_values[3,])
table(munsell_col)


#####
# Statistical Feature Extraction

SFE <- function(seg_im) {
  # Step 1: Extract pixel values from the segmented area
  # Convert the segmented image to a grayscale image (if it's not already)
  segmented_gray <- grayscale(seg_im)
  
  # Extract pixel values as a vector (flattened matrix) only for non-zero pixels (segmented area)
  # We exclude background pixels (which are zero after masking) from the calculations
  segmented_pixels <- as.numeric(segmented_gray[segmented_gray > 0])
  
  # Step 2: Calculate statistical features
  
  # Mean
  mean_value <- mean(segmented_pixels)
  
  # Median
  median_value <- median(segmented_pixels)
  
  # Variance
  variance_value <- var(segmented_pixels)
  
  # Skewness
  skewness_value <- skewness(segmented_pixels)
  
  # Kurtosis
  kurtosis_value <- kurtosis(segmented_pixels)
  
  # Step 3: Display the calculated statistical features
  cat("Statistical Features of Segmented Image:\n")
  cat("Mean: ", mean_value, "\n")
  cat("Median: ", median_value, "\n")
  cat("Variance: ", variance_value, "\n")
  cat("Skewness: ", skewness_value, "\n")
  cat("Kurtosis: ", kurtosis_value, "\n")
  return(c(mean_value,median_value, variance_value, skewness_value,kurtosis_value))
}

s1<-SFE(segmented_image1)
scon<-SFE(control_im)
s2<-SFE(segmented_image2) #light soil

dist(rbind(s1,s2), method = "euclidean")
dist(rbind(s1,scon), method = "euclidean")

#####
# Entropy

Entropy_im <- function(seg_im) {
  # Step 1: Convert the segmented image to grayscale (if it's not already in grayscale)
  segmented_gray <- grayscale(seg_im)
  
  # Step 2: Extract pixel values for the segmented area (excluding background pixels)
  segmented_pixels <- as.numeric(segmented_gray[segmented_gray > 0])
  
  # Step 3: Calculate the histogram of pixel intensities for the segmented area
  # Define the number of bins for the histogram (e.g., 256 for 8-bit grayscale images)
  num_bins <- 256
  
  # Calculate histogram and normalize it to obtain probability distribution
  pixel_hist <- hist(segmented_pixels, breaks = num_bins, plot = FALSE)
  pixel_prob <- pixel_hist$counts / sum(pixel_hist$counts)  # Normalize to get probabilities
  
  # Step 4: Calculate Entropy
  # Using entropy package to compute entropy based on the pixel probability distribution
  entropy_value <- entropy(pixel_prob, unit = "log2")  # Use log base 2 for bits
  
  # Step 5: Display the calculated entropy
  cat("Entropy of the Segmented Image: ", entropy_value, " bits\n")
}

Entropy_im(segmented_image1)
Entropy_im(segmented_image2)  #light soil
Entropy_im(control_im)


#####
# Color hist


color_hist <- function(seg_im) {
  # Step 1: Split the segmented image into R, G, and B channels
  # 'segmented_image' should be in color format; if not, load it in color.
  R_channel <- R(seg_im)
  G_channel <- G(seg_im)
  B_channel <- B(seg_im)
  
  # Step 2: Extract pixel values for each channel for the segmented area (excluding background pixels)
  # Convert to vectors and filter out background pixels (assuming 0 represents background)
  R_pixels <- as.numeric(R_channel[R_channel > 0])
  G_pixels <- as.numeric(G_channel[G_channel > 0])
  B_pixels <- as.numeric(B_channel[B_channel > 0])
  
  # Step 3: Calculate color histograms for each channel
  # Define the number of bins for the histogram (e.g., 256 for 8-bit images)
  num_bins <- 25
  breaks <- seq(0, 1, length.out = num_bins + 1)  # Bins range from 0 to 1
  
  # Calculate histograms for each color channel
  R_hist <- hist(R_pixels, breaks = breaks, plot = FALSE)
  G_hist <- hist(G_pixels, breaks = breaks, plot = FALSE)
  B_hist <- hist(B_pixels, breaks = breaks, plot = FALSE)
  
  # Normalize the histograms to get probability distributions
  R_prob <- R_hist$counts / sum(R_hist$counts)
  G_prob <- G_hist$counts / sum(G_hist$counts)
  B_prob <- B_hist$counts / sum(B_hist$counts)
  
  # Step 4: Plot the color histograms for each channel
  par(mfrow = c(1, 3))
  barplot(R_prob, col = "red", main = "Red Channel Histogram", xlab = "Intensity", ylab = "Frequency")
  barplot(G_prob, col = "green", main = "Green Channel Histogram", xlab = "Intensity", ylab = "Frequency")
  barplot(B_prob, col = "blue", main = "Blue Channel Histogram", xlab = "Intensity", ylab = "Frequency")
  
  # Step 5: Combine the histograms into a feature vector (if needed for analysis)
  color_histogram_vector <- c(R_prob, G_prob, B_prob)
  
  # Display the combined feature vector
  cat("Color Histogram Feature Vector (R, G, B):\n", color_histogram_vector)
  return(color_histogram_vector)
}

chv1<-color_hist(segmented_image1)
chv2<-color_hist(segmented_image2) #light soil
chvcon <-color_hist(control_im)

dist(rbind(chv1,chv2), method = "euclidean") # dist dependant on bins
dist(rbind(chv1,chvcon), method = "euclidean")


#####
# Color Moments


color_moments <- function(seg_im) {
  # Step 1: Split the segmented image into R, G, and B channels
  # 'segmented_image' should be in color format; if not, load it in color.
  R_channel <- R(seg_im)
  G_channel <- G(seg_im)
  B_channel <- B(seg_im)
  
  # Step 2: Extract pixel values for each channel for the segmented area (excluding background pixels)
  # Convert to vectors and filter out background pixels (assuming 0 represents background)
  R_pixels <- as.numeric(R_channel[R_channel > 0])
  G_pixels <- as.numeric(G_channel[G_channel > 0])
  B_pixels <- as.numeric(B_channel[B_channel > 0])
  
  # Step 3: Calculate color moments for each channel
  
  # Red Channel Moments
  R_mean <- mean(R_pixels)
  R_variance <- var(R_pixels)
  R_skewness <- skewness(R_pixels)
  
  # Green Channel Moments
  G_mean <- mean(G_pixels)
  G_variance <- var(G_pixels)
  G_skewness <- skewness(G_pixels)
  
  # Blue Channel Moments
  B_mean <- mean(B_pixels)
  B_variance <- var(B_pixels)
  B_skewness <- skewness(B_pixels)
  
  # Step 4: Combine moments into a feature vector
  color_moments_vector <- c(R_mean, R_variance, R_skewness,
                            G_mean, G_variance, G_skewness,
                            B_mean, B_variance, B_skewness)
  
  # Step 5: Display the calculated color moments
  cat("Color Moments Feature Vector (Mean, Variance, Skewness for R, G, B):\n")
  print(color_moments_vector)
}

c1<-color_moments(segmented_image1)
c2<-color_moments(segmented_image2)  # light soil
c3 <-color_moments(control_im)

dist(rbind(c1,c2), method = "euclidean")  # "minkowski"  "manhattan"
dist(rbind(c1,c3), method = "euclidean")



#####
# Color Auto-correologram

col_autocorr <- function(seg_im, max_dist = 5, num_bins = 16) {
  # Step 1: Split the color channels into R, G, and B components
  R_channel <- R(seg_im)
  G_channel <- G(seg_im)
  B_channel <- B(seg_im)
  
  # Step 2: Define a function to calculate the auto-correlogram for each channel
  calculate_autocorrelogram <- function(channel, max_dist =max_dist , num_bins = num_bins) {
    # Get the dimensions of the channel
    channel_dim <- dim(channel)
    height <- channel_dim[1]
    width <- channel_dim[2]
    
    # Quantize the channel into 'num_bins' bins
    quantized <- as.integer(floor(channel * num_bins))
    
    # Reshape the quantized array into a 2D matrix
    quantized <- matrix(quantized, nrow = height, ncol = width)
    
    # Initialize the correlogram vector for storing results
    correlogram <- numeric(num_bins)
    
    # Loop over each pixel and calculate color correlation at various distances
    for (d in 1:max_dist) {
      for (i in 1:height) {
        for (j in 1:width) {
          current_value <- quantized[i, j]
          
          # Check the neighbors at distance 'd'
          if (i + d <= height && quantized[i + d, j] == current_value) correlogram[current_value + 1] <- correlogram[current_value + 1] + 1
          if (j + d <= width && quantized[i, j + d] == current_value) correlogram[current_value + 1] <- correlogram[current_value + 1] + 1
          if (i - d > 0 && quantized[i - d, j] == current_value) correlogram[current_value + 1] <- correlogram[current_value + 1] + 1
          if (j - d > 0 && quantized[i, j - d] == current_value) correlogram[current_value + 1] <- correlogram[current_value + 1] + 1
        }
      }
    }
    
    # Normalize the correlogram
    correlogram <- correlogram / sum(correlogram)
    return(correlogram)
  }
  
  # Step 3: Calculate the auto-correlogram for each color channel
  R_correlogram <- calculate_autocorrelogram(R_channel)
  G_correlogram <- calculate_autocorrelogram(G_channel)
  B_correlogram <- calculate_autocorrelogram(B_channel)
  
  # Step 4: Combine the auto-correlograms into a single vector
  color_correlogram_vector <- c(R_correlogram, G_correlogram, B_correlogram)
  
  # Display the result
  cat("Color Auto-Correlogram Vector:\n")
  print(color_correlogram_vector)
}

auto1<-col_autocorr(segmented_image1)*1000
auto2<-col_autocorr(segmented_image2)*1000  # light soil
auto3<-col_autocorr(control_im)*1000

dist(rbind(auto1,auto2), method = "euclidean")  # "minkowski"  "manhattan"
dist(rbind(auto1,auto3), method = "euclidean")



#####
# Hu Moments

## image shape dependant
## Grayscale is possible

hu_moments <- function(seg_im) {
  # Step 1: Convert the segmented image to grayscale (if it's not already in grayscale)
  segmented_gray <- grayscale(seg_im)
  
  # Step 2: Binarize the segmented image to create a binary mask (1 for object, 0 for background)
  #binary_image <- segmented_gray > 0.5  # Adjust threshold as needed
  
  binary_matrix <- as.matrix(segmented_gray[,,1,1])
  
  # Step 3: Calculate Image Moments
  calc_moments <- function(binary_image) {
    # Image dimensions
    h <- dim(binary_image)[1]
    w <- dim(binary_image)[2]
    
    # Initialize moments
    m00 <- sum(binary_image)  # Zeroth moment (area)
    m10 <- sum(row(binary_image) * binary_image)
    m01 <- sum(col(binary_image) * binary_image)
    
    # Centroid
    x_bar <- m10 / m00
    y_bar <- m01 / m00
    
    # Central moments
    mu20 <- sum(((row(binary_image) - x_bar)^2) * binary_image)
    mu02 <- sum(((col(binary_image) - y_bar)^2) * binary_image)
    mu11 <- sum((((binary_image) - x_bar) * (col(binary_image) - y_bar)) * binary_image)
    mu30 <- sum(((row(binary_image) - x_bar)^3) * binary_image)
    mu03 <- sum(((col(binary_image) - y_bar)^3) * binary_image)
    mu21 <- sum(((row(binary_image) - x_bar)^2 * (col(binary_image) - y_bar)) * binary_image)
    mu12 <- sum((((binary_image) - x_bar) * (col(binary_image) - y_bar)^2) * binary_image)
    
    # Return moments
    list(mu20 = mu20, mu02 = mu02, mu11 = mu11, mu30 = mu30, mu03 = mu03, mu21 = mu21, mu12 = mu12)
  }
  
  
  # Step 4: Calculate Hu Moments based on central moments
  calc_hu_moments <- function(moments) {
    # Extract moments
    mu20 <- moments$mu20
    mu02 <- moments$mu02
    mu11 <- moments$mu11
    mu30 <- moments$mu30
    mu03 <- moments$mu03
    mu21 <- moments$mu21
    mu12 <- moments$mu12
    
    # Calculate Hu moments using the formulas
    hu1 <- mu20 + mu02
    hu2 <- (mu20 - mu02)^2 + 4 * mu11^2
    hu3 <- (mu30 - 3 * mu12)^2 + (3 * mu21 - mu03)^2
    hu4 <- (mu30 + mu12)^2 + (mu21 + mu03)^2
    hu5 <- (mu30 - 3 * mu12) * (mu30 + mu12) * ((mu30 + mu12)^2 - 3 * (mu21 + mu03)^2) + 
      (3 * mu21 - mu03) * (mu21 + mu03) * (3 * (mu30 + mu12)^2 - (mu21 + mu03)^2)
    hu6 <- (mu20 - mu02) * ((mu30 + mu12)^2 - (mu21 + mu03)^2) + 4 * mu11 * (mu30 + mu12) * (mu21 + mu03)
    hu7 <- (3 * mu21 - mu03) * (mu30 + mu12) * ((mu30 + mu12)^2 - 3 * (mu21 + mu03)^2) - 
      (mu30 - 3 * mu12) * (mu21 + mu03) * (3 * (mu30 + mu12)^2 - (mu21 + mu03)^2)
    
    # Return Hu moments
    c(hu1, hu2, hu3, hu4, hu5, hu6, hu7)
  }
  
  # Step 5: Calculate Moments and Hu Moments for the segmented image
  moments <- calc_moments(binary_matrix)
  hu_moments <- calc_hu_moments(moments)
  
  # Step 6: Display the Hu Moments
  cat("Hu Moments for the Segmented Image:\n")
  print(hu_moments)
  
  
  
  # Apply log transformation to Hu Moments for better interpretation
  normalize_hu_moments <- function(hu_moments) {
    sapply(hu_moments, function(x) {
      if (x != 0) {
        -sign(x) * log10(abs(x))
      } else {
        0
      }
    })
  }
  
  # Normalize the Hu Moments
  normalized_hu_moments <- normalize_hu_moments(hu_moments)
  
  # Display the normalized Hu Moments
  cat("Normalized Hu Moments for the Segmented Image:\n")
  print(normalized_hu_moments)
  return(normalized_hu_moments)
}

h1 <- hu_moments(segmented_image1)
h2 <- hu_moments(segmented_image2)   # light soil
hcon <- hu_moments(control_im)

dist(rbind(h1,h2), method = "euclidean")  #"minkowski"   "manhattan"
dist(rbind(h1,hcon), method = "euclidean") # dist dependant on shape


#####
# Aspect Ratio and Extent

## Binary image is used

# Step 1: Convert the segmented image to grayscale (if it's not already in grayscale)
segmented_gray <- grayscale(segmented_image)

# Step 2: Binarize the segmented image to create a binary mask (1 for object, 0 for background)
binary_image <- segmented_gray > 0.5  # Adjust threshold as needed

# Step 3: Calculate the bounding box of the segmented area
# Get the non-zero pixels' coordinates (object region)
coords <- which(binary_image == TRUE, arr.ind = TRUE)

# Calculate bounding box dimensions
x_min <- min(coords[, 1])
x_max <- max(coords[, 1])
y_min <- min(coords[, 2])
y_max <- max(coords[, 2])

# Width and Height of the Bounding Box
bounding_box_width <- x_max - x_min + 1
bounding_box_height <- y_max - y_min + 1

# Step 4: Calculate Aspect Ratio
aspect_ratio <- bounding_box_width / bounding_box_height

# Step 5: Calculate Extent
# Area of the object (number of non-zero pixels)
object_area <- sum(binary_image)

# Area of the bounding box
bounding_box_area <- bounding_box_width * bounding_box_height

# Calculate Extent
extent <- object_area / bounding_box_area

# Step 6: Display the Aspect Ratio and Extent
cat("Aspect Ratio of the Segmented Object: ", aspect_ratio, "\n")
cat("Extent of the Segmented Object: ", extent, "\n")



######
# Compactness and Roundness

## Binary image is used

# Step 1: Convert the segmented image to grayscale (if it's not already in grayscale)
segmented_gray <- grayscale(segmented_image)

# Step 2: Binarize the segmented image to create a binary mask (1 for object, 0 for background)
binary_image <- segmented_gray > 0.5  # Adjust threshold as needed

# Step 3: Calculate the area of the segmented object
# Area is the number of non-zero pixels
object_area <- sum(binary_image)


# Step 4: Calculate the perimeter of the segmented object
# Use the contours function to extract the contour of the segmented object
contour_list <- contours(binary_image)
contour_points <- contour_list[[1]]

# Calculate perimeter using the Euclidean distance between consecutive contour points
object_perimeter <- sum(sqrt(diff(contour_points$x)^2 + diff(contour_points$y)^2))

# Step 5: Calculate Compactness
compactness <- (object_perimeter^2) / (4 * pi * object_area)

# Step 6: Calculate Roundness using the convex hull to find the major axis length
# Compute convex hull to find the minimum bounding box
library(grDevices)

# Extract unique coordinates of the segmented region
coords <- which(binary_image == TRUE, arr.ind = TRUE)

# Compute the convex hull
convex_hull <- chull(coords)
hull_coords <- coords[convex_hull, ]

# Calculate the major axis length as the maximum Euclidean distance between any two points on the convex hull
major_axis_length <- max(dist(hull_coords))

# Step 7: Calculate Roundness
roundness <- (4 * object_area) / (pi * major_axis_length^2)

# Step 8: Display the Compactness and Roundness
cat("Compactness of the Segmented Object: ", compactness, "\n")
cat("Roundness of the Segmented Object: ", roundness, "\n")



######
# Zernike Moments

## robust against shape

## Grayscale is enough


Zernike_moments <- function(seg_im) {
  # Step 1: Convert the segmented image to grayscale (if it's not already in grayscale)
  segmented_gray <- grayscale(seg_im)
  
  # Step 2: Binarize the segmented image to create a binary mask (1 for object, 0 for background)
  #binary_image <- segmented_gray > 0.5  # Adjust threshold as needed
  
  binary_matrix <- as.matrix(segmented_gray[,,1,1])
  
  # Step 3: Define functions to compute Zernike polynomials and moments
  
  # Zernike Polynomial Calculation
  zernike_polynomial <- function(n, m, rho, theta) {
    if ((n - abs(m)) %% 2 != 0) return(0)  # Condition for Zernike polynomial
    
    radial_poly <- numeric(length(rho))
    for (s in 0:((n - abs(m)) / 2)) {
      c <- (-1)^s * factorial(n - s) /
        (factorial(s) * factorial((n + abs(m)) / 2 - s) * factorial((n - abs(m)) / 2 - s))
      radial_poly <- radial_poly + c * rho^(n - 2 * s)
    }
    
    return(radial_poly * exp(1i * m * theta))
  }
  
  # Zernike Moment Calculation
  zernike_moment <- function(n, m, image) {
    #image <- as.cimg(image)
    # Convert image to polar coordinates centered at the image center
    cx <- dim(image)[1] / 2
    cy <- dim(image)[2] / 2
    radius <- min(cx, cy)  # Unit circle radius
    
    # Prepare the polar coordinate grid
    x <- row(image) - cx
    y <- col(image) - cy
    rho <- sqrt(x^2 + y^2) / radius
    theta <- atan2(y, x)
    
    # Consider only points within the unit circle
    mask <- rho <= 1
    rho <- rho[mask]
    theta <- theta[mask]
    image_values <- as.numeric(image)[mask]
    
    # Compute the Zernike polynomial
    Znm <- zernike_polynomial(n, m, rho, theta)
    
    # Compute the Zernike moment
    Vnm <- sum(image_values * Znm) * (n + 1) / pi
    
    return(Vnm)
  }
  
  # Step 4: Calculate Zernike Moments for the segmented image
  # Example for Zernike moments of order 0 <= n <= 4
  zernike_moments <- list()
  for (n in 0:4) {
    for (m in seq(-n, n, 2)) {
      moment <- zernike_moment(n, m, binary_matrix)
      zernike_moments[[paste0("Z_", n, "_", m)]] <- moment
    }
  }
  
  # Step 5: Display the Zernike Moments
  cat("Zernike Moments for the Segmented Image:\n")
  print(zernike_moments)
  
  
  # Calculate magnitudes of Zernike moments
  zernike_magnitudes <- sapply(zernike_moments, function(z) Mod(z))
  
  # Display magnitudes of Zernike moments
  cat("Magnitudes of Zernike Moments:\n")
  print(zernike_magnitudes)
  return(zernike_magnitudes)
}

z1 <- Zernike_moments(segmented_image1)
z2 <- Zernike_moments(segmented_image2)  #light soil
zcon <- Zernike_moments(control_im) #imrotate(control_im,90)

dist(rbind(z1,z2), method = "euclidean")
dist(rbind(z1,zcon), method = "euclidean") 



######
# Local Binary Patterns (LBP)

LBP <- function(seg_im) {
  # Assuming 'segmented_image' is already obtained and is in grayscale
  
  # Step 1: Convert the image to grayscale and then to a 2D matrix
  grayscale_image <- grayscale(seg_im)
  
  # Convert the grayscale image to a 2D matrix
  grayscale_matrix <- as.matrix(grayscale_image[,,1,1])
  
  # Step 2: Define a function to calculate LBP for a given pixel neighborhood
  calculate_lbp <- function(matrix, x, y) {
    center_value <- matrix[x, y]  # Get the center pixel value
    # Extract the 3x3 neighborhood around the pixel (ignoring boundary pixels)
    neighborhood <- matrix[(x - 1):(x + 1), (y - 1):(y + 1)]  # 3x3 neighborhood
    
    # Compare the center pixel with its neighbors
    binary_pattern <- as.integer(neighborhood >= center_value)
    
    # Ensure binary_pattern is treated as a matrix
    binary_pattern <- matrix(binary_pattern, nrow = 3, ncol = 3)
    
    # Remove the center pixel from the comparison (set it to 0)
    binary_pattern[2, 2] <- 0
    
    # Convert the binary pattern to a vector in clockwise order starting from the top-left
    binary_pattern_vector <- c(binary_pattern[1, 1], binary_pattern[1, 2], binary_pattern[1, 3],
                               binary_pattern[2, 3], binary_pattern[3, 3], binary_pattern[3, 2],
                               binary_pattern[3, 1], binary_pattern[2, 1])
    
    # Calculate the LBP code by converting binary to decimal
    lbp_code <- sum(2^(0:7) * binary_pattern_vector)
    return(lbp_code)
  }
  
  # Step 3: Apply the LBP function to each pixel in the image
  # Initialize an empty matrix to store the LBP image
  lbp_image <- matrix(0, nrow = nrow(grayscale_matrix), ncol = ncol(grayscale_matrix))
  
  # Loop through each pixel, ignoring the border pixels
  for (x in 2:(nrow(grayscale_matrix) - 1)) {
    for (y in 2:(ncol(grayscale_matrix) - 1)) {
      lbp_image[x, y] <- calculate_lbp(grayscale_matrix, x, y)
    }
  }
  
  
  # Step 4: Convert the LBP matrix to an image format for visualization
  lbp_image_cimg <- as.cimg(lbp_image)
  
  # Step 5: Display the original and LBP images for comparison
  par(mfrow = c(1, 2))
  plot(grayscale_image, main = "Original Grayscale Image")
  plot(lbp_image_cimg, main = "LBP Image")
  
  # Step 6: Calculate the histogram of LBP codes for texture classification
  lbp_histogram <- hist(lbp_image, breaks = seq(0, 256, by = 1), plot = FALSE)$counts
  
  # Display the LBP histogram
  cat("LBP Histogram:\n")
  print(lbp_histogram) # vector of 256 numbers
  return(list(lbp_histogram, lbp_image))
}

lbp1<-LBP(segmented_image1[,,,-4, drop = FALSE])
lbp2<-LBP(original_image2)  #light soil
lbp_con <- LBP(control_im)

dist(rbind(lbp1,lbp2), method = "euclidean")
dist(rbind(lbp1,lbp_con), method = "euclidean")


######
## mean shift segmentation
library(reshape2)
library(ggplot2)
library(dbscan)

preprocess_lbp_image <- function(lbp_image, spatial_weight = 1.0, texture_weight = 2.0) {
  # Ensure the input is a matrix
  if (!is.matrix(lbp_image)) {
    stop("lbp_image must be a 2D matrix.")
  }
  
  # Get the dimensions of the image
  nrow_img <- nrow(lbp_image)
  ncol_img <- ncol(lbp_image)
  
  # Generate pixel coordinates:
  # x: column indices (repeated for each row)
  # y: row indices (repeated for each column)
  x_coords <- rep(1:ncol_img, each = nrow_img)
  y_coords <- rep(1:nrow_img, times = ncol_img)
  
  # Flatten the LBP image to a vector
  lbp_values <- as.vector(lbp_image)
  
  # Normalize the spatial coordinates to the range [0, 1]
  x_norm <- x_coords / ncol_img
  y_norm <- y_coords / nrow_img
  
  # Apply the weights
  x_weighted <- x_norm * spatial_weight
  y_weighted <- y_norm * spatial_weight
  lbp_weighted <- lbp_values * texture_weight
  
  # Combine into a feature matrix: each row is (x, y, lbp)
  data_mat <- cbind(x = x_weighted, y = y_weighted, lbp = lbp_weighted)
  
  return(data_mat)
}

# Example usage:
# Suppose lbp_image is your LBP image matrix, you can create data_mat as follows:
# data_mat <- preprocess_lbp_image(lbp_image, spatial_weight = 1.0, texture_weight = 2.0)
data_mat <- preprocess_lbp_image(lbp1[[2]])

# --- 3. Implement a Basic Mean Shift Clustering Function ---
mean_shift_cluster_fast <- function(data, bandwidth, max_iter = 10, tol = 1e-3) {
  n <- nrow(data)
  shifted <- data
  # We define a search radius (e.g., 3 times the bandwidth) to restrict neighbor search
  search_radius <- 3 * bandwidth
  
  for (i in 1:n) {
    point <- data[i, ]
    for (iter in 1:max_iter) {
      # Find neighbors within the search radius using dbscan's fast radius search
      nn <- frNN(x = data, eps = search_radius, query = matrix(point, nrow = 1))
      idx <- nn$id[[1]]
      if (length(idx) == 0) break  # no neighbors found
      
      # Compute differences and squared distances for only the nearby points
      diffs <- data[idx, , drop = FALSE] - matrix(point, nrow = length(idx), ncol = ncol(data), byrow = TRUE)
      distances_sq <- rowSums(diffs^2)
      weights <- exp(-distances_sq / (2 * bandwidth^2))
      
      # Compute the weighted mean
      new_point <- colSums(data[idx, , drop = FALSE] * weights) / sum(weights)
      if (sqrt(sum((new_point - point)^2)) < tol) {
        point <- new_point
        break
      }
      point <- new_point
    }
    shifted[i, ] <- point
  }
  
  # --- Group modes that are close ---
  cluster_assignment <- rep(0, n)
  cluster_id <- 0
  for (i in 1:n) {
    if (cluster_assignment[i] == 0) {
      cluster_id <- cluster_id + 1
      cluster_assignment[i] <- cluster_id
      for (j in (i + 1):n) {
        if (cluster_assignment[j] == 0 && sqrt(sum((shifted[i, ] - shifted[j, ])^2)) < (bandwidth / 2)) {
          cluster_assignment[j] <- cluster_id
        }
      }
    }
  }
  
  list(assignment = cluster_assignment, modes = shifted)
}

# Example usage:
# Suppose data_mat is the feature matrix (from your pre-processing function)
# bandwidth <- 0.1  # (adjust based on your data scale)
# ms_result <- mean_shift_cluster_fast(data = data_mat, bandwidth = bandwidth)
# assignments <- ms_result$assignment


# --- 5. Apply Mean Shift Clustering ---
# Choose an appropriate bandwidth; you may need to experiment with this value.
bandwidth <- 0.1

#ms_result <- mean_shift_cluster(data = data_mat, bandwidth = bandwidth)
#assignments <- ms_result$assignment

mean_shift_cluster_subsample <- function(data, bandwidth, subsample_fraction = 0.1, max_iter = 10, tol = 1e-3) {
  n <- nrow(data)
  # Subsample indices
  sample_indices <- sample(1:n, size = floor(n * subsample_fraction))
  data_sample <- data[sample_indices, , drop = FALSE]
  
  # Run mean shift on the subsample using the fast version
  shifted_sample <- data_sample
  search_radius <- 3 * bandwidth
  for (i in 1:nrow(data_sample)) {
    point <- data_sample[i, ]
    for (iter in 1:max_iter) {
      nn <- frNN(x = data_sample, eps = search_radius, query = matrix(point, nrow = 1))
      idx <- nn$id[[1]]
      if (length(idx) == 0) break
      diffs <- data_sample[idx, , drop = FALSE] - matrix(point, nrow = length(idx), ncol = ncol(data), byrow = TRUE)
      distances_sq <- rowSums(diffs^2)
      weights <- exp(-distances_sq / (2 * bandwidth^2))
      new_point <- colSums(data_sample[idx, , drop = FALSE] * weights) / sum(weights)
      if (sqrt(sum((new_point - point)^2)) < tol) {
        point <- new_point
        break
      }
      point <- new_point
    }
    shifted_sample[i, ] <- point
  }
  
  # Group similar modes in the subsample
  cluster_assignment_sample <- rep(0, nrow(data_sample))
  cluster_id <- 0
  for (i in 1:nrow(data_sample)) {
    if (cluster_assignment_sample[i] == 0) {
      cluster_id <- cluster_id + 1
      cluster_assignment_sample[i] <- cluster_id
      for (j in (i + 1):nrow(data_sample)) {
        if (cluster_assignment_sample[j] == 0 &&
            sqrt(sum((shifted_sample[i, ] - shifted_sample[j, ])^2)) < (bandwidth / 2)) {
          cluster_assignment_sample[j] <- cluster_id
        }
      }
    }
  }
  
  # Assign each point in the original data to the nearest mode from the subsample
  assignments <- rep(0, n)
  for (i in 1:n) {
    dists <- sqrt(rowSums((shifted_sample - matrix(data[i, ], nrow = nrow(shifted_sample), ncol = ncol(data), byrow = TRUE))^2))
    nearest <- which.min(dists)
    assignments[i] <- cluster_assignment_sample[nearest]
  }
  
  list(assignment = assignments,
       sample_modes = shifted_sample,
       sample_indices = sample_indices,
       sample_clusters = cluster_assignment_sample)
}

# Example usage:
ms_result_sub <- mean_shift_cluster_subsample(data = data_mat, bandwidth = 5.5, subsample_fraction = 0.2)
assignments <- ms_result_sub$assignment

ncol_img = dim(lbp1[[2]])[2]
nrow_img = dim(lbp1[[2]])[1]

# --- 6. Reshape and Visualize the Segmented Image ---
segmented_image <- matrix(assignments, nrow = nrow_img, ncol = ncol_img, byrow = FALSE)

# Option 1: Basic plot using base R
image(segmented_image, col = rainbow(max(assignments) + 1), 
      main = "Segmented Image (Mean Shift)", xlab = "X", ylab = "Y")


# K-means
segment_kmeans <- function(img, k = 3, iter.max = 100, nstart = 5) {
  # Vectorize the image: each pixel becomes an observation.
  # If your image has more channels (or additional features), you can modify this step.
  img_vec <- as.vector(img)
  
  # Run k-means clustering on the pixel intensities.
  # For grayscale or LBP images, the intensity (or LBP value) is the feature.
  km <- kmeans(img_vec, centers = k, iter.max = iter.max, nstart = nstart)
  
  # Reshape the cluster assignments back into the image dimensions.
  segmentation <- matrix(km$cluster, nrow = nrow(img), ncol = ncol(img))
  
  return(segmentation)
}

K <- 2  # For example, we want to segment the image into 3 regions.
segmented_image <- segment_kmeans(lbp1[[2]], k = K, iter.max = 10)

image(segmented_image, col = rainbow(K), main = "K-Means Segmentation", 
      xlab = "X", ylab = "Y")

#####
# Gray-Level Co-Occurrence Matrix (GLCM) Features

## Robust with shapes

## Grayscale , and uses mean

library(glcm)

GLCM <- function(seg_im) {
  # Assuming 'segmented_image' is already obtained and is in grayscale
  
  # Step 1: Convert the image to grayscale if it's not already
  grayscale_image <- grayscale(seg_im)
  
  # Convert the grayscale image to a matrix
  grayscale_matrix <- as.matrix(grayscale_image[,,1,1])
  
  # Step 2: Calculate the GLCM for different directions and distances
  # Using the glcm() function from the 'glcm' package
  
  glcm_result <- glcm(grayscale_matrix, 
                      window = c(5, 5),  # GLCM computed in windows of size 5x5 (change as needed)
                      shift = c(0, 1))   # Shift specifies direction and distance (here, 0° and distance 1)
  
  
  # Step 3: Extract GLCM Features MEAN
  glcm_features_summary <- list(
    mean = mean(glcm_result[,, "glcm_mean"], na.rm = TRUE),
    variance = mean(glcm_result[,, "glcm_variance"], na.rm = TRUE),
    homogeneity = mean(glcm_result[,, "glcm_homogeneity"], na.rm = TRUE),
    contrast = mean(glcm_result[,, "glcm_contrast"], na.rm = TRUE),
    dissimilarity = mean(glcm_result[,, "glcm_dissimilarity"], na.rm = TRUE),
    entropy = mean(glcm_result[,, "glcm_entropy"], na.rm = TRUE),
    second_moment = mean(glcm_result[,, "glcm_second_moment"], na.rm = TRUE),
    correlation = mean(glcm_result[,, "glcm_correlation"], na.rm = TRUE)
  )
  
  
  # Display GLCM features
  cat("GLCM Features Summary (Mean Values):\n")
  print(glcm_features_summary)
  
  
  
  glcm_mean_image <- as.cimg(glcm_result[, , "glcm_mean"])
  glcm_variance_image <- as.cimg(glcm_result[, , "glcm_variance"])
  glcm_homogeneity_image <- as.cimg(glcm_result[, , "glcm_homogeneity"])
  glcm_contrast_image <- as.cimg(glcm_result[, , "glcm_contrast"])
  glcm_dissimilarity_image <- as.cimg(glcm_result[, , "glcm_dissimilarity"])
  glcm_entropy_image <- as.cimg(glcm_result[, , "glcm_entropy"])
  glcm_second_moment_image <- as.cimg(glcm_result[, , "glcm_second_moment"])
  glcm_correlation_image <- as.cimg(glcm_result[, , "glcm_correlation"])
  
  # Step 4: Plot the original grayscale image and the GLCM feature images
  par(mfrow = c(3, 3))  # 3x3 layout for plots
  
  plot(grayscale_image, main = "Original Grayscale Image")
  plot(glcm_mean_image, main = "GLCM Mean")
  plot(glcm_variance_image, main = "GLCM Variance")
  plot(glcm_homogeneity_image, main = "GLCM Homogeneity")
  plot(glcm_contrast_image, main = "GLCM Contrast")
  plot(glcm_dissimilarity_image, main = "GLCM Dissimilarity")
  plot(glcm_entropy_image, main = "GLCM Entropy")
  plot(glcm_second_moment_image, main = "GLCM Second Moment")
  plot(glcm_correlation_image, main = "GLCM Correlation")
  
  return(list(glcm_result, glcm_features_summary))
}

g1 <- GLCM(segmented_image1)
g2 <- GLCM(segmented_image2)  # light soil
gcon <- GLCM(control_im ) # imrotate(control_im,90)  

dist(rbind(g1[2][[1]],g2[2][[1]]), method = "euclidean")
dist(rbind(g1[2][[1]],gcon[2][[1]]), method = "euclidean")


#####
# Haralick Features

## robust against shape

# Step 3: Define a function to compute Haralick features from GLCM
compute_haralick_features <- function(glcm) {
  # Initialize a list to store the Haralick features
  haralick_features <- list()
  
  # Angular Second Moment (ASM) or Energy
  haralick_features$ASM <- mean(glcm[, , "glcm_second_moment"], na.rm = TRUE)
  
  # Contrast
  haralick_features$Contrast <- mean(glcm[, , "glcm_contrast"], na.rm = TRUE)
  
  # Correlation
  haralick_features$Correlation <- mean(glcm[, , "glcm_correlation"], na.rm = TRUE)
  
  # Variance
  haralick_features$Variance <- mean(glcm[, , "glcm_variance"], na.rm = TRUE)
  
  # Inverse Difference Moment (Homogeneity)
  haralick_features$Homogeneity <- mean(glcm[, , "glcm_homogeneity"], na.rm = TRUE)
  
  # Sum Average (sum of averages)
  sum_avg <- glcm[, , "glcm_mean"] + glcm[, , "glcm_mean"]
  haralick_features$SumAverage <- mean(sum_avg, na.rm = TRUE)
  
  # Sum Variance
  sum_var <- apply(glcm[, , "glcm_variance"], c(1, 2), function(x) var(x, na.rm = TRUE))
  haralick_features$SumVariance <- mean(sum_var, na.rm = TRUE)
  
  # Sum Entropy
  sum_entropy <- -sum(glcm[, , "glcm_entropy"] * log(glcm[, , "glcm_entropy"] + 1e-10), na.rm = TRUE)
  haralick_features$SumEntropy <- mean(sum_entropy, na.rm = TRUE)
  
  # Entropy
  haralick_features$Entropy <- mean(glcm[, , "glcm_entropy"], na.rm = TRUE)
  
  # Difference Variance
  diff_var <- apply(glcm[, , "glcm_contrast"], c(1, 2), function(x) var(x, na.rm = TRUE))
  haralick_features$DifferenceVariance <- mean(diff_var, na.rm = TRUE)
  
  # Difference Entropy
  diff_entropy <- -sum(glcm[, , "glcm_entropy"] * log(glcm[, , "glcm_entropy"] + 1e-10), na.rm = TRUE)
  haralick_features$DifferenceEntropy <- mean(diff_entropy, na.rm = TRUE)
  
  return(haralick_features)
}

# Step 4: Compute Haralick features for the GLCM result
(haralick_features1 <- compute_haralick_features(g1[1][[1]]))
(haralick_features2 <- compute_haralick_features(g2[1][[1]]))
(haralick_features3 <- compute_haralick_features(gcon[1][[1]]))

dist(rbind(haralick_features1,haralick_features2), method = "euclidean")
dist(rbind(haralick_features1,haralick_features3), method = "euclidean")

