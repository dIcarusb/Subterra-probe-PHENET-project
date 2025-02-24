# Build new experiment data file
# This script extract the data from raw h5 files, and adds the data into a 
# new h5 file, which has the following hierarchy: 

# Experiment <- file name
# Description (metadata) [attr]
#   Plot name
#     Insertion point name
#     White and dark reference
#     Wavelength names
#     Date [attr]
#       GPS [attr]
#       Depth
#       Force [data]
#       Spectra1 [data] (raw)
#       Spectra2 [data] (raw)
#       Water content [data]
#       Bulk [data]
#       Other [data]

# The whole script has several steps which should be followed precisely.
# However, some steps can be re-used or re-run several times, at exception
# of the "Step 0" which is the creation of the new file


##### 
# Step 0
# Create new h5 file, where to store the experiment's data

library(hdf5r)
setwd("D:/Soil_spectra_h5files") # set directory

filename <- "AGES_lab.h5"
file.h5 <- H5File$new(filename, mode = "w") # "w" the first time, next time use "r+"

#####
# Add metadata description to the file
description <- "This file contains experimental data for the AGES project. It structures in the following way: 	
Plot  > Insertion point, GPS, white, dark, wavelength names, Date > Depth, Force data, Spectra data (raw), Water content data, Bulk data, Other data ."
file.h5$create_attr(attr_name = "Description", robj = description)
## To check the it:
h5attributes(file.h5)
## If wrong, for deletion use this line:
file.h5$attr_delete("Description")


# Add plots
file.h5$create_group("RT_I")
file.h5$create_group("NT_I")
file.h5$create_group("CT_I")
file.h5$create_group("NT_II")
file.h5$create_group("RT_II")
file.h5$create_group("CT_II")
file.h5$create_group("RT_III")
file.h5$create_group("CT_III")
file.h5$create_group("NT_III")
file.h5$create_group("grassland")
## To check it:
file.h5
list.groups(file.h5)
## To delete groups:
file.h5$link_delete("plotA")

## to delete attributes:
file.h5$attr_delete("bla bla")


#####
# Step 1
# Read and extract the data from raw h5 files
f <- H5File$new("lab_ages.h5", mode = "r") #ages_gl
f

all_groups <- list.groups(f)
(spectrometer1_groups <- grep("^session\\d{3}/cal\\d{3}/ins\\d{3}/spectrometer1$", all_groups, value = TRUE))
(spectrometer2_groups <- grep("^session\\d{3}/cal\\d{3}/ins\\d{3}/spectrometer2$", all_groups, value = TRUE))

# select session
(spectrometer1_groups<-spectrometer1_groups[1:16])
(spectrometer2_groups<-spectrometer2_groups[63:77])

# example
f[[spectrometer1_groups[7]]]
f[[spectrometer1_groups[7]]][["spectrum_depths"]][]
f[[spectrometer1_groups[7]]][["spectrum_forces"]][]

#extract data, be careful with session and call numbers!
white <- h5attr(f[["session001"]][["cal001"]],"white_spectrum")
dark <- h5attr(f[["session001"]][["cal001"]],"dark_spectrum")
white2 <- h5attr(f[["session001"]][["cal001"]],"white_spectrum2")
dark2 <- h5attr(f[["session001"]][["cal001"]],"dark_spectrum2")

wavelengths_name <- h5attr(f[["session001"]][["cal001"]],"spec1_wavelengths_vector")
wavelengths_name2 <- h5attr(f[["session001"]][["cal001"]],"spec2_wavelengths_vector")


# How it works, not necesary to do
# ins<-1 # insertion 1
# # depths
# depths<-f[[spectrometer1_groups[ins]]][["spectrum_depths"]][]
# # forces
# forces<-f[[spectrometer1_groups[ins]]][["spectrum_forces"]][]
# # spectra
# # first spectra
# VISNIR_spectra <- t(f[[spectrometer1_groups[ins]]][["spectra"]][,])
# colnames(VISNIR_spectra) <- wavelengths_name
# 
# # second spectra
# NIRSWIR_spectra2 <- t(f[[spectrometer2_groups[ins]]][["spectra"]][,])
# colnames(NIRSWIR_spectra2) <- wavelengths_name2


#####
# Step 2
# Write new data into the new h5 experiment file

# plot information
# references
file.h5[["grassland"]]$create_attr(attr_name = "white_ref", robj = white)
file.h5[["grassland"]]$create_attr(attr_name = "white_ref2", robj = white2)
file.h5[["grassland"]]$create_attr(attr_name = "dark_ref", robj = dark)
file.h5[["grassland"]]$create_attr(attr_name = "dark_ref2", robj = dark2)
# wavelengths names
file.h5[["grassland"]]$create_attr(attr_name = "wavelengths_names", robj = wavelengths_name)
file.h5[["grassland"]]$create_attr(attr_name = "wavelengths_names2", robj = wavelengths_name2)
# date dd-mm-yyyy
file.h5[["grassland"]]$create_attr(attr_name = "date", robj = "2024-05-07") 

# check
h5attr_names(file.h5[["grassland"]])
h5attributes(file.h5[["RT_III"]])

# Insertion 1 , shows how it works, not necesary to run
# path <- "grassland/ins001"
# file.h5$create_group(path)
# file.h5[[path]]$create_dataset(name="depths", robj=depths)
# file.h5[[path]]$create_dataset(name="forces", robj=forces)
# file.h5[[path]]$create_dataset(name="VISNIR_spectra", 
#                                              robj=VISNIR_spectra)
# file.h5[[path]]$create_dataset(name="NIRSWIR_spectra2", 
#                                              robj=NIRSWIR_spectra2)
# WATER CONTENT !!!!!!
# file.h5[[path]]$create_dataset(name="water_content", robj=water_content)
# SOC !!!!!!!!
# file.h5[[path]]$create_dataset(name="SOC", robj=SOC)
# GPS !!!!!!!!!!
# file.h5[[path]]$create_attr(attr_name = "GPS", robj = GPS)
# file.h5[[path]]


# load the functions
source("Useful_functions.R")
# show insertions depths
sess01 <- spectrometer1_groups
depths_table <- get_insertion_depths(sess01, f)
depths_table


#filter insertions
(cleaned_groups<-clean_insertions(spectrometer1_groups, spectrometer2_groups,c(3)))
# save and keep the changes within the same variables
spectrometer1_groups <- cleaned_groups$spectrometer1_groups
spectrometer2_groups <- cleaned_groups$spectrometer2_groups

## Here is where it runs in a loop for every insertion
#loop


process_insertions(file.h5, seq(1,17,1), "grassland/", spectrometer1_groups, 
                   spectrometer2_groups, wavelengths_name, 
                   wavelengths_name2, order=TRUE, lab = TRUE) # order if false is from 15 to 1, if true is from 1 to 15

# check the addings
file.h5[["grassland"]]

#####
# Step 3
# Close h5 file, after everything is done
h5close(file.h5)


