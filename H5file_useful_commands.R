# Check the structure of the h5 file
library(hdf5r)
setwd("D:/Soil_spectra_h5files") # your directory
filename <- "AT_lab.h5" #AGES_exp
file.h5 <- H5File$new(filename, mode = "r")

# print the file, it shows the attributes and groups that it contains
file.h5

# list all the groups
list.groups(file.h5)
# list all the attributes and what they contain
h5attributes(file.h5)


# Example with plot CT_II
file.h5[["session002"]] # print everything
list.groups(file.h5[["session007"]]) # list groups inside the plot
h5attr_names(file.h5[["session011"]]) # list attributes names
h5attributes(file.h5[["session001"]]) # list attributes and what they contain
h5attr(f[["session001"]], "session_date") # read what a particular attribute contains
# Now with the insertion:
file.h5[["session001/cal001/ins001"]]
list.groups(file.h5[["session001/cal001/ins001"]]) # NO groups
h5attr_names(file.h5[["session001/cal001/ins001"]]) # NO attributes
# Insertions contain data, you can access it like this:
depths <- file.h5[["CT_II/ins001/depths"]][] # for 1D data
spectra <- file.h5[["CT_II/ins001/NIRSWIR_spectra2"]][,] # for 2D data




