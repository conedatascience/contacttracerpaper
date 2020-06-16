# Purpose: Generate Simulated Data to Mimic EMR Record Data


# packages ----------------------------------------------------------------
#install.packages("wakefield")
library(wakefield)
library(tidyverse)

# generate fake data ------------------------------------------------------
# Using Wakefield, generate fake patient population
df <- r_data_frame(
	n = 500,
	x = rnorm,
	y = rnorm,
	race,
	id,
	race,
	age,
	sex,
	language,
	date_stamp
)

# clean data --------------------------------------------------------------
# 

