# Purpose: Generate Simulated Data to Mimic EMR Record Data


# packages ----------------------------------------------------------------
#install.packages("wakefield")
library(wakefield)
library(tidyverse)

# setup -------------------------------------------------------------------
set.seed(336)

# generate fake data ------------------------------------------------------
# Using Wakefield, generate fake patient population
df <- r_data_frame(
	n = 5500,
	x = rnorm,
	y = rnorm,
	race,
	id,
	age,
	sex,
	language,
	date_stamp
)

# clean data --------------------------------------------------------------
# 
df_2 <- df %>% 
	mutate(location = paste0(round(x,2),"-",round(y,2))) %>% 
	count(location, sort = T)
				 