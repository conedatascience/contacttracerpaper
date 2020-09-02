# Purpose: Generate Simulated Data to Mimic EMR Record Data
# packages ----------------------------------------------------------------

library(wakefield)
library(dplyr)
library(purrr)
library()
library(data.table)

# setup -------------------------------------------------------------------
set.seed(336)

# generate fake data ------------------------------------------------------
# Using Wakefield, generate fake patient population
df <- r_data_frame(
	n = 10500,
	x = rnorm,
	y = rnorm,
	race,
	id,
	age,
	sex,
	language,
)

# Add generic test dates
df$patient_id <- 1:nrow(df)
df$date <- sample(seq.Date(as.Date("2020-06-01"), 
																as.Date("2020-07-01"), by = 1),
														size = 10500, replace = T)

# clean data --------------------------------------------------------------
# 
df_2 <- df %>% 
	mutate(location = paste0(round(x,2),"-",round(y,2))) %>% 
	count(location, sort = T)


# connect cases -----------------------------------------------------------

#' Connect Probable Cases
#'
#' A helper function that can join cases together
#'
#' @param dat the dataframe of likely connect cases
#' @param weights_in the weights to use for the serial interval if available
#' @param exposure_link how these cases are connected
#'
#' @importFrom data.table .SD `:=`
#' @export

connect_probable_cases <- function(dat, weights_in = NULL, exposure_link = NULL){
	
	if(is.null(weights_in)){
		SI_param = epitrix::gamma_mucv2shapescale(4.7, 2.9/4.7)
		SI_distribution <- distcrete::distcrete("gamma", interval = 1,
																						shape = SI_param$shape,
																						scale = SI_param$scale, w = 0.5)
		
		w <- SI_distribution$d(1:21)
	}
	
	dat <- dat[order(dat$date),]
	
	dat$id <- 1:nrow(dat)
	dates <- dplyr::pull(dat, date)
	
	# If all the dates are equivalent, the make sure that
	if(length(unique(dates))==1){
		dates[1] <- dates[1]-1
	}
	
	cluster <- outbreaker2::outbreaker_data(dates = dates, w_dens = w)
	
	res <- outbreaker2::outbreaker(cluster)
	
	out <- summary(res)
	
	tree <- data.table::setDT(out$tree)
	
	# If non of the cases can be linked, will assign the oldest case as the index
	# case. This isn't super elegent, but it keeps the cluster together.
	if(sum(!is.na(tree$from))==0){
		tree$from[2:nrow(tree)] <-1L
	}
	
	
	tree <- tree[time<=30, .(from, to)]
	
	tree <- merge(tree, dat[ ,c("id", "patient_id", "date")], by.x = "from", by.y="id", all.x = TRUE)
	
	tree <-tree[,.(patient_id,to, date)]
	
	names(tree) <- c("from", "to", "from_date")
	
	tree <- merge(tree, dat[ ,c("id", "patient_id", "date")], by.x = "to", by.y="id", all.x = TRUE)
	
	tree <- tree[,.(from, patient_id, from_date, date)]
	
	names(tree) <- c("from", "to", "from_date", "to_date")
	
	tree <- tree[,diff_dt_onset := abs(as.numeric(to_date-from_date))]
	tree <- tree[!is.na(from), .(from,to, diff_dt_onset)]
	
	if(!is.null(exposure_link)){
		tree <- cbind(tree, exposure_link = exposure_link)
	}
	# Returning as DF to avoid conflicts with purrr
	as.data.frame(tree)
	
}

contact_matrix <- df %>%
	mutate(location = paste0(round(x,2),"-",round(y,2))) %>% 
	group_by(location) %>%
	mutate(ct = n()) %>%
	filter(ct>1) %>%
	arrange(-ct) %>% 
	dplyr::select(-ct) %>%
	nest() %>%
	.[1:10,] %>% 
	mutate(connections = map(data, connect_probable_cases, exposure_link = location)) %>%
	dplyr::select(connections) %>%
	unnest(connections)

library(epicontacts)

x <- epicontacts::make_epicontacts(linelist = df,
																	 contacts = contact_matrix)
plot(x)
