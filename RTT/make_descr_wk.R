library(dplyr)
library(lubridate)
library(ggplot2)
library(COVID19)


# Load in Sequence Description Data

desc <- read.csv("sequences_descr.csv")
desc$date <- ymd(desc$CollectionDate)

# The dates that failed to parse only have YYYY-MM, no day
# For our purposes, we can ignore these.
#desc[is.na(desc$date), c("CollectionDate", "date")]
desc <- filter(desc, !is.na(date), tolower(Host) == "homo sapiens")

desc$yr <- year(desc$date)
desc$wk <- week(desc$date)

# Load in COVID19 Cases

covid <- covid19(verbose = FALSE)
covid$yr <- as.numeric(year(covid$date))
covid$wk <- as.numeric(week(covid$date))

covid$yw <- covid$yr + covid$wk/52

covid2 <- covid %>%
	group_by(yr, wk) %>% 
	summarise(new = sum(confirmed, na.rm = TRUE), 
		recovered = sum(recovered, na.rm = TRUE),
		.groups = "drop") %>% 
	mutate(yw = yr + wk/52)

if (FALSE) {
	ggplot(covid2) + 
		aes(
			x = as.numeric(yw), 
			y = new - recovered
			) +
		geom_col()
}

desc2 <- left_join(desc, covid2, by = c("yr", "wk"))

#head(desc2)

# Save this modification (will be useful)
write.csv(desc2, file = "sequences_descr_wk.csv", row.names = FALSE)










