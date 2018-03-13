#packages
library(tidyverse) 
library(lubridate)
library(sqldf)
library(ggthemes)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(sp)
library(spdplyr)
library(rgdal)
library(rgeos)
library(spdep)
library(INLA)
library(coefINLA) #devtools::install_github("hesscl/coefINLA")
library(zoo)
library(xts)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


DB <- dbConnect(SQLite(), dbname="../data/cl/craigslistDB.sqlite")
cl <- tbl(DB, "clean") #clean listing table

#read in tract shapefile for king county
kc_shp <- readOGR(dsn = "./input/KCTract2010/KCTract2010.shp",
                  layer = "KCTract2010",
                  GDAL1_integer64_policy = TRUE,
                  stringsAsFactors = F)

#compute tract aggregates for CL listing count
kc_cl <- cl %>%
  filter(!is.na(GISJOIN), !is.na(cleanBeds), !is.na(cleanRent), !is.na(cleanSqft)) %>% #only listings with valid Bed/Rent
  filter(GISJOIN %in% kc_shp@data$GISJOIN) %>% #filter to KC only (db has metro area)
  dplyr::select(listingMoYr, listingDate, GISJOIN, seattle, matchAddress, matchAddress2, matchType, cleanBeds, cleanRent) %>% #SELECT these columns
  collect %>% #bring db query into memory
  mutate(listingMoYr = as.Date(listingMoYr),
         listingDate = as.Date(listingDate)) %>%
  filter(listingMoYr < "2018-03-12") %>% #everything up till march 12, 2018  (i.e cut off new data)
  filter(!grepl("Google", matchType)) %>% #no Google geocodes, only Smartystreets (precise to Zip9)
  distinct(matchAddress, matchAddress2, cleanBeds, cleanRent, cleanSqft, .keep_all = T)
dbDisconnect(DB)


glimpse(kc_cl)
length(unique(kc_cl$listingDate))

computeCal <- function(x){
  y <- x %>%
    group_by(listingDate) %>%
    summarize(n = n())
  
  fullDates <- seq(from = min(y$listingDate), to = max(y$listingDate), by = 1)
  index <- seq(from = 1, length.out = length(fullDates), by = 1)
  dateTable <- cbind.data.frame(fullDates, index)
  
  n <- NULL
  for(i in 1:length(fullDates)){
    if(fullDates[i] %in% y$listingDate){
      n[i] <- y$n[y$listingDate == fullDates[i]]
    } else{
      n[i] <- NA
    }
  }
  n <- cbind.data.frame(n, fullDates)
}


cal <- computeCal(kc_cl) %>%
  mutate(weekdayf = weekdays(fullDates) %>% as.factor,
         monthf = month(fullDates, label = T),
         yearmon = as.yearmon(fullDates),
         year = year(fullDates),
         week = as.numeric(format(fullDates,"%W")))

cal$weekdayf <- factor(cal$weekdayf, levels = rev(levels(cal$weekdayf)[c(2,6,7,5,1,3,4)]))

#credit to https://rpubs.com/haj3/calheatmap for code and idea
cal<-plyr::ddply(cal,plyr::.(yearmon),transform,monthweek=1+week-min(week))


ggplot(cal, aes(x = monthweek, y = weekdayf, fill = n)) +
  scale_fill_viridis_c(na.value = "grey80", direction = 1) +
  scale_x_continuous(breaks = c(1, 2, 3, 4)) +
  scale_size(range=c(5,20))+
  geom_tile(color = "white") +
  facet_grid(year~monthf) +
  theme_minimal() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        panel.spacing = unit(.1, "lines"),
        axis.text = element_text(size = 6),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        strip.text = element_text(size = 8),
        panel.grid = element_blank()) +
  xlab("\nWeek of Month") +
  ylab("Day of Week\n") +
  labs(fill = "N Listings") +
  ggsave(filename = "./output/graphics/timeseries.pdf",
         width = 6, height = 3)

mean(cal$n, na.rm=T)
median(cal$n, na.rm=T)
