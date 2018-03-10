#### Incidence Rates of Craigslist Rental Listing -----------------------------
#### GBH Spatial Epidemiology Project

#### Preamble -----------------------------------------------------------------

#packages
library(tidyverse)
library(sqldf)
library(ggthemes)
library(RColorBrewer)
library(viridis)
library(sp)
library(spdplyr)
library(rgdal)
library(rgeos)
library(INLA)

#setwd to location of file (REQUIRES RSTUDIO)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load 2012-2016 ACS data (2016 tracts same as 2010 for KC)
census <- read_csv(file = "./input/nhgis0051_ds225_20165_2016_tract.csv")

#read in tract shapefile for king county
kc_shp <- readOGR(dsn = "./input/KCTract2010/KCTract2010.shp",
                  layer = "KCTract2010",
                  GDAL1_integer64_policy = TRUE,
                  stringsAsFactors = F)

#determine config (dev or github)
if(file.exists("../data/cl/craigslistDB.sqlite")){
  DB <- dbConnect(SQLite(), dbname="../data/cl/craigslistDB.sqlite")
  cl <- tbl(DB, "clean") #clean listing table
  
  #compute tract aggregates for CL listing count
  tractCl <- cl %>%
    filter(!is.na(GISJOIN), !is.na(cleanBeds), !is.na(cleanRent), !is.na(cleanSqft)) %>% #only listings with valid Bed/Rent
    filter(GISJOIN %in% kc_shp@data$GISJOIN) %>% #filter to KC only (db has metro area)
    dplyr::select(listingMoYr, GISJOIN, seattle, matchAddress, matchAddress2, matchType, cleanBeds, cleanRent) %>% #SELECT these columns
    collect %>% #bring db query into memory
    filter(!grepl("Google", matchType)) %>% #no Google geocodes, only Smartystreets (precise to Zip9)
    distinct(matchAddress, matchAddress2, cleanBeds, cleanRent, cleanSqft, .keep_all = T) %>% #dedupe unique address-bed-rent combos
    group_by(GISJOIN) %>% #group listings by tract
    summarize(nListings = n(),
              medCL = median(cleanRent),
              med0B = median(cleanRent[cleanBeds==0]),
              med1B = median(cleanRent[cleanBeds==1]),
              med2B = median(cleanRent[cleanBeds==2]),
              med3Plus = median(cleanRent[cleanBeds>2]),
              seattle = max(seattle))
  dbDisconnect(DB)
  write_csv(tractCl, "./input/tractCl.csv")
  
} else{
  tract <- read_csv("./input/tractCl.csv")
}



#check listing var
summary(tractCl$nListings)


#### A. Prep for `test` table for analyses -----------------------------------------

#merge tables together, allow missing values for tracts with no listings (i.e. no value in tract)
test <- left_join(census, tractCl)

#compute some new vars/mutate sensibly named ACS vars
test <- test %>%
  mutate(nListings = ifelse(is.na(nListings), 0, nListings), #if NA then 0, since no listings were observed
         nHU = AF7NE001,
         GISJOIN = as.character(GISJOIN)) %>%
  mutate(tpop = AF2LE001,
         tnhw = AF2UE003,
         tnhb = AF2UE004,
         tnha = AF2UE006,
         thsp = AF2UE012,
         tnho =  AF2UE005+AF2UE007+AF2UE008+AF2UE009,
         nvachu = AF7OE003,
         nforrent = AF7ZM002,
         nocc =  AF7OE002,
         ppov = (AF43E002+AF43E003)/AF43E001,
         nownocc = AF7PE002,
         medGRentACS = AF89M001,
         medStrYr = AF8IE001,
         medVal = AF9LE001/1000,
         medHHInc = AF49M001/1000) %>%
  mutate(pnhw = tnhw/tpop,
         pnhb = tnhb/tpop,
         pnho = tnho/tpop,
         pnha = tnha/tpop,
         phsp = thsp/tpop,
         pvac = nvachu/nHU,
         pforrent = ifelse(nvachu == 0, 0, AF7ZM002/nvachu),
         pownocc = nownocc/nocc,
         pblt14lat = AF8HE002/nHU) %>%
  filter(tpop > 0)


#### B. Descriptive Analyses --------------------------------------------------

pal <- brewer.pal(7, "Purples")

#craigslist listing count density
ggplot(test, aes(x = nListings)) +
  geom_density(fill = pal[7]) +
  theme_minimal() +
  xlab("\nNumber of Craigslist Listings") +
  ylab("Density\n") +
  labs(title = "Density of Tract N CL Listings") +
  ggsave(filename = "./output/graphics/nListingDensity.pdf",
         width = 8, height = 6, dpi = 300)

#nHU ACS density
ggplot(test, aes(x = nHU)) +
  geom_density(fill = pal[7]) +
  theme_minimal() +
  xlab("\nNumber of Housing Units 2012-2016") +
  ylab("Density\n") +
  labs(title = "Density of Total Housing Units in Tract") +
  ggsave(filename = "./output/graphics/nHUDensity.pdf",
         width = 8, height = 6, dpi = 300)

#median HH income
ggplot(test, aes(x = medHHInc)) +
  geom_density(fill = pal[7]) +
  theme_minimal() +
  xlab("\nMedian Household Income") +
  ylab("Density\n") +
  labs(title = "Density of Median HH Income") +
  ggsave(filename = "./output/graphics/medHHIncDensity.pdf",
         width = 8, height = 6, dpi = 300)

#racial and ethnic composition
blk <- test %>% 
  mutate(race = "Black",
         percent = pnhb) %>%
  dplyr::select(race, percent, GISJOIN)

hsp <- test %>% 
  mutate(race = "Hispanic",
         percent = phsp) %>%
  dplyr::select(race, percent, GISJOIN)

asi <- test %>% 
  mutate(race = "Asian",
         percent = pnha) %>%
  dplyr::select(race, percent, GISJOIN)

wht <- test %>% 
  mutate(race = "White",
         percent = pnhw) %>%
  dplyr::select(race, percent, GISJOIN)

oth <- test %>%
  mutate(race = "Other",
         percent = pnho) %>%
  dplyr::select(race, percent, GISJOIN)

race_dens <- bind_rows(asi, blk, hsp, oth, wht)
race_dens$race <- as.factor(race_dens$race)

ggplot(race_dens, aes(x = race, y = percent, group = race)) + 
  geom_boxplot(color = pal[7]) +
  theme_minimal() +
  ggsave(filename = "./output/graphics/raceBoxplot.pdf",
         width = 8, height = 6, dpi = 300)

#### D. Spatial Visualizations ------------------------------------------------

options(scipen = 99)

kc_shp@data$id <- rownames(kc_shp@data)
kc_shp@data <- left_join(kc_shp@data, test)

kc_f <- fortify(kc_shp)
kc_f <- inner_join(kc_f, kc_shp@data, "id")

ggplot(kc_f, aes(x = long, y = lat, group = group, fill = nHU)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="N HUs", palette = "Blues", 
                       breaks = scales::pretty_breaks(n = 5), direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/ACS_nHU.pdf",
         dpi = 300)

ggplot(kc_f, aes(x = long, y = lat, group = group, fill = nforrent)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="N HU For Rent", palette = "Blues", 
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/ACS_forrent.pdf",
         dpi = 300)


ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(nListings))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="N Listings", palette = "Blues", 
                       breaks = c(0, 1, 2, 3),
                       labels = c("1", "10", "100", "1000"),
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/log10nListings.pdf",
         dpi = 300)


# make grid graphic for different compositions
blk_f <- kc_f %>% 
  mutate(race = "Non-Hispanic Black",
         percent = pnhb) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

hsp_f <- kc_f %>% 
  mutate(race = "Hispanic",
         percent = phsp) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

asi_f <- kc_f %>% 
  mutate(race = "Non-Hispanic Asian",
         percent = pnha) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

wht_f <- kc_f %>% 
  mutate(race = "Non-Hispanic White",
         percent = pnhw) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

oth_f <- kc_f %>%
  mutate(race = "Non-Hispanic Other",
         percent = pnho) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

demog_f <- bind_rows(asi_f, blk_f, hsp_f, oth_f, wht_f)
demog_f$race <- as.factor(demog_f$race)
demog_f$race <- factor(demog_f$race, levels = levels(demog_f$race)[c(2, 3, 1, 4, 5)])

ggplot(demog_f, aes(x = long, y = lat, group = group, fill = percent)) +
  facet_wrap(~ race) +
  geom_polygon(color = "white", lwd = .1) +
  scale_fill_viridis_c(labels = scales::percent) +
  theme_map() +
  coord_map() +
  theme(strip.background = element_rect(fill = "white", color = "white")) +
  theme(strip.text = element_text(colour = "Black", size = 12)) +  
  theme(legend.position = c(.85, .1), legend.justification = c(1, 0)) +
  labs(fill = "Percent of Tract Population") +
  ggsave(filename = "./output/maps/ACS_ethrace.pdf",
         width = 12, height = 7, dpi = 300)

ggplot(wht_f, aes(x = long, y = lat, group = group, fill = percent)) +
  geom_polygon(color = "white", lwd = .1) +
  scale_fill_viridis_c(labels = scales::percent) +
  theme_map() +
  coord_map() +
  theme(strip.background = element_rect(fill = "white", color = "white")) +
  theme(strip.text = element_text(colour = "Black", size = 12)) +  
  theme(legend.position = c(1.05, .1), legend.justification = c(1, 0)) +
  labs(fill = "Percent Non-Hispanic White") +
  ggsave(filename = "./output/maps/ACS_nhwMap.pdf",
         width = 12, height = 7, dpi = 300)

ggplot(test, aes(x = pnhb, y = nListings)) + geom_smooth() + geom_point() +
  theme_minimal() +
  xlab("\n% Non-Hispanic White") +
  ylab("N Listings\n")
  ggsave(filename = "./output/graphics/nListingNHBscatter.pdf",
         dpi = 300)



#### E. INLA models -----------------------------------------------------------

#save(test, file = "R:/Project/seattle_rental_market/scripts/spatial_epi/working.RData")
#load("R:/Project/seattle_rental_market/scripts/spatial_epi/working.RData")

library(INLA)
  

  

#King County listing incidence ------------------------------------------------------------
kc_shp <- readOGR(dsn = "./input/KCTract2010/KCTract2010.shp",
                  layer = "KCTract2010",
                  GDAL1_integer64_policy = TRUE,
                  stringsAsFactors = F, verbose = F)

test <- test[match(kc_shp$GISJOIN, test$GISJOIN),]
test$GISJOIN == kc_shp$GISJOIN #dbl check order
rownames(test) <- row.names(kc_shp)

#spCbind the data frame onto shapefile
kc_shp <- maptools::spCbind(kc_shp, test) #both are complete tables of KC tracts
kc_df <- as.tbl(kc_shp@data) #save table (now in the order of the shapefile)
kc_adj <- poly2nb(kc_shp)
nb2INLA("./output/graphINLA/kctract.graph", kc_adj)


#create a couple numeric ID columns to use later for INLA 
test$idx <- 1:nrow(test) 
test$idxx <- test$idx

kc_df <- as.data.frame(kc_shp@data)

form2 <- nListings ~ 1 + log(tpop) + nHU + pvac + pownocc + pblt14lat +
  pnhb + phsp + pnha + pnho + medHHInc + seattle +
  f(idx, model = "iid")

m2 <- inla(form2, 
           family = "poisson", 
           data = test,
           control.predictor = list(compute = TRUE))
summary(m2)
plot(m2)

form3 <- nListings ~ 1 + log(tpop) + nHU + pvac + pownocc + pblt14lat +
  pnhb + phsp + pnha + pnho + medHHInc + seattle +
  f(idx, model = "iid") +
  f(idxx, model = "besag", graph = "./output/graphINLA/kctract.graph")

m3 <- inla(form3, 
           family = "poisson", 
           E = nHU, 
           data = kc_df,
           control.predictor = list(compute = TRUE))
summary(m3)
plot(m3)

kc_df$fittedListings <- m1$summary.fitted[, "0.5quant"]

kc_shp <- readOGR(dsn = "./input/KCTract2010/KCTract2010.shp",
                  layer = "KCTract2010",
                  GDAL1_integer64_policy = TRUE,
                  stringsAsFactors = F, verbose = F)
kc_shp <- kc_shp %>%
  inner_join(kc_df)

pdf(width = 8, height = 6, file = "./maps/smoothListing.pdf")
spplot(kc_shp, c("nListings", "fittedListings"),
       col.regions = brewer.pal(9, "Blues"), cuts = 8)
dev.off()








#Seattle listing incidence ----------------------------------------------------------

sea_shp <- readOGR(dsn = "./sea_tract_2010/sea_tract_2010.shp",
                   layer = "sea_tract_2010",
                   GDAL1_integer64_policy = TRUE,
                   stringsAsFactors = F)
sea_adj <- poly2nb(sea_shp)
nb2INLA("./graphINLA/seatract.graph", sea_adj)

sea_df <- as.data.frame(test) %>% filter(seattle == 1)
sea_df$idx <- 1:nrow(sea_df)
sea_df$idxx <- sea_df$idx

form3 <- nListings ~ offset(log(nHU)) +  tpop + pvac + medVal + pownocc + pblt14lat +
  pnhb + phsp + pnha + pnho + medHHInc +
  f(idx, model = "iid") +
  f(idxx, model = "besag", graph = "./graphINLA/seatract.graph")

sea_df <- as.data.frame(test)

m3 <- inla(form3, 
           family = "poisson", 
           data = sea_df,
           verbose = T,
           control.predictor = list(compute = TRUE))
summary(m3)

sea_df$logListings <- log(sea_df$nListings)
sea_df$smoothedCL <- log(m3$summary.fitted[, "mean"])

sea_shp <- sea_shp %>%
  inner_join(sea_df)

pdf(width = 8, height = 6, file = "R:/Project/seattle_rental_market/report/spatial_epi/maps/smoothSeaListing.pdf")
spplot(sea_shp, c("logListings", "smoothedCL"),
       col.regions = brewer.pal(9, "Blues"), cuts = 8)
dev.off()






#### J. Close out db Connection -----------------------------------------------

dbDisconnect(DB)












