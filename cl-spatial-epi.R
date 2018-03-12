#### Frequency of Craigslist Rental Listings ----------------------------------
#### GBH Spatial Epidemiology Project

#### Preamble -----------------------------------------------------------------

#packages
library(tidyverse) 
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

#setwd to location of file (REQUIRES RSTUDIO)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load 2012-2016 ACS data (2016 tracts same as 2010 for KC)
census <- read_csv(file = "./input/nhgis0052_ds225_20165_2016_tract.csv")

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
         nforborn = AF95E005,
         nvachu = AF7OE003,
         nforrent = AF7ZM002,
         nocc =  AF7OE002,
         ppov = ((AF43E002+AF43E003)/AF43E001),
         nownocc = AF7PE002,
         medGRentACS = AF89M001,
         medStrYr = AF8IE001,
         medVal = AF9LE001,
         medHHInc = AF49M001) %>%
  mutate(pnhw = (tnhw/tpop),
         pnhb = (tnhb/tpop),
         pnho = (tnho/tpop),
         pnha = (tnha/tpop),
         phsp = (thsp/tpop),
         pforborn = (nforborn/tpop),
         pvac = (nvachu/nHU),
         pforrent = ifelse(nvachu == 0, 0, (AF7ZM002/nvachu)),
         pownocc = (nownocc/nocc),
         pblt14lat = (AF8HE002/nHU)) %>%
  filter(tpop > 0)


#check listing var
summary(test$nListings) #overdispersed
sum(test$nListings == 0) #not zero-inflated


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

#median HH income
ggplot(test, aes(x = pforborn)) +
  geom_histogram(fill = pal[7], bins = 20) +
  scale_x_continuous(labels = scales::percent) +
  theme_minimal() +
  xlab("\nPercent Foreign-Born") +
  ylab("Frequency\n") +
  labs(title = "Histogram of Percent Foreign Born") +
  ggsave(filename = "./output/graphics/forbornHist.pdf",
         width = 8, height = 6, dpi = 300)

#racial and ethnic composition
blk <- test %>% 
  mutate(race = "Black",
         percent = pnhb) %>%
  dplyr::select(race, percent, nListings, GISJOIN)

hsp <- test %>% 
  mutate(race = "Hispanic",
         percent = phsp) %>%
  dplyr::select(race, percent, nListings, GISJOIN)

asi <- test %>% 
  mutate(race = "Asian",
         percent = pnha) %>%
  dplyr::select(race, percent, nListings, GISJOIN)

wht <- test %>% 
  mutate(race = "White",
         percent = pnhw) %>%
  dplyr::select(race, percent, nListings, GISJOIN)

oth <- test %>%
  mutate(race = "Other",
         percent = pnho) %>%
  dplyr::select(race, percent, nListings, GISJOIN)

race_dens <- bind_rows(asi, blk, hsp, oth, wht)
race_dens$race <- as.factor(race_dens$race)

ggplot(race_dens, aes(x = race, y = percent, group = race)) + 
  geom_boxplot(color = pal[7]) +
  theme_minimal() +
  ggsave(filename = "./output/graphics/raceBoxplot.pdf",
         width = 8, height = 6, dpi = 300)

ggplot(race_dens, aes(x = percent, y = nListings, group = race, color = race)) + 
  geom_smooth(se = F) + 
  scale_x_continuous(labels = scales::percent) +
  theme_minimal() +
  xlab("\n% of Neighborhood Population") +
  ylab("N Listings\n") +
  ggsave(filename = "./output/graphics/nListingEthRace.pdf",
       width = 11, height = 8.5, dpi = 300)

ggplot(test, aes(x = medHHInc, y = nListings)) + 
  geom_point(alpha = .25) +
  geom_smooth(se = F) + 
  theme_minimal() +
  xlab("\nMedian HH Income") +
  ylab("N Listings\n") +
  ggsave(filename = "./output/graphics/nListingAvgInc.pdf",
         width = 8, height = 6, dpi = 300)

ggplot(test, aes(x = ppov, y = nListings)) +
  geom_point(alpha = .25) +
  geom_smooth(se = F) + 
  theme_minimal() +
  xlab("\n% of Population Below Poverty Line") +
  ylab("N Listings\n") +
  ggsave(filename = "./output/graphics/nListingPov.pdf",
         width = 8, height = 6, dpi = 300)

ggplot(test, aes(x = pforborn, y = nListings)) +
  geom_point(alpha = .25) +
  geom_smooth(se = F) + 
  theme_minimal() +
  xlab("\n% of Population Foreign-Born") +
  ylab("N Listings\n") +
  ggsave(filename = "./output/graphics/nListingForborn.pdf",
         width = 8, height = 6, dpi = 300)


#### C. Spatial Visualizations ------------------------------------------------

options(scipen = 99)

kc_shp@data$id <- rownames(kc_shp@data)
kc_shp@data <- left_join(kc_shp@data, test)

kc_f <- fortify(kc_shp)
kc_f <- inner_join(kc_f, kc_shp@data, "id")

ggplot(kc_f, aes(x = long, y = lat, group = group, fill = nHU)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="N HUs", palette = "Purples", 
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
  scale_fill_distiller(name="N HU For Rent", palette = "Purples", 
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
  scale_fill_distiller(name="N Listings", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
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

ggplot(kc_f, aes(x = long, y = lat, group = group, fill = medHHInc)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="Median HH Income", palette = "Purples", direction = 1) + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/medHHInc.pdf",
         dpi = 300)

ggplot(kc_f, aes(x = long, y = lat, group = group, fill = ppov)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="% of Persons below Poverty Line",
                       labels = scales::percent, palette = "Purples", direction = 1) + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/ppov.pdf",
         dpi = 300)

ggplot(kc_f, aes(x = long, y = lat, group = group, fill = pforborn)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="% of Persons Foreign-Born",
                       labels = scales::percent, palette = "Purples", direction = 1) + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/pforborn.pdf",
         dpi = 300)

ggplot(kc_f, aes(x = long, y = lat, group = group, fill = pforborn)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name="% Persons Foreign-Born",
                       labels = scales::percent, palette = "Purples", direction = 1) + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/pforborn.pdf",
         dpi = 300)

# make grid graphic for different compositions
blk_f <- kc_f %>% 
  mutate(race = "Black",
         percent = pnhb) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

hsp_f <- kc_f %>% 
  mutate(race = "Hispanic",
         percent = phsp) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

asi_f <- kc_f %>% 
  mutate(race = "Asian",
         percent = pnha) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

wht_f <- kc_f %>% 
  mutate(race = "White",
         percent = pnhw) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

oth_f <- kc_f %>%
  mutate(race = "Other",
         percent = pnho) %>%
  dplyr::select(long, lat, group, race, percent, nListings)

demog_f <- bind_rows(asi_f, blk_f, hsp_f, oth_f)
demog_f$race <- as.factor(demog_f$race)

ggplot(demog_f, aes(x = long, y = lat, group = group, fill = percent)) +
  facet_wrap(~ race) +
  geom_polygon(color = "white", lwd = .1) +
  scale_fill_viridis_c(labels = scales::percent) +
  theme_map() +
  coord_map() +
  theme(strip.background = element_rect(fill = "white", color = "white")) +
  theme(strip.text = element_text(colour = "Black", size = 12)) +  
  theme(legend.position = c(1.1, .1), legend.justification = c(1, 0)) +
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
  labs(fill = "Percent of Neighborhood Population") +
  ggsave(filename = "./output/maps/ACS_nhw.pdf",
         width = 12, height = 7, dpi = 300)


#### D. INLA models -----------------------------------------------------------

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
idx <- 1:nrow(kc_df) 
idxx <- idx

#INLA seems to be finicky and want either option?
kc_df$idx <- 1:nrow(kc_df) 
kc_df$idxx <- kc_df$idx

#extract the dataframe to hand to INLA
kc_df <- as.data.frame(kc_shp@data)

### Model 0: Negative Binomial (no random effects)
form0 <- nListings ~ 1 + log(tpop) + nHU + pforrent + pownocc + pblt14lat +
  pnhb + phsp + pnha + pnho + medHHInc + ppov + pforborn + seattle

m0 <- inla(form0, 
           family = "nbinomial", 
           data = kc_df,
           control.predictor = list(compute = TRUE),
           control.compute = list(dic = TRUE, waic = TRUE))
summary(m0)
plot(m0)

#predicted values
m0median <- exp(m0$summary.linear.predictor["0.5quant"])
m0lower <- exp(m0$summary.linear.predictor["0.025quant"])
m0upper <- exp(m0$summary.linear.predictor["0.975quant"])
kc_shp@data$m0median <- m0median[,1]
kc_shp@data$m0lower <- m0lower[,1]
kc_shp@data$m0upper <- m0upper[,1]
kc_shp@data$m0post95wid <- kc_shp@data$m0upper - kc_shp@data$m0lower


### Model 1: Lognormal Non-Spatial Negative Binomial Model
form1 <- nListings ~ 1 + log(tpop) + nHU + pforrent + pownocc + pblt14lat +
  pnhb + phsp + pnha + pnho + medHHInc + ppov + pforborn + seattle +
  f(idx, model = "iid")

m1 <- inla(form1, 
           family = "nbinomial", 
           data = kc_df,
           control.predictor = list(compute = TRUE),
           control.compute = list(dic = TRUE, waic = TRUE))
summary(m1)
plot(m1)

#lognormal random effect medians
m1RE <- exp(m1$summary.random$idx[5])
kc_shp@data$m1RE <- m1RE[,1]

#predicted values
m1median <- exp(m1$summary.linear.predictor["0.5quant"])
m1lower <- exp(m1$summary.linear.predictor["0.025quant"])
m1upper <- exp(m1$summary.linear.predictor["0.975quant"])
kc_shp@data$m1median <- m1median[,1]
kc_shp@data$m1lower <- m1lower[,1]
kc_shp@data$m1upper <- m1upper[,1]
kc_shp@data$m1post95wid <- kc_shp@data$m1upper - kc_shp@data$m1lower



### Model 2: Lognormal Spatial Negative Binomial Model
form2 <- nListings ~ 1 + log(tpop) + nHU + pforrent + pownocc + pblt14lat +
  pnhb + phsp + pnha + pnho + medHHInc + ppov + pforborn + seattle +
  f(idx, model = "iid") +
  f(idxx, model = "besag", graph = "./output/graphINLA/kctract.graph")

m2 <- inla(form2, 
           family = "nbinomial", 
           data = kc_df,
           control.predictor = list(compute = TRUE),
           control.compute = list(dic = TRUE, waic = TRUE))
summary(m2)
plot(m2)

#lognormal random effect medians
m2RE <- exp(m2$summary.random$idx[5])
kc_shp@data$m2RE <- m2RE[,1]

#spatial effect median
m2SE <- exp(m2$summary.random$idxx[5])
kc_shp@data$m2SE <- m2SE[,1]

#predicted values
m2median <- exp(m2$summary.linear.predictor["0.5quant"])
m2lower <- exp(m2$summary.linear.predictor["0.025quant"])
m2upper <- exp(m2$summary.linear.predictor["0.975quant"])
kc_shp@data$m2median <- m2median[,1]
kc_shp@data$m2lower <- m2lower[,1]
kc_shp@data$m2upper <- m2upper[,1]
kc_shp@data$m2post95wid <- kc_shp@data$m2upper - kc_shp@data$m2lower

kc_shp@data$idx <- idx

#### E. Model Diagnostics -----------------------------------------------------

#gof
rmse <- function(error)
{
  sqrt(mean(error^2))
}

#marginal likelihood
m0mlik <- m0$mlik[1]
m1mlik <- m1$mlik[1]
m2mlik <- m2$mlik[1]

#DIC
m0DIC <- m0$dic$dic
m1DIC <- m1$dic$dic
m2DIC <- m2$dic$dic

#WAIC
m0WAIC <- m0$waic$waic
m1WAIC <- m1$waic$waic
m2WAIC <- m2$waic$waic

#RMSE
m0RMSE <- rmse(kc_shp@data$nListings-kc_shp@data$m0median)
m1RMSE <- rmse(kc_shp@data$nListings-kc_shp@data$m1median)
m2RMSE <- rmse(kc_shp@data$nListings-kc_shp@data$m2median)

#compile df
mfit <- data.frame(
  Model = c("GLM", "Non-Spatial RE", "Spatial RE"),
  MLik = c(m0mlik, m1mlik, m2mlik),
  RMSE = c(m0RMSE, m1RMSE, m2RMSE),
  DIC = c(m0DIC, m1DIC, m2DIC),
  WAIC = c(m0WAIC, m1WAIC, m2WAIC)
)

#save mfit to files
print(xtable::xtable(mfit), type = "html", "./output/modelFit.html")
print(xtable::xtable(mfit), type = "latex", "./output/modelFit.tex")

#for maps
kc_shp@data$id <- rownames(kc_shp@data)
kc_f <- fortify(kc_shp)
kc_f <- inner_join(kc_f, kc_shp@data, "id")

#labels for coefINLA
var_labs <- c(
  `(Intercept)`="Intercept",
  `log(tpop)`="log(Total Pop)",
  nHU="N HU",
  pforrent="Prp HU For Rent",
  pblt14lat="Prp HU Blt >2014",
  pownocc="Prp HU Own-Occ",
  pnhb="Prp Black",
  phsp="Prp Hispanic",
  pnha="Prp Asian",
  pnho="Prp Other",
  medHHInc="Med. HH Income",
  ppov="Poverty Rate",
  pforborn="Prp Foreign Born",
  seattle="Seattle"
)

#labeller for coefINLA
var_labeller <- labeller(
     var = var_labs
)

#coefINLA ropeladder-like plots
coefINLA(m0, exp = T, labeller = var_labeller) + 
  labs(title = "GLM Posterior Distributions",
       subtitle = "with 95% interval") +
  ggsave(filename = "./output/graphics/coefM0.pdf",
         width = 8.5, height = 11)

coefINLA(m1, exp = T, labeller = var_labeller) + 
  labs(title = "Non-Spatial RE Posterior Distributions",
       subtitle = "with 95% interval") +
ggsave(filename = "./output/graphics/coefM1.pdf",
       width = 8.5, height = 11)

coefINLA(m2, exp = T, labeller = var_labeller) + 
  labs(title = "Spatial RE Posterior Distributions",
       subtitle = "with 95% interval") +
ggsave(filename = "./output/graphics/coefM2.pdf",
       width = 8.5, height = 11)



### Spatial visualizations for models

#idx labels
ggplot(kc_f, aes(x = long, y = lat, group = group, label = idx)) +
  geom_polygon(fill = "grey90", color = "grey10") +
  geom_text(aes(x = type.convert(INTPTLON10), y = type.convert(INTPTLAT10)), size = 1) +
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/idxLabeled.pdf",
         width = 11, height = 8.5, dpi = 300)

#median fitted value
ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m0median))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Posterior Prediction\nN Listings (Non-Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/fittedM0.pdf",
         dpi = 300)

#median fitted value (for faceted plot)
m0_plot1 <- ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m0median))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Posterior Prediction\nN Listings (Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(angle = 45)) +
  xlab("") +
  ylab("")

#95% credible interval width
m0_plot2 <- ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m0post95wid))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "95% Credible Int. Width\n(Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") +  
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(angle = 45)) +
  xlab("") +
  ylab("")

ggsave("./output/maps/model0.pdf", arrangeGrob(m0_plot1, m0_plot2,
                                               nrow = 1),
       width = 12, height  = 6, dpi = 300)


#median fitted value
ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m1median))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Posterior Prediction\nN Listings (Non-Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/fittedM1.pdf",
         dpi = 300)

#median fitted value (for faceted plot)
m1_plot1 <- ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m1median))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Posterior Prediction\nN Listings (Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(angle = 45)) +
  xlab("") +
  ylab("")

#95% credible interval width
m1_plot2 <- ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m1post95wid))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "95% Credible Int. Width\n(Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") +  
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(angle = 45)) +
  xlab("") +
  ylab("")

ggsave("./output/maps/model1.pdf", arrangeGrob(m1_plot1, m1_plot2,
                                               nrow = 1),
       width = 12, height  = 6, dpi = 300)

#residual map
ggplot(kc_f, aes(x = long, y = lat, group = group, fill = nListings - m1median)) +
  geom_polygon(color = "grey80", lwd = .15) +
  scale_fill_gradient2(name = "Non-Spatial Model Residual") +
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/fittedM1_res.pdf",
         dpi = 300)

#median random effect
ggplot(kc_f, aes(x = long, y = lat, group = group, fill = m1RE)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Non-Spatial Effect", palette = "Purples", 
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/fittedM1RE.pdf",
         dpi = 300)



#### Model 2

#median fitted value
ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m2median))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Posterior Prediction\nN Listings (Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/fittedM2.pdf",
         dpi = 300)

#median fitted value (for faceted plot)
m2_plot1 <- ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m2median))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Posterior Prediction\nN Listings (Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(angle = 45)) +
  xlab("") +
  ylab("")

#95% credible interval width
m2_plot2 <- ggplot(kc_f, aes(x = long, y = lat, group = group, fill = log10(m2post95wid))) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "95% Credible Int. Width\n(Spatial)", palette = "Purples", 
                       limits = c(0, 4),
                       breaks = c(1, 2, 3, 4),
                       labels = c("10", "100", "1000", "10000"),
                       direction = 1, na.value = "grey80") +  
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(angle = 45)) +
  xlab("") +
  ylab("")

ggsave("./output/maps/model2.pdf", arrangeGrob(m2_plot1, m2_plot2,
                                               nrow = 1),
       width = 12, height  = 6, dpi = 300)

#median spatial effect
ggplot(kc_f, aes(x = long, y = lat, group = group, fill = m2SE)) +
  geom_polygon(color = "white", lwd = .15) +
  scale_fill_distiller(name = "Median Spatial Effect", palette = "Purples", 
                       direction = 1, na.value = "grey80") + 
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/fittedM2SE.pdf",
         dpi = 300)

#model 2 residual
ggplot(kc_f, aes(x = long, y = lat, group = group, fill = nListings - m2median)) +
  geom_polygon(color = "grey80", lwd = .15) +
  scale_fill_gradient2(name = "Spatial Model Residual") +
  coord_map() +
  theme_minimal() +
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(), # get rid of x ticks/text
        axis.ticks.x = element_blank(),axis.text.x = element_blank(), # get rid of y ticks/text
        panel.background = element_blank(),
        panel.grid = element_blank()) +
  xlab("") +
  ylab("") +
  ggsave(filename = "./output/maps/fittedM2_res.pdf",
         dpi = 300)









