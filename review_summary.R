library(dagitty)
library(ggdag)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(tidygraph)
library(dplyr)
library(tidyr)
library(ggplot2)
# library(ComplexHeatmap)
library(data.table)
library(terra)
library(tidyverse)
library(tidygeocoder)
library(maps)
library(mapdata)
library(ggnewscale)
library(rnaturalearth)
library(ggrepel)
library(ggpattern)
library(stringr)
################################################################################
##### counting reasons for exclusion in each category. 
################################################################################

## 2026 screen
excluded_all <- fread("./excluded_reasons/articles.csv")

## pull out what is after
extractedreasons <- sub(".*RAYYAN-EXCLUSION-REASONS:\\s*", "", excluded_all$notes)

table(grepl("isolation in air|cigarette|host|Host|No meteorological|no meteorological|no environmental|no air|proxy|quantified|monthly", extractedreasons))
extractedreasons[grepl("isolation in air|cigarette|host|Host|No meteorological|no meteorological|no environmental|no air|proxy|quantified|monthly", extractedreasons)]

## article type
table(grepl("wrong publication|peer|abstract|comment|thesis|english|poster", extractedreasons))
extractedreasons[grepl("wrong publication|peer|abstract|comment|thesis|english|poster", extractedreasons)]

## outcome type
table(grepl("pneumo dynamics|pneumoniae|pneumococcal|pneumococci|wrong outcome", extractedreasons))
extractedreasons[grepl("pneumo dynamics|pneumoniae|pneumococcal|pneumococci|wrong outcome", extractedreasons)]

## other (duplicates)
table(!grepl("wrong outcome|isolation in air|cigarette|pneumo dynamics|pneumoniae|pneumococcal|pneumococci|wrong publication|peer|abstract|comment|thesis|english|poster|host|Host|No meteorological|no meteorological|no environmental|no air|proxy|quantified|monthly", extractedreasons))
extractedreasons[!grepl("wrong outcome|isolation in air|cigarette|pneumo dynamics|pneumoniae|pneumococcal|pneumococci|wrong publication|peer|abstract|comment|thesis|english|poster|host|Host|No meteorological|no meteorological|no environmental|no air|proxy|quantified|monthly", extractedreasons)]

#### clean up authorship ###
    excluded_all$authors[which(excluded_all$title=="Burden of lower respiratory infections in the Eastern Mediterranean Region between 1990 and 2015: findings from the Global Burden of Disease 2015 study")] <- "Mokdad"
    # Clean curly brackets
    clean_authors <- gsub("[{}]", "", excluded_all$authors)
    # Extract first author (everything before first " and ")
    first_author_full <- sub(" and .*", "", clean_authors)
    # If format is "Surname, First" → take before comma
    first_author <- ifelse(
      grepl(",", first_author_full),
      sub(",.*", "", first_author_full),
      # Otherwise remove initials (e.g., "Samuelson D.R." → "Samuelson")
      sub("\\s+[A-Z.]+$", "", first_author_full)
    )
    # Trim whitespace
    first_author <- trimws(first_author)
    # Extract year (if you already have a year column use that instead)
    year <- excluded_all$year
    # Avoid double "et al"
    author_year <- ifelse(
      grepl("et al", excluded_all$authors),
      paste(first_author, year),
      paste(first_author, "et al", year)
    )

    author_year[14:17]

##### count number in each
env_flag <- grepl("isolation in air|cigarette|host|Host|No meteorological|no meteorological|no environmental|no air|proxy|quantified|monthly", extractedreasons, ignore.case = TRUE)
article_flag <- grepl("wrong publication|peer|abstract|comment|thesis|english|poster", extractedreasons, ignore.case = TRUE)
outcome_flag <- grepl("pneumo dynamics|pneumoniae|pneumococcal|pneumococci|wrong outcome", extractedreasons, ignore.case = TRUE)
other_flag <- !grepl("wrong outcome|isolation in air|cigarette|pneumo dynamics|pneumoniae|pneumococcal|pneumococci|wrong publication|peer|abstract|comment|thesis|english|poster|host|Host|No meteorological|no meteorological|no environmental|no air|proxy|quantified|monthly", extractedreasons, ignore.case = TRUE)
reason_flags <- data.frame(
  citation = author_year,
  # reason = extractedreasons,
  env = env_flag,
  article = article_flag,
  outcome = outcome_flag,
  other = other_flag
)
reason_flags$match_count <- rowSums(reason_flags[, c("env", "article", "outcome", "other")])
multi_matches <- reason_flags[reason_flags$match_count > 1, ]
colnames(reason_flags) <- c("citation", "wrong/unquantified exposure","wrong article type", "wrong outcome", "duplicate")
write.table(reason_flags, file = "./excluded_reasons/full_text_exclusions_2026.csv", sep = ",", quote = FALSE , col.names = TRUE, row.names = FALSE)




################################################################################
##### summarize experimental and epidemiological results with effect direction
################################################################################
expr <- fread("./extracted_results/experimental_results_2026.csv")
# expr$PM2.5 <- gsub(".","",expr$PM2.5)
expr$type <- "experimental"
exposures <- colnames(expr)[10:22]
tmp <- expr[, exposures, with = FALSE]
tmp[tmp == ""] <- NA
expr[, exposures] <- tmp

expr$Studyb <- sub(" -.*", "", expr$Studya)

summexp <- expr %>%
  group_by(simple_effect,simple_outcome) %>%
  summarise(unique_studies = n_distinct(Studyb)) %>%
  arrange(simple_effect)
summexp <- expr %>%
  group_by(simple_effect) %>%
  summarise(unique_studies = n_distinct(Studyb)) %>%
  arrange(simple_effect)
sum(summexp[grep("cig|air",summexp$simple_effect),"unique_studies"])
######## read in epi data
epir <- fread("./extracted_results/epidemiological_results_2026.csv")# 
epir$type <- "epidemiological"
exposures <- colnames(epir)[15:36]
tmp <- epir[, exposures, with = FALSE]
tmp[tmp == ""] <- NA
tmp[is.na(tmp)] <- NA
epir[, exposures] <- tmp

# studs <- sort(unique(c(epir$Study, expr$Study)))

summepi <- epir %>%
  group_by(simple_effect, simple_outcome, outcome) %>%
  summarise(unique_studies = n_distinct(Study)) %>%
  arrange(simple_effect, simple_outcome)
summepi <- epir %>%
  group_by(simple_effect, simple_outcome) %>%
  summarise(unique_studies = n_distinct(Study)) %>%
  arrange(simple_effect, simple_outcome)

#### EXTRACT LATITUDE AND LONGITUDE location ####
epir$Location[which(epir$Location == "Spain (15 provinces)")] <- "Spain"
epir$Location[which(epir$Location == "ABC (7states), USA")] <- "USA"
epir$Location[which(epir$Location == "North Jutland Couny, Denmark")] <- "Denmark"
epir$Location[which(epir$Location == "Canada, Australia, USA")] <- "Global"
epir$Location[which(epir$Location == "Barcelona, Spain")] 


geo_data <- epir %>%
  distinct(Location) %>%
  mutate(search_location = ifelse(Location == "Barcelona, Spain", "Barcelona City, Spain", Location)) %>%
  geocode(search_location, method = "osm", full_results = TRUE)  # OpenStreetMap free geocoding

# Join back to the main dataset
epir_geo <- epir %>%
  left_join(geo_data, by = "Location")

### save long lat data set
saveRDS(epir_geo, file = "./extracted_results/epidemiological_results_latlon_2026.rds")

################################################################################
####### ADD KOPPEN CODES AND PLOT LOCATIONS
################################################################################
epir_geo <- readRDS( file = "./extracted_results/epidemiological_results_latlon_2026.rds")
epir_geo$simple_effect <- factor(epir_geo$simple_effect , levels = c("air_quality","meteorology + air quality", "meteorology"),
                                 labels = c("air pollutants","meteorology + air pollutants", "meteorology"))
# Get world map data
world_map <- map_data("world")
world_map <- world_map[which(world_map$region!="Antarctica"),]
# Plot
# mapepi <- ggplot() +
#   geom_polygon(data = world_map,aes(x = long, y = lat, group = group),fill = "gray90", color = "gray70", size = 0.2) +
#   geom_point(data = epir_geo, aes(x = long, y = lat, color = as.factor(Year),
#   shape = simple_effect, group = interaction(Year, simple_effect, simple_outcome)),size = 5, alpha = 0.9, fill = "white") +
#   scale_fill_gradient(low = "grey85",high = "grey20",name = "Year") +
#   theme_bw() +
#   xlim(-180,195)+
#   labs(title = "",x = "Longitude",y = "Latitude",
#     color = "Year",shape = "Exposure") +
#   theme(panel.grid = element_blank(),axis.text = element_text(size = 10)) 
mapepi <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey70", linewidth = 0.2) +
  geom_point(data = epir_geo, aes(x = long, y = lat, color = as.numeric(Year), shape = simple_effect), size = 3, stroke = 1, alpha = 0.9) +
  # scale_color_viridis_c(option = "plasma", direction = -1) +
  scale_color_gradient(low = "#D73027", high = "#0000CD",name = "Year") +
  scale_shape_manual(values = c("air pollutants" = 21, "meteorology + air pollutants" = 22, "meteorology" = 23)) +
  coord_quickmap(xlim = c(-160, 170), ylim = c(-55, 70)) +
  labs(x = "Longitude", y = "Latitude", color = "Year", shape = "Exposure") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "right")

  print(mapepi)
  ggsave(file = "./manuscript/figures/epidemiological_location_map_2026.png", width = 8, height = 6)
  pdf(file = "./manuscript/figures/epidemiological_location_map_2026.pdf", width = 8, height = 6)
  print(mapepi)
  dev.off()
  
  

### zoom in on China and Western Europe
  cities <- ne_download(scale = 50, type = "populated_places", category = "cultural", returnclass = "sf")
  cities_we <- subset(cities, LONGITUDE > -12 & LONGITUDE < 20 & LATITUDE > 35 & LATITUDE < 60 & SCALERANK <= 2)
  china_cities <- data.frame(
    city = c("Beijing", "Shanghai", "Guangzhou", "Shenzhen", "Wuhan", "Chengdu", "Hangzhou", "Tianshui"),
    long = c(116.40, 121.47, 113.27, 114.06, 114.31, 104.06, 120.16, 105.72),
    lat  = c(39.90, 31.23, 23.13, 22.54, 30.59, 30.67, 30.27, 34.58)
  )
  china_cities <- data.frame(
    city = c("Beijing", "Shanghai", "Guangzhou", "Shenzhen", "Wuhan", "Chengdu", "Hangzhou", "Tianshui", "Maela, Thailand"),
    long = c(116.40, 121.47, 113.27, 114.06, 114.31, 104.06, 120.16, 105.72, 98.57),
    lat  = c(39.90, 31.23, 23.13, 22.54, 30.59, 30.67, 30.27, 34.58, 16.73)
  )
  
  cities_cn <- subset(cities, LONGITUDE > 73 & LONGITUDE < 135 & LATITUDE > 18 & LATITUDE < 54 & SCALERANK <= 2)

  
  mapepi_westeurope <- ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey70", linewidth = 0.2) +
    geom_point(data = epir_geo, aes(x = long, y = lat, color = as.numeric(Year), shape = simple_effect), size = 3, stroke = 1, alpha = 0.9) +
    geom_text_repel(data = cities_we, aes(x = LONGITUDE, y = LATITUDE, label = NAME), size = 3, color = "grey20", max.overlaps = 15) +
    # geom_text_repel(data = countries_we, aes(x = st_coordinates(geometry)[,1], y = st_coordinates(geometry)[,2], label = name), size = 3, color = "grey30", max.overlaps = 20) +
    scale_color_gradient(low = "#D73027", high = "#0000CD", name = "Year") +
    scale_shape_manual(values = c("air pollutants" = 21, "meteorology + air pollutants" = 22, "meteorology" = 23)) +
    coord_quickmap(xlim = c(-12, 20), ylim = c(35, 60)) +
    labs(x = "Longitude", y = "Latitude", color = "Year", shape = "Exposure") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none")
  

  mapepi_china <- ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey70", linewidth = 0.2) +
    geom_point(data = epir_geo, aes(x = long, y = lat, color = as.numeric(Year), shape = simple_effect), size = 3, stroke = 1, alpha = 0.9) +
    # geom_text_repel(data = cities_cn, aes(x = LONGITUDE, y = LATITUDE, label = NAME), size = 3, color = "grey20", max.overlaps = 15) +
    # geom_point(data = china_cities, aes(x = long, y = lat), shape = 21, size = 2.8, fill = "white", color = "black", stroke = 0.8) +
    geom_text_repel(data = china_cities, aes(x = long, y = lat, label = city), size = 3, color = "grey20") +
    scale_color_gradient(low = "#D73027", high = "#0000CD", name = "Year") +
    scale_shape_manual(values = c("air pollutants" = 21, "meteorology + air pollutants" = 22, "meteorology" = 23)) +
    coord_quickmap(xlim = c(90, 132), ylim = c(15, 43)) +
    labs(x = "Longitude", y = "Latitude", color = "Year", shape = "Exposure") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none")
  
  print(mapepi_china)
  
  pdf(file = "./manuscript/figures/epidemiological_WesternEurope_map_2026.pdf", width = 3.8, height = 3.8)
  print(mapepi_westeurope)
  dev.off()
  pdf(file = "./manuscript/figures/epidemiological_China_map_2026.pdf", width = 3.8, height = 3.8)
  print(mapepi_china)
  dev.off()
  
  
  
  
########## read in climate zones and harmonize ########
epir_geo <- readRDS( file = "./extracted_results/epidemiological_results_latlon_2026.rds")
colnames(epir_geo)[which(colnames(epir_geo)=="long")] <- "lon"
pts <- vect(epir_geo, geom = c("lon", "lat"), crs = "EPSG:4326")
kg <- rast("./koppen/1991_2020/koppen_geiger_0p00833333.tif")
kglookup <- fread("./koppen/koppen_lookup.csv")
kglookup$Description <- gsub(", ","_",kglookup$Description)
kglookup$Description <- gsub(",","",kglookup$Description)
## extract
epir_geo$koppen_code <-  terra::extract(kg, pts, method="nearest")[,2]
ext_vals <- terra::extract(kg, pts, search_radius = 3000)
epir_geo$koppen_code <- ext_vals[, 2]
### find those which are missing
epir_geo$Location[which(is.na(epir_geo$koppen_code))]
na_idx <-which(is.na(epir_geo$koppen_code))
epir_geo[grepl("Spain",epir_geo$Country),]

#################### CHECK KOPPEN RASTER ############
# library(terra)
# library(ggplot2)
# library(tidyterra)
# library(dplyr)
# 
# # 1. Crop raster to Barcelona bounding box
# bcn_extent <- ext(2.0, 2.3, 41.25, 41.5)
# kg_bcn <- crop(kg, bcn_extent)
# 
# # 2. Convert point data to sf/SpatVector
# bcn_pt <- data.frame(lon = 2.17727, lat = 41.38232, label = "Barcelona Study Point")
# 
# # 3. Build color palette from your kglookup table
# # (Ensures Class 7 = BSk and Class 8 = Csa match their standard colors)
# kglookup_sub <- kglookup %>%
#   mutate(hex = rgb(Red, Green, Blue, maxColorValue = 255))
# 
# color_map <- setNames(kglookup_sub$hex, kglookup_sub$koppen_code)
# 
# # 4. Plot using ggplot2
# ggplot() +
#   geom_spatraster(data = as.factor(kg_bcn)) +
#   scale_fill_manual(
#     values = color_map,
#     na.value = "lightblue", # Color ocean cells blue
#     name = "Köppen Code"
#   ) +
#   geom_point(data = bcn_pt, aes(x = lon, y = lat), color = "black", size = 4) +
#   geom_point(data = bcn_pt, aes(x = lon, y = lat), color = "red", size = 2.5) +
#   coord_sf(xlim = c(2.0, 2.3), ylim = c(41.25, 41.5), expand = FALSE) +
#   labs(
#     title = "Köppen-Geiger Grid Cells around Barcelona",
#     subtitle = "Red dot shows the geocoded city point relative to Csa vs BSk pixels",
#     x = "Longitude",
#     y = "Latitude"
#   ) +
#   theme_minimal()
################################################################################

### Brisbane Australia fix ## Temperate_nodryseason_hotsummer
epir_geo$Location[na_idx[2]]
epir_geo$koppen_code[na_idx[2]] <- 14

### North Jutland Denmark ## Temperate_nodryseason_warmsummer
epir_geo$Location[na_idx[1]]
epir_geo$koppen_code[na_idx[1]] <- c(15)


epir_geo <- epir_geo %>%
  left_join(kglookup, by = c("koppen_code" = "Class"))
epir_geo$koppen_code[which(epir_geo$Country=="Global")] <- NA
epir_geo$Description[which(epir_geo$Country=="Global")] <- "Global"

sub <- epir_geo[,c("Study","Location","koppen_code","Description")]
saveRDS(epir_geo, file = "./extracted_results/epidemiological_results_latlon_koppen_2026.rds")


################################################################################
############SUMMARIZE DATA
################################################################################
epir_geo <- readRDS( file = "./extracted_results/epidemiological_results_latlon_koppen_2026.rds")
epir_geo$Study[grepl("Australia",epir_geo$Country)]

epir_geoair <- epir_geo[grepl("air",epir_geo$simple_effect),]
################################################################################
############ PLOT HEATMAPS WITH EPIDEMIOLOGICAL RESULTS  #################
################################################################################
epir_geo <- readRDS( file = "./extracted_results/epidemiological_results_latlon_koppen_2026.rds")
exposures <- colnames(epir_geo)[15:37]
cols_interest <- c("Study","StudyA","simple_outcome","outcome",exposures)
epirhm <- epir_geo[,..cols_interest]
epirhm$BC[is.na(epirhm$BC)] <- ""
epirhm$Diesel[is.na(epirhm$Diesel)] <- ""
epirhm <- left_join(epirhm, unique(epir_geo[,c("Study","Description","Major")]), by = "Study")
## those with multiple results fix manually
epirhm$Study[is.na(epirhm$Major)]
# had to fix before but now are fine
# epir_geo[which(epir_geo$Study=="Kim et al. 1996"),] ### Major - 3, Temperate_nodryseason_hotsummer, Symbol Cfa
# epir_geo[which(epir_geo$Study=="Sahuquillo-Arce et al. 2017"),] ### Major - 2, Arid_steppe_hot, Symbol BSh
# epir_geo[which(epir_geo$Study=="Numminen et al. 2015"),] ### Major - 1, Tropical_savannah, Symbol Aw
# epirhm$Major[grepl("Kim",epirhm$Study)] <- 3
# epirhm$Major[grepl("Sahuquillo-Arce",epirhm$Study)] <- 2
# epirhm$Major[grepl("Numminen",epirhm$Study)] <- 1

epirhm$Major <- factor(epirhm$Major, levels = c(1,2,3,4), labels = c("Tropical","Arid","Temperate","Continental"))

df_long <- epirhm %>%
  pivot_longer(cols = c(PM2.5, PM10, CO, NO, NO2, NOx, O3, SO2, VOC, HM, BC, Diesel, Dust,
                        Vis, Tmin, Tmax, Tmean, Rf, AHmean, RHmean, VPmean, SunD, WS), names_to = "Exposure",values_to = "Effect")

df_long <- df_long %>%
  mutate(Facet_Group = case_when(
    Exposure %in% c("PM2.5", "PM10", "CO", "NO", "NO2","NOx", "O3", "SO2", "VOC","HM", "BC", "Diesel", "Dust") ~ "Air Quality",
    # Exposure == "CS" ~ "PCS Exposure",
    TRUE ~ "Meteorological" # Catch all other (Temp, RH, WS, etc.)
  ))

# df_long$StudyOutcome <- factor(df_long$StudyOutcome, levels = df_long$StudyOutcome, labels = df_long$Study)
#### create labels and colors
df_long$Effect[which(df_long$Effect=="")]<-"_"
df_long$Effect[which(df_long$Effect=="inc")]<-"_"

df_long$Effect <- factor(df_long$Effect, levels = c("pos","neg","ns"),
                           labels = c("positive","negative","not significant"))

effect_colors <- c("positive" = "#4575B4","negative" = "#D73027","not significant"= "grey20","_" = "white")

outcome_colors <- c("carriage"= "#4d1b7b","disease" = "#944654","resistance"= "darkred", "clearance" = "#c79dd7")
koppen_colors <- c("Tropical"="#1a472a","Arid"="#E8D5AA","Temperate"="#d0881b","Continental"="#6A5ACD")

# drop the column that doesn't include multiple outcomes
# df_long <- df_long[,-which(colnames(df_long)=="StudyA")]


df_plot <- df_long %>%
  mutate(Effect_clean = ifelse(Effect %in% c("", "_", NA), NA, Effect)) %>%
  group_by(Facet_Group, Exposure) %>%
  filter(any(!is.na(Effect_clean))) %>%
  ungroup()

#### add annotation columns to go alongside
ann_outcome <- epirhm %>%
  distinct(Study, simple_outcome) %>%
  mutate(Exposure = "Outcome",Facet_Group = "Annotations",Effect = NA, Effect_clean = NA)
ann_outcome$Facet_Group <- factor("Annotations",levels = c("Annotations", "Air Quality", "Meteorological"),ordered = TRUE)
# colnames(ann_outcome)[which(colnames(ann_outcome)=="StudyA")] <- "Study"

ann_climate <- epirhm %>%
  distinct(Study, Major) %>%
  mutate(Exposure = "Climate",Facet_Group = "Annotations",Effect = NA,Effect_clean = NA)
ann_climate$Facet_Group <- factor("Annotations",levels = c("Annotations", "Air Quality", "Meteorological"),ordered = TRUE)
# colnames(ann_climate)[which(colnames(ann_climate)=="StudyA")] <- "Study"

df_plot2 <- bind_rows(ann_outcome,ann_climate,df_plot)
df_plot2$Exposure <- factor(
  df_plot2$Exposure,
  levels = c("Climate","Outcome",
             "VOC","HM","SO2","PM2.5","PM10","O3","NO2","NO","NOx","Dust","CO","Diesel","BC",
             "WS","VPmean","Vis","Tmin","Tmean","Tmax", "Temperature","SunD","RHmean","AHmean","Rf"),
  labels = c("Climate","Outcome",
             "VOC","HM","SO2","PM2.5","PM10","O3","NO2","NO","NOx","Dust","CO","Diesel","Black Carbon",
             "Wind Speed","Vapor Pressure","Visibility","Temp. min","Temp. mean","Temp. max","Temperature","Sunlight","Relative Humidity","Absolute Humidity","Rainfall"))
df_plot2$Facet_Group <- factor(
  df_plot2$Facet_Group,
  levels = c("Annotations", "Air Quality", "Meteorological")
)

## order study by outcome
outcome_order <- c("carriage", "clearance", "disease", "resistance")
# outcome_order <- c("Arid","Continental","Temperate","Tropical")
df_plot2 <- df_plot2 %>%
  mutate(simple_outcome = factor(simple_outcome, levels = outcome_order)) %>%
  arrange(simple_outcome, Study) %>%
  mutate(Study = factor(Study, levels = unique(Study)))

epi_heatmap <- ggplot(df_plot2, aes(x = Study, y = Exposure, fill = Effect, group = Facet_Group)) +
  geom_tile(color = "grey85") +
  scale_fill_manual(name = "Effect", values = effect_colors, na.value = "white")+
  new_scale_fill()+
  # Outcome annotation bar
  geom_tile(data = ann_outcome,aes(x = Study, y = Exposure, fill = simple_outcome, )) +
  scale_fill_manual(name = "Outcome", values = outcome_colors) +
  # Climate region annotation bar
  new_scale_fill() +
  geom_tile(data = ann_climate, aes(x = Study, y = Exposure, fill = Major)) +
  scale_fill_manual(name = "Climate", values = koppen_colors) +
  # Put study labels at the top
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = TRUE) +
  # Expand the space so the annotation bars appear above exposures
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0),axis.title.x = element_blank(),
        axis.title.y = element_blank(),legend.position = "right", strip.text.y = element_text(angle = 0))+
  facet_grid(fct_relevel(Facet_Group,'Annotations','Air Quality','Meteorological')~.,scales = "free_y", space = "free_y")

print(epi_heatmap)
ggsave(file = "./manuscript/figures/epidemiological_heatmap_2026.png", width = 8, height = 7)
cairo_pdf(file = "./manuscript/figures/epidemiological_heatmap_2026.pdf", width = 8, height = 7)
print(epi_heatmap)
dev.off()


#### reorder the annotations so they are clustered by outcome and climate order
outcome_order  <- c("carriage", "clearance", "disease", "resistance")
climate_order  <- c("Tropical", "Arid", "Temperate", "Continental")
study_order <- unique(
  merge(
    ann_outcome[, .(Study, simple_outcome)],
    ann_climate[, .(Study, Major)],
    by = "Study",
    all = TRUE
  )
)
study_order[, simple_outcome := factor(simple_outcome, levels = outcome_order)]
study_order[, Major := factor(Major, levels = climate_order)]
study_order[, order_key := interaction(simple_outcome, Major, drop = TRUE, lex.order = TRUE)]
study_order <- study_order[order(simple_outcome, Major)]

ordered_studies <- study_order$Study

df_plot2[, Study := factor(Study, levels = ordered_studies)]
ann_outcome[, Study := factor(Study, levels = ordered_studies)]
ann_climate[, Study := factor(Study, levels = ordered_studies)]
df_plot2[, Exposure := fct_relevel(Exposure, "Outcome", "Climate")]
ann_outcome[, Exposure := "Outcome"]
ann_climate[, Exposure := "Climate"]

outcome_colors <- c("carriage"= "#4d1b7b","disease" = "#3d5b59","resistance"= "darkred", "clearance" = "#c79dd7")
koppen_colors <- c("Tropical"="#025607","Arid"="#E8D5AA","Temperate"="#d0881b","Continental"="#007bbe")
########## vertical epidemiological heatmap rather than horizontal
epi_heatmap_vertical <- ggplot(df_plot2, aes(y = Study, x = Exposure, fill = Effect)) +
  geom_tile(color = "grey85") +
  scale_fill_manual(name = "Effect", values = effect_colors, na.value = "white", drop = TRUE) +
  new_scale_fill() +
  geom_tile(data = ann_outcome,
            aes(y = Study, x = Exposure, fill = simple_outcome)) +
  scale_fill_manual(name = "Outcome", values = outcome_colors) +
  new_scale_fill() +
  geom_tile(data = ann_climate,
            aes(y = Study, x = Exposure, fill = Major)) +
  scale_fill_manual(name = "Climate", values = koppen_colors) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = TRUE, limits =rev) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0),axis.title.x = element_blank(),axis.title.y = element_blank(),
    legend.position = "right",strip.text.x = element_text(angle = 0), strip.text.x.top = element_blank()) +
  facet_grid(. ~ fct_relevel(Facet_Group, 'Annotations','Air Quality','Meteorological'),scales = "free_x",space = "free_x")

print(epi_heatmap_vertical)
ggsave(file = "./manuscript/figures/epidemiological_heatmap_vertical_2026.png", width = 7, height = 8)
pdf(file = "./manuscript/figures/epidemiological_heatmap_vertical_2026.pdf", width = 7, height = 8)
print(epi_heatmap_vertical)
dev.off()

################################################################################
#### PLOT number of studies in each region by climate ####
################################################################################
library(countrycode)
epir_geo <- readRDS( file = "./extracted_results/epidemiological_results_latlon_koppen_2026.rds")
epir_geo$Country_clean <- recode(
  epir_geo$Country,
  "The_Gambia"   = "Gambia",
  "USA"          = "United States",
  "USA_Canada"   = "United States",
  "Global"       = NA_character_,
  "Alaska"       = "United States",
  "Canada, Australia, USA" = NA_character_
)
epir_geo$Continent <- countrycode(
  sourcevar = epir_geo$Country_clean,
  origin    = "country.name",
  destination = "continent"
)
epir_geo$Continent[is.na(epir_geo$Continent)] <- "Global"
# koppen_colors <- c("Tropical"="green","Arid"="#E8D5AA","Temperate"="#d0881b","Continental"="#6A5ACD")
koppen_colors <- c("Tropical"="#025607","Arid"="#E8D5AA","Temperate"="#d0881b","Continental"="#007bbe")


studies_over_time <- epir_geo[, .(n_studies = .N), by = .(Major, Continent)]
studies_over_time$Major <- factor(studies_over_time$Major, levels = c(1,2,3,4), labels = c("Tropical","Arid","Temperate","Continental"))
studies_over_time <- subset(studies_over_time, studies_over_time$Continent!="Global")

### number per region
climcount <- ggplot(studies_over_time, aes(x=Continent, y = n_studies, fill = Major)) + 
  geom_bar(size = 1.2, stat = "identity") + 
  scale_fill_manual(values = koppen_colors) +
  labs(x = "Continent", y = "Number of Studies", fill = "Climate Zone") + theme_minimal() + 
  theme(text = element_text(size = 12), legend.position = "right") +
  theme(axis.text.x = element_text(angle=45,hjust=1))
print(climcount)
ggsave(file = "./manuscript/figures/clim_count_2026.png", width = 3, height = 3)
pdf(file = "./manuscript/figures/clim_count_2026.pdf", width = 3, height = 3)
print(climcount)
dev.off()

######################### climate zone/continent/temperature rainfall and humidity
library(data.table)
library(ggplot2)

# 1. Filter for disease outcomes
epir_disease <- epir_geo[simple_outcome == "disease"]

# 2. Convert Major climate zone codes to descriptive factors
epir_disease[, Major := factor(Major, 
                               levels = c(1, 2, 3, 4), 
                               labels = c("Tropical", "Arid", "Temperate", "Continental"))]

# Helper function to check for valid exposure values
is_active <- function(x) !is.na(x) & x != 0 & x != "" & x != FALSE

# 3. Create combined indicator flags for Temperature, Humidity, and Rainfall
epir_disease[, `:=`(
  Temperature = is_active(Tmin) | is_active(Tmax) | is_active(Tmean),
  Humidity    = is_active(AHmean) | is_active(RHmean),
  Rainfall    = is_active(Rf)
)]

# 4. Reshape data: Keep Exposure and Climate Zone (Major)
long_disease <- melt(
  epir_disease, 
  id.vars = c("Study", "Major"),
  measure.vars = c("Temperature", "Humidity", "Rainfall"),
  variable.name = "Exposure", 
  value.name = "active"
)

# Filter for rows where the exposure was measured
long_disease <- long_disease[active == TRUE]

# 5. Aggregate study counts by Exposure and Climate Zone
# uniqueN(Study) ensures each study is counted once per exposure category
studies_exposure <- long_disease[, .(n_studies = uniqueN(Study)), 
                                 by = .(Exposure, Major)]

# 6. Plot: Exposure on X-axis, Climate Zone as Fill
climcount_exposure <- ggplot(studies_exposure, aes(x = Exposure, y = n_studies, fill = Major)) + 
  geom_col(width = 0.6) + 
  scale_fill_manual(values = koppen_colors) +
  labs(
    x = "Meteorological Exposure", 
    y = "Number of Studies", 
    fill = "Climate Zone"
  ) + 
  theme_minimal() + 
  theme(
    text = element_text(size = 12), 
    legend.position = "right"
  )

print(climcount_exposure)


################################################################################
#### PLOT FREQUENCY OF EVIDENCE PER EXPOSURE FROM EPIDEMIOLOGICAL DATA
################################################################################
## df plot is from the epi heatmap section
#### epi frequency plot
df_plot3 <- df_plot %>%
  filter(
    !simple_outcome %in% c("resistance", "clearance"),
         !is.na(Effect))
df_plot4 <- df_plot3 %>%
  mutate(
    Exposure = case_when(
      Exposure == "Dust" & Effect %in% c("positive", "negative", "not significant") ~ "PM10",
      Exposure %in% c("Diesel", "Black Carbon", "HM") & Effect %in% c("positive", "negative", "not significant") ~ "PM2.5",
      Exposure %in% c("RHmean", "AHmean") & Effect %in% c("positive", "negative", "not significant") ~ "Humidity",
      TRUE ~ as.character(Exposure)
    )
  ) %>%
  filter(!Exposure %in% c("Vis", "Dust")) %>%
  # Remove phantom rows that were left empty when recoding exposures
  filter(!is.na(Effect)) %>%
  mutate(Exposure = factor(Exposure))


# Create Temperature rows
temperature_rows <- df_plot4 %>%
  filter(Exposure %in% c("Tmin", "Tmean", "Tmax")) %>%
  group_by(Study, simple_outcome) %>%
  summarise(
    Exposure = "Temperature",
    Effect = if (all(is.na(Effect))) {
      NA_character_
    } else {
      first(Effect[!is.na(Effect)])
    },
    .groups = "drop"
  )

# Remove original temperature metrics
df_plot4 <- df_plot4 %>%
  filter(!Exposure %in% c("Tmin", "Tmean", "Tmax")) %>%
  bind_rows(temperature_rows)

#### counts
df_counts <- df_plot4 %>%
  group_by(Exposure, Effect, simple_outcome, Facet_Group) %>%
  summarise(n = n(), .groups = "drop")
df_counts <- df_counts %>%
  group_by(Exposure) %>%
  mutate(total_n = sum(n)) %>%
  ungroup()

#############################
exposure_order <- c("VOC","HM", "Dust","CO", "NO", "NO2","NOx","O3", "PM2.5", "PM10", "SO2",
   "Vis","VPmean", "WS", "SunD", "Humidity","Rf", "Tmin","Tmax", "Tmean", "Temperature")
df_counts$Exposure <- factor(df_counts$Exposure, levels = exposure_order,
                             labels = c("VOC","Heavy Metals", "Dust","CO", "NO", "NO2", "NOx", "O3", "PM2.5", "PM10", "SO2",
                                        "Visibility","Vapor Pressure", "Wind Speed", "Sunlight", "Humidity","Rainfall", "Temp. min","Temp. max", "Temp. mean", "Temperature"))
effect_colors <- c("positive" = "#4575B4", "negative" = "#D73027", "not significant" = "grey20")
freqev <- ggplot(df_counts, 
       aes(x= Exposure, y = n, fill = Effect,group=interaction(Facet_Group,simple_outcome)))+
  geom_bar(stat="identity", color="black", linewidth=0.2) +
  labs(x = "Exposure",y = "Evidence (n)",fill = "Effect") +
  scale_fill_manual(values = effect_colors) + 
  theme_bw(base_size=12) +
  theme(axis.text.x = element_text(angle=45, hjust=1),panel.grid.minor = element_blank()) +
  facet_grid(. ~ simple_outcome)

print(freqev)
ggsave(file = "./manuscript/figures/epidemiological_evidence_frequency_2026.png", width = 8, height = 4)
pdf(file = "./manuscript/figures/epidemiological_evidence_frequency_2026.pdf", width = 8, height = 4)
print(freqev)
dev.off()


################################################################################
############# FIND STUDIES WITH POSITIVE VERSUS NEGATIVE RAINFALL EFFECTS
################################################################################
library(tidyverse)
library(maps)

# 1. Reshape and filter
epir_geo$Tmin[which(epir_geo$Tmin=="inc")] <- NA
epir_geo$Tmax[which(epir_geo$Tmax=="inc")] <- NA
epir_geo$RHmean[which(epir_geo$RHmean=="inc")] <- NA
epir_geo$Tmean[which(epir_geo$Tmean=="inc")] <- NA
epir_geo$AHmean[which(epir_geo$AHmean=="inc")] <- NA

epir_filtered <- epir_geo %>%
  pivot_longer(cols = c(Tmin, Tmax, Tmean, AHmean, RHmean, Rf), names_to = "Exposure_Var", values_to = "Effect_Direction") %>%
  filter(!is.na(Effect_Direction) & Effect_Direction != "") %>%
  mutate(
    Exposure_Group = case_when(Exposure_Var %in% c("Tmin", "Tmax", "Tmean") ~ "Temperature", Exposure_Var %in% c("AHmean", "RHmean") ~ "Humidity", Exposure_Var == "Rf" ~ "Rainfall"),
    Effect_Direction = case_when(tolower(Effect_Direction) %in% c("pos", "positive") ~ "positive", tolower(Effect_Direction) %in% c("neg", "negative") ~ "negative", tolower(Effect_Direction) %in% c("ns", "not significant") ~ "not significant", TRUE ~ tolower(Effect_Direction))
  ) %>%
  filter(!is.na(Exposure_Group))

epir_filtered$Major <- factor(epir_filtered$Major, levels = c(1, 2, 3, 4), labels = c("Tropical", "Arid", "Temperate", "Continental"))

# 2. Plot Map
world_map <- map_data("world")

mapepiclim <- ggplot() +
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), fill = "grey90", color = "grey70", linewidth = 0.2) +
  geom_point(data = epir_filtered, aes(x = lon, y = lat, fill = Effect_Direction, shape = Major), size = 5, stroke = 0.5, alpha = 0.6) +
  scale_fill_manual(values = c("positive" = "#4575B4", "negative" = "#D73027", "not significant" = "purple"), name = "Effect") +
  scale_shape_manual(values = c("Tropical" = 21, "Arid" = 24, "Temperate" = 22, "Continental" = 23), name = "Major Climate", na.translate = FALSE) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_quickmap(xlim = c(-160, 170), ylim = c(-55, 70)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), legend.position = "right") +
  facet_wrap(~ Exposure_Group, ncol = 1)

print(mapepiclim)
ggsave(file = "./manuscript/figures/epidemiological_evidence_HumRFTemp_2026.png", width = 8, height = 8)
pdf(file = "./manuscript/figures/epidemiological_evidence_HumRFTemp_2026.pdf", width = 8, height = 8)
print(mapepiclim)
dev.off()
################################################################################
############ PLOT HEATMAP OF EXPERIMENTAL EXPOSURES, OUTCOMES, AND EFFECTS #####
################################################################################
expr <- fread("./extracted_results/experimental_results_2026.csv")
expr_long <- expr %>%
  mutate(study_row = row_number(),
  study_label = Studya) %>% 
  pivot_longer(
    cols = c(PM2.5, PM10, BC, Wood_Smoke, Diesel, Dust, HM, Hot, RT, Cold, Rh, Rf, Sunlight),
    names_to = "exposure",
    values_to = "effect"
  ) %>%
  filter(effect != "") 
# expr_long$simple_outcome[which(expr_long$simple_outcome=="clearance_phagocyosis")] <- "clearance_phagocytosis"
expr_long$exposure <- factor(expr_long$exposure, levels = c("PM2.5","PM10","Wood_Smoke","Dust","Diesel","BC", "HM",
                                                            "Cold","RT","Hot","Rh","Rf","Sunlight"),
                             labels = c("PM2.5","PM10", "Wood Smoke","Dust","Diesel","Black Carbon", "Heavy Metals (Fe, Ni, Cadmium)",
                                        "Cold","Room Temperature","Hot","Relative Humidity (50%)","Rainfall", "Sunlight"))
expr_long$simple_outcome <- factor(expr_long$simple_outcome, levels = c("shedding","colonization","epithelial_adhesion", "inflammation",
                                                                        "bacterial_burden","biofilm","lung_translocation","lung_injury","antigen_presentation",
                                                                        "clearance_macrophagebinding","clearance_phagocytosis","clearance","pneumococcal_survival","host_mortality"))
expr_long$effect <- factor(expr_long$effect, levels = c("pos","neg","ns"),
                           labels = c("positive","negative","not significant"))

effect_colors <- c("positive" = "#4575B4","negative" = "#D73027","not significant"= "grey20","_" = "white")
# effect_colors <- c("positive" = "#0072B2","negative" = "#D55E00","not significant" = "grey20")
ann_method <- expr_long %>%
  distinct(study_label, simple_outcome, method) %>%
  mutate(x = "Experimental Type")
method_colors <- c(human = "#E8D5AA",mouse = "#054213",rat= "#acd1af",cells = "#c79dd7",sampling = "#007bbe")

### previous versions
# ann_serotype <- expr_long %>%
#   distinct(study_label, simple_outcome, Serotype) %>%
#   mutate(x = "Seroype")
# ann_serotype <- expr_long %>%
#   distinct(study_label, simple_outcome, Serotype) %>%
#   mutate(Serotype_grouped = case_when(
#       str_detect(Serotype, ",| and ") ~ "Multiple serotypes",
#       TRUE ~ Serotype),x = "Serotype")

ann_serotype <- expr_long %>%
  distinct(study_row, study_label, simple_outcome, Serotype) %>%
  mutate(Serotype_grouped = case_when(
    Serotype == "19F (ATCC 49169), 6, 7, 9, 19, 20 (clinical)" ~ "Serotypes: 19F, 6, 7, 9, 19, 20",
    Serotype == "2 (D39), 3 (ATCC 6303), 19A and 19F (ATCC 49619)" ~ "Serotypes: 2, 3, 19A, 19F",
    Serotype %in% c("1 (ATCC 33400)", "1 (unspecified)") ~ "1",
    str_detect(Serotype, regex(",|\\band\\b", ignore_case = TRUE)) ~ "Multiple",
    str_detect(Serotype, regex("^before serotypes$", ignore_case = TRUE)) ~ "Pre-serotype",
    str_detect(Serotype, regex("^unspecified$", ignore_case = TRUE)) ~ "Unspecified",
    str_detect(Serotype, regex("19F", ignore_case = TRUE)) ~ "19F",
    str_detect(Serotype, regex("6B", ignore_case = TRUE)) ~ "6B",
    str_detect(Serotype, regex("TIGR4|\\b4\\b", ignore_case = TRUE)) ~ "4",
    str_detect(Serotype, regex("D39|^2", ignore_case = TRUE)) ~ "2",
    str_detect(Serotype, regex("ATCC 6303|A66\\.1|^3", ignore_case = TRUE)) ~ "3",
    str_detect(Serotype, regex("^type 1$", ignore_case = TRUE)) ~ "Type 1",
    TRUE ~ "Other"),
    x = "Serotype")

serotype_colors <- c("1" = "#C44E52","2" = "#8172B2","3" = "#55A868","4" = "#CCB974",
  "6B" = "#64B5CD","19F" = "#4C72B0","Serotypes: 19F, 6, 7, 9, 19, 20" = "#041D56","Serotypes: 2, 3, 19A, 19F" = "#3c1361",
  "Pre-serotype" = "grey20","Type 1"= "grey","Unspecified" = "white")


exp_plot <- ggplot() +
  geom_tile(data = ann_method, aes(y = study_label, x = x, fill = method), width = 0.9) +
  scale_fill_manual(name = "Experimental Type",values = method_colors) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = ann_serotype, aes(y = study_label, x = x, fill = Serotype_grouped), color = "grey", width = 0.9) +
  scale_fill_manual(name = "Serotype",values = serotype_colors) +
  guides(fill = guide_legend(ncol = 2))+
  ggnewscale::new_scale_fill() +
  geom_tile(data = expr_long, aes(y = study_label, x = exposure, fill = effect), color = "white") +
  scale_fill_manual(values = effect_colors) +
  facet_grid(simple_outcome ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.05),
        strip.text.y = element_text(size = 8, angle = 0),
        axis.text.y = element_text(size = 7),
        panel.spacing = unit(0.8, "lines"),
        strip.background = element_rect(fill = "white"))+
  scale_x_discrete(position = "top") +
  labs(x = "",y = "",fill = "Effect") 

print(exp_plot)
ggsave(file = "./manuscript/figures/experimental_heatmap_2026.png", width = 10, height = 9)
pdf(file = "./manuscript/figures/experimental_heatmap_2026.pdf", width = 10, height = 9)
print(exp_plot)
dev.off()
# 
# geom_tile_pattern(data = ann_serotype, aes(y = study_label, x = x, fill = Serotype_grouped, pattern = Serotype_grouped),
#                   width = 0.9, color = "grey", pattern_fill = NA, pattern_density = 0.1, pattern_colour = "white", pattern_spacing = 0.025, pattern_angle = 45, pattern_alpha = 1) +
#   scale_pattern_manual(values = c("1" = "none", "2" = "none", "3" = "none", "4" = "none", "6B" = "none", "19F" = "none", "Multiple: 19F" = "stripe", 
#                                   "Multiple: 2" = "stripe", "Pre-serotype" = "none", "Type 1" = "none",  "Unspecified" = "none"), guide = "none") +

################################################################################
################## aggregate heatmap PM2.5 and PM10s ###########
################################################################################
#### SUBSET DFPLOT2 TO CLUSTER BLACK CARBON AND DIESEL WITH PM2.5 AND DUST WITH PM10
############## EPIDEMIOLOGICAL REORDER AND CLUSTER###############################


### subset and aggregate, ensure that dropping the NA row for Relative humidity when merging Absolute and Relative Humidity into just humidity
df_plot4 <- df_plot2 %>%
  mutate(Exposure = case_when(
      Exposure == "Dust" & Effect %in% c("positive", "negative", "not significant") ~ "PM10",
      Exposure %in% c("Diesel", "Black Carbon", "HM", "Wood Smoke") & Effect %in% c("positive", "negative", "not significant", NA) ~ "PM2.5",
      Exposure %in% c("Relative Humidity", "Absolute Humidity") ~ "Humidity",
      TRUE ~ as.character(Exposure))) %>%
  filter(!Exposure %in% c("Visibility", "Dust")) %>%
  group_by(Study, Exposure, Facet_Group, simple_outcome) %>%
  summarise(
    Effect = if (all(is.na(Effect))) {
      NA_character_
    } else {
      Effect[which(!is.na(Effect))[1]]
    },
    .groups = "drop"
  ) %>%
  mutate(Exposure = factor(Exposure))

# Create Temperature rows
temperature_rows <- df_plot4 %>%
  filter(Exposure %in% c("Temp. max", "Temp. mean", "Temp. min")) %>%
  group_by(Study, simple_outcome) %>%
  summarise(
    Exposure = "Temperature",
    Effect = if (all(is.na(Effect))) {
      NA_character_
    } else {
      first(Effect[!is.na(Effect)])
    },
    .groups = "drop"
  )

# Remove original temperature metrics
df_plot4 <- df_plot4 %>%
  filter(!Exposure %in% c("Temp. max", "Temp. mean", "Temp. min")) %>%
  bind_rows(temperature_rows)
df_plot4$Facet_Group[which(df_plot4$Exposure=="Temperature")] <- "Meteorological"

#### reorder the annotations so they are clustered by outcome and climate order
outcome_order  <- c("carriage", "clearance", "disease", "resistance")
climate_order  <- c("Tropical", "Arid", "Temperate", "Continental")
study_order <- unique(
  merge(
    ann_outcome[, .(Study, simple_outcome)],
    ann_climate[, .(Study, Major)],
    by = "Study",
    all = TRUE
  )
)
study_order[, simple_outcome := factor(simple_outcome, levels = outcome_order)]
study_order[, Major := factor(Major, levels = climate_order)]
study_order[, order_key := interaction(simple_outcome, Major, drop = TRUE, lex.order = TRUE)]
study_order <- study_order[order(simple_outcome, Major)]

ordered_studies <- study_order$Study
setDT(df_plot4)
df_plot4[, Study := factor(Study, levels = ordered_studies)]
ann_outcome[, Study := factor(Study, levels = ordered_studies)]
ann_climate[, Study := factor(Study, levels = ordered_studies)]
df_plot4[, Exposure := fct_relevel(Exposure, "Outcome", "Climate")]
ann_outcome[, Exposure := "Outcome"]
ann_climate[, Exposure := "Climate"]


outcome_colors <- c("carriage"= "#4d1b7b","disease" = "#3d5b59","resistance"= "darkred", "clearance" = "#c79dd7")
koppen_colors <- c("Tropical"="#025607","Arid"="#E8D5AA","Temperate"="#d0881b","Continental"="#007bbe")

df_plot4$Effect[which(is.na(df_plot4$Effect))] <- NA
df_plot4$Effect <- factor(df_plot4$Effect, levels = c("positive","negative","not significant"))

df_plot4[which(grepl("Berry",df_plot4$Study)&!is.na(df_plot4$Effect)),]
ann_outcome[which(grepl("Berry",ann_outcome$Study)),]
ann_climate[which(grepl("Berry",ann_climate$Study)),]

epi_heatmap_vertical <- ggplot(df_plot4, aes(y = Study, x = Exposure, fill = Effect)) +
  geom_tile(color = "grey85") +
  scale_fill_manual(name = "Effect", values = effect_colors, na.value = "white", drop = TRUE) +
  new_scale_fill() +
  geom_tile(data = ann_outcome,
            aes(y = Study, x = Exposure, fill = simple_outcome)) +
  scale_fill_manual(name = "Outcome", values = outcome_colors) +
  new_scale_fill() +
  geom_tile(data = ann_climate,
            aes(y = Study, x = Exposure, fill = Major)) +
  scale_fill_manual(name = "Climate", values = koppen_colors) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(drop = TRUE, limits =rev) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0),axis.title.x = element_blank(),axis.title.y = element_blank(),
        legend.position = "right",strip.text.x = element_text(angle = 0), strip.text.x.top = element_blank()) +
  facet_grid(. ~ fct_relevel(Facet_Group, 'Annotations','Air Quality','Meteorological'),scales = "free_x",space = "free_x")

print(epi_heatmap_vertical)
ggsave(file = "./manuscript/figures/epidemiological_heatmap_vertical_aggregatePM_2026.png", width = 7, height = 8)
pdf(file = "./manuscript/figures/epidemiological_heatmap_vertical_aggregatePM_2026.pdf", width = 7, height = 8)
print(epi_heatmap_vertical)
dev.off()


############## EXPERIMENTAL REORDER AND CLUSTER ################################
expr <- fread("./extracted_results/experimental_results_2026.csv")
expr_long <- expr %>%
  mutate(study_row = row_number(),
         study_label = Studya) %>% 
  pivot_longer(
    cols = c(PM2.5, PM10, BC, Wood_Smoke, Diesel, Dust, HM, Hot, RT, Cold, Rh, Rf, Sunlight),
    names_to = "exposure",
    values_to = "effect"
  ) %>%
  filter(effect != "")  

expr_longagg <- expr_long %>%
  mutate(
    exposure_agg = case_when(
      exposure == "Dust" ~ "PM10",
      exposure %in% c("Diesel", "BC", "HM", "Wood_Smoke") ~ "PM2.5",
      TRUE ~ as.character(exposure)
    )
  ) %>%
  ## prioritize significant results over non-significant
  group_by(study_label, outcome, exposure_agg) %>%
  filter(!(effect == "ns" & any(effect %in% c("pos", "neg")))) %>%
  ungroup() %>%
  mutate(exposure = factor(exposure_agg)) %>%
  select(-exposure_agg)

## COMBINE DIESEL, BC, HM INTO PM2.5 AND DUST INTO PM10
expr_longagg$exposure <- factor(expr_longagg$exposure, levels = c("PM2.5","PM10",
                                                            "Cold","RT","Hot","Rh","Rf","Sunlight"),
                             labels = c("PM2.5", "PM10", 
                                        "Cold","Room Temperature","Hot","Relative Humidity (50%)","Rainfall", "Sunlight"))
expr_longagg$simple_outcome <- factor(expr_longagg$simple_outcome, levels = c("shedding","colonization","epithelial_adhesion", "inflammation",
                                                                        "bacterial_burden","biofilm","lung_translocation","lung_injury","antigen_presentation",
                                                                        "clearance_macrophagebinding","clearance_phagocytosis","clearance","pneumococcal_survival","host_mortality"))
expr_longagg$effect <- factor(expr_longagg$effect, levels = c("pos","neg","ns"),
                           labels = c("positive","negative","not significant"))

effect_colors <- c("positive" = "#4575B4","negative" = "#D73027","not significant"= "grey20","_" = "white")
# effect_colors <- c("positive" = "#0072B2","negative" = "#D55E00","not significant" = "grey20")
ann_method <- expr_longagg %>%
  distinct(study_label, simple_outcome, method) %>%
  mutate(x = "Experimental Type")
method_colors <- c(human = "#E8D5AA",mouse = "#054213",rat= "#acd1af",cells = "#c79dd7",sampling = "#007bbe")


ann_serotype <- expr_longagg %>%
  distinct(study_row, study_label, simple_outcome, Serotype) %>%
  mutate(Serotype_grouped = case_when(
    Serotype == "19F (ATCC 49169), 6, 7, 9, 19, 20 (clinical)" ~ "Serotypes: 19F, 6, 7, 9, 19, 20",
    Serotype == "2 (D39), 3 (ATCC 6303), 19A and 19F (ATCC 49619)" ~ "Serotypes: 2, 3, 19A, 19F",
    Serotype %in% c("1 (ATCC 33400)", "1 (unspecified)") ~ "1",
    str_detect(Serotype, regex(",|\\band\\b", ignore_case = TRUE)) ~ "Multiple",
    str_detect(Serotype, regex("^before serotypes$", ignore_case = TRUE)) ~ "Pre-serotype",
    str_detect(Serotype, regex("^unspecified$", ignore_case = TRUE)) ~ "Unspecified",
    str_detect(Serotype, regex("19F", ignore_case = TRUE)) ~ "19F",
    str_detect(Serotype, regex("6B", ignore_case = TRUE)) ~ "6B",
    str_detect(Serotype, regex("TIGR4|\\b4\\b", ignore_case = TRUE)) ~ "4",
    str_detect(Serotype, regex("D39|^2", ignore_case = TRUE)) ~ "2",
    str_detect(Serotype, regex("ATCC 6303|A66\\.1|^3", ignore_case = TRUE)) ~ "3",
    str_detect(Serotype, regex("^type 1$", ignore_case = TRUE)) ~ "Type 1",
    TRUE ~ "Other"),
    x = "Serotype")

serotype_colors <- c("1" = "#C44E52","2" = "#8172B2","3" = "#55A868","4" = "#CCB974",
                     "6B" = "#64B5CD","19F" = "#4C72B0","Serotypes: 19F, 6, 7, 9, 19, 20" = "#041D56","Serotypes: 2, 3, 19A, 19F" = "#3c1361",
                     "Pre-serotype" = "grey20","Type 1"= "grey","Unspecified" = "white")


exp_plot <- ggplot() +
  geom_tile(data = ann_method, aes(y = study_label, x = x, fill = method), width = 0.9) +
  scale_fill_manual(name = "Experimental Type",values = method_colors) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = ann_serotype, aes(y = study_label, x = x, fill = Serotype_grouped), color = "grey", width = 0.9) +
  scale_fill_manual(name = "Serotype",values = serotype_colors) +
  guides(fill = guide_legend(ncol = 2))+
  ggnewscale::new_scale_fill() +
  geom_tile(data = expr_longagg, aes(y = study_label, x = exposure, fill = effect), color = "white") +
  scale_fill_manual(values = effect_colors) +
  facet_grid(simple_outcome ~ ., scales = "free_y", space = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = -0.05),
        strip.text.y = element_text(size = 8, angle = 0),
        axis.text.y = element_text(size = 7),
        panel.spacing = unit(0.8, "lines"),
        strip.background = element_rect(fill = "white"))+
  scale_x_discrete(position = "top") +
  labs(x = "",y = "",fill = "Effect") 
print(exp_plot)
ggsave(file = "./manuscript/figures/experimental_heatmap_aggregatePM_2026.png", width = 10, height = 9)
pdf(file = "./manuscript/figures/experimental_heatmap_aggregatePM_2026.pdf", width = 10, height = 9)
print(exp_plot)
dev.off()


################################################################################
############# PLOT THE STRENGTH OF EVIDENCE IN A HEATMEAP ##################
################################################################################
library(ggplot2)
library(dplyr)
library(data.table)
library(tidyr)
library(textshaping)
## load in experimental data and combine the exposures of interest
expr <- fread("./extracted_results/experimental_results_2026.csv")
expr_long <- expr %>%
  mutate(study_row = row_number(),study_label = Studya) %>% 
  pivot_longer(cols = c(PM2.5, PM10, BC, Wood_Smoke, Diesel, Dust, HM, Hot, RT, Cold, Rh, Rf, Sunlight),
    names_to = "exposure",values_to = "effect") %>%
  filter(effect != "")  

expr_longagg <- expr_long %>%
  mutate(
    exposure_agg = case_when(
      exposure == "Dust" ~ "PM10",
      exposure %in% c("Diesel", "BC", "HM", "Wood_Smoke") ~ "PM2.5",
      TRUE ~ as.character(exposure)
    )
  ) %>%
  group_by(study_label, outcome, exposure_agg) %>%
  filter(!(effect == "ns" & any(effect %in% c("pos", "neg")))) %>%
  ungroup() %>%
  mutate(exposure = factor(exposure_agg)) %>%
  select(-exposure_agg)

ftable(expr_longagg$exposure,expr_longagg$simple_outcome,expr_longagg$effect)
## save the columsn I need.
exp_evidence <- expr_longagg[,c("exposure","simple_outcome","effect")]

#### load in the epidemiological data
epir_geo <- readRDS( file = "./extracted_results/epidemiological_results_latlon_koppen_2026.rds")
df_long <- epir_geo %>%
  pivot_longer(cols = c(PM2.5, PM10, CO, NO, NO2, NOx, O3, SO2, VOC, HM, BC, Diesel, Dust,
                        Vis, Tmin, Tmax, Tmean, Rf, AHmean, RHmean, VPmean, SunD, WS), names_to = "Exposure",values_to = "Effect")
df_plot4 <- df_long %>%
  mutate(Exposure = case_when(
    Exposure == "Dust" & Effect %in% c("positive", "negative", "not significant") ~ "PM10",
    Exposure %in% c("Diesel", "BC", "HM", "WS") & Effect %in% c("positive", "negative", "not significant", NA) ~ "PM2.5",
    Exposure %in% c("RHmean", "AHmean") ~ "Humidity",
    TRUE ~ as.character(Exposure))) %>%
  filter(!Exposure %in% c("Visibility", "Dust")) %>%
  group_by(Study, Exposure,outcome) %>%
  summarise(Effect = if (all(is.na(Effect))) {
      NA_character_} else {
      Effect[which(!is.na(Effect))[1]]},.groups = "drop")


# Create Temperature rows combining maximum mean and minimum
temperature_rows <- df_plot4 %>%
  filter(Exposure %in% c("Tmax", "Tmean", "Tmin")) %>%
  group_by(Study, outcome) %>%
  summarise(
    Exposure = "Temperature",
    Effect = if (all(is.na(Effect))) {
      NA_character_
    } else {first(Effect[!is.na(Effect)])},
    .groups = "drop" )
# Remove original temperature metrics
df_plot4 <- df_plot4 %>%
  filter(!Exposure %in% c("Tmax", "Tmean", "Tmin")) %>%
  bind_rows(temperature_rows)

### clean up outcomes
df_plot4$outcome_clean <- ifelse(grepl("carriage|acquisition",df_plot4$outcome), "Carriage", 
                                 ifelse(grepl("resist",df_plot4$outcome), "Resistance", 
                                        ifelse(grepl("blood|mening|invasive|bacteremia",df_plot4$outcome), "Invasive_Disease",
                                                     ifelse(grepl("infection|otitis|pneumonia|acute",df_plot4$outcome), "Noninvasive_Disease", "Clearance"))))
## save the columns I need
epi_evidence <- df_plot4[,c("Exposure","outcome_clean","Effect")]
colnames(epi_evidence) <- c("Exposure","Outcome","Effect")
epi_evidence$type <- "epidemiological"
colnames(exp_evidence) <-  c("Exposure","Outcome","Effect")
exp_evidence$type <- "experimental"

expr_longagg[which(expr_longagg$Study=="Cheliotis et al 2022"),]
#### save all together
all_evidence <- rbind(epi_evidence,exp_evidence)
all_evidence$Effect[which(all_evidence$Effect=="inc")] <- NA
all_evidence$Effect[which(all_evidence$Effect=="")] <- NA
all_evidence <- all_evidence[!is.na(all_evidence$Effect),]
#### prepare heatmap - aggregate data and calculate net agreement score
df_summary <- all_evidence %>%
  count(type, Exposure, Outcome, Effect) %>%
  pivot_wider(names_from = Effect, values_from = n, values_fill = 0) %>%
  # Ensure all necessary effect columns exist
  { if (!"neg" %in% names(.)) mutate(., neg = 0) else . } %>%
  { if (!"ns"  %in% names(.)) mutate(., ns = 0)  else . } %>%
  { if (!"pos" %in% names(.)) mutate(., pos = 0) else . } %>%
  mutate(
    total_studies = pos + neg + ns,
    # Score ranges from -1 (100% negative) to +1 (100% positive)
    net_score = ifelse(total_studies > 0, (pos - neg) / total_studies, 0)
  )

# Create wide matrix for clustering calculation
score_mat <- df_summary %>%
  select(Outcome, Exposure, net_score) %>%
  pivot_wider(names_from = Exposure, values_from = net_score, values_fill = 0) %>%
  column_to_rownames("Outcome") %>%
  as.matrix()

# Calculate distance and dendrogram order for rows and columns
row_order <- hclust(dist(score_mat))$order
col_order <- hclust(dist(t(score_mat)))$order

# Re-factor levels based on dendrogram clustering
df_summary <- df_summary %>%
  mutate(Outcome = factor(Outcome, levels = rownames(score_mat)[row_order]),
    Exposure       = factor(Exposure,       levels = colnames(score_mat)[col_order]))

### drop those we don't want
df_summary2 <- df_summary
df_summary2 <- subset(df_summary, df_summary$Exposure %in% c("PM2.5","PM10","SO2","NO2","CO","Temperature","Humidity","Rf","Cold","Hot","Sunlight"))
df_summary2$Exposure <- factor(df_summary2$Exposure, levels = c("PM2.5","PM10","SO2","NO2","CO","Temperature","Humidity","Rf","Cold","Hot","Sunlight"),
                               labels = c("PM2.5","PM10","SO2","NO2","CO","Temperature","Humidity","Rainfall","Cold","Hot","Sunlight"))

df_summary2$Outcome <- factor(df_summary2$Outcome, levels = rev(c("Carriage","Noninvasive_Disease","Invasive_Disease","Clearance","Resistance",
                                                              "colonization","epithelial_adhesion","inflammation","bacterial_burden","biofilm",
                                                              "lung_translocation","lung_injury","antigen_presentation","clearance_macrophagebinding",
                                                              "clearance_phagocytosis","pneumococcal_survival","host_mortality")),
                                                  labels =  rev(c("Carriage","Non Invasive Disease","Invasive Disease","Clearance","Resistance",
                                                                "Colonization","Epithelial Adhesion","Inflammation","Bacterial Burden","Biofilm",
                                                                "Lung Translocation","Lung Injury","Antigen Presentation","Macrophage Binding",
                                                                "Phagocytosis","Pneumococcal Survival","Host Mortality")))

support <- ggplot(df_summary2 %>% filter(total_studies > 0), aes(x = Exposure, y = Outcome)) +
  # Heatmap tiles mapped to net effect score and study volume
  geom_tile(aes(fill = net_score), color = "white", linewidth = 0.6) +
  # Numerical study count labels inside active tiles
  geom_text(aes(label = total_studies), color = "black", size = 3.5, fontface = "bold") +
  # Color scale: Red (Negative) -> White (Neutral) -> Blue (Positive)
  scale_fill_gradient2(
    low = "#D73027",
    mid = "grey",
    high = "#4575B4",
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    labels = c("-1.0 Negative", "0", "1.0 Positive"),
    name = "Agreement Score") +
  facet_grid( type ~ ., scales = "free_y", space = "free_y") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
        strip.background = element_rect(fill = "white"),
    axis.title = element_blank(),plot.title = element_text( size = 14) ) 


print(support)
ggsave(file = "./manuscript/figures/supportingevidence_heatmap_2026.png", width = 7, height = 5)
pdf(file = "./manuscript/figures/supportingevidence_heatmap_2026.pdf", width = 7, height = 5)
print(support)
dev.off()



############# METEOROLOGICAL FACTORS DISAGGREGATION FOR DISEASE BY CLIMATE ZONE AND COUNTRY OF STUDY
epir_geo <- readRDS( file = "./extracted_results/epidemiological_results_latlon_koppen_2026.rds")
epir_geo$Country_clean <- recode(
  epir_geo$Country,
  "The_Gambia"   = "Gambia",
  "USA"          = "United States",
  "USA_Canada"   = "United States & Canada",
  "Global"       = NA_character_,
  "Alaska"       = "United States",
  "New_Zealand"       = "New Zealand",
  "Canada, Australia, USA" = NA_character_
)
clean <- function(x) ifelse(as.character(x) %in% c("inc", "", "NA"), NA_character_, as.character(x))

plot_datacountrystrat <- epir_geo %>%
  filter(simple_outcome %in% c("disease")) %>%
  mutate(
    Major = factor(Major, levels = 1:4, labels = c("Tropical", "Arid", "Temperate", "Continental")),
    Temperature = coalesce(clean(Tmean), clean(Tmax), clean(Tmin)),
    Humidity    = coalesce(clean(AHmean), clean(RHmean)),
    Rainfall    = clean(Rf)
  ) %>%
  pivot_longer(cols = c(Temperature, Humidity, Rainfall), names_to = "Exposure", values_to = "Effect_Direction") %>%
  filter(!is.na(Effect_Direction) & !is.na(Country_clean) & !is.na(Major)) %>%
  mutate(Effect_Direction = factor(Effect_Direction, levels = c("pos", "neg", "ns"))) %>%
  # Retain Major alongside Country_clean to group by Climate Zone
  summarize(n_studies = n_distinct(Study), .by = c(Major, Country_clean, Exposure, Effect_Direction)) 
plot_datacountrystrat$Exposure <- factor(plot_datacountrystrat$Exposure, levels = c("Temperature", "Humidity", "Rainfall"))
climstratmet <- ggplot(plot_datacountrystrat, aes(x = Country_clean, y = n_studies, fill = Effect_Direction)) +
  geom_col(width = 0.6) +
  # Grid with Exposure on rows, Climate Zone on columns
  facet_grid(Exposure ~ Major, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c("pos" = "#4575B4", "neg" = "#D73027", "ns" = "#999999"),
                    labels = c("Positive", "Negative", "Non-significant")) +
  labs(x = "Country", y = "Number of Studies (disease)", fill = "Effect Direction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank())

print(climstratmet)
ggsave(file = "./manuscript/figures/supportingevidence_meteorologicalstratification_2026.png", width = 6, height = 3.5)
pdf(file = "./manuscript/figures/supportingevidence_meteorologicalstratification_2026.pdf", width = 6, height = 3.5)
print(climstratmet)
dev.off()



################################################################################
##### SUMMARY OF EXCLUSIONS FOR THE APPENDIX
################################################################################

library(dplyr)
library(purrr)
library(stringr)

cigex <- fread("./excluded_reasons/cigarette_experiments/articles.csv")
cigex$type <- "cigarette_experiment"
pas<- fread("./excluded_reasons/passive_cigarettesmoke/articles.csv")
pas$type <- "passive_smoking"
air <- fread("./excluded_reasons/isolationinair/articles.csv")
air$type <- "isolation_in_air"
# 1. Helper function to extract the first surname (for sorting)
get_first_surname <- function(author_str) {
  if (is.na(author_str) || author_str == "") return("ZZZ")
  authors <- unlist(str_split(author_str, "\\s+and\\s+"))
  str_trim(str_split(authors[1], ",")[[1]][1])
}

# 2. Helper function to format the citation string
format_citation <- function(author_str, yr) {
  if (is.na(author_str) || author_str == "") return(paste0("Unknown, ", yr))
  
  authors <- unlist(str_split(author_str, "\\s+and\\s+"))
  first_surname <- str_trim(str_split(authors[1], ",")[[1]][1])
  num_authors <- length(authors)
  
  if (num_authors == 1) {
    paste0(first_surname, ", ", yr)
  } else if (num_authors == 2) {
    second_surname <- str_trim(str_split(authors[2], ",")[[1]][1])
    paste0(first_surname, " & ", second_surname, ", ", yr)
  } else {
    paste0(first_surname, " et al., ", yr)
  }
}

# 3. Processing function: cleans braces, formats, sorts A-Z, and drops DOI
clean_and_sort <- function(df) {
  df %>%
    mutate(
      # Remove all '{' and '}' characters from the author column
      authors_clean = str_remove_all(authors, "[{}]"),
      
      # Extract surname for sorting
      sort_surname = map_chr(authors_clean, get_first_surname),
      
      # Build intact "Author et al., Year" string
      `Author et al., Year` = map2_chr(authors_clean, year, format_citation)
    ) %>%
    arrange(tolower(sort_surname)) %>% # Sort A to Z
    select(
      `Author et al., Year`,
      Journal = journal,
      `Article Title` = title
    )
}

# 4. Process all three datasets
pas_clean   <- clean_and_sort(pas)
cigex_clean <- clean_and_sort(cigex)
air_clean   <- clean_and_sort(air)
# Save as three separate CSV files
write.table(pas_clean,   "./excluded_reasons/pas_clean.csv",   row.names = FALSE, quote = FALSE, sep = "\t")
write.table(cigex_clean, "./excluded_reasons/cigex_clean.csv", row.names = FALSE, quote = FALSE, sep = "\t")
write.table(air_clean,   "./excluded_reasons/air_clean.csv",   row.names = FALSE, quote = FALSE, sep = "\t")
################################################################################
# 
# ########### STUDIES IN BOXES FROM EVIDENCE SUMMARY
# library(dplyr)
# 
# # Pull study-level details for each box
# studies_in_boxes <- expr_longagg %>%
#   select(exposure, simple_outcome, study_label, effect) %>% # replace study_id with your actual study ID column name
#   arrange(simple_outcome, exposure, study_label)
# 
# # Pull study-level details for each box
# studies_in_boxes <- df_plot4 %>%
#   select(Exposure, outcome_clean, Study, Effect) %>% # replace study_id with your actual study ID column name
#   arrange(outcome_clean, Exposure, Study)
# studies_in_boxesepi <- studies_in_boxes[which(!is.na(studies_in_boxes$Effect)),]

