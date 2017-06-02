library(tidyverse)
library(stringr)
library(rvest)
library(treemap)
library(highcharter)
library(viridis)

path <- function(x) paste0("https://raw.githubusercontent.com/phalt/pokeapi/master/data/v2/csv/", x)

dfpkmn <- read_csv(path("pokemon.csv")) %>% 
  select(-order, -is_default) %>% 
  rename(pokemon = identifier)
# Parsed with column specification:
#   cols(
#     id = col_integer(),
#     identifier = col_character(),
#     species_id = col_integer(),
#     height = col_integer(),
#     weight = col_integer(),
#     base_experience = col_integer(),
#     order = col_integer(),
#     is_default = col_integer()
#   )
dfstat <- read_csv(path("stats.csv")) %>% 
  rename(stat_id = id) %>% 
  right_join(read_csv(path("pokemon_stats.csv")),
             by = "stat_id") %>% 
  mutate(identifier = str_replace(identifier, "-", "_")) %>% 
  select(pokemon_id, identifier, base_stat) %>% 
  spread(identifier, base_stat) %>% 
  rename(id = pokemon_id)
# Parsed with column specification:
#   cols(
#     id = col_integer(),
#     damage_class_id = col_integer(),
#     identifier = col_character(),
#     is_battle_only = col_integer(),
#     game_index = col_integer()
#   )
# Parsed with column specification:
#   cols(
#     pokemon_id = col_integer(),
#     stat_id = col_integer(),
#     base_stat = col_integer(),
#     effort = col_integer()
#   )
dftype <- read_csv(path("types.csv")) %>% 
  rename(type_id = id) %>% 
  right_join(read_csv(path("pokemon_types.csv")), by = "type_id") %>% 
  select(pokemon_id, identifier, slot) %>% 
  mutate(slot = paste0("type_", slot)) %>% 
  spread(slot, identifier) %>% 
  rename(id = pokemon_id)
# Parsed with column specification:
#   cols(
#     id = col_integer(),
#     identifier = col_character(),
#     generation_id = col_integer(),
#     damage_class_id = col_integer()
#   )
# Parsed with column specification:
#   cols(
#     pokemon_id = col_integer(),
#     type_id = col_integer(),
#     slot = col_integer()
#   )
dfegg <- read_csv(path("egg_groups.csv")) %>% 
  rename(egg_group_id = id) %>% 
  right_join(read_csv(path("pokemon_egg_groups.csv")), by = "egg_group_id") %>% 
  group_by(species_id) %>% 
  mutate(ranking = row_number(),
         ranking = paste0("egg_group_", ranking)) %>% 
  select(species_id, ranking, identifier) %>% 
  spread(ranking, identifier) 

dfimg <- "https://github.com/phalt/pokeapi/tree/master/data/Pokemon_XY_Sprites" %>% 
  read_html() %>% 
  html_nodes("tr.js-navigation-item > .content > .css-truncate a") %>% 
  map_df(function(x){
    url <- x %>% html_attr("href")
    data_frame(
      id = str_extract(basename(url), "\\d+"),
      url_image = basename(url)
    )
  }) %>%
  mutate(id = as.numeric(id))

url_bulbapedia_list <- "http://bulbapedia.bulbagarden.net/wiki/List_of_Pok%C3%A9mon_by_base_stats_(Generation_VI-present)" 

id <- url_bulbapedia_list %>% 
  read_html(encoding = "UTF-8") %>% 
  html_node("table.sortable") %>% 
  html_table() %>% 
  .[[1]] %>% 
  as.numeric()

url_icon <-  url_bulbapedia_list %>% 
  read_html() %>%
  html_nodes("table.sortable img") %>% 
  html_attr("src")

dficon <- data_frame(id, url_icon) %>% 
  filter(!is.na(id)) %>% 
  distinct(id)

dfcolor <- map_df(na.omit(unique(c(dftype$type_1, dftype$type_2))), function(t){
  # t <- "bug"
  col <- "http://pokemon-uranium.wikia.com/wiki/Template:%s_color" %>% 
    sprintf(t) %>%
    read_html() %>% 
    html_nodes("span > b") %>% 
    html_text()
  data_frame(type = t, color = paste0("#", col))
})

#type,color 18*2 

dfcolorf <- expand.grid(color_1 = dfcolor$color, color_2 = dfcolor$color,
                        stringsAsFactors = FALSE) %>% 
  tbl_df() %>% 
  group_by(color_1, color_2) %>% 
  do({
    n = 100;p = 0.25
    data_frame(color_f = colorRampPalette(c(.$color_1, .$color_2))(n)[round(n*p)])
  })

#color_1 color_2 color_f #324*3

# THE join
df <- dfpkmn %>% 
  left_join(dftype, by = "id") %>% 
  left_join(dfstat, by = "id") %>% 
  left_join(dfcolor %>% rename(type_1 = type, color_1 = color), by = "type_1") %>% 
  left_join(dfcolor %>% rename(type_2 = type, color_2 = color), by = "type_2") %>% 
  left_join(dfcolorf, by =  c("color_1", "color_2")) %>% 
  #left_join(dfegg, by = "species_id") %>% 
  #left_join(dfimg, by = "id") %>% 
  #left_join(dficon, by = "id")
  
#   > length(unique(df$color_1))
# [1] 18
# > length(unique(df$color_2))
# [1] 19
# > length(unique(df$color_f))
# [1] 137

rm(dftype, dfstat, dfcolor, dfcolorf, dfegg, dfimg, dficon)
rm(id, url_bulbapedia_list, url_icon)

df <- df %>% 
  mutate(color_f = ifelse(is.na(color_f), color_1, color_f)) %>% 
  filter(!is.na(url_image)) 

head(df)

dstype <- df %>% 
  count(type_1, color_1) %>% 
  ungroup() %>% 
  arrange(desc(n)) %>% 
  mutate(x = row_number()) %>% 
  rename(name = type_1,
         color = color_1,
         y = n) %>% 
  select(y, name, color) %>% 
  list_parse()
  #list.parse3()

hcbar <- highchart() %>% 
  hc_xAxis(categories = unlist(pluck(dstype, i = 2))) %>% 
  hc_yAxis(title = NULL) %>% 
  hc_add_series(data = dstype, type = "bar", showInLegend = FALSE,
                name = "Number of species")



tm <- df %>% 
  mutate(type_2 = ifelse(is.na(type_2), paste("only", type_1), type_2),
         type_1 = type_1) %>% 
  group_by(type_1, type_2) %>%
  summarise(n = n()) %>% 
  ungroup() %>% 
  treemap::treemap(index = c("type_1", "type_2"),
                   vSize = "n", vColor = "type_1")

# > summary(tm)
# Length Class      Mode     
# tm        12     data.frame list     
# type       1     -none-     character
# vSize      1     -none-     character
# vColor     1     -none-     logical  
# stdErr     1     -none-     character
# algorithm  1     -none-     character
# vpCoorX    2     -none-     numeric  
# vpCoorY    2     -none-     numeric  
# aspRatio   1     -none-     numeric  
# range      1     -none-     logical  
# mapping    3     -none-     logical  
# draw       1     -none-     logical 

tm$tm <- tm$tm %>%
  tbl_df() %>% 
  left_join(df %>% select(type_1, type_2, color_f) %>% distinct(), by = c("type_1", "type_2")) %>%
  left_join(df %>% select(type_1, color_1) %>% distinct(), by = c("type_1")) %>% 
  mutate(type_1 = paste0("Main ", type_1),
         color = ifelse(is.na(color_f), color_1, color_f))

hctm <- highchart() %>% 
  hc_add_series_treemap(tm, allowDrillToNode = TRUE,
                        layoutAlgorithm = "squarified") #hc_add_series_treemap
