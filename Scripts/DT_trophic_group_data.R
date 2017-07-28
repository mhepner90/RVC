RVCdata_DT_rds = 'big_csv/RVCdata_DT.rds'

RVCdata_DT <- getRvcData(1999:2016, "DRY TORT", server = 'https://www.sefsc.noaa.gov/rvc_analysis20/')
write_rds(RVCdata_DT, 'big_csv/RVCdata_DT.rds')

RVCdata_DT = read_rds('big_csv/RVCdata_DT.rds')


forage_dt_abun_csv     = "trophic_groups/abundance/forage_dt_abun.csv"
forage_dt_den_csv      = "trophic_groups/density/forage_dt_den.csv"
forage_dt_bio_csv      = "trophic_groups/density/forage_dt_den.csv"
forage_dt_diversity_csv= "trophic_groups/diversity/forage_dt_diversity.csv"

forage_dt_abun = getDomainAbundance(RVCdata_DT, forage_fish_spp_list, merge_protected = F)
write_csv(forage_dt_abun, "big_csv/trophic_groups/abundance/forage_dt_abun.csv")

forage_dt_den = getDomainDensity(RVCdata_DT, forage_fish_spp_list, merge_protected = F)
write_csv(forage_dt_den, "big_csv/trophic_groups/density/forage_dt_den.csv")

# dt diversity 
forage_dt_abundance = forage_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)
forage_dt_diversity = forage_dt_abundance %>%
  select(YEAR, protected_status, richness, shannon, simpson) 
write_csv(forage_dt_diversity, "trophic_groups/diversity/forage_dt_diversity.csv")


grouper_dt_abun_csv     = "trophic_groups/abundance/grouper_dt_abun.csv"
grouper_dt_den_csv      = "trophic_groups/density/grouper_dt_den.csv"
grouper_dt_diversity_csv= "trophic_groups/diversity/grouper_dt_diversity.csv"

grouper_dt_abun = getDomainAbundance(RVCdata_DT, grouper_spp_list, merge_protected = F)
write_csv(grouper_dt_abun, "big_csv/trophic_groups/abundance/grouper_dt_abun.csv")
grouper_dt_den = getDomainDensity(RVCdata_DT, grouper_spp_list)
write_csv(grouper_dt_den, "big_csv/trophic_groups/density/grouper_dt_den.csv")

#diversity metrics 
grouper_dt_abundance = grouper_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill =0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

grouper_dt_diversity = grouper_dt_abundance %>%
  select(YEAR, protected_status, richness, shannon, simpson) 
write_csv(grouper_dt_diversity, "big_csv/trophic_groups/diversity/grouper_dt_diversity.csv")

grunt_dt_abun_csv     = "trophic_groups/abundance/grunt_dt_abun.csv"
grunt_dt_den_csv      = "trophic_groups/density/grunt_dt_den.csv"
grunt_dt_diversity_csv= "trophic_groups/diversity/grunt_dt_diversity.csv"

grunt_dt_abun = getDomainAbundance(RVCdata_DT, grunt_spp_list, merge_protected = F)
write_csv(grunt_dt_abun, "trophic_groups/abundance/grunt_dt_abun.csv")

grunt_dt_den = getDomainDensity(RVCdata_DT, grunt_spp_list, merge_protected = F)
write_csv(grunt_dt_den, "trophic_groups/density/grunt_dt_den.csv")

#dt diversity
grunt_dt_abundance = grunt_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

grunt_dt_diversity = grunt_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(grunt_dt_diversity, "trophic_groups/diversity/grunt_dt_diversity.csv")

higher_reef_dt_abun_csv     = "trophic_groups/abundance/higher_reef_dt_abun.csv"
higher_reef_dt_den_csv      = "trophic_groups/density/higher_reef_dt_den.csv"
higher_reef_dt_diversity_csv= "trophic_groups/diversity/higher_reef_dt_diversity.csv"

higher_reef_dt_abun = getDomainAbundance(RVCdata_DT, higher_reef_spp_list, merge_protected = F)
write_csv(higher_reef_dt_abun, "trophic_groups/abundance/higher_reef_dt_abun.csv")
higher_reef_dt_den = getDomainDensity(RVCdata_DT, higher_reef_spp_list, merge_protected = F)
write_csv(higher_reef_dt_den, "trophic_groups/density/higher_reef_dt_den.csv")

#dt diversity 
higher_reef_dt_abundance = higher_reef_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

higher_reef_diversity = higher_reef_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(higher_reef_diversity, "trophic_groups/diversity/higher_reef_dt_diversity.csv")

hogfish_dt_den_csv      = "trophic_groups/density/hogfish_dt_den.csv"
hogfish_fk_diversity_csv= "trophic_groups/diversity/hogfish_fk_diversity.csv"
hogfish_dt_diversity_csv= "trophic_groups/diversity/hogfish_dt_diversity.csv"

hogfish_dt_abun = getDomainAbundance(RVCdata_DT, hogfish_spp_list, merge_protected = F)
write_csv(hogfish_dt_abun, "big_csv/trophic_groups/abundance/hogfish_dt_abun.csv")

hogfish_dt_den = getDomainDensity(RVCdata_DT, hogfish_spp_list, merge_protected = F)
write_csv(hogfish_dt_den, "big_csv/trophic_groups/density/hogfish_dt_den.csv")

#DT diversity
hogfish_dt_abundance = hogfish_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

hogfish_dt_diversity = hogfish_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(hogfish_dt_diversity, "trophic_groups/diversity/hogfish_dt_diversity.csv")

lionfish_dt_abun_csv     = "trophic_groups/abundance/lionfish_dt_abun.csv"
lionfish_dt_den_csv      = "trophic_groups/density/lionfish_dt_den.csv"
lionfish_dt_abun = getDomainAbundance(RVCdata_DT, "PTE VOLI", merge_protected = F)
write_csv(lionfish_dt_abun, "big_csv/trophic_groups/abundance/lionfish_dt_abun.csv")
lionfish_dt_den = getDomainDensity(RVCdata_DT, "PTE VOLI", merge_protected = F)
write_csv(lionfish_dt_den, "big_csv/trophic_groups/density/lionfish_dt_den.csv")

opportunists_dt_abun_csv     = "trophic_groups/abundance/opportunists_dt_abun.csv"
opportunists_dt_den_csv      = "trophic_groups/density/opportunists_dt_den.csv"
opportunists_dt_diversity_csv= "trophic_groups/diversity/opportunists_dt_diversity.csv"
opportunists_dt_den = getDomainDensity(RVCdata_DT, opportunists_spp_list, merge_protected = F)
write_csv(opportunists_dt_den, "trophic_groups/density/opportunists_dt_den.csv")
opportunists_dt_abun = getDomainAbundance(RVCdata_DT, opportunists_spp_list, merge_protected = F)
write_csv(opportunists_dt_abun, "trophic_groups/abundance/opportunists_dt_abun.csv")
#dt diversity
opportunists_dt_abundance = opportunists_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

opportunists_dt_diversity = opportunists_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(opportunists_dt_diversity, "big_csv/trophic_groups/diversity/opportunists_dt_diversity.csv")

parrotfish_dt_diversity_csv = "trophic_groups/diversity/parrotfish_dt_diversity.csv"
parrotfish_dt_abun_csv      = "trophic_groups/abundance/parrotfish_dt_abun.csv"
parrotfish_dt_den_csv       = "trophic_groups/density/parrotfish_dt_den.csv"
parrotfish_dt_abun = getDomainAbundance(RVCdata_DT, parrotfish_spp_list, merge_protected = F)
write_csv(parrotfish_dt_abun, "big_csv/trophic_groups/abundance/parrotfish_dt_abun.csv")
parrotfish_dt_den = getDomainDensity(RVCdata_DT, parrotfish_spp_list, merge_protected = F)
write_csv(parrotfish_dt_den, "big_csv/trophic_groups/density/parrotfish_dt_den.csv")
#dt diversity
parrotfish_dt_abundance = parrotfish_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

parrotfish_dt_diversity = parrotfish_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(parrotfish_dt_diversity, "big_csv/trophic_groups/diversity/parrotfish_dt_diversity.csv")

rays_dt_abun_csv     = "trophic_groups/abundance/rays_dt_abun.csv"
rays_dt_den_csv      = "trophic_groups/density/rays_dt_den.csv"
rays_dt_diversity_csv= "trophic_groups/diversity/rays_dt_diversity.csv"
#diversity
rays_fk_abundance = rays_fk_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

rays_fk_diversity = rays_fk_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(rays_fk_diversity, "big_csv/trophic_groups/diversity/rays_fk_diversity.csv")

#dt diversity
rays_dt_abundance = rays_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

rays_dt_diversity = rays_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(rays_dt_diversity, "big_csv/trophic_groups/diversity/rays_dt_diversity.csv")
rays_dt_den = getDomainDensity(RVCdata_DT, rays_spp_list, merge_protected = F)
write_csv(rays_dt_den, "big_csv/trophic_groups/density/rays_dt_den.csv")
rays_dt_abun = getDomainAbundance(RVCdata_DT, rays_spp_list, merge_protected = F)
write_csv(rays_dt_abun, "big_csv/trophic_groups/abundance/rays_dt_abun.csv")

seabass_dt_abun_csv     = "trophic_groups/abundance/seabass_dt_abun.csv"
seabass_dt_den_csv      = "trophic_groups/density/seabass_dt_den.csv"
seabass_dt_diversity_csv= "trophic_groups/diversity/seabass_dt_diversity.csv"

seabass_dt_abun = getDomainAbundance(RVCdata_DT, sea_bass_hamlet_spp_list, merge_protected = F)
write_csv(seabass_dt_abun, "big_csv/trophic_groups/abundance/seabass_dt_abun.csv")
seabass_dt_den = getDomainDensity(RVCdata_DT, sea_bass_hamlet_spp_list, merge_protected = F)
write_csv(seabass_dt_den, "big_csv/trophic_groups/density/seabass_dt_den.csv")
#dt diversity
seabass_dt_abundance = seabass_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

seabass_dt_diversity = seabass_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(seabass_dt_diversity, "big_csv/trophic_groups/diversity/seabass_dt_diversity.csv")

sharks_dt_abun_csv      = "trophic_groups/abundance/sharks_dt_abun.csv"
sharks_dt_den_csv       = "trophic_groups/density/sharks_dt_den.csv"
sharks_fk_diversity_csv = "trophic_groups/diversity/sharks_fk_diversity.csv"
sharks_dt_diversity_csv = "trophic_groups/diversity/sharks_dt_diversity.csv"

sharks_dt_abun = getDomainAbundance(RVCdata_DT, sharks_spp_list, merge_protected = F)
write_csv(sharks_dt_abun, "big_csv/trophic_groups/abundance/sharks_dt_abun.csv")
sharks_dt_den = getDomainDensity(RVCdata_DT, sharks_spp_list, merge_protected = F)
write_csv(sharks_dt_den, "big_csv/trophic_groups/density/sharks_dt_den.csv")
#diversity
sharks_fk_abundance = sharks_fk_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

sharks_fk_diversity = sharks_fk_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(sharks_fk_diversity, "big_csv/trophic_groups/diversity/sharks_fk_diversity.csv")

#dt diversity
sharks_dt_abundance = sharks_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

sharks_dt_diversity = sharks_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(sharks_dt_diversity, "big_csv/trophic_groups/diversity/sharks_dt_diversity.csv")

small_reef_dt_abun_csv      = "trophic_groups/abundance/small_reef_dt_abun.csv"
small_reef_dt_den_csv       = "trophic_groups/density/small_reef_dt_den.csv"
small_reef_dt_diversity_csv = "trophic_groups/diversity/small_reef_dt_diversity.csv"

small_reef_dt_abun = getDomainAbundance(RVCdata_DT, small_reef_fish_spp_list, merge_protected = F)
write_csv(small_reef_dt_abun, "big_csv/trophic_groups/abundance/small_reef_dt_abun.csv")
small_reef_dt_den = getDomainDensity(RVCdata_DT, small_reef_fish_spp_list, merge_protected = F)
write_csv(small_reef_dt_den, "big_csv/trophic_groups/density/small_reef_dt_den.csv")
#dt diversity
small_reef_dt_abundance = small_reef_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

small_reef_dt_diversity = small_reef_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(small_reef_dt_diversity, "big_csv/trophic_groups/diversity/small_reef_dt_diversity.csv")

snapper_dt_abun_csv      = "trophic_groups/abundance/snapper_dt_abun.csv"
snapper_dt_den_csv       = "trophic_groups/density/snapper_dt_den.csv"
snapper_dt_diversity_csv = "trophic_groups/diversity/snapper_dt_diversity.csv"

snapper_dt_abun = getDomainAbundance(RVCdata_DT, snapper_spp_list, merge_protected = F)
write_csv(snapper_dt_abun, "big_csv/trophic_groups/abundance/snapper_dt_abun.csv")
snapper_dt_den = getDomainDensity(RVCdata_DT, snapper_spp_list, merge_protected = F)
write_csv(snapper_dt_den, "big_csv/trophic_groups/density/snapper_dt_den.csv")
#dt diversity
snapper_dt_abundance = snapper_dt_abun %>% 
  select(YEAR, protected_status, SPECIES_CD, abundance) %>%  
  group_by(YEAR, protected_status) %>% #48 groups 
  nest(-YEAR) %>%  
  mutate(
    data_wide = map(data, ~ spread(data=.x, SPECIES_CD, abundance, fill = 0)), 
    #species richness 
    richness = map(
      data_wide, 
      function(x) specnumber(x)),
    #simpson diversity as effective number of species 
    simpson = map( 
      data_wide,
      function(x) 1/(1 - diversity(x, index = 'simpson'))),
    #shannon diversity as effective number of species 
    shannon = map(
      data_wide,
      function(x) exp(diversity(x, index = 'shannon')))) %>%
  unnest(richness, simpson, shannon)

snapper_dt_diversity = snapper_dt_abundance %>%
  select(YEAR, protected_status, richness, simpson, shannon) 
write_csv(snapper_dt_diversity, "big_csv/trophic_groups/diversity/snapper_dt_diversity.csv")

parrotfish_fk_abun = read_csv('trophic_groups/abundance/parrotfish_fk_abun.csv')

parrotfish_fk_abundance = parrotfish_fk_abun %>%
  select(YEAR, protected_status, SPECIES_CD, abundance) %>% 
  group_by(YEAR, protected_status) %>%
  summarise(
    ntot = sum(abundance)) %>%
  filter(protected_status != "0") %>%
  filter(protected_status != "1") %>% 
  mutate(
    protected_status = recode(
      protected_status, '0'='unprotected', '1'='protected')) %>%
  spread(protected_status, ntot)

dygraph(parrotfish_fk_abundance, main = 'Parrotfish Abundance') %>% 
  dyAxis("y", label = "Abundance") %>%
  dyAxis("x", label = "Year") %>%
  dyOptions(stackedGraph = T)

#mean density 
snapper_fk_mean_density = snapper_fk_den %>%
  select(YEAR, protected_status, SPECIES_CD, density) %>% 
  group_by(YEAR, protected_status) %>%
  summarise(
    meanden = mean(density)) %>%
  filter(protected_status != "0") %>%
  filter(protected_status != "1") %>% 
  mutate(
    protected_status = recode(
      protected_status, '0'='unprotected', '1'='protected')) %>%
  spread(protected_status, meanden)

dygraph(snapper_fk_mean_density, main = 'Snapper Density') %>% 
  dyAxis("y", label = "Density") %>%
  dyAxis("x", label = "Year") %>%
  dyOptions(stackedGraph = T, fillAlpha = 0.6, axisLineWidth = 2)

exploited_dt_abun_csv = "exploited_species/exploited_dt_abun.csv"
exploited_dt_den_csv  = "exploited_species/exploited_dt_den.csv"
exploited_dt_abun <- getDomainAbundance(RVCdata_DT, species = exploited_spp_list, merge_protected = F)
write_csv(exploited_dt_abun,"big_csv/exploited/exploited_dt_abun.csv")
exploited_dt_den <- getDomainDensity(RVCdata_DT, species = exploited_spp_list, merge_protected = F)
write_csv(exploited_dt_den,"big_csv/exploited/exploited_dt_den.csv")

domain_dt_abun_diversity = read_csv("big_csv/abundance_domain/domain_dt_abun_diversity.csv")
strata_dt_abun_diversity = read_csv("big_csv/abundance_strata/strata_dt_abun_diversity.csv")
psu_dt_abun_diversity = read_csv("big_csv/abundance_psu/psu_dt_abun_diversity.csv")