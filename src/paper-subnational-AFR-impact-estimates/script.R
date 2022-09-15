## Load all function files
R.utils::sourceDirectory("R", modifiedOnly = FALSE)


# 1. INPUT DATA ----
## 1.1 Coverage data -----

# MCV1 coverage
MCV1_cov <- readRDS("./dat/MCV1_vax_cov_2000_2019_Africa_gadm36.RDS")
MCV1_cov_samples <- readRDS("./dat/mcv1_cov_draws_vimc_countries.RDS")

# YF coverage
YF_cov <- readRDS("./dat/YF_vax_cov_2000_2019_Africa_gadm36.RDS")

# Get DTP3 coverage
DTP13 <- readRDS("./dat/DTP1_3_vax_cov_2000_2016_Africa_gadm36.RDS")

# 1.2. Population data (subnational values) ----
pop_dat <- readRDS("dat/Pop_ADMIN1_Africa_gadm36.RDS")

# Filters
pop_dat$country = pop_dat$iso3
countries_to_include <- unique(pop_dat$country)
period = 2000:2030

## Read live births data
live_births <- readRDS("./dat/live_birth_under_allcase_under5_mort.rds")

# 1.3. Impact estimates cohort view  ----
im_cohort <- readRDS("./dat/impact-ratios-cohort-act-strat.rds")

# 1.4. Admin1 maps and general shapefiles ----
adm1_map <- readRDS("./dat/AFR_map_smooth.RDS")
AF_map <- readRDS("./dat/AF_shp.RDS")

## get standarise names 
adm1_names_var <- read.csv("dat/adm1_names_variations.csv")

### Include country names
country_names <- 
  tbl(con, "country") %>% 
  distinct(id,name) %>% 
  collect() %>% 
  dplyr::rename(iso3 = id) %>% 
  mutate(name = ifelse(iso3 == "COD", "Congo, DR", name)) %>% 
  mutate(name = ifelse(iso3 == "TZA", "Tanzania", name)) 


# 2. STANDARISED DATA ----
## 2.1 Standarising names for population data ----
pop_dat_age0 <-
  pop_dat %>% 
  ungroup() %>% 
  filter(age_from == 0 & age_to == 0) %>% 
  select(iso3, ADM1_NAME = region, year, total_pop, unadjusted) %>% 
  stand_name() %>% 
  group_by(iso3, ADM1_NAME, year) %>% 
  summarise(total_pop = sum(total_pop, na.rm = TRUE), 
            unadjusted = sum(unadjusted, na.rm = TRUE))

pop_dat_age0 <- 
  pop_dat_age0 %>% 
  ungroup() %>% 
  group_by(iso3, year) %>% 
  mutate(region_code = paste(iso3, 1:length(ADM1_NAME), sep = "_")) %>% 
  rename(ADM1_NAME_VIMC = ADM1_NAME) %>% 
  rename(ADM1_NAME = region_code)


pop_dat_total <-
  pop_dat %>% 
  ungroup() %>% 
  select(iso3, ADM1_NAME = region, year, total_pop, unadjusted) %>% 
  stand_name() %>% 
  group_by(iso3, ADM1_NAME, year) %>% 
  summarise(total_pop = sum(total_pop, na.rm = TRUE), 
            unadjusted = sum(unadjusted, na.rm = TRUE))

pop_dat_total <- 
  pop_dat_total %>% 
  ungroup() %>% 
  group_by(iso3, year) %>% 
  mutate(region_code = paste(iso3, 1:length(ADM1_NAME), sep = "_")) %>% 
  rename(ADM1_NAME_VIMC = ADM1_NAME) %>% 
  rename(ADM1_NAME = region_code)


## Adjusted live births
pop_live_births <- 
  pop_dat_age0 %>% 
  left_join(
    live_births %>% 
      select(iso3, year = cohort, wpp_live_b = pop)) %>% 
  group_by(iso3, year) %>% 
  mutate(total_pop_age0 = sum(total_pop, na.rm = TRUE)) %>% 
  #normalise live births
  mutate(weighted_pop = (total_pop/total_pop_age0)*wpp_live_b) %>% 
  select(iso3, ADM1_NAME, ADM1_NAME_VIMC, year, adjusted_liveB = weighted_pop)


pop_dat <- list(age0 = pop_dat_age0, total = pop_dat_total, live_births = pop_live_births)
saveRDS(pop_dat, "processed_pop_dat.RDS")

## 2.2 Standarising names for cov data ----
## MCV1 coverage (fix some adm1 names)
MCV1_cov2 <- 
  MCV1_cov %>% 
  mutate(ADM1_NAME = stringi::stri_trans_general(ADM1_NAME, "latin-ascii")) %>% 
  mutate(ADM1_NAME = tolower(ADM1_NAME)) %>% 
  mutate(ADM1_NAME = gsub("-", " ", ADM1_NAME)) %>% 
  mutate(ADM1_NAME = gsub("&", "and", ADM1_NAME)) %>% 
  group_by(iso3, year) %>% arrange(ADM1_NAME) %>% 
  mutate(region_code = paste(iso3, 1:length(ADM1_NAME), sep = "_")) %>% 
  filter(iso3%in%countries_to_include) %>% 
  rename(ADM1_NAME_IHME = ADM1_NAME) %>% 
  rename(ADM1_NAME = region_code)


## coverage samples 
cov_names <-
MCV1_cov2 %>% 
  ungroup() %>% 
  distinct(ADM1_CODE, iso3, ADM1_NAME)

MCV1_cov_samples_unad <- 
  MCV1_cov_samples %>% 
  mutate(ADM1_NAME_orig = ADM1_NAME) %>% 
  mutate(ADM1_CODE = as.character(ADM1_CODE)) %>% 
  select(-ADM1_NAME) %>% 
  left_join(cov_names) %>% 
  pivot_longer(names_to = "sample", cols = starts_with("V"), values_to = "value")
  
# 3. GET FVPS ----

## National fvps from 201910gavi coverage
fvps_VIMC <- vimpact::extract_vaccination_history(con, 
                                                  touchstone_cov = "201910gavi", 
                                                  year_min = 1990, 
                                                  year_max = 2030,
                                                  gavi_support_levels = c("with"), 
                                                  countries_to_extract = countries_to_include, 
                                                  disease_to_extract = c("Measles","YF","HepB"))


fvps_VIMC_mcv1 <- 
  fvps_VIMC %>% 
  filter(vaccine == "MCV1") %>% 
  select(iso3 = country, year, nat_cov = coverage_adjusted) %>% 
  filter(iso3%in%countries_to_include)


# 4. ARRANGE NAT COV ----

## compare national coverage data from IHME and VIMC
IHME_cov <- 
  MCV1_cov2 %>% 
  group_by(iso3, year) %>% 
  summarise(IHME_cov = round(mean(mean_cov, na.rm = TRUE),2))

cov_comparison <- 
  fvps_VIMC_mcv1 %>% 
  filter(year%in%2000:2019) %>% 
  rename(VIMC_cov = nat_cov) %>% 
  left_join(IHME_cov) %>% 
  pivot_longer(cols= contains("cov"), names_to = "source", values_to = "cov")

saveRDS(cov_comparison, "national_coverage_IHME_vs_201910gavi.RDS")


# 5. ADJUST SUBNAT COV ----
MCV1_cov2_new <- normalise_fvps_MCV1(cov_dat = MCV1_cov2,
                                     fvps_dat = fvps_VIMC,
                                     pop_dat = pop_dat_age0)

MCV1_cov2_new <- 
  MCV1_cov2_new  %>%  
  mutate(mean_cov  = ifelse(mean_cov > 1, 1, mean_cov)) %>% 
  filter(!is.na(year)) %>% 
  filter(!is.na(mean_cov)) %>%
  filter(iso3!= "STP") %>% 
  left_join(MCV1_cov2 %>% 
              select(iso3, ADM1_NAME, year, unadjusted_cov = mean_cov))


## include total population
MCV1_cov2_new <- MCV1_cov2_new %>% left_join(pop_dat_age0 %>% select(iso3, ADM1_NAME, year, total_pop))

saveRDS(MCV1_cov2_new, "normalised_fvps_coverage.RDS")

# Adjust coverage samples
MCV1_cov_samples_adjus <- normalise_fvps_MCV1(cov_dat = MCV1_cov_samples_unad,
                                              fvps_dat = fvps_VIMC,
                                              pop_dat = pop_dat_age0, 
                                              has_samples = TRUE)

MCV1_cov_samples_adjus <- 
MCV1_cov_samples_adjus %>% 
  group_by(iso3, ADM1_NAME, ADM1_NAME_orig, year) %>% 
  summarise(mean_cov = round(mean(adj_cov, na.rm = TRUE), 3), 
            lo_cov = round(quantile(adj_cov, 0.025, na.rm = TRUE), 3), 
            hi_cov = round(quantile(adj_cov, 0.975, na.rm = TRUE), 3))


saveRDS(MCV1_cov_samples_adjus, "normalised_coverage_uncertainty.RDS")

# 6. GET IMPACT RATIOS ----

## Stratified by cohort and activity
stratifications <- c("un-stratified", "activity", "cohort", "both")
fvps_VIMC_tmp <- 
  fvps_VIMC %>% 
  filter(vaccine == "MCV1") %>% 
  mutate(fvps = fvps_adjusted, 
         cohort = year - age)

## Compare the impact ratios
da_rate_all <- 
  foreach(i =  seq(stratifications), .combine = rbind)%do%{
    
    dat_tmp <- prepare_impact_ratio(d = im_cohort, fvps = fvps_VIMC_tmp, stratifications[i], period)
    dat_tmp$stratification = stratifications[i]
    dat_tmp
  }

## Using impact ratios stratified by cohort and activity type 
da_rate_activity_cohort <- 
  da_rate_all %>% 
  filter(stratification == "both") %>% 
  filter(vaccine == "MCV1") %>% 
  filter(country%in%countries_to_include) %>% 
  mutate(iso3 = country) %>% 
  mutate(deaths_averted_rate = impact_ratio)


## implement the new subnational impact ratios approach
national_impact_ratios <- 
  da_rate_activity_cohort %>% 
  select(boots_id, year, iso3, impact_ratio, impact, vaccine)

all_dat <- 
  MCV1_cov2_new %>% 
  select(-ADM0_NAME, -ADM1_CODE) %>% ## subnational coverage
  left_join(national_impact_ratios) %>% ## national impact ratios
  left_join(fvps_VIMC_mcv1) ## include national coverage

da_rate_subnational <- 
  all_dat %>% 
  group_by(boots_id, iso3) %>% 
  mutate(max_propn_change = max(impact_ratio, na.rm = TRUE)-min(impact_ratio,  na.rm = TRUE)) %>% 
  mutate(dis_impact_ratio = impact_ratio*(1 - ((mean_cov - nat_cov)/nat_cov)*max_propn_change)) 

da_rate_subnational_tmp <- 
  da_rate_subnational %>% 
  distinct(vaccine, iso3, ADM1_NAME, year, dis_impact_ratio, boots_id) %>% 
  left_join(da_rate_activity_cohort) %>% 
  mutate(deaths_averted_rate = dis_impact_ratio)

da_rate_activity_cohort = da_rate_subnational_tmp

vimc_impact_ratios_dist <- 
  da_rate_subnational %>% 
  filter(iso3!="IND") %>%
  filter(year %in% c(2000:2019)) %>% 
  group_by(iso3, year, vaccine) %>% 
  summarise(max_propn_change_mean = mean(max_propn_change, na.rm = TRUE), 
            max_propn_change_lo =  quantile(max_propn_change, 0.025, na.rm = TRUE), 
            max_propn_change_hi = quantile(max_propn_change, 0.975, na.rm = TRUE)) %>% 
  left_join(country_names)

saveRDS(vimc_impact_ratios_dist, file = "deaths_averted_rate_ranges.RDS")

## get table of impact ratios for the paper
vimc_impact_ratios <- 
  da_rate_activity_cohort %>% 
  filter(iso3!="IND") %>%
  filter(year %in% c(2000:2019)) %>% 
  group_by(iso3, year, vaccine) %>% 
  summarise(deaths_averted_rate_mean = mean(deaths_averted_rate, na.rm = TRUE), 
            deaths_averted_rate_lo =  quantile(deaths_averted_rate, 0.025, na.rm = TRUE), 
            deaths_averted_rate_hi = quantile(deaths_averted_rate, 0.975, na.rm = TRUE)) %>% 
  left_join(country_names) %>% 
  select(Country = name, 
         iso_code = iso3, 
         Year = year, 
         deaths_averted_rate_mean, deaths_averted_rate_lo, deaths_averted_rate_hi)

write.csv(vimc_impact_ratios, file = "deaths_averted_rate_2000_2019.csv", row.names = FALSE)

### Impact ratios per period
vimc_impact_ratios_tab <- 
  vimc_impact_ratios %>% 
  filter(Year %in% c(2000, 2010, 2019)) %>% 
  mutate_at(vars(contains("deaths_averted")), ~round(., 3)) %>% 
  mutate(values = paste0(deaths_averted_rate_mean, 
                         " (", deaths_averted_rate_lo, "–", 
                         deaths_averted_rate_hi, ")")) %>% 
  select(Country, iso_code, Year, values) %>% 
  pivot_wider(names_from = Year, values_from = values) %>% 
  arrange(Country)

write.csv(vimc_impact_ratios_tab, file = "deaths_averted_rate_period.csv", row.names = FALSE)

print(xtable(vimc_impact_ratios_tab, type = "latex"), 
      file = "deaths_averted_rate_period.tex", 
      include.rownames=FALSE)

# 7. GET FVPS UNDER DIFF SCENARIOS ----

# Using the standarised subnational coverage, get the fvps under the different scenarios
options_l <- c("sub_cov", "at_nat_cov", "max_cov", "at_gvap_cov")

## Get the fvps under the differente coverage scenarios
fvps_scen_df <- foreach(i =  seq(options_l), .combine = rbind)%do%{
  
  print(paste("Get estimates for", options_l[i]))
  
  fvps_H <- get_fvps_under_H(sub_cov = MCV1_cov2_new, 
                             vimc_nat_cov = fvps_VIMC_mcv1, 
                             sub_pop = pop_dat_age0, 
                             DTP3_cov = DTP3_cov, 
                             option_cov = options_l[i])
  
  
  fvps_H %>% 
    filter(!is.na(mean_cov))
  
}

# Save coverage
cov_summary <- 
  MCV1_cov2_new %>% 
  filter(year == 2019) %>% 
  group_by(iso3) %>% 
  summarise(max_cov_y = max(mean_cov, na.rm = TRUE)) %>% 
  mutate(max_cov_y = ifelse(max_cov_y > 1, 1, max_cov_y)) %>% 
  left_join(fvps_VIMC_mcv1 %>% filter(year == 2019))

cov_scenarios <- 
  fvps_scen_df %>%
  filter(year == 2019) %>%
  group_by(iso3, year, scenario) %>%
  summarise(total_fvps = sum(mean_fvps, na.rm = TRUE),
            tot_pop = sum(total_pop, na.rm = TRUE)) %>%
  mutate(new_cov = round(total_fvps/tot_pop, 2)) %>% 
  select(iso3, year, scenario, new_cov)

cov_scenarios_global <- 
  fvps_scen_df %>%
  filter(year == 2019) %>%
  group_by(year, scenario) %>%
  summarise(total_fvps = sum(mean_fvps, na.rm = TRUE),
            tot_pop = sum(total_pop, na.rm = TRUE)) %>%
  mutate(new_cov = round(total_fvps/tot_pop, 2)) %>% 
  select(year, scenario, new_cov)

cov_all <- list(cov_summary = cov_summary,
                cov_scenarios = cov_scenarios, 
                cov_scenarios_global = cov_scenarios_global)

saveRDS(cov_all, "coverage_under_scenarios.RDS")

# 8. CALCULATE IMPACT ----

## Calculate impact per country per year per scenario of inequality
impact_df_cohort <- IU_impact_calculation(new_fvps = fvps_scen_df, 
                                          impact_rates = da_rate_activity_cohort, 
                                          period = NULL)  
impact_df_cohort <- 
  impact_df_cohort %>% 
  filter(iso3!="IND") %>%
  filter(iso3!="STP")

## Calculate average impact across boostrap samples
vimc_impact <- 
  da_rate_activity_cohort %>% 
  group_by(iso3, year, vaccine) %>% 
  summarise(deaths_averted_new_mean = mean(impact, na.rm = TRUE), 
            deaths_averted_new_q1 =  quantile(impact, 0.025, na.rm = TRUE), 
            deaths_averted_new_q3 = quantile(impact, 0.975, na.rm = TRUE)) %>% 
  mutate(scenario = "VIMC")

# Include live births into the impact data
impact_df_cohort <- 
  rbind(impact_df_cohort, vimc_impact) %>% 
  left_join(live_births %>% mutate(year = cohort))

saveRDS(impact_df_cohort, "impact_inequality_scenarios.RDS")



# 9. RELATIVE CHANGES IN IMPACT ----

## 9.1 National estimates ----
## Impact and relative change (2000-2019 national estimates)
rel_ch_df_both <- rel_change_impact(new_fvps = fvps_scen_df,
                                    impact_rates = da_rate_activity_cohort,
                                    period = 2000:2019, 
                                    ref_scenario = "sub_cov")

## Impact and relative change (only 2019 national estimates)
rel_ch_2019 <- rel_change_impact(new_fvps = fvps_scen_df,
                                 impact_rates = da_rate_activity_cohort,
                                 period = 2019, 
                                 ref_scenario = "sub_cov")


## calculate impact per periods
periods_l <- list(p0 = 2019, p1 = 2000:2009, p2 = 2010:2019)

## Impact periods national estimates
all_dat_periods <- 
  
  foreach(i = seq_along(periods_l), .combine = rbind)%do%{
    
    impact_tmp <- IU_impact_calculation(new_fvps = fvps_scen_df, 
                                        impact_rates = da_rate_activity_cohort,
                                        period = periods_l[[i]])
    
    pop_tmp <- 
      live_births %>% 
      mutate(year = cohort) %>% 
      filter(year %in% periods_l[[i]]) %>% 
      group_by(iso3) %>% 
      summarise(total_pop = sum(pop, na.rm = TRUE))
    
    impact_tmp %>% 
      left_join(pop_tmp)
  }


## 9.2 subnational estimates ----
## Impact and relative change (2000-2019 subnational (adm1) estimates)
rel_ch_adm1_both <- rel_change_impact(new_fvps = fvps_scen_df,
                                      impact_rates = da_rate_activity_cohort,
                                      ref_scenario = "sub_cov", 
                                      period = 2000:2019, 
                                      per_adm1 = TRUE)

## Impact and relative change (only 2019 subnational (adm1) estimates)
rel_ch_adm1_2019_both <- rel_change_impact(new_fvps = fvps_scen_df,
                                           impact_rates = da_rate_activity_cohort,
                                           ref_scenario = "sub_cov", 
                                           period = 2019, 
                                           per_adm1 = TRUE)


## Impact periods subnational (adm1) estimates
all_dat_periods_adm1 <- 
  
  foreach(i = seq_along(periods_l), .combine = rbind)%do%{
    
    impact_tmp <- IU_impact_calculation(new_fvps = fvps_scen_df, 
                                        impact_rates = da_rate_activity_cohort,
                                        period = periods_l[[i]], 
                                        per_adm1 = TRUE)
    
    pop_tmp <- 
      pop_live_births %>% 
      filter(year %in% periods_l[[i]]) %>% 
      group_by(iso3, ADM1_NAME, ADM1_NAME_VIMC) %>% 
      summarise(total_pop = sum(adjusted_liveB, na.rm = TRUE))
    
    impact_tmp %>% 
      left_join(pop_tmp)
  }


## 9.3 Global estimates ----
all_dat_periods_global <- 
  
  foreach(i = seq_along(periods_l), .combine = rbind)%do%{
    
    impact_tmp <- IU_impact_calculation(new_fvps = fvps_scen_df %>% filter(iso3!= "IND"), 
                                        impact_rates = da_rate_activity_cohort %>% filter(iso3!= "IND"),
                                        period = periods_l[[i]], 
                                        global = TRUE)
    
    pop_tmp <- 
      pop_live_births %>% 
      filter(year %in% periods_l[[i]]) %>% 
      mutate(vaccine = "MCV1") %>% 
      group_by(vaccine) %>% 
      summarise(total_pop = sum(adjusted_liveB, na.rm = TRUE))
    
    impact_tmp %>% 
      left_join(pop_tmp)
  }



glob_rel_ch <- list(p_split_2000_2019 = all_dat_periods_global)


nat_rel_ch <- list(p_2000_2019 = rel_ch_df_both, 
                   p_2019 = rel_ch_2019, 
                   p_split_2000_2019 = all_dat_periods)

sub_rel_ch <- list(p_2000_2019 = rel_ch_adm1_both, 
                   p_2019 = rel_ch_adm1_2019_both, 
                   p_split_2000_2019 = all_dat_periods_adm1)


saveRDS(glob_rel_ch, "glob_rel_ch.RDS")
saveRDS(nat_rel_ch, "nat_rel_ch.RDS")
saveRDS(sub_rel_ch, "sub_rel_ch.RDS")


## 9.3 Impact estimates per year and per adm1 
tab_df_admin1_py <- IU_impact_calculation(new_fvps = fvps_scen_df %>% filter(scenario == "sub_cov"), 
                                          impact_rates = da_rate_activity_cohort,
                                          per_adm1 = TRUE)


tab_df_admin1_py_df <- 
  tab_df_admin1_py %>% 
  left_join(all_dat_periods_adm1 %>% distinct(iso3, ADM1_NAME, ADM1_NAME_VIMC)) %>% 
  left_join(country_names) %>% 
  ungroup() %>% 
  select(country = name,
         iso3,
         ADM1_NAME = ADM1_NAME_VIMC, 
         ADM1_CODE = ADM1_NAME, 
         vaccine, year, contains("deaths_averted")) %>% 
  mutate(ADM1_NAME = str_to_title(ADM1_NAME)) %>% 
  arrange(iso3, ADM1_NAME) %>% 
  ## include population
  left_join(pop_live_births %>% rename(total_pop = adjusted_liveB, ADM1_CODE = ADM1_NAME) %>% select(-ADM1_NAME_VIMC))

write.csv(tab_df_admin1_py_df, "impact_ADM1_per_year.csv", row.names = FALSE)

# Live births pop data
pop_tmp_2019 <- 
  pop_live_births %>% 
  filter(year %in% 2019) %>% 
  group_by(iso3) %>% 
  summarise(total_pop = sum(adjusted_liveB, na.rm = TRUE))

ch_df_both <- 
  rel_ch_2019$rel_ch %>% 
  filter(iso3!="IND") %>% 
  filter(!is.na(vaccine)) %>% 
  filter(rel_ch != Inf) %>% 
  left_join(pop_tmp_2019) %>% 
  left_join(country_names) 

saveRDS(ch_df_both, "Change_impact_national_2019.RDS")

# 10. DISSIMILARITY INDEX -----
## Get indices of dissimilarity
diss_dat <- get_dissimilarity_index(new_fvps = fvps_scen_df,
                                    impact_rates = da_rate_activity_cohort,
                                    pop_dat = pop_live_births %>% rename(total_pop = adjusted_liveB))


CI_dat <-
  diss_dat %>%
  filter(iso3 != "IND") %>%
  group_by(iso3, year) %>%
  summarise_at(vars(matches("WID|INS")), 
               list(mean = ~ mean(., na.rm = TRUE), 
                    lo = ~ quantile(., 0.025, na.rm = TRUE), 
                    hi = ~ quantile(., 0.975, na.rm = TRUE))) 

saveRDS(CI_dat, "dissimilarity_2000_2019.RDS")

# 11. DIF IMPACT FROM SAMPLES ---- 
diff_tb_2019 <- dif_impact_from_samples(new_fvps = fvps_scen_df, 
                                        impact_rates = da_rate_activity_cohort, 
                                        period = 2019, 
                                        get_boots = TRUE)

diff_tb_2000_2019 <- 
  lapply(list(2000:2009, 2010:2019), function(x)
    dif_impact_from_samples(new_fvps = fvps_scen_df, 
                            impact_rates = da_rate_activity_cohort, 
                            period = x, 
                            get_boots = TRUE)
  )

diff_tb_2000_2019 <- bind_rows(diff_tb_2000_2019)


## Global changes
diff_tb_2019_global <- dif_impact_from_samples(new_fvps = fvps_scen_df %>%  filter(iso3!="IND"), 
                                               impact_rates = da_rate_activity_cohort %>%  filter(iso3!="IND"), 
                                               period = 2019, 
                                               get_boots = TRUE, 
                                               global = TRUE)


ch_impact_samples <- list(diff_tb_2019 = diff_tb_2019,
                          diff_tb_2000_2019 = diff_tb_2000_2019, 
                          diff_tb_2019_global = diff_tb_2019_global)

saveRDS(ch_impact_samples, "ch_impact_samples_2019.RDS")

# 12. GET MAPS ----

AF_adm0 <- readRDS("./dat/AF_shp.RDS")

## Malawi district division https://commons.wikimedia.org/wiki/File:Malawi_district_map_2020.svg
## https://data.humdata.org/dataset/malawi-administrative-level-0-3-boundaries
mwi_map <- readRDS("./dat/mwi_adm1.RDS")

mwi_map <- 
  mwi_map %>%
  mutate(GID_0 = "MWI") %>% 
  rename(NAME_0 = ADM0_EN,
         GID_1 = ADM1_PCODE,
         NAME_1= ADM1_EN) %>% 
  select(GID_0, GID_1, NAME_0, NAME_1, geometry)


adm1_map2 <- 
  adm1_map %>% 
  filter(GID_0%in% c("IND", "MWI") == FALSE) %>% 
  select(GID_0, GID_1, NAME_0, NAME_1, geometry) %>% 
  bind_rows(mwi_map) %>% 
  mutate(NAME_1 = stringi::stri_trans_general(NAME_1, "latin-ascii")) %>% 
  mutate(NAME_1 = tolower(NAME_1)) %>% 
  mutate(NAME_1 = gsub("-", " ", NAME_1)) %>% 
  mutate(NAME_1 = gsub("&", "and", NAME_1))

adm1_map2 <- 
  adm1_map2 %>% 
  group_by(GID_0) %>% 
  mutate(ADM1_NAME = paste(GID_0, 1:length(NAME_1), sep = "_"))


AF_maps <- list(adm0 = AF_adm0, adm1 = adm1_map2)
saveRDS(AF_maps, "AFR_maps.RDS")


## Print markd
rmarkdown::render("report.Rmd", "html_document", output_file = "report.html")