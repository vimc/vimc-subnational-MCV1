### utility functions
save_csv <- function(dat, file_name) {
  write.csv(dat, file_name, row.names = FALSE)
}


not_is_finite <- function(x) {
  is.numeric(x) & !is.finite(x)
}

set_as_na <- function(x) {
  x <-NA
}

### output x with n significant numbers (ceiling)
signif.ceiling <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- ceiling(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}

### output x with n significant numbers (floor)
signif.floor <- function(x, n){
  pow <- floor( log10( abs(x) ) ) + 1 - n
  y <- floor(x / 10 ^ pow) * 10^pow
  # handle the x = 0 case
  y[x==0] <- 0
  y
}

mean_ <- function(x){
  mean(x, na.rm = TRUE)
}

quantile_ <- function(x){
  quantile(x, na.rm = TRUE)
}



## Get country metadata
add_country_metadata <- function(dat){
  
  country_metadata <- 
    tbl(con, "country_metadata") %>% 
    left_join(tbl(con, "country"), by = c("country" = "id")) %>%
    filter(grepl("201910", touchstone)) %>% 
    select(country_name = name, country, gavi73, who_region) %>%
    collect()
  
  dat <- dat %>% left_join(country_metadata, by = "country")
  
  return(dat)
}


### Functions to query data from annex

## Get mean and quantile ranges
get_stats <- function(d, group, get_all_quantiles = FALSE){
  
  if(get_all_quantiles == TRUE){
    
    b1 <- 
      d %>%
      group_by_at(group) %>%
      summarise_at(vars(matches("deaths|dalys|cases")), 
                   list(mid = ~ mean(., na.rm = TRUE), 
                        med =  ~ quantile(., 0.50, na.rm = TRUE),
                        `95_lo` = ~ quantile(., 0.025, na.rm = TRUE), 
                        `95_hi` = ~ quantile(., 0.975, na.rm = TRUE), 
                        `80_lo` = ~ quantile(., 0.1, na.rm = TRUE), 
                        `80_hi` = ~ quantile(., 0.9, na.rm = TRUE), 
                        `50_lo` = ~ quantile(., 0.25, na.rm = TRUE), 
                        `50_hi` = ~ quantile(., 0.75, na.rm = TRUE))) 
    
  }else{
    
    b1 <- 
      d %>%
      group_by_at(group) %>%
      summarise_at(vars(matches("deaths|dalys|cases")), 
                   list(mean = ~ mean(., na.rm = TRUE), 
                        q1 = ~ quantile(., 0.025, na.rm = TRUE), 
                        q3 = ~ quantile(., 0.975, na.rm = TRUE))) 
    
  }
  
  return(b1)
}


## global function to get estimates with different aggregations set in the group argument
aggregate_stochastics_global <- function(annex, 
                                         table_name = "cross_all_2019", 
                                         period = c(2000:2019), 
                                         disease_name = NULL, 
                                         country_name = NULL, 
                                         country_group = NULL, 
                                         group = c("disease", "modelling_group", "who_region"), 
                                         get_all_quantiles = FALSE, 
                                         is_test = FALSE, 
                                         calculate_NNV = FALSE, 
                                         touchstone = NULL){
  
  
  
  cohort = ifelse(grepl("cohort", table_name), TRUE, FALSE)
  
  under5 = ifelse(grepl("under5", table_name), TRUE, FALSE)
  
  touchstone_name = ifelse(grepl("2019", table_name), "201910gavi", 
                           ifelse(grepl("2021", table_name), "202110gavi", "201710gavi"))
  
  boots_table = ifelse(grepl("2019", table_name), "bootstrap_2019", 
                       ifelse(grepl("2021", table_name), "bootstrap_2021",
                              ifelse(grepl("ia2030", table_name), "bootstrap_ia2030", "bootstrap")))
  
  who_reg <- tbl(con, "country_metadata") %>% 
    left_join(tbl(con, "country"), by = c("country" = "id")) %>%
    filter(grepl("201910", touchstone)) %>% 
    select(country, who_region) %>%
    collect()
  
  
  if(is.null(disease_name)){
    
    disease_name <- DBI::dbGetQuery(con, "SELECT DISTINCT id FROM disease")[["id"]]
  }
  
  
  if(is.null(country_name) & is.null(country_group)){
    
    country_name <-  DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata")[["country"]]
    
  }else if(!is.null(country_name) & !is.null(country_group)){
    
    stop("Provide either country_name OR country_group")
  }
  
  
  ## Main country filters ----
  ## The idea with these filters is to be able to get estimates directly for the different country groups we used in the VIMC
  if(is.null(country_name) && country_group == "gavi73"){
    
    ## Get GAVI countries
    country_name <- DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE gavi73")[["country"]]
    
  }else if(is.null(country_name) && country_group == "vimc98"){
    
    ## get VIMC98 countries
    country_name <- DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE touchstone = '201710gavi-1'")[["country"]]
    
  }else if(is.null(country_name) && (country_group == "vimc112" | "who_region" %in% group) ){
    
    ## get VIMC112 countries
    country_name <- DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE touchstone = '201910gavi-1'")[["country"]]
    
  }else if(is.null(country_name) && country_group == "dove94"){
    
    ## get DOVE94 countries
    country_name <- DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE dove94")[["country"]]
    
  }else if(is.null(country_name) &&  !(country_group%in%c("gavi73", "vimc98", "vimc112", "dove94"))){
    
    stop("Only gavi73, vimc98, vimc112 and dove94 are allowed")
  }
  
  if(is_test){
    
    tmp <- 
      dplyr::tbl(annex, table_name) %>% 
      filter(run_id == 1) %>% 
      filter(disease%in%disease_name) %>% 
      filter(country%in%c("PAK", "IND", "NGA", "ETH")) %>% 
      filter(year%in%period)
    
    
  }else{
  
    tmp <- 
      dplyr::tbl(annex, table_name) %>% 
      filter(disease%in%disease_name) %>% 
      filter(country%in%country_name) %>% 
      filter(year%in%period)
    
  }
  
  ## add WHO region
  if("who_region" %in% group){
    tmp <- tmp %>% 
      left_join(who_reg, by = "country", copy = TRUE) 
  }
  
  
  if(grepl("2021", table_name)){
    
    var_to_get1 <- c("deaths_averted", "dalys_averted", "cases_averted", "deaths_averted_rate")
    var_to_get2 <- c("deaths_impact", "deaths_default", "deaths_novac", 
                     "dalys_impact", "dalys_default", "dalys_novac",
                     "cases_impact", "cases_default", "cases_novac")
  }else{
    
    var_to_get1 <- c("deaths_averted", "dalys_averted", "deaths_averted_rate")
    var_to_get2 <- c("deaths_impact", "deaths_default", "deaths_novac", 
                     "dalys_impact", "dalys_default", "dalys_novac")
    
  }
    
    
    ### include cases into the tables
    if(!"stochastic_file_id" %in% colnames(tmp) & !"deaths_impact" %in% colnames(tmp)){
      
      ## Get aggregated data per disease and modeling group
      dat <-  
        ## Join boostrap with the table
        tmp %>% 
        left_join(
          tbl(annex, boots_table) %>% 
            filter(disease%in%disease_name) %>% 
            select(-stochastic_file_id) %>% 
            distinct()
        ) %>% 
        filter(!is.na(boots_id)) %>% 
        group_by_at(c("boots_id", group)) %>% 
        summarise_at(var_to_get1, sum, na.rm = TRUE) %>% 
        collect()
      
      names(dat) <- gsub("_averted", "_impact", names(dat))
      
    }else{
      
      ## To get the correct modelling group names for Rubella cohort data
      if(disease_name == "Rubella" & cohort){
        
        cohort = FALSE
        
      }
      
      ## Sum across all years from the period and collect samples
      dat <- 
        tmp %>% 
        group_by_at(c(group[group!="modelling_group"], "run_id", "stochastic_file_id")) %>% 
        summarise_at(var_to_get2, sum, na.rm = TRUE) %>%
        
        mutate(boots_id = (run_id * 1000 + stochastic_file_id)) %>% 
        ungroup() %>% 
        ## Include modelling group name
        left_join(
          
          tbl(annex, "stochastic_file") %>% 
            filter(disease%in%disease_name & grepl(touchstone_name, touchstone) & is_cohort == cohort & is_under5 == under5) %>% 
            select(modelling_group, id), by = c("stochastic_file_id" = "id")
        ) %>% 
        select(-stochastic_file_id, -run_id) %>% 
        collect()
      
      
    }
    
  if(!grepl("intervention", table_name)){
    
    dat <- 
      dat %>% 
      mutate(proportion_deaths_averted = deaths_impact/deaths_novac)
  }
  
  if(calculate_NNV & grepl("cohort", table_name)){
    
    fvps <- vimpact::extract_vaccination_history(con, 
                                                 touchstone_cov = touchstone,  
                                                 vaccine_to_ignore = c("DTP3", "HepB_BD_home"),
                                                 gavi_support_levels = c("with"))
    
    fvps_df <- 
      fvps %>% 
      filter(disease%in%disease_name) %>% 
      filter(country%in%country_name) %>% 
      add_country_metadata() %>%
      mutate(cohort = year - age) %>% 
      select(-year) %>% 
      mutate(year = cohort) %>% 
      filter(year%in%period) %>% 
      ### sum fvps over the groups
      group_by_at(c(group[group!="modelling_group"])) %>%
      summarise(fvps = sum(fvps_adjusted, na.rm = TRUE))
    
    if(grepl("2021", table_name)){
      
      dat <- 
        dat %>% 
        left_join(fvps_df) %>% 
        filter(!is.na(fvps)) %>% 
        mutate(nnv_deaths = fvps/deaths_impact, 
               nnv_dalys = fvps/dalys_impact, 
               nnv_cases = fvps/cases_impact, 
               IR_deaths = deaths_impact/fvps,
               IR_dalys = dalys_impact/fvps,
               IR_cases = cases_impact/fvps)
      
    }else{
      
      dat <- 
        dat %>% 
        left_join(fvps_df) %>% 
        filter(!is.na(fvps)) %>% 
        mutate(nnv_deaths = fvps/deaths_impact, 
               nnv_dalys = fvps/dalys_impact)
      
    }
}
  
  
  if(!("year"%in%group)){ ### When year isn't within the aggregation group then create a time label. 
    
    ## Create label for period
    p_label <- ifelse(length(period) == 1, 
                      as.character(period), 
                      paste(min(period), max(period), sep = "-"))
    
    d <- 
      dat %>%
      mutate(time = p_label) %>% 
      get_stats(group = c(group, "time"), get_all_quantiles = get_all_quantiles)
    
    
  }else{
    
    d <- 
      dat %>%
      get_stats(group = group, get_all_quantiles = get_all_quantiles)
    
  }
  
  if(!is.null(country_group)){
    
    d$country_group <- country_group
  }
  
  
  return(d)
}



###############################################################################
### extract estimates with double counting ####################################
###############################################################################
## get aggregated data
aggregate_per_disease <- function(annex, table, period = c(2000:2019)){
  
  ## Get GAVI countries
  GAVI73 <- DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE gavi73")[["country"]]
  
  ## get VIMC98 countries
  VIMC98 <- DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE touchstone = '201710gavi-1'")[["country"]]
  
  boots_table = ifelse(grepl("2019", table), "bootstrap_2019", 
                       ifelse(grepl("2021", table), "bootstrap_2021",
                              ifelse(grepl("ia2030", table), "bootstrap_ia2030", "bootstrap")))
  
  tmp <- 
    dplyr::tbl(annex, table) %>% 
    mutate(gavi73 = country %in% GAVI73) %>%
    mutate(vimc98 = country %in% VIMC98) %>%
    filter(year%in%period)
  
  if(!"stochastic_file_id" %in% colnames(tmp) & !"deaths_impact" %in% colnames(tmp)){
    
    ## Join boostrap with the table
    dat <- 
      tbl(annex, boots_table) %>% 
      select(-stochastic_file_id) %>% 
      distinct() %>% 
      left_join(tmp) %>% 
      group_by(boots_id, disease, country) %>% 
      summarise(deaths_impact = sum(deaths_averted, na.rm = TRUE), 
                dalys_impact = sum(dalys_averted, na.rm = TRUE),
                cases_impact = sum(cases_averted, na.rm = TRUE)) %>% 
      collect() %>% 
      mutate(deaths_default = NA, 
             deaths_novac = NA,
             dalys_default = NA,
             dalys_novac = NA,
             cases_default = NA,
             cases_novac = NA)
    
  }else{
    
    dat <- tmp %>% group_by(disease, country, run_id, stochastic_file_id) %>% 
      summarise(deaths_impact = sum(deaths_impact, na.rm = TRUE), 
                deaths_default = sum(deaths_default, na.rm = TRUE),
                deaths_novac = sum(deaths_novac, na.rm = TRUE), 
                dalys_impact = sum(dalys_impact, na.rm = TRUE),
                dalys_default = sum(dalys_default, na.rm = TRUE),
                dalys_novac = sum(dalys_novac, na.rm = TRUE),
                cases_impact = sum(cases_impact, na.rm = TRUE),
                cases_default = sum(cases_default, na.rm = TRUE),
                cases_novac = sum(cases_novac, na.rm = TRUE)) %>% 
      mutate(boots_id = (run_id * 1000 + stochastic_file_id)) %>% 
      ungroup() %>% 
      select(disease, boots_id, country, contains("death"), contains("dalys"), contains("cases")) %>% 
      collect()
    
  }
  
  ## Include GAVI column and prop deaths and dalys averted
  ## Create label for period
  p_label <- ifelse(length(period) == 1, 
                    as.character(period), 
                    paste(min(period), max(period), sep = "-"))
  
  d <- 
    dat %>%
    mutate(gavi73 = country %in% GAVI73) %>%
    mutate(vimc98 = country %in% VIMC98) %>%
    mutate_if(is.numeric, ~ set_as_zero(.)) %>%
    mutate(time = p_label)
  
  return(d)
}


get_gavi_73 <- function() {
  DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE gavi73")[["country"]]
}

get_vimc_98 <- function() {
  DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE touchstone = '201710gavi-1'")[["country"]]
}

format_output <- function(data) {
  
  if ("start_time" %in% colnames(data)) {
    p_label <- ifelse(data$start_time == data$end_time, 
                      as.character(data$start_time), 
                      paste(data$start_time, data$end_time, sep = "-"))
    data <- data %>%
      mutate(time = p_label) %>%
      select(-c(end_time, start_time))
  }
  
  if ("country" %in% colnames(data)) {
    GAVI73 <- get_gavi_73()
    VIMC98 <- get_vimc_98()
    data <- data %>%
      mutate(gavi73 = country %in% GAVI73) %>%
      mutate(vimc98 = country %in% VIMC98)
  }
  data
}


## get global aggregated 

## Global aggregation without correcting for double counting. 
global_withDB <- function(annex, table, periods = list(c(2000:2019)), countries_filter = NULL, 
                          get_all_quantiles = FALSE){
  
  
  if(countries_filter == "gavi73"){
    
    countries_included = DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE gavi73")[["country"]]
    countries_group = "gavi73"
    
  }else if(countries_filter == "vimc98"){
    
    countries_included = DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE touchstone = '201710gavi-1'")[["country"]]
    countries_group = "vimc98"
    
  }else if(countries_filter == "dove94"){
    
    countries_included = DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE dove94")[["country"]]
    countries_group = "dove94"
    
  }else{
    print("Using the 112 VIMC countries")
    
    countries_included = DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE touchstone = '201910gavi-1'")[["country"]]
    countries_group = "vimc112"
  }
  
  boots_table = ifelse(grepl("2019", table), "bootstrap_2019", 
                       ifelse(grepl("2021", table), "bootstrap_2021",
                              ifelse(grepl("ia2030", table), "bootstrap_ia2030", "bootstrap")))
  
  
  print(paste("extracting table", table))
  
  d <- lapply(periods, function(x){
    ## Get aggregated data per disease and modeling group
    dat <-  
      ## Join boostrap with the table
      tbl(annex, boots_table) %>% 
      select(-stochastic_file_id) %>% 
      distinct() %>% 
      left_join(tbl(annex, table)) %>% 
      filter(year %in% x) %>% 
      filter(country%in%countries_included) %>% 
      group_by(boots_id, disease, modelling_group) %>% 
      summarise(deaths_impact = sum(deaths_averted, na.rm = TRUE), 
                dalys_impact = sum(dalys_averted, na.rm = TRUE), 
                cases_impact = sum(dalys_averted, na.rm = TRUE)) %>% 
      collect()
    
    ## Create label for period
    p_label <- ifelse(length(x) == 1, 
                      as.character(x), 
                      paste(min(x), max(x), sep = "-"))
    
    ## Aggregate across all countries
    dat <- 
      dat %>%  
      ### Calculate mean values per disease
      group_by(boots_id, disease) %>% 
      summarise(deaths_impact = mean(deaths_impact, na.rm = TRUE), 
                dalys_impact = mean(dalys_impact, na.rm = TRUE),
                cases_impact = mean(cases_impact, na.rm = TRUE)) %>%
      ungroup() %>% 
      ## sum by boots_id 
      group_by(boots_id) %>% 
      summarise(deaths_impact = sum(deaths_impact, na.rm = TRUE), 
                dalys_impact = sum(dalys_impact, na.rm = TRUE)) %>% 
      mutate(time = p_label)
    
    
    dat
  })
  
  d <- bind_rows(d)
  
  d <- 
    get_stats(d, group = c("time"), get_all_quantiles = get_all_quantiles) %>% 
    mutate(country_group = countries_group)
  
  return(d)
}

## get aggregated data per disease 
global_per_disease <- function(annex, table, periods = list(c(2000:2019)), get_all_quantiles = FALSE, is_test = FALSE){
  print(paste("global per disease extracting table", table))
  
  d <- lapply(periods, function(x){
    ## Get aggregated data per country
    aggregate_stochastics_global(annex, table, 
                                 period = x, 
                                 group = c("disease"), 
                                 get_all_quantiles = get_all_quantiles, 
                                 is_test = is_test)
  })
  
  d <- bind_rows(d)
  d <- format_output(d)
  
  return(d)
}

## get aggregated data per country 
iso_per_disease <- function(annex, table, 
                            periods = list(c(2000:2019)), 
                            per_year = NULL,
                            get_all_quantiles = NULL, is_test = FALSE){
  
  print(paste("iso per disease extracting table", table))
  
  d <- lapply(periods, function(x){
    ## Get aggregated data per country
    aggregate_stochastics_global(annex, table, 
                                 period = x, 
                                 group = c("country", "disease"), 
                                 get_all_quantiles = get_all_quantiles,
                                 is_test = is_test)
  })
  
  d <- bind_rows(d)
  
  if(is.null(per_year)){
    
    d <- format_output(d)
    
  }else{
    
    d_tmp <- aggregate_stochastics_global(annex, table, 
                                          period = per_year, 
                                          group = c("country", "disease", "year"), 
                                          get_all_quantiles = get_all_quantiles,
                                          is_test = is_test)
    
    d <- 
      d_tmp %>% 
      mutate(time = as.character(year)) %>% 
      select(-year) %>% 
      bind_rows(d)
    
    d <- format_output(d)
    
  }
  
  
  return(d)
  
}


## get aggregated data per country 
who_region_per_disease <- function(annex, table, periods = list(c(2000:2019)),
                                   get_all_quantiles = FALSE, is_test = FALSE, per_year = NULL, touchstone_name = "202110gavi-3"){
  print(paste("who_region per disease extracting table", table))
  
  d <- lapply(periods, function(x){
    ## Get aggregated data per country
    aggregate_stochastics_global(annex, table, 
                                 period = x, 
                                 group = c("who_region", "disease"), 
                                 get_all_quantiles = get_all_quantiles, 
                                 is_test = is_test, 
                                 calculate_NNV = TRUE, 
                                 touchstone = touchstone_name)
  })
  
  d <- bind_rows(d)
  
  if(is.null(per_year)){
    
    d <- format_output(d)
    
  }else{
    
    d_tmp <- aggregate_stochastics_global(annex, table, 
                                          period = per_year, 
                                          group = c("who_region", "disease", "year"), 
                                          get_all_quantiles = get_all_quantiles,
                                          is_test = is_test,
                                          calculate_NNV = TRUE, 
                                          touchstone = touchstone_name)
    
    d <- 
      d_tmp %>% 
      mutate(time = as.character(year)) %>% 
      select(-year) %>% 
      bind_rows(d)
    
    d <- format_output(d)
    
  }
  
  return(d)
  
}

## get aggregated data gavi
gavi_per_disease <- function(annex, table, periods = list(c(2000:2019)), get_all_quantiles = FALSE, is_test = FALSE){
  print(paste("gavi per disease extracting table", table))
  
  d <- lapply(periods, function(x){
    ## Get aggregated data per country
    aggregate_stochastics_global(annex, table, 
                                 period = x, 
                                 group = c("disease"), 
                                 get_all_quantiles = get_all_quantiles, 
                                 country_group = "gavi73", 
                                 is_test = is_test)
  })
  
  d <- bind_rows(d)
  d <- format_output(d)
  
  return(d)
}

## get aggregated data vimc98
vimc98_per_disease <- function(annex, table, periods = list(c(2000:2019)), get_all_quantiles = FALSE, is_test = FALSE){
  print(paste("vimc98 per disease extracting table", table))
  
  d <- lapply(periods, function(x){
    ## Get aggregated data per country
    aggregate_stochastics_global(annex, table, 
                                 period = x, 
                                 group = c("disease"), 
                                 get_all_quantiles = get_all_quantiles, 
                                 country_group = "vimc98", 
                                 is_test = is_test)
  })
  
  d <- bind_rows(d)
  d <- format_output(d)
  
  return(d)
}

###############################################################################
### extract estimates without double counting #################################
###############################################################################
## Get aggregate data
aggregate_annex <- function(annex, 
                            table = "impact_cluster_cross_all_2021", 
                            view = "global", 
                            period = 2000:2030){
  if(view == "iso"){
    group = c("boots_id", "country")
  }else if(view == "who_region"){
    group = c("boots_id", "who_region")
  } else{
    group = "boots_id"
  }
  
  ## Sum across all years per boot_id and country
  if(view == "gavi"){
    
    d <- 
      tbl(annex, table) %>% 
      filter(year%in%period & gavi73 == TRUE) %>%
      mutate(deaths_default = deaths_impact - deaths_novac) %>% 
      mutate(dalys_default = dalys_impact - dalys_novac) %>%
      group_by_at(group) %>% 
      summarise(deaths_impact =  sum(deaths_impact, na.rm = TRUE),
                deaths_default = sum(deaths_default, na.rm = TRUE),
                deaths_novac =   sum(deaths_novac, na.rm = TRUE), 
                dalys_impact =   sum(dalys_impact, na.rm = TRUE),
                dalys_default =  sum(dalys_default, na.rm = TRUE),
                dalys_novac =    sum(dalys_novac, na.rm = TRUE)) %>% 
      collect()
    
  }else if(view == "vimc98"){
    
    ## get VIMC98 countries
    VIMC98 <- DBI::dbGetQuery(con, "SELECT DISTINCT country FROM country_metadata WHERE touchstone = '201710gavi-1'")[["country"]]
    
    d <- 
      tbl(annex, table) %>% 
      filter(year%in%period & country %in% VIMC98) %>%
      mutate(deaths_default = deaths_impact - deaths_novac) %>% 
      mutate(dalys_default = dalys_impact - dalys_novac) %>%
      group_by_at(group) %>% 
      summarise(deaths_impact = sum(deaths_impact, na.rm = TRUE),
                deaths_default = sum(deaths_default, na.rm = TRUE),
                deaths_novac = sum(deaths_novac, na.rm = TRUE), 
                dalys_impact = sum(dalys_impact, na.rm = TRUE),
                dalys_default = sum(dalys_default, na.rm = TRUE),
                dalys_novac = sum(dalys_novac, na.rm = TRUE)) %>% 
      collect()
  } else if(view == "who_region"){
    d <- tbl(annex, table) %>%
      
      left_join(
        tbl(con, "country_metadata") %>% 
          left_join(tbl(con, "country"), by = c("country" = "id")) %>%
          filter(grepl("201910", touchstone)) %>% 
          select(country, who_region) %>%
          collect(),
        by = "country", copy = TRUE) %>%
      
      filter(year%in%period) %>%
      mutate(deaths_default = deaths_novac - deaths_impact) %>% 
      mutate(dalys_default = dalys_novac - dalys_impact) %>%
      group_by_at(group) %>% 
      summarise(deaths_impact = sum(deaths_impact, na.rm = TRUE),
                deaths_default = sum(deaths_default, na.rm = TRUE),
                deaths_novac = sum(deaths_novac, na.rm = TRUE), 
                dalys_impact = sum(dalys_impact, na.rm = TRUE),
                dalys_default = sum(dalys_default, na.rm = TRUE),
                dalys_novac = sum(dalys_novac, na.rm = TRUE)) %>% 
      collect()
  } else{
    d <- 
      tbl(annex, table) %>% 
      filter(year%in%period) %>%
      mutate(deaths_default = deaths_novac - deaths_impact) %>% 
      mutate(dalys_default = dalys_novac - dalys_impact) %>%
      group_by_at(group) %>% 
      summarise(deaths_impact = sum(deaths_impact, na.rm = TRUE),
                deaths_default = sum(deaths_default, na.rm = TRUE),
                deaths_novac = sum(deaths_novac, na.rm = TRUE), 
                dalys_impact = sum(dalys_impact, na.rm = TRUE),
                dalys_default = sum(dalys_default, na.rm = TRUE),
                dalys_novac = sum(dalys_novac, na.rm = TRUE)) %>% 
      collect()
    
  }
  
  
  ## Include period label
  p_label <- ifelse(length(period) == 1, 
                    as.character(period), 
                    paste(min(period), max(period), sep = "-"))
  
  ## Calculate prop deaths and dalys averted
  d <- 
    d %>% 
    mutate(
      proportion_deaths_averted = deaths_impact / deaths_novac,
      proportion_dalys_averted = dalys_impact / dalys_novac,
      time = p_label)
  
  return(d)
}

### per country all diseases
get_cluster_agg <- function(annex, table, periods = list(c(2000:2019)), view = "global", 
                            get_all_quantiles = FALSE){
  print(paste("extracting table", table))
  
  d <- 
    foreach(i = seq_along(periods), .combine = rbind)%do%{
      aggregate_annex(annex, table, period = periods[[i]], view = view)
    }
  
  ## Aggregate depending on the view
  if(view == "iso"){
    
    d <- get_stats(d, group = c("time", "country"), get_all_quantiles = get_all_quantiles)
    
  }else if(view%in%c("gavi", "global", "vimc98")){
    
    d <- get_stats(d, group = "time", get_all_quantiles = get_all_quantiles)
    
  }else if(view == "who_region"){
    d <- get_stats(d, group = c("time", "who_region"), get_all_quantiles = get_all_quantiles)
  } else {
    stop("Invalid view option")
  }
  
  return(d)
  
}

set_as_zero <- function(x){
  i <- !is.na(x) & is.finite(x)
  i <- !i
  if(any(i)){
    x[i] <- 0
  }
  return(x)
}
