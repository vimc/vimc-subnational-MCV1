firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

match_adm1_names<- function(vax_df, adm1_names_ref){
  
  ## Fix some admin1 mismatches
  adm1_names_ref$vimc <- as.character(adm1_names_ref$vimc)
  ind1 <- which(vax_df$ADM1_NAME%in%adm1_names_ref$who)
  vax_df$ADM1_NAME[ind1] <- adm1_names_ref$vimc[match(vax_df$ADM1_NAME, adm1_names_ref$who)][ind1]
  
  return(vax_df)
}


stand_name <- function(df){
  
  df_f <- 
    df %>% 
    mutate(ADM1_NAME = stringi::stri_trans_general(ADM1_NAME, "latin-ascii")) %>% 
    mutate(ADM1_NAME = tolower(ADM1_NAME)) %>% 
    ## remove unicode-replacement-character
    mutate(ADM1_NAME = gsub("\uFFFD", "", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = gsub(" region", "", ADM1_NAME)) %>%
    mutate(ADM1_NAME = gsub(" - ", " ", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = gsub("   ", " ", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = gsub("-", " ", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = gsub("&", "and", ADM1_NAME)) %>% 
    ## fix Angola district names
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "bi", "bie", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "hula", "huila", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "uge", "uige", ADM1_NAME)) %>% 
    
    ## Fix UGANDA classification
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "kyotera", "kyotara", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "kalaki", "kaberamaido", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "karenga", "kaabong", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "kazo", "kiruhura", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "kitagwenda", "kamwenge", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "madi okollo", "arua", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "obongi", "moyo", ADM1_NAME)) %>% 
    mutate(ADM1_NAME = ifelse(ADM1_NAME == "rwampara", "mbarara", ADM1_NAME))
    
    
  return(df_f)
}

get_HepB3_cov <- function(raw_cov, adm1_names_ref, countries_to_include, pop_ref){
  
  HepB3_cov <-
    raw_cov %>% 
    filter(Vaccine == "HepB3") %>% 
    mutate(Admin1 = ifelse(iso == "KEN", Admin2, Admin1)) %>% 
    mutate(Admin1 = ifelse(iso == "STP", Admin2, Admin1)) %>% 
    filter(iso%in%countries_to_include) %>% 
    mutate(Numerator = ifelse(Numerator == -4444, NA, Numerator)) %>% 
    mutate(Numerator = ifelse(Numerator == -2222, NA, Numerator))
  
  ## Getting estimates for admin1
  HepB3_cov2 <- 
    HepB3_cov %>% 
    group_by(iso, CountryName, Year, Admin1, Vaccine) %>% 
    summarise(Numerator = sum(Numerator, na.rm = TRUE)) %>% 
    mutate(Admin1 = tolower(Admin1)) %>% 
    mutate(Admin1 = gsub("-", " ", Admin1)) %>% 
    mutate(Admin1 = gsub("&", "and", Admin1))
  
  HepB3_cov2$Admin1 <- stringi::stri_trans_general(HepB3_cov2$Admin1, "latin-ascii")
  
  ## Fix some admin1 mismatches
  ind1 <- which(HepB3_cov2$Admin1%in%adm1_names_ref$who)
  HepB3_cov2$Admin1[ind1] <- adm1_names_ref$vimc[match(HepB3_cov2$Admin1, adm1_names_ref$who)][ind1]
  
  
  pop_ref_2018 <- 
    pop_ref %>% 
    filter(year == 2018)
  
  HepB3_cov_pop <- 
    HepB3_cov2 %>% 
    select(iso3 = iso, ADM1_NAME = Admin1, year = Year, Numerator) %>% 
    left_join(pop_ref_2018) %>% 
    filter(!is.na(total_pop)) %>% 
    mutate(mean_cov_ref = round(Numerator/total_pop,1)) %>% 
    mutate(mean_cov = ifelse(mean_cov_ref > 1, 1, mean_cov_ref)) %>% 
    mutate(lo_cov = NA, 
           hi_cov = NA)
  
  return(HepB3_cov_pop)
}


prepare_impact_ratio <- function(d, fvps, stratification, period, get_all_dat = FALSE){
  
  ## in the case of MCV1 data, cohort and year are the same
  d$year = d$cohort
  d$country = d$iso3
  
  ## transformation for calculating impact ratio
  ## see the readme for details
  if(stratification %in% c("un-stratified", "activity")){
    tot_impact_rout <- d %>%
      filter(year %in% period & activity_type == "routine") %>%
      group_by(boots_id, disease, vaccine, activity_type, country) %>%
      summarise(tot_impact = round(sum(deaths_averted_MCV1, na.rm = TRUE)), 
                tot_deaths_novac = round(sum(deaths_novac_MCV1, na.rm = TRUE))) %>%
      as.data.frame() 
    
    tot_impact_sia <- d %>%
      filter(activity_type == "campaign") %>%
      group_by(boots_id, disease, vaccine, activity_type, country) %>%
      summarise(tot_impact = round(sum(deaths_averted_MCV1, na.rm = TRUE)),
                tot_deaths_novac = round(sum(deaths_novac_MCV1, na.rm = TRUE))) %>%
      as.data.frame() 
    
    tot_impact <- rbind(tot_impact_rout, tot_impact_sia)
    
    if(stratification == "un-stratified"){
      tot_impact <- tot_impact %>%
        group_by(disease,  boots_id, country) %>%
        summarise(tot_impact = round(sum(tot_impact, na.rm = TRUE)), 
                  tot_deaths_novac = round(sum(tot_deaths_novac, na.rm = TRUE))) %>%
        as.data.frame()
      
      tot_fvps <- fvps %>%
        filter(year %in% period) %>%
        group_by(disease, country) %>%
        summarise(tot_fvps = round(sum(fvps, na.rm = TRUE)))
    } else {
      tot_fvps <- fvps %>%
        filter(year %in% period) %>%
        group_by(disease, country, vaccine, activity_type) %>%
        summarise(tot_fvps = round(sum(fvps, na.rm = TRUE)))
    }
  } else if(stratification == "cohort"){
    tot_impact <- d %>%
      filter(cohort <= max(period)) %>%
      group_by(disease, boots_id, country, cohort) %>%
      summarise(tot_impact = round(sum(deaths_averted_MCV1, na.rm = TRUE)), 
                tot_deaths_novac = round(sum(deaths_novac_MCV1, na.rm = TRUE))) %>%
      as.data.frame()
    
    tot_fvps <- fvps %>%
      filter(cohort <= max(period)) %>%
      group_by(disease, country, cohort) %>%
      summarise(tot_fvps = round(sum(fvps, na.rm = TRUE)))
    
  }else if(stratification == "both"){
    tot_impact <- d %>%
      filter(cohort <= max(period)) %>%
      group_by(disease, vaccine, activity_type, boots_id, country, cohort) %>%
      summarise(tot_impact = round(sum(deaths_averted_MCV1, na.rm = TRUE)), 
                tot_deaths_novac = round(sum(deaths_novac_MCV1, na.rm = TRUE))) %>%
      as.data.frame()
    
    tot_fvps <- fvps %>%
      filter(cohort <= max(period)) %>%
      group_by(disease, country, vaccine, activity_type, cohort) %>%
      summarise(tot_fvps = round(sum(fvps, na.rm = TRUE)))
  } else {
    stop(print("unspecified stratification level"))
  }
  
  impact_ratio <- tot_impact %>%
    left_join(tot_fvps) %>%
    mutate(impact_ratio = tot_impact / tot_fvps) %>%
    select(-tot_impact, -tot_fvps)
  
  dat <- impact_ratio %>%
    left_join(fvps) %>% 
    filter(!is.na(impact_ratio)) %>%
    mutate(impact = impact_ratio * fvps)  %>% 
    select(disease, boots_id, country, impact_ratio, vaccine, activity_type, year, age, fvps_adjusted, fvps, cohort, impact, tot_deaths_novac)
  
  if(get_all_dat){
    
    return(list(tot_impact = tot_impact, tot_fvps = tot_fvps, impact_ratio = impact_ratio, final_impact = dat))
    
  }else{
    return(dat)
  }
  
}

normalise_fvps_MCV1 <- function(cov_dat, fvps_dat, pop_dat, has_samples = FALSE){
  
  fvps_tmp <- 
    fvps_dat %>% 
    filter(vaccine == "MCV1") %>% 
    select(iso3 = country, year, nat_cov = coverage_adjusted, fvps_adjusted) %>% 
    filter(iso3%in%countries_to_include)
  
  if(has_samples){
    
    new_cov <- 
      cov_dat %>% 
      left_join(pop_dat) %>% 
      left_join(fvps_tmp) %>% 
      mutate(tmp_fvps_m = value*total_pop) %>% 
      group_by(iso3, year, sample) %>% 
      mutate(t_fvps_m = sum(tmp_fvps_m, na.rm=TRUE)) %>% 
      ungroup() %>% 
      mutate(mean_fvps = tmp_fvps_m/t_fvps_m*fvps_adjusted) %>% 
      select(-contains("t_fvps"), -nat_cov) %>% 
      mutate(adj_cov = mean_fvps/total_pop) %>% 
      select(sample, ADM0_NAME, iso3, ADM1_NAME, ADM1_CODE, ADM1_NAME_orig, year, adj_cov, adjusted_fvps = mean_fvps, unadjusted_fvps = tmp_fvps_m)
    
    
  }else{
  
    new_cov <- 
      cov_dat %>% 
      left_join(pop_dat) %>% 
      left_join(fvps_tmp) %>% 
      mutate(tmp_fvps_m = mean_cov*total_pop) %>% 
      group_by(iso3, year) %>% 
      mutate(t_fvps_m = sum(tmp_fvps_m, na.rm=TRUE)) %>% 
      ungroup() %>% 
      mutate(mean_fvps = tmp_fvps_m/t_fvps_m*fvps_adjusted) %>% 
      select(-contains("t_fvps"), -nat_cov) %>% 
      mutate(mean_cov = mean_fvps/total_pop) %>% 
      select(vaccine, ADM0_NAME, iso3, ADM1_NAME, ADM1_CODE, year, mean_cov, adjusted_fvps = mean_fvps, unadjusted_fvps = tmp_fvps_m)
    
  }
  
  return(new_cov)
}

get_fvps_under_H <- function(sub_cov, vimc_nat_cov, sub_pop, DTP3_cov, option_cov = "sub_cov"){
  
  options_l <- c("sub_cov", "nat_cov", "at_nat_cov", "at_DTP3_cov", "max_cov", "at_gvap_cov")
  
  if(option_cov%in%options_l == FALSE){
  
    stop("Provide one of the following options 'sub_cov', 'nat_cov', 'at_nat_cov', 'at_DTP3_cov', 'at_gvap_cov' OR 'max_cov'")
    
  }
  
  ## Get new fvps: fvps at subnational level using pop data and reference subnational coverage
  if(option_cov == "sub_cov"){
    
    ## option sub_cov:  subnational cov
    sub_fvps_new <- 
      sub_cov %>% 
      left_join(sub_pop) %>% 
      mutate(mean_fvps = mean_cov*total_pop) %>% 
      select(-vaccine, -ADM1_CODE)
    
  }else if(option_cov == "nat_cov"){
    
    ## option nat_cov:  national coverage
    sub_fvps_new <- 
      sub_cov %>% 
      left_join(sub_pop) %>% 
      left_join(vimc_nat_cov) %>% 
      mutate(mean_fvps = nat_cov*total_pop, 
             mean_cov = nat_cov) %>% 
      select(-vaccine, -ADM1_CODE, -nat_cov)
    
  }else if(option_cov == "at_nat_cov"){
    
    ## option at_nat_cov:  at least national coverage
    sub_fvps_new <- 
      sub_cov %>% 
      left_join(sub_pop) %>% 
      left_join(vimc_nat_cov) %>% 
      mutate(mean_cov = ifelse(mean_cov < nat_cov, nat_cov, mean_cov)) %>% 
      mutate(mean_fvps = mean_cov*total_pop) %>% 
      select(-vaccine, -ADM1_CODE, -nat_cov)
    
  }else if(option_cov == "at_gvap_cov"){
    
    ## option at_gvap_cov:  at least the gvap target (80%)
    sub_fvps_new <- 
      sub_cov %>% 
      left_join(sub_pop) %>% 
      left_join(vimc_nat_cov) %>% 
      mutate(mean_cov = ifelse(mean_cov < 0.8, 0.8, mean_cov)) %>% 
      mutate(mean_fvps = mean_cov*total_pop) %>% 
      select(-vaccine, -ADM1_CODE, -nat_cov)
    
  }else if(option_cov == "at_DTP3_cov"){
    
    ## option at_dtp3_cov:  DTP3 coverage
    # assuming at least DTP3 coverage
    DTP3_cov_tmp <- 
      DTP3_cov %>% 
      select(ADM0_NAME, iso3, ADM1_NAME, mean_cov_dtp3 = mean_cov, year)
    
    n_dis_DTP3 <- 
      DTP3_cov_tmp %>% 
      group_by(iso3) %>% 
      summarise(n = n_distinct(ADM1_NAME))
    
    n_dis <- 
      sub_pop %>% 
      group_by(iso3) %>% 
      summarise(n_pop = n_distinct(ADM1_NAME)) %>% 
      left_join(n_dis_DTP3)
      
    cou_to_drop <- n_dis$iso3[which(n_dis$n_pop!=n_dis$n & n_dis$iso3!="MWI")]
    
    sub_fvps_new <- 
      sub_cov %>% 
      filter(iso3%in%cou_to_drop == FALSE) %>% 
      left_join(sub_pop) %>% 
      left_join(DTP3_cov_tmp) %>% 
      mutate(mean_cov = ifelse(mean_cov < mean_cov_dtp3, mean_cov_dtp3, mean_cov)) %>% 
      mutate(mean_fvps = mean_cov*total_pop) %>% 
      select(-vaccine, -ADM1_CODE, -mean_cov_dtp3)
    
  }else if(option_cov == "max_cov"){
    
    ## option max_cov:  max coverage
    ## get fvps if we assume all districts have achieved the max coverage in a specific year
    get_max_cov_regions <- 
      sub_cov %>% 
      group_by(iso3, year) %>% 
      summarise(max_cov_mid = max(mean_cov, na.rm = TRUE))
    
    sub_fvps_new <- 
      sub_cov %>% 
      left_join(sub_pop) %>% 
      left_join(get_max_cov_regions) %>% 
      mutate(mean_fvps = max_cov_mid*total_pop, 
             mean_cov = max_cov_mid) %>% 
      select(-vaccine, -ADM1_CODE, -max_cov_mid)
    
  }
  
  sub_fvps_new$scenario = option_cov
  
  return(sub_fvps_new)
}

IU_impact_period <- function(new_fvps, 
                             impact_rates,
                             period = 2000:2019,
                             per_adm1 = FALSE, 
                             get_boots = FALSE, 
                             global = FALSE){
  
  if(per_adm1){
    
    vars1 <- c("iso3", "boots_id", "scenario", "ADM1_NAME", "vaccine")
    vars2 <- c("iso3", "scenario", "ADM1_NAME", "vaccine")
    
  }else if(global){
    
    vars1 <- c("boots_id", "scenario", "vaccine")
    vars2 <- c("scenario", "vaccine")
    
  }else{
    vars1 <- c("iso3", "boots_id", "scenario", "vaccine")
    vars2 <- c("iso3", "scenario", "vaccine")
    
  }
  
  new_impact <- 
    new_fvps %>% 
    filter(year%in%period) %>% 
    left_join(impact_rates) %>%
    ## Calculate sub-national impact (using MCV1 coverage)
    mutate(deaths_averted_new = deaths_averted_rate*mean_fvps)
  
  
  impact_df <- 
    new_impact %>% 
    ungroup() %>% 
    group_by_at(vars1) %>% 
    summarise(deaths_averted_new = sum(deaths_averted_new, na.rm = TRUE))
  
  ## get mean estimates
  impact_df_mean <- 
    impact_df %>% 
    group_by_at(vars2) %>% 
    summarise(deaths_averted_new_mean = mean(deaths_averted_new, na.rm = TRUE), 
              deaths_averted_new_q1 =  quantile(deaths_averted_new, 0.025, na.rm = TRUE), 
              deaths_averted_new_q3 = quantile(deaths_averted_new, 0.975, na.rm = TRUE))
  
  
  ## Create label for period
  p_label <- ifelse(length(period) == 1, 
                    as.character(period), 
                    paste(min(period), max(period), sep = "-"))
  
  
  impact_df_mean$period = p_label
  impact_df$period = p_label
  
  if(get_boots == TRUE){
    return(impact_df)
  }else{
    return(impact_df_mean)
  }
}

IU_impact_calculation <- function(new_fvps, 
                                  impact_rates, 
                                  period = NULL, 
                                  per_adm1 = FALSE, 
                                  get_boots = FALSE, 
                                  global = FALSE){
  
  if(per_adm1){
    
    vars1 <- c("iso3", "boots_id", "scenario", "ADM1_NAME", "vaccine")
    vars2 <- c("iso3", "scenario", "ADM1_NAME", "vaccine")
    
  }else if(global){
    
    vars1 <- c("boots_id", "scenario", "vaccine")
    vars2 <- c("scenario", "vaccine")

  }else{
  
    vars1 <- c("iso3", "boots_id", "scenario", "vaccine")
    vars2 <- c("iso3", "scenario", "vaccine")
    
  }
  
  
  if(!is.null(period)){
  
    if(global){
      
      impact_df_mean <- 
        new_fvps %>% 
        IU_impact_period(new_fvps =., 
                         impact_rates = impact_rates, 
                         period = period, 
                         per_adm1 = per_adm1, 
                         get_boots = get_boots, 
                         global = TRUE)
      
      impact_df = impact_df_mean
      
    }else{
    
      ## Due to memory issues I need to the the table joinings and impact calculation per country
      country_list <- unique(new_fvps$iso3)
      
      impact_df_mean <- 
        foreach(i = seq_along(country_list), .combine = rbind)%do%{
          
          new_fvps %>% 
            filter(iso3 == country_list[i]) %>% 
            IU_impact_period(new_fvps =., 
                             impact_rates = impact_rates, 
                             period = period, 
                             per_adm1 = per_adm1, 
                             get_boots = get_boots)
          
          
        }
      
      impact_df = impact_df_mean
      
    } 
    
  }else{
    
    vars1 <- c(vars1, "year")
    vars2 <- c(vars2, "year")
    
    
    new_impact <- 
      new_fvps %>% 
      left_join(impact_rates) %>%
      ## Calculate sub-national impact (using MCV1 coverage)
      mutate(deaths_averted_new = deaths_averted_rate*mean_fvps) %>% 
      mutate(da_ratio = deaths_averted_new/impact) %>% 
      mutate(deaths_novac_new = da_ratio*tot_deaths_novac)
    
    impact_df <- 
      new_impact %>% 
      ungroup() %>% 
      group_by_at(vars1) %>% 
      summarise(deaths_averted_new = sum(deaths_averted_new, na.rm = TRUE), 
                tot_deaths_novac = sum(deaths_novac_new, na.rm = TRUE))
    
    ## get mean estimates
    impact_df_mean <- 
      impact_df %>% 
      group_by_at(vars2) %>% 
      summarise(deaths_averted_new_mean = mean(deaths_averted_new, na.rm = TRUE), 
                deaths_averted_new_q1 =  quantile(deaths_averted_new, 0.025, na.rm = TRUE), 
                deaths_averted_new_q3 = quantile(deaths_averted_new, 0.975, na.rm = TRUE), 
                deaths_novac_mean = mean(tot_deaths_novac, na.rm = TRUE), 
                deaths_novac_q1 =  quantile(tot_deaths_novac, 0.025, na.rm = TRUE), 
                deaths_novac_q3 = quantile(tot_deaths_novac, 0.975, na.rm = TRUE))
    
  }
  
  if(get_boots == TRUE){
  
    return(impact_df)
    
  }else{
  
    return(impact_df_mean)
    
  }
    
  
}


get_dissimilarity_index <- function(new_fvps, impact_rates, pop_dat, scale_pop = 1e5){
  
  ##  WID: Calculate the weighted index of disparity (WID)
  ## as in https://dx.doi.org/10.2471%2FBLT.20.279232 and  https://sciwheel.com/fulltext/doi/10.1093/phr/117.3.273
  
  ## INS: Index of dissimilarity 
  ## as in 10.1186/s12939-016-0307-y  
  
  impact_all_y_uncer <- 
    fvps_scen_df %>% filter(scenario == "sub_cov") %>% 
    left_join(da_rate_activity_cohort) %>%
    ## Calculate sub-national impact (using MCV1 coverage)
    mutate(deaths_averted_new = deaths_averted_rate*mean_fvps) %>% 
    select(iso3, ADM0_NAME, ADM1_NAME, year, boots_id, deaths_averted_rate, deaths_averted_new, scenario, vaccine)
  
  
  tmp_all_ind <- 
    impact_all_y_uncer %>% 
    left_join(pop_dat) %>% 
    
    group_by(iso3, year, boots_id) %>% 
    mutate(deaths_averted_pop = deaths_averted_new/total_pop*scale_pop) %>% 
    mutate(deaths_averted_pop_nal = mean(deaths_averted_pop, na.rm = TRUE)) %>% 
    mutate(total_pop_nal = sum(total_pop, na.rm = TRUE)) %>% 
    mutate(shared_pop = total_pop/total_pop_nal) %>% 
    
    mutate(WID = (sum(shared_pop*abs(deaths_averted_pop-deaths_averted_pop_nal))/deaths_averted_pop_nal)*100) %>% 
    mutate(INS = 0.5*(sum(
      abs(deaths_averted_pop/deaths_averted_pop_nal * total_pop/total_pop_nal - total_pop/total_pop_nal)
    ))*100) %>% 
    mutate(WID = round(WID, 1), 
           INS = round(INS, 1)) %>% 
    distinct(iso3, year, boots_id, WID, INS)
  
  return(tmp_all_ind)
}

get_numbers_text <- function(mid, lo, hi){
  
  df <- paste0(mid, " (", lo, " - ", hi, ")")
}



violin_plots_vaccines <- function(df, years = seq(2000, 2019, 3), 
                                  col_pal = c("#00AFBB", "#E7B800", "#FC4E07")){
  
  dodge <- position_dodge(width = 0.6)
  
  df %>% 
    filter(year%in%years) %>% 
    mutate(year_f = as.factor(year)) %>% 
    group_by(iso3) %>% 
    mutate(n = n_distinct(ADM1_NAME), 
           n_state = paste0(iso3, " (", n, ")")) %>% 
    mutate(cov_per = mean_cov * 100) %>% 
    ggplot(aes(x = year_f, y = cov_per, fill = vaccine))+
    geom_violin(position = dodge, aes(colour = vaccine))+
    geom_boxplot(width=.1, position = dodge)+
    ylab("Coverage %") +
    xlab("Year")+
    facet_wrap(.~n_state)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          legend.position="none")+
    scale_fill_manual(values = col_pal)+
    scale_color_manual(values = col_pal)
  
}



plot_rel_changes <- function(rel_ch_df, 
                             impact_df,
                             scenario_im = c("0_VIMC", "max_cov"),
                             col_pal = c("grey", "#009ADE") ){
  dat_tmp <- 
    rel_ch_df %>% 
    filter(scenario%in%scenario_im) %>% 
    filter(scenario!=unique(ref_scenario)) %>% 
    filter(!is.na(vaccine) & !is.na(rel_ch)) %>% 
    mutate(dif_pop = abs_dif/total_pop*1e5)
  
  dat_tmp$iso3 = with(dat_tmp, reorder(iso3, dif_pop, max))
  
  
  p1 <- 
    dat_tmp %>% 
    ggplot(aes(x = dif_pop, y = iso3))+
    geom_col(aes(fill = scenario))+
    #geom_text(aes(label = paste0(rel_ch*100, "%")))+
    theme_cowplot()+
    xlab("Differences in deaths averted (per 100K)")+
    ylab("")+
    scale_fill_manual(values = col_pal[2])+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    
  
  
  impact_df <- impact_df %>% filter(iso3%in%unique(dat_tmp$iso3))
  impact_df$iso3 <- factor(impact_df$iso3, levels = levels(dat_tmp$iso3))
  
  
  ref_scale = 1e5
  
  p2 <- 
    impact_df %>% 
    filter(!is.na(vaccine)) %>% 
    filter(deaths_averted_new_mean!=0) %>% 
    filter(scenario%in%scenario_im) %>% 
    ggplot(aes(y = iso3, x = deaths_averted_new_mean/total_pop*ref_scale)) +
    #geom_point(aes(color = scenario), position = position_dodge(0.4))+
    geom_col(position = position_dodge(0.4), aes(fill = scenario), alpha = 0.8)+
    geom_errorbar(aes(xmin = deaths_averted_new_q1/total_pop*ref_scale, 
                      xmax = deaths_averted_new_q3/total_pop*ref_scale, colour = scenario),
                  position = position_dodge(0.4), width = 0.2)+
    theme_cowplot()+
    facet_wrap(period~.)+
    labs(y = "", x = "Deaths averted (per 100K)")+
    theme(axis.text.x = element_text(angle = 90))+
    scale_fill_manual(values = col_pal)+
    scale_color_manual(values = col_pal)+
    #theme(legend.position = "none")+
    theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  plot_grid(p1, p2, align = "h", axis = "tb", rel_widths = c(0.5,1))  
  
}

## Taken from https://github.com/grssnbchr/bivariate-maps-ggplot2-sf/blob/master/analysis/index.Rmd
theme_map <- function(...) {
  
  
  # some constants
  default_font_color <- "#4e4d47"
  default_background_color <- "#f5f5f2"
  default_font_family <- "Ubuntu Regular"
  
  
  theme_minimal() +
    theme(
      text = element_text(family = default_font_family,
                          color = default_font_color),
      # remove all axes
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      # add a subtle grid
      panel.grid.major = element_line(color = "#dbdbd9", size = 0.2),
      panel.grid.minor = element_blank(),
      
      # borders and margins
      plot.margin = unit(c(.5, .5, .2, .5), "cm"),
      panel.border = element_blank(),
      panel.spacing = unit(c(-.1, 0.2, .2, 0.2), "cm"),
      # titles
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 9, hjust = 0,
                                 color = default_font_color),
      plot.title = element_text(size = 15, hjust = 0.5,
                                color = default_font_color),
      plot.subtitle = element_text(size = 10, hjust = 0.5,
                                   color = default_font_color,
                                   margin = margin(b = -0.1,
                                                   t = -0.1,
                                                   l = 2,
                                                   unit = "cm"),
                                   debug = F),
      # captions
      plot.caption = element_text(size = 7,
                                  hjust = .5,
                                  margin = margin(t = 0.2,
                                                  b = 0,
                                                  unit = "cm"),
                                  color = "#939184"),
      ...
    )
}


rel_change_impact <- function(new_fvps,
                              impact_rates,
                              ref_round = 2, 
                              ref_scenario = "VIMC",
                              period = 2000:2019, 
                              per_year = FALSE, 
                              per_adm1 = FALSE){
  
  
  impact_all_y <- IU_impact_calculation(new_fvps = new_fvps, 
                                        impact_rates = impact_rates,
                                        period = period, 
                                        per_adm1 = per_adm1)
  
  ## Include VIMC estimates
  vimc_impact <- 
    impact_rates %>% 
    filter(year%in%period) %>% 
    group_by(iso3, vaccine) %>% 
    summarise(deaths_averted_new_mean = mean(impact, na.rm = TRUE), 
              deaths_averted_new_q1 =  quantile(impact, 0.025, na.rm = TRUE), 
              deaths_averted_new_q3 = quantile(impact, 0.975, na.rm = TRUE)) %>% 
    mutate(scenario = "VIMC") %>% 
    mutate(period = unique(impact_all_y$period))
  
  impact_all_y <- rbind(impact_all_y, vimc_impact)
  
  df_with_ref_scen <- 
    impact_all_y %>% 
    select(-deaths_averted_new_q1, -deaths_averted_new_q3) %>% 
    pivot_wider(names_from = scenario, values_from = deaths_averted_new_mean) %>% 
    select(iso3, ref_scenario_val = contains(ref_scenario), period) 
  
  
  df_with_rel_ch <- 
    df_with_ref_scen %>% 
    left_join(impact_all_y) %>% 
    mutate(rel_ch = round((deaths_averted_new_mean - ref_scenario_val)/ref_scenario_val, ref_round)) %>%
    mutate(abs_dif = round(deaths_averted_new_mean - ref_scenario_val)) %>% 
    mutate(ref_scenario = ref_scenario) %>% 
    select(iso3, period, vaccine, scenario, ref_scenario, rel_ch, abs_dif)
  
  return(list(impact_est = impact_all_y, rel_ch = df_with_rel_ch))
}


dif_impact_from_samples <- function(new_fvps, 
                                    impact_rates,
                                    period = 2019, 
                                    per_adm1 = FALSE, 
                                    get_boots = TRUE, 
                                    global = FALSE){
  
  impact_all_y <- IU_impact_calculation(new_fvps = new_fvps, 
                                        impact_rates = impact_rates,
                                        period = period, 
                                        per_adm1 = per_adm1, 
                                        get_boots = get_boots, 
                                        global = global)
  
  if(global){
  
    pop_tmp <- 
      live_births %>% 
      mutate(year = cohort) %>% 
      filter(year %in% period) %>% 
      mutate(vaccine = "MCV1") %>% 
      group_by(vaccine) %>% 
      summarise(total_pop = sum(pop, na.rm = TRUE))
    
    impact_all_y <- 
      impact_all_y %>% 
      left_join(pop_tmp) %>% 
      mutate(impact_pop = deaths_averted_new/total_pop*1e5)
    
    ref_sc <- 
      impact_all_y %>% 
      ungroup() %>% 
      filter(scenario == "sub_cov") %>% 
      select(boots_id, period, sub_cov = impact_pop, sub_cov_da = deaths_averted_new) 
    
    estimates_df <- 
      impact_all_y %>% 
      left_join(ref_sc) %>% 
      mutate(diff_impact = impact_pop - sub_cov) %>% 
      mutate(diff_da = deaths_averted_new - sub_cov_da) %>% 
      group_by(scenario, period) %>% 
      summarise(diff_impact_mean = mean(diff_impact, na.rm = TRUE), 
                diff_impact_q1 =  quantile(diff_impact, 0.025, na.rm = TRUE), 
                diff_impact_q3 = quantile(diff_impact, 0.975, na.rm = TRUE), 
                diff_da_mean = mean(diff_da, na.rm = TRUE), 
                diff_da_q1 =  quantile(diff_da, 0.025, na.rm = TRUE), 
                diff_da_q3 = quantile(diff_da, 0.975, na.rm = TRUE))
    
    estimates_df_1 <- 
      estimates_df %>% 
      filter(scenario%in%c("sub_cov","at_nat_cov", "max_cov", "at_gvap_cov")) %>% 
      mutate_at(vars(contains("diff_impact")), ~round(.)) %>%
      mutate_at(vars(contains("diff_impact")), ~prettyNum(.,big.mark=",")) %>% 
      mutate(values = paste0(diff_impact_mean, 
                             " (", diff_impact_q1, "–", 
                             diff_impact_q3, ")")) %>% 
      select(scenario, period, values) %>% 
      pivot_wider(names_from = scenario, values_from = values) %>% 
      select(period, at_nat_cov, max_cov, at_gvap_cov) %>% 
      mutate(estimates = "per 100,000 live births")
    
    
    estimates_df_2 <- 
      estimates_df %>% 
      filter(scenario%in%c("sub_cov","at_nat_cov", "max_cov", "at_gvap_cov")) %>% 
      mutate_at(vars(contains("diff_da")), ~round(.)) %>%
      mutate_at(vars(contains("diff_da")), ~prettyNum(.,big.mark=",")) %>% 
      mutate(values = paste0(diff_da_mean, 
                             " (", diff_da_q1, "–", 
                             diff_da_q3, ")")) %>% 
      select(scenario, period, values) %>% 
      pivot_wider(names_from = scenario, values_from = values) %>% 
      select(period, at_nat_cov, max_cov, at_gvap_cov) %>% 
      mutate(estimates = "absolute")
    
    estimates_df_f <- bind_rows(estimates_df_1, estimates_df_2)
    
  }else{
  
    pop_tmp <- 
      live_births %>% 
      mutate(year = cohort) %>% 
      filter(year %in% period) %>% 
      group_by(iso3) %>% 
      summarise(total_pop = sum(pop, na.rm = TRUE))
    
    impact_all_y <- 
      impact_all_y %>% 
      left_join(pop_tmp) %>% 
      mutate(impact_pop = deaths_averted_new/total_pop*1e5)
    
    ref_sc <- 
      impact_all_y %>% 
      ungroup() %>% 
      filter(scenario == "sub_cov") %>% 
      select(iso3, boots_id, period, sub_cov = impact_pop) 
    
    estimates_df <- 
      impact_all_y %>% 
      left_join(ref_sc) %>% 
      mutate(diff_impact = impact_pop - sub_cov) %>% 
      group_by(iso3, scenario, period) %>% 
      summarise(diff_impact_mean = mean(diff_impact, na.rm = TRUE), 
                diff_impact_q1 =  quantile(diff_impact, 0.025, na.rm = TRUE), 
                diff_impact_q3 = quantile(diff_impact, 0.975, na.rm = TRUE))
    
    estimates_df_f <- 
      estimates_df %>% 
      filter(iso3!="IND") %>% 
      filter(scenario%in%c("sub_cov","at_nat_cov", "max_cov", "at_gvap_cov")) %>% 
      mutate_at(vars(contains("diff_impact")), ~round(.)) %>%
      mutate_at(vars(contains("diff_impact")), ~prettyNum(.,big.mark=",")) %>% 
      mutate(values = paste0(diff_impact_mean, 
                             " (", diff_impact_q1, "–", 
                             diff_impact_q3, ")")) %>% 
      select(iso3, scenario, period, values) %>% 
      pivot_wider(names_from = scenario, values_from = values) %>% 
      select(iso3, period, at_nat_cov, max_cov, at_gvap_cov)
    
  }
  
  return(estimates_df_f)  
}    


