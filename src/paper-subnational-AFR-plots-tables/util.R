plot_ribbons <- function(all_dat, 
                         col_pal = c("grey", "#009ADE", "#00CD6C", "#AF58BA", "#FFC61E"), 
                         scen_to_plot = NULL){
  
  if(is.null(scen_to_plot)){
    scen_to_plot = unique(all_dat$scenario)
  }
  
  scale_ref <- 1e3
  
  all_dat %>% 
    filter(scenario%in%scen_to_plot) %>% 
    ggplot(aes(y = deaths_averted_new_mean/scale_ref, x = year, fill = scenario))+
    geom_ribbon(aes(ymin = deaths_averted_new_q1/scale_ref, ymax = deaths_averted_new_q3/scale_ref, x = year), 
                alpha = 0.6)+
    geom_line(aes(linetype = scenario))+
    facet_wrap(.~iso3, scales = "free_y")+
    xlab("Year")+
    ylab("Deaths averted (thousands)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values = col_pal)+
    scale_color_manual(values = col_pal)
  
}


plot_impact_ratios <- function(df, 
                         col_pal = c("#009ADE", "#00CD6C", "#AF58BA", "#FFC61E"), 
                         scen_to_plot = NULL){
  
 
  df %>% 
    mutate(scenario = "VIMC") %>% 
    ggplot(aes(y = deaths_averted_rate_mean, x = Year, fill = scenario))+
    geom_ribbon(aes(ymin = deaths_averted_rate_lo, ymax = deaths_averted_rate_hi, x = Year), 
                alpha = 0.6)+
    geom_line(aes(linetype = scenario))+
    facet_wrap(.~iso_code, scales = "free_y")+
    xlab("Year")+
    ylab("Deaths averted rates")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), 
          legend.position="none")+
    geom_text(
      mapping = aes(x = -Inf, y = -Inf, label = values),
      hjust   = -0.8,
      vjust   = -8, 
      size = 3.2, 
      family = "serif"
    )+
    scale_fill_manual(values = col_pal)+
    scale_color_manual(values = col_pal)
  
}

bar_plot_combine <- function(rel_ch_df, ref_order, scenario_im, ref_scenario, col_pal){
  
  dat_tmp <- 
    rel_ch_df %>% 
    filter(scenario%in%scenario_im) %>% 
    filter(scenario!=unique(ref_scenario)) %>% 
    filter(!is.na(vaccine) & !is.na(rel_ch)) %>% 
    mutate(dif_pop = abs_dif/total_pop*1e5) %>% 
    left_join(ref_order)
  
  
  ## REf countries
  ref_co <-
    dat_tmp %>% 
    filter(scenario == "at_nat_cov") %>% 
    arrange(desc(INS_mean))
  
  dat_tmp$name = with(dat_tmp, reorder(name, dif_pop, max))
  
  
  p1 <- 
    dat_tmp %>% 
    filter(name%in%ref_co$name[1:10]) %>% 
    mutate(scenario = fct_relevel(scenario, "at_nat_cov", "max_cov", "at_gvap_cov")) %>% 
    ggplot(aes(x = dif_pop, y = name))+
    geom_col(position = position_dodge(0.6), aes(fill = scenario))+
    theme_cowplot()+
    xlab("Differences in deaths averted (per 100K)")+
    ylab("")+
    scale_fill_manual(values = col_pal[-1])+
    theme(legend.position = "none")+
    theme(axis.text.x = element_text(angle = 90))+
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
  
  p1
}

get_map_dat <- function(dat_ch, scenario_name, pop_dat, period, ref_scale = 1e4){
  
  dat_map <- 
    dat_ch %>% 
    filter(iso3!= "IND") %>% 
    filter(scenario%in%scenario_name) %>%
    filter(!is.na(vaccine)) %>%
    left_join(
      pop_dat %>% 
        filter(year%in%period) %>% 
        group_by(ADM1_NAME, iso3) %>% 
        summarise(tot_pop = sum(total_pop, na.rm = TRUE))) %>% 
    mutate(dif_wpop = abs_dif/tot_pop*ref_scale)
  
  return(dat_map)
}

plot_map_all <- function(df, adm1_map, adm0_map, var_name, name_plot = "", countries_hl = NULL, legend_name = ""){
  
  
  theme_for_map <- 
    theme(
      axis.ticks = element_blank(),
      axis.text= element_blank(),
      axis.line = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_line(color='transparent'),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent",color='transparent')
    )
  
  sf::sf_use_s2(FALSE)
  
  maps_df <- 
    adm1_map %>% 
    filter(GID_0%in%unique(df$iso3)) %>% 
    left_join(df %>% 
                select(GID_0 = iso3, ADM1_NAME, contains(var_name), scenario), 
              by = c("GID_0" = "GID_0", "ADM1_NAME" = "ADM1_NAME")) %>% 
    filter(!is.na(scenario))
  
  maps_df$to_plot <- maps_df[[var_name]]
  
  maps_df <- st_crop(maps_df, ymin = -35, ymax = 37,  xmin = -20, xmax = 51)
  AF_map <- st_crop(adm0_map, ymin = -35, ymax = 37, xmin = -20, xmax = 51)
  
  if(!is.null(countries_hl)){
    
    suppressMessages(
      ggplot()+
        geom_sf(data = maps_df, aes(fill = to_plot), color = "grey90", size = 0.3)+
        geom_sf(fill = "transparent", color = "#4e4d47", size = 0.5, 
                data = AF_map %>% group_by(adm0_a3) %>% summarise()) +
        geom_sf(fill = "transparent", color = "black", size = 0.8, 
                data =  maps_df %>% filter(adm0_a3 %in%countries_hl) %>% 
                  group_by(adm0_a3) %>% summarise())+
        scale_fill_viridis(option = "inferno",
                           name = legend_name,
                           #alpha = 0.7, # make fill a bit brighter
                           begin = 0.1, # this option seems to be new (compared to 2016):
                           # with this we can truncate the
                           # color scale, so that extreme colors (very dark and very bright) are not
                           # used, which makes the map a bit more aesthetic
                           #end = 0.9,
                           na.value = "#4e4d47") +
        #theme_for_map+
        theme_void(18)+
        labs(fill = name_plot)+
        facet_grid(cols = vars(scenario)) +
        #theme_bw()+
        theme(legend.position = "bottom", 
              legend.key.height = unit(0.5, "cm"),
              legend.key.width = unit(2, "cm"), 
              plot.margin = margin(-2, -2, -2, -2, "cm"))
    ) 
    
  }else{
    
    suppressMessages(
      ggplot()+
        geom_sf(data = maps_df, aes(fill = to_plot), color = "gray60", size = 0.2)+
        geom_sf(fill = "transparent", color = "black", size = 0.4, 
                data = AF_map %>% group_by(adm0_a3) %>% summarise()) +
        scale_fill_viridis(option = "inferno",
                           name = legend_name,
                           #alpha = 0.7, # make fill a bit brighter
                           begin = 0.1, # this option seems to be new (compared to 2016):
                           # with this we can truncate the
                           # color scale, so that extreme colors (very dark and very bright) are not
                           # used, which makes the map a bit more aesthetic
                           #end = 0.9,
                           na.value = "#4e4d47") +
        #theme_for_map+
        theme_void(16)+
        labs(fill = name_plot)+
        facet_grid(cols = vars(scenario)) +
        #theme_bw()+
        theme(legend.position = "bottom", 
              #legend.key.height = unit(0.5, "cm"),
              legend.key.width = unit(1.5, "cm"))
    )
  }
  
}


plot_map_unique<- function(df, adm1_map, adm0_map, var_name, name_plot = "", legend_name = ""){
  
  sf::sf_use_s2(FALSE)
  
  maps_df <- 
    adm1_map %>% 
    filter(GID_0%in%unique(df$iso3)) %>% 
    left_join(df %>% 
                select(GID_0 = iso3, ADM1_NAME, contains(var_name)), 
              by = c("GID_0" = "GID_0", "ADM1_NAME" = "ADM1_NAME"))
  
 
  maps_df$to_plot <- maps_df[[var_name]]
  
  maps_df <- st_crop(maps_df, ymin = -35, ymax = 37,  xmin = -20, xmax = 51)
  AF_map <- st_crop(adm0_map, ymin = -35, ymax = 37, xmin = -20, xmax = 51)
  
  suppressMessages(
    ggplot()+
      geom_sf(data = maps_df, aes(fill = to_plot), color = "gray60", size = 0.2)+
      geom_sf(fill = "transparent", color = "black", size = 0.4, 
              data = AF_map %>% group_by(adm0_a3) %>% summarise()) +
      scale_fill_viridis(option = "viridis",
                         name = legend_name,
                         na.value = "#4e4d47") +
      #theme_for_map+
      theme_void(18)+
      labs(fill = name_plot) +
      theme(#legend.position = "bottom", 
        #           legend.key.height = unit(0.5, "cm"),
        #           legend.key.width = unit(1.5, "cm"), 
        plot.margin = margin(-3, -1, -3, -1, "cm"))
  )
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
    ylab("Adjusted coverage %") +
    facet_wrap(vars(n_state))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values = col_pal)+
    scale_color_manual(values = col_pal)
  
}