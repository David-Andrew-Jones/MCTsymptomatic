library(tidyverse)
library(readxl)

# Useful - https://freerangestats.info/blog/2019/08/25/fitting-bins

#file <- read.csv(file = "S2_CountsPerDay10day.csv", na.string = "")
file <- readxl::read_xlsx(path = "SYMPLIFY2_CancerCountsPerDay10Intervals.xlsx", sheet = "HTA age45+ - 12m cancer")

CPRD_waitingtimes <- file %>%
  mutate(across(countstageOverall:countstageMissing, as.numeric)) %>%
  dplyr::select(-window) %>%
  mutate(across(everything(), ~replace(., is.na(.), 0))) %>%
  pivot_longer(cols = countstageOverall:countstageMissing, names_to = "stage", values_to = "count") %>%
  mutate(stage = case_when(stage == "countstageOverall" ~ "Overall",
                           stage == "countstageMissing" ~ "Missing",
                           stage == "countstage1" ~ "I",
                           stage == "countstage2" ~ "II",
                           stage == "countstage3" ~ "III",
                           stage == "countstage4" ~ "IV",
                           .default = NA
                           )) %>%
  uncount(count) %>%
  mutate(left = as.numeric(forstringr::str_extract_part(windowcat, before = TRUE, pattern = "-"))) %>%
  #mutate(left = case_when(left == 0 ~1, .default = left)) %>%
  mutate(right = as.numeric(forstringr::str_extract_part(windowcat, before = FALSE, pattern = "-"))) %>%
  mutate(day_midpoint = case_when(windowcat == "0-10" ~ 5,
                                  windowcat =="11-20"~ 15,
                                  windowcat =="21-30"~ 25,
                                  windowcat =="31-40"~ 35,
                                  windowcat =="41-50"~ 45,
                                  windowcat =="51-60"~ 55,
                                  windowcat =="61-70"~ 65,
                                  windowcat == "71-80"~ 75,
                                  windowcat =="81-90"~ 85,
                                  windowcat =="91-100"~ 95,
                                  windowcat =="101-110"~ 105,
                                  windowcat =="111-120"~ 115,
                                  windowcat =="121-130"~ 125,
                                  windowcat =="131-140"~ 135,
                                  windowcat =="141-150"~ 145,
                                  windowcat =="151-160"~ 155,
                                  windowcat =="161-170"~ 165,
                                  windowcat =="171-180"~ 175,
                                  windowcat =="181-190"~ 185,
                                  .default = NA))


# Get best fitting parameterisatins of log-normal, Gamma and weibull
f_fits <- function( site, Stage) {

  temp <- as.data.frame(CPRD_waitingtimes %>%
                           filter(cancersite == site, stage == Stage) %>%
                           dplyr::select(left, right) )

  #vog normal
  out <-  tibble(`log norm mean` = as.numeric(fitdistrplus::fitdistcens(temp, "lnorm")$estimate[[1]]),
                 `log norm SD` = as.numeric(fitdistrplus::fitdistcens(temp, "lnorm")$estimate[[2]]),
                 `Gamma shape` = as.numeric(fitdistrplus::fitdistcens(temp, "gamma")$estimate[[1]]),
                 `Gamma rate` = as.numeric(fitdistrplus::fitdistcens(temp, "gamma")$estimate[[2]]),
                 `Weibull shape` = as.numeric(fitdistrplus::fitdistcens(temp, "weibull")$estimate[[1]]),
                 `Weibull scale` = as.numeric(fitdistrplus::fitdistcens(temp, "weibull")$estimate[[2]])
         )

  return(out)

}

df_fits <- expand.grid(`Cancer site` =  c("Overall", "Bowel", "Lung" , "Gastro-Oesophageal" , "Pancreas" , "Hepatobiliary", "Lymphoma"),
                       Stage = levels(as.factor(CPRD_waitingtimes$stage)))

df_fits[,3:8] <- map2_df(as.character(df_fits$`Cancer site`), as.character(df_fits$`Stage`), f_fits)

f_fitplots <- function(site,
                       Stage,
                       lnorm_mean,
                       lnorm_sd,
                       gamma_shape,
                       gamma_rate,
                       weibull_shape,
                       weibull_scale){


  plot <- CPRD_waitingtimes %>%
    rowwise() %>%
    mutate(rsample = runif(1, min = left , max = right)) %>%
    ggplot() +
    geom_bar(aes(x = rsample, y=after_stat(density)), stat="bin", binwidth = 10) +
    theme_bw(base_size = 14) +

    stat_function(aes(colour = paste0("Log normal:", " \U003BC", " = ", round(lnorm_mean,2),
                                      ", \u03C3"," = ", round(lnorm_sd, 2) )),
                  fun = dlnorm,
                  args = list(meanlog = lnorm_mean, sdlog = lnorm_sd)) +

    stat_function(aes(colour = paste0("Gamma:", " \U03B1", " = ", round(gamma_shape,2),
                                      ", \U03BB"," = ", round(gamma_rate, 2) )),
                  fun = dgamma,
                  args = list(shape = gamma_shape, rate = gamma_rate)) +

    stat_function(aes(colour = paste0("Weibull:", " \U03B1", " = ", round(weibull_shape,2),
                                      ", \U03B8"," = ", round(weibull_scale, 2) )),
                  fun = dweibull,
                  args = list(shape = weibull_shape, scale = weibull_scale)) +
    xlim(0, 400) +
    labs(title = paste0("Cancer site: ", site , ", Stage: ", Stage),colour = "Distribution", x = "Time from symptom to diagnosis (days)") +
    theme(legend.position = "bottom") +
    guides(colour=guide_legend(nrow=3))
    ggsci::scale_color_jco()

  ggsave( paste0("C:/Users/davidj/Documents/S2_stageshift_predict/Waiting_time_dist/","waiting_time_dist_",
                 site, "_", Stage,".png"), plot ,device = "png",
          width = 200,
          height = 120,
          units = "mm")


  return(plot)
}


l_plots <- pmap(list(df_fits$`Cancer site`,
                     df_fits$Stage,
                     df_fits$`log norm mean`,
                     df_fits$`log norm SD`,
                     df_fits$`Gamma shape`,
                     df_fits$`Gamma rate`,
                     df_fits$`Weibull shape`,
                     df_fits$`Weibull scale`),
                f_fitplots)


Overall <- patchwork::wrap_plots(l_plots[[36]],
                      l_plots[[1]],
                      l_plots[[8]],
                      l_plots[[15]],
                      l_plots[[22]],
                      l_plots[[29]]
                       )


ggsave( paste0("C:/Users/davidj/Documents/S2_stageshift_predict/Waiting_time_dist/","waiting_time_dist_",
               "all_overall",".png"), Overall ,device = "png",
        width = 350,
        height = 200,
        units = "mm")

#write.csv(df_fits, file = "C:/Users/davidj/Documents/S2_stageshift_predict/Waiting_time_dist/parametric_fits.csv")

#' _____________________________________________________________________________






