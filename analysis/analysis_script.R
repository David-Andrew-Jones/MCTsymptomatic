#' _____________________________________________________________________________
#' Create output table workbook for results
excel_output <- createWorkbook()
#' Create cancer type and stage vectors
cancer_types <- dwelltime_bystypestage$cancer_type
stages <- dwelltime_bystypestage$stage

#' _____________________________________________________________________________
#' Perform checks that all the input data aligns:
#' Symptom cluster are the same in sens_bystypestage, CPRD_revised_any, waiting_bystage_symp, & symp_cluster_prev

bind_cols(sort(as.character(colnames(sens_bystypestage %>% select(-cancer_type,-stage ))) ),
          sort(as.character(colnames(CPRD_revised_any %>% select(-cancer_type,-stage ))) ),
          sort(as.character(levels(as.factor(waiting_bystage_symp$symp_cluster)))),
          sort(as.character(levels(as.factor(symp_cluster_prev$symp_cluster)))) )

#' _____________________________________________________________________________
#' Derive stage shift percentages by cancer type and stage
#' requires dwelltime_bystypestage, sens_bystypestage, CPRD_revised_any, waiting_bystage_symp

df_all_cohort_combinations <- expand.grid(symp_cluster = c("HTA","madeup"),
                                          # Limit to just one dwell now!
                                         dwell = c("dwell1"),
                                         timepoint_of_impact = c(14, 28, 42),
                                         reduction_after = c(50, 70, 90),
                                         # Set higher sim number on main run
                                         nsim = 10000)

l_out <- pmap(df_all_cohort_combinations, f_sym_stage_shift) %>%
  set_names(., pull(unite(data = df_all_cohort_combinations %>% mutate(across(where(is.numeric), as.character)), col = "all",1:ncol(df_all_cohort_combinations),sep="_")) ) %>%
  bind_rows(., .id = "ID") %>%
  separate_wider_delim(ID, delim = "_", names = c("symp_cluster", "dwell", "timepoint_of_impact", "reduction_after", "nsim")) %>%
  select(-nsim)

#' _____________________________________________________________________________
#' Apply the above calculated stage shift percentages to the symptom cohorts.
#' Assume a population of 100,000 to get a per 100,000 results
#' Creates large table with each row being a scenario, cancer type, and stage.

df_byscen_type_stage <- tibble(CPRD_revised_any,
               population = 100000,
               ) %>%
  # Convert to longer format which will be grouped_by later
  pivot_longer(!(c(cancer_type,stage, population)), values_to = "proportion", names_to = "symp_cluster" ) %>%
  # Join on prevalence for the symptom clusters
  left_join(symp_cluster_prev,by = c("symp_cluster")) %>%
  # adjust population prevalence
  mutate(pop_prev = (prevalence/100) * population) %>%
  mutate(pop_UC_bytypestage = pop_prev* (proportion/100)) %>%
  # check total cancer - summarise(sum = sum(pop_UC_bytypestage), .by = symp_cluster )
  # combine above dataframe to the percentage stage shift outcome
  right_join(l_out %>%
               filter(cancer_type %in% levels(as.factor(CPRD_revised_any$cancer_type))) %>%
               pivot_longer(I:IV, names_to = "origin_stage", values_to = "percent_shift"),
             by = c("symp_cluster","cancer_type", "stage" = "origin_stage")) %>%
  mutate(across(percent_shift, ~replace_na(.x, 0))) %>%
  # reorder
  select(symp_cluster, dwell, timepoint_of_impact, reduction_after,
         cancer_type, stage, proportion, pop_UC_bytypestage,
         percent_shift       ) %>%
  # Coombine the combinations into one scenario column
  unite(scenario, symp_cluster, dwell, timepoint_of_impact, reduction_after) %>%
  # check total cancer - summarise(sum = sum(pop_UC_bytypestage), .by = c(scenario)     )
  mutate(num_MCT_shifted_from = pop_UC_bytypestage*(percent_shift/100)) %>%
  arrange(scenario, cancer_type, stage ) %>%
  group_by(scenario ,cancer_type) %>%
  # create number shift to column, 0 for stage IV
  mutate(num_MCT_shifted_to = case_when( stage %in% c("I", "II", "III") ~ lead(num_MCT_shifted_from),
                                         .default = 0) ) %>%
  ungroup() %>%
  mutate(pop_MCT_bytypestage = pop_UC_bytypestage - num_MCT_shifted_from + num_MCT_shifted_to) %>%
  # check summarise(sum = sum(pop_UC_bytypestage), sum2 = sum(pop_MCT_bytypestage),.by = c(scenario)     )
  select(scenario, cancer_type, stage, pop_UC_bytypestage, pop_MCT_bytypestage, percent_shift, num_MCT_shifted_from  ) %>%
  rename(`Cancer type` = cancer_type,
         `Stage` = stage,
         `Number of cancers under usual care` = pop_UC_bytypestage ,
         `Number of cancers with MCT` = pop_MCT_bytypestage,
         `Percentage of cancers shifted` = percent_shift,
         `Number of cancers shifted` = num_MCT_shifted_from) %>%
  # check ssummarise(total_control = sum(`Number of cancers under usual care`), total_MCT =  sum(`Number of cancers with MCT`), .by = c(scenario)  )
  left_join(netsurv_bystage, by = c("Cancer type" = "cancer_type", "Stage" = "stage")) %>%
  mutate(netsurv_5year =  netsurv_5year/100) %>%
  rename(`Five year net survival` = netsurv_5year) %>%
  mutate(`Number of cancer deaths at 5 years under usual care` = (1 - `Five year net survival`) *  `Number of cancers under usual care`,
         `Number of cancer deaths at 5 years with MCT` = (1 - `Five year net survival`) *  `Number of cancers with MCT`) %>%
  group_by(scenario, `Cancer type`) %>%
  mutate(`Total expected 5-year cancer-related deaths in control arm` = round(sum(`Number of cancer deaths at 5 years under usual care`), 3),
         `Total expected 5-year cancer-related deaths with MCT` = round(sum(`Number of cancer deaths at 5 years with MCT`), 3)) %>%
  ungroup()

addWorksheet(excel_output, sheetName = "1. Stage-shift percent")
writeData(excel_output, sheet = "1. Stage-shift percent", df_byscen_type_stage)


#' _____________________________________________________________________________
#' Summary tables, compressing the above table to by scenario and stage only


df_shift_byscen_stage <- df_byscen_type_stage %>%
  summarise(`Number of cancers under usual care` = round(sum(`Number of cancers under usual care`), 0),
            `Number of cancers stage-shifted from` = round(sum(`Number of cancers shifted`), 0),
            `Number of cancers with MCT` = round(sum(`Number of cancers with MCT`), 0),.by = c(scenario, Stage)) %>%
  mutate(`Percentage of cancer stage-shifted from` = paste0(round((`Number of cancers stage-shifted from`/`Number of cancers under usual care`) *100 , 1), "%" ), .before = `Number of cancers with MCT`) %>%
  # make wide and use GT table
  select(scenario, Stage, `Number of cancers under usual care`, `Number of cancers stage-shifted from`, `Percentage of cancer stage-shifted from`) %>%
  filter(Stage != "I") %>%
  pivot_wider(names_from = "Stage", values_from = c("Number of cancers under usual care", "Number of cancers stage-shifted from", "Percentage of cancer stage-shifted from"))

# Save
addWorksheet(excel_output, sheetName = "2. Stage reduction")
writeData(excel_output, sheet = "2. Stage reduction", df_shift_byscen_stage)

# Save gt format as word doc to copy in to manuscrupt - cant save into excel in gt format
# Focuses on stage-shift nyumber and percetange by symptom (rows), stage (cols), and scenario (cols)
gt(df_shift_byscen_stage %>%
     mutate(`Number (%) of cancers shifted to stage_II` = paste0(`Number of cancers stage-shifted from_II`, " (", `Percentage of cancer stage-shifted from_II`,")" ),
            `Number (%) of cancers shifted to stage_III` = paste0(`Number of cancers stage-shifted from_III`, " (", `Percentage of cancer stage-shifted from_III`,")" ),
            `Number (%) of cancers shifted to stage_IV` = paste0(`Number of cancers stage-shifted from_IV`, " (", `Percentage of cancer stage-shifted from_IV` ,")")) %>%
     select(scenario, `Number of cancers under usual care_II`:`Number of cancers under usual care_IV`, `Number (%) of cancers shifted to stage_II`:`Number (%) of cancers shifted to stage_IV`) %>%
     separate_wider_delim(scenario, delim = "_", names = c("Symptom cluster", "Dwell set", "Time point of impact", "Reduction % impact")) %>%
     # Filter down to the three impact scenarios in paper
     filter( (grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) |
               (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) |
               (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`))) %>%
     mutate(`MCT impact` = case_when((grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) ~ "High",
                                     (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) ~ "Medium",
                                     (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`)) ~ "Low",
                                     .default = NA), .before = `Dwell set`) %>%
     select(`Symptom cluster`, `MCT impact`, ends_with("_I"), ends_with("_II"), ends_with("_III"), ends_with("_IV"))) %>%
  tab_spanner_delim(
    delim = "_",
    columns = (`Number of cancers under usual care_II`:`Number (%) of cancers shifted to stage_IV`),
    split = c("first","last"),
    limit = NULL,
    reverse = TRUE
  ) %>%
  gtsave(filename = "table_2.docx", path =  "C:/Users/davidj/Documents/MCTsymptomatic/analysis/raw_results/")


#' mortality by scenario  and stage
df_mort_byscen <- df_byscen_type_stage %>%
  summarise(`Number of cancer deaths at 5 years under usual care` = round(sum(`Number of cancer deaths at 5 years under usual care`), 0),
            `Number of cancer deaths at 5 years with MCT` = round(sum(`Number of cancer deaths at 5 years with MCT`), 0),.by = c(scenario)) %>%
  mutate(`Number of cancer deaths at 5 years averted` = `Number of cancer deaths at 5 years under usual care` - `Number of cancer deaths at 5 years with MCT`,
         `Percentage averted` = paste0(round((`Number of cancer deaths at 5 years averted`/`Number of cancer deaths at 5 years under usual care`) *100 , 1), "%" ),
         .before = `Number of cancer deaths at 5 years with MCT`)
# Save
addWorksheet(excel_output, sheetName = "3. Mortality reduction")
writeData(excel_output, sheet = "3. Mortality reduction", df_mort_byscen)

# Save gt format as word doc to copy in to manuscrupt - cant save into excel in gt format
# Note gtsave can't handle ~, therefore I've had to specify my c drive, feel free to change.
gt(df_mort_byscen %>%
     mutate(`Number (%) of cancer deaths at 5 years averted` = paste0(`Number of cancer deaths at 5 years averted`, " (", `Percentage averted`,")" )) %>%
     select(scenario,  `Number of cancer deaths at 5 years under usual care`, `Number (%) of cancer deaths at 5 years averted`) %>%
     separate_wider_delim(scenario, delim = "_", names = c("Symptom cluster", "Dwell set", "Time point of impact", "Reduction % impact")) %>%
     filter( (grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) |
               (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) |
               (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`))) %>%
     mutate(`MCT impact` = case_when((grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) ~ "High MCT impact",
                                     (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) ~ "Medium MCT impact",
                                     (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`)) ~ "Low MCT impact",
                                     .default = NA), .before = `Dwell set`) %>%
     select(-(`Dwell set`:`Reduction % impact`)) %>%
     pivot_wider(names_from = "MCT impact", values_from = c("Number (%) of cancer deaths at 5 years averted"))) %>%
  tab_spanner(., label = "Number (%) of cancer deaths at 5 years averted", columns = `High MCT impact`:`Low MCT impact`) %>%
  gtsave(filename = "table_3.docx", path =  "C:/Users/davidj/Documents/MCTsymptomatic/analysis/raw_results/")


saveWorkbook(excel_output, "~/MCTsymptomatic/analysis/raw_results/all_raw_results.xlsx",
             overwrite = TRUE)


#' _____________________________________________________________________________
#' Figures

#' ____________________________________
#' Manuscript figure 3 - radar plot of deaths by cancer type for each symptom group

fig_3 <- df_byscen_type_stage %>%
  select(scenario, `Cancer type`, `Total expected 5-year cancer-related deaths in control arm`, `Total expected 5-year cancer-related deaths with MCT`) %>%
  #mutate(`Cancer type` = case_when(`Cancer type` == "Head and neck" ~ "Head\nand\nneck", .default = `Cancer type`)) %>%
  slice(1, .by = c(scenario , `Cancer type`)) %>%
  mutate(`Deaths averted` = `Total expected 5-year cancer-related deaths in control arm` - `Total expected 5-year cancer-related deaths with MCT`) %>%
  mutate(`Total deaths averted in scenario` = sum(`Deaths averted`), .by = scenario) %>%
  mutate(`Proportion averted by cancer type` = round((`Deaths averted`/ `Total deaths averted in scenario`) *100, 4)) %>%
  # check summarise(check = sum(`Proportion averted by cancer type`), .by = scenario )
  separate_wider_delim(scenario, delim = "_", names = c("Symptom cluster", "Dwell set", "Time point of impact", "Reduction % impact")) %>%
 # filter(`Symptom cluster` == symp_cluster) %>%
  filter( (grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) |
            (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) |
            (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`))) %>%
  mutate(`MCT impact` = case_when((grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) ~ "High",
                                  (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) ~ "Medium",
                                  (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`)) ~ "Low",
                                  .default = NA), .before = `Dwell set`) %>%
  # Proportions should be similar across impacts, so limit to high
  filter(`MCT impact` == "High") %>%
  select(-(`Dwell set`:`Reduction % impact`))%>%
  ggplot() +
  geom_bar(aes(y = `Proportion averted by cancer type`, x = fct_rev(`Cancer type`)),
           stat = "identity") +
  facet_wrap(`Symptom cluster` ~. , ncol = 2) +
  theme_bw() +
  coord_flip() +
  labs(x = "Cancer type",
       y = "Proportion of total cancer deaths averted") +
  # Scale y axis so bars don't start in the center
  scale_y_continuous(
    limits = c(-5, 40), expand = c(0, 0),  breaks = c(0, 10,20,30 ), labels = scales::label_percent(scale = 1)) +
  theme(legend.position="none",
        text = element_text(size = 16, face = "bold"),
        strip.background = element_blank(),
        strip.text = element_text(size = 18, face = "bold"))

ggsave( "~/MCTsymptomatic/analysis/raw_results/figure_2.png", fig_3 ,device = "png",
        width =300,
        height = 200,
        units = "mm")






#' ________________________________________________________________________
#' OLD
#' ____________________________________
#' Figure - bar plot of deaths by symtom and impact - maybee just the table of this?
#' calculate scalar transformation for the percentage change
# v_mort_percent_transform <- df_mort_byscen %>%
#   mutate(percent = (`Number of cancer deaths at 5 years averted`/`Number of cancer deaths at 5 years under usual care`) *100) %>%
#   mutate(transform = `Number of cancer deaths at 5 years averted`/ percent ) %>%
#   slice(1) %>%
#   pull(transform)
#
#
# fig_3 <- df_mort_byscen %>%
#   separate_wider_delim(scenario, delim = "_", names = c("Symptom cluster", "Dwell set", "Time point of impact", "Reduction % impact")) %>%
#   filter( (grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) |
#             (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) |
#             (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`))) %>%
#   mutate(`MCT impact` = case_when((grepl(c("14"), `Time point of impact`) & grepl(c("90"),`Reduction % impact`)) ~ "High",
#                                   (grepl(c("28"), `Time point of impact`) & grepl(c("70"),`Reduction % impact`)) ~ "Medium",
#                                   (grepl(c("42"), `Time point of impact`) & grepl(c("50"),`Reduction % impact`)) ~ "Low",
#                                   .default = NA), .before = `Dwell set`) %>%
#   mutate(`MCT impact` = factor(`MCT impact`, levels = c("High", "Medium", "Low"))) %>%
#   select(-(`Dwell set`:`Reduction % impact`))%>%
#   ggplot() +
#   geom_bar(aes(y = `Number of cancer deaths at 5 years averted`, x = `MCT impact`, fill = `MCT impact` ),
#            stat = "identity",
#            position = 'dodge') +
#   facet_grid(`Symptom cluster` ~  ., switch = "y") +
#   labs(fill = "Percentage reduction in time to diagnosis",
#        x = "Time point (days since test) of impact",
#        y = "Number of 5-year cancer deaths averted") +
#   scale_y_continuous(
#
#     sec.axis = sec_axis(~ . / v_mort_percent_transform, labels = scales::label_percent(scale = 1), name = "Reduction in mortality")
#   ) +
#   theme_bw() +
#   theme(legend.position="none",
#         axis.title = element_text(face = "bold", size = 16),
#         strip.background = element_rect(fill="white"),
#         text = element_text(size = 14),
#         axis.text.x = element_text(face = "bold"),
#         axis.text.y = element_text(face = "bold"),
#         strip.text.y.left = element_text(face = "bold",angle = 0),
#         strip.text.x = element_text(face = "bold", angle = 0)) +
#   scale_fill_discrete(labels=c("50%", "70%", "90%"),
#                       type=c("#002147","#D9D8D6", "#426A5A"))
#
#
# ggsave( "~/MCTsymptomatic/analysis/raw_results/figure_3.png", fig_3 ,device = "png",
#         width =220,
#         height = 130,
#         units = "mm")



#'
#' #' Figure - radar plot of cancer type and stage distributions by cancer type and stage - maybe for supplement?
#' fig_2 <- CPRD_revised_any %>%
#'   pivot_longer(HTA:madeup, names_to = "Symptom cluster", values_to = "Percent") %>%
#'   mutate(stage = factor(stage, levels = c("I", "II", "III", "IV"))) %>%
#'   ggplot() +
#'   geom_bar(aes(y = `Percent`, x = `cancer_type`, fill = forcats::fct_rev(`stage` )),
#'            stat = "identity",
#'            position = 'stack') +
#'   facet_wrap(`Symptom cluster` ~ ., ncol =2) +
#'   coord_polar(, clip = FALSE) +
#'   labs(fill = "Stage",
#'        x = "",
#'        y = "") +
#'   theme_minimal() +
#'   # Make custom panel grid
#'   geom_hline(
#'     aes(yintercept = y),
#'     data.frame(y = c(0,5,10,15,20, 25 )),
#'     color = "lightgrey"
#'   ) +
#'   # Annotate custom scale inside plot
#'   annotate(x = "Lymphoma", y = 10, label = "10%", geom = "text", color = "gray12") +
#'   annotate(x = "Lymphoma", y = 20, label = "20%", geom = "text", color = "gray12") +
#'   # Scale y axis so bars don't start in the center
#'   scale_y_continuous(
#'     limits = c(-5, 25),
#'     expand = c(0, 0),
#'     breaks = c(0, 5, 10, 15,20, 25)
#'   ) +
#'   theme(legend.position="bottom",
#'         axis.title =element_blank(),
#'         text = element_text(size = 14),
#'         axis.text.y = element_blank(),
#'         strip.text.x = element_text(size = 18, face = "bold", angle = 0),
#'         panel.background = element_rect(fill = "white", color = "white"),
#'         panel.border = element_rect(fill = NA, color = NA),
#'         panel.grid = element_blank(),
#'         panel.grid.major.x = element_blank(),
#'         strip.background = element_blank(),
#'         panel.spacing = unit(50, "pt")) +
#'   scale_fill_discrete(type=c("#D9D8D6", "#426A5A", "#776885", "#002147"))
#'
#' ggsave( "~/MCTsymptomatic/analysis/raw_results/figure_2.png", fig_2 ,device = "png",
#'         width =220,
#'         height = 130,
#'         units = "mm")
