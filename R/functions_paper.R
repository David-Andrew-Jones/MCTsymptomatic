#' Runs the model using stated inputs
#'
#' Runs the model using stated inputs
#'
#' @param symp_cluster # A character vector corresponding to a column in waiting_bystage_symp and sens_bystypestage
#' @param dwell # A character vector corresponding to a column in dwelltime_bystypestage
#' @param timepoint_of_impact A scalar e.g. 42, giving the number of days after which the MCT has an effect
#' @param reduction_after A scalar, e.g. 50, giving the percentage reduction in waiting time due to the MCT
#' @param nsim number of Monte Carlo simulations
#'
#' @return dataframe with stageshift percentages for each cancer type and stage combination
#' @examples
#' f_sym_stage_shift()
f_sym_stage_shift <- function( symp_cluster,
                               dwell,
                               timepoint_of_impact,
                               reduction_after,
                               nsim){

  #' Check args - these will be updated when the dwell sets and symptom groups are fully decided upon
  #' Current thinking is one dwell set
  #' And to limit results to three timepoint_of_impact X reduction_after combinations
  if(!(dwell %in% c("dwell1"))){

    stop("dwell argument must be either dwell1, dwell2, dwell3")
  }

  if(!(symp_cluster %in% c("HTA","madeup"))){

    stop("symp_cluster argument must be either HTA or madeup")

  }

  if(timepoint_of_impact < 0 | reduction_after < 0){

    stop("timepoint_of_impact and reduction_after arguments must be zero or greater")

  }

  #' An internal function which calculates the percentage stage shift for a specific cancer type and stage, flows through the above function arguments
  f_internal <- function(cancer_type,
                         clinical_stage,
                         dwell_2 = dwell,
                         timepoint_of_impact2 = timepoint_of_impact,
                         reduction_after2 = reduction_after,
                         symp_cluster_2 = symp_cluster,
                         nsim2 = nsim){

    #' First, get the cancer type and stage specific dwell time for the usual care stage
    get_dwell <- dwelltime_bystypestage %>%
      # N.b. !!enquo() used here do that the character string dwell_2 flows through to dplyr's select filter and pull
      select(cancer_type, stage, !!enquo(dwell_2)) %>%
      filter(cancer_type == !!enquo(cancer_type)) %>%
      filter(stage == !!enquo(clinical_stage)) %>%
      pull(!!enquo(dwell_2))

    #' Next, get cancer symptom and stage specific waiting time distribution parameters (mean log and sd log)
    get_wait <- waiting_bystage_symp %>%
      select(symp_cluster, param, !!enquo(clinical_stage)) %>%
      filter(symp_cluster == !!enquo(symp_cluster_2))

    #' Then, create a datrame the length of the number of simulations
    #' and determine for each simulated patient whether they are all shifted or not
    stage_shift_prop <- data.frame(
      # For each simulated patient sample their indidivudal dwell
      rand_sample_dwell = rexp(nsim2, 1/get_dwell),
      # To speed up the run time, the place in the dwell is sampled for each patient here
      # in the interval 0-1, this is then scaled down to the sampled dwell later
      rand_unif = runif(nsim2, min = 0, max = 1),
      # Sample usual care waiting time for each individual
      rand_waiting_time = rlnorm(nsim2, meanlog = get_wait[1, 3], sdlog = get_wait[2, 3])) %>%
      # Generate the MCT waiting time counterfactual
      mutate(MCT_waiting_time = case_when(
        # If the usual care waiting time is greater or equal to the time point of impact but the percentage reduction induces a waiting time reduction
        # greater than the specified time point, then the MCT counterfactual is capped as the time point of impact
        rand_waiting_time >= timepoint_of_impact2 & (rand_waiting_time * ((100-reduction_after2)/100)) <=    timepoint_of_impact2  ~ timepoint_of_impact2,
        # If it does not reach the time point of impact then perform the percentage reduction.
                                          rand_waiting_time * ((100-reduction_after2)/100) >   timepoint_of_impact2 ~ rand_waiting_time * ((100-reduction_after2)/100),
        # Else (i.e. the sampled waiting time is less than the time point of impact), the MCT does not change anything
                                          .default = rand_waiting_time)) %>%
      # Find the difference between usual care and and MCT waiting time
      mutate(reduction = rand_waiting_time - MCT_waiting_time) %>%
      # If there is a reduction and it's greater than the dwell then the patient shifts
      # First, adjust dwell for point of diagnosis
      mutate(rand_diag_in_dwell = rand_unif * rand_sample_dwell) %>%
      # Find difference between time point in dwell diagnosed and the reduction - simulations with negative difference_to_dwell would have been found at earlier stage
      mutate(difference_to_dwell = rand_diag_in_dwell - (reduction/365.25)) %>%
      # Sum up the proportion with negative  difference_to_dwell
      summarise(perc_shifted = (sum(difference_to_dwell <0) / nsim2) * 100  ) %>%
      mutate(cancer_type = cancer_type,
             Stage = clinical_stage, .before = perc_shifted)

    return(stage_shift_prop)
  }


  #' _____
  #' This internal function is then applied to each posible cancer type and stage combination

  cancer_types <- dwelltime_bystypestage$cancer_type
  stages <- dwelltime_bystypestage$stage

  #' Run for all cancer type and stage combinations
  all_perc <- purrr::map2_df(cancer_types, stages, f_internal) %>%
    #' Adjust percentage stage shift for test sensitivity
    left_join( sens_bystypestage %>% select(cancer_type, stage, !!enquo(symp_cluster)),
               by = c("cancer_type", "Stage" = "stage")) %>%
    rename(sensitivity = !!enquo(symp_cluster)) %>%
    group_by(cancer_type) %>%
    rename(Stage_from = Stage) %>%
    mutate(Stage_to = lag(Stage_from), .after = "Stage_from") %>%
    mutate(stage_before_sens = lag(sensitivity)) %>%
    ungroup() %>%
    #' tidy up
    mutate(perc_shifted = round((perc_shifted * stage_before_sens), 2)) %>%
    select(-Stage_to, -sensitivity, -stage_before_sens) %>%
    pivot_wider(names_from = Stage_from, values_from = perc_shifted)

  return(all_perc)

}
