

# parallel_subsets <- function(subset, original_data, C_candidates, rc, rcj, S_cand, k_cand, kg_cand, maxit = 30) {
#
#
# print(paste("subset:", subset))
# temp <- make_subsamples(original_data, subset)
# Y <- temp[[1]]
# X <- temp[[2]]
# g_true <- temp[[3]]
#
# df_results <- initialise_df_results(use_robust = TRUE) #df_results will contain results for each configuration (for example adjusted randindex, ...)
# df_pic <- data.frame(t(matrix(rep(NA, length(C_candidates))))) #df_pic will contain PIC for each C and each subset
#
# index_configuration <- 0
# for(S in S_cand) { #loop over the possible number of groups
#   print(paste("S:", S))
#   for(k in k_cand) { #loop over the possible number of common factors
#     #define set of possible combinations of number of group factors
#     kg_candidates <- define_kg_candidates(S, min(kg_cand), max(kg_cand), TRUE)
#
#     for(w in 1:nrow(kg_candidates)) {
#
#       index_configuration <- index_configuration + 1
#       kg <- unlist(kg_candidates[w,])
#
#
#       ########## initialisation
#
#       iteration <- 0 #number of the iteration; 0 indicates being in the initialisation phase
#       beta_est <- initialise_beta(use_robust = TRUE, Y, X, S)
#       #initial grouping
#       g <- initialise_clustering(use_robust = TRUE, Y, g, beta_est, S, k, kg, comfactor, max_percent_outliers_tkmeans = 0)
#       #initial common factorstructure
#       temp <- initialise_commonfactorstructure_macropca(use_robust = TRUE, Y, X, beta_est, g, NA, k, kg)
#       comfactor <- temp[[1]]
#       lambda <- temp[[2]]
#       #initial group specific factorstructure
#       factor_group <- estimate_factor_group(use_robust = TRUE, Y, X, beta_est, g, NA, NA, NA, S, k, kg, initialise = TRUE)
#
#       lambda_group <- calculate_lambda_group(use_robust = TRUE, Y, X, beta_est, factor_group, g, NA, NA, S, k, kg, initialise = TRUE)
#
#
#       ######### estimations
#       obj_funct_values = c()
#       speed = 999999 #convergence speed: set to initial high value
#       while(iteration < maxit & !check_stopping_rules(iteration, speed, obj_funct_values, verbose = TRUE)) {
#         temp <- iterate(use_robust = TRUE, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor, S, k, kg, verbose = FALSE)
#         beta_est <- temp[[1]]
#         g <- temp[[2]]
#         comfactor <- temp[[3]]
#         lambda <- temp[[4]]
#         factor_group <- temp[[5]]
#         lambda_group <- temp[[6]]
#         value <- temp[[7]]
#         if(!exists("i")) i = 1
#         print(paste(c(i,"   subset:", subset, "/", max(indices_subset), "config:", S,"-", k, "-", kg[1:S], "-", "iteration", iteration, ":", round(value)), collapse = " "))
#         obj_funct_values <- c(obj_funct_values, value)
#         iteration <- iteration + 1
#         speed <- get_convergence_speed(iteration, obj_funct_values / nrow(Y) / ncol(Y))
#       }
#       print(show_gtable(g, g_true))
#
#       # calculate the estimation errors
#       pic_e2 <- calculate_error_term(Y, X, beta_est, g,
#                                      factor_group, lambda_group,
#                                      comfactor, lambda,
#                                      S, k, kg)
#
#       # calculate PIC for all C's and collect in a dataframe df_pic
#       df_pic <- add_pic(df_pic, index_configuration, use_robust = TRUE, Y, beta_est, g, S, k, kg,
#                         pic_e2, C_candidates)
#
#
#       # add results of this configuration to df_results
#       pic_sigma2 <- calculate_sigma2(pic_e2, nrow(Y), ncol(Y))
#       df_results <- df_results %>%
#         add_configuration(S, k, kg) %>%
#         add_metrics(index_configuration, pic_sigma2, g, g_true, ncol(Y), iteration)
#
#       #df_results <- add_reduced_pic4(df_results, ncol(Y)) #only for test purposes
#     }
#   }
# }
# #(MOVED TO FUNCTION)df_results <- df_results[-1,] #delete empty row
#
# #Ensure that all PIC-values use the same sigma2 by replacing the previously used configuration-dependent "sigma2" by "sigma2_max_model" (this is the sigma2 of the last configuration).
# df_pic <- adapt_pic_with_sigma2maxmodel(df_pic, df_results, pic_sigma2)
#
#
# #calculate for each candidate value for C the best S, k and kg
# all_best_values <- calculate_best_config(df_results, df_pic, C_candidates)
# return(all_best_values)
# }


#' Wrapper around the non-parallel algorithm.
#'
#' The function estimates beta, group membership and the common and group specific factorstructures for one configuration.
#' @param configs dataframe where each row contains one configuration of groups and factors
#' @param i index of the current configuration
#' @inheritParams calculate_VCsquared
#' @inheritParams estimate_beta
#' @param maxit maximum limit for the number of iterations
#' @export
run_parallel_config <- function(configs, i, C_candidates, Y, X, indices_subset, maxit = 30) {
  config <- configs[i,]
  print("start:")
  S <- config %>%
    dplyr::select(S) %>%
    dplyr::pull()
  k <- config %>%
    dplyr::select(k) %>%
    dplyr::pull()
  kg <- unlist(config %>% dplyr::select(.data$k1:.data$k20)) # must be a vector (class "integer")


  print("initialise:")
  ########## initialisation

  iteration <- 0 # number of the iteration; 0 indicates being in the initialisation phase
  beta_est <- initialise_beta(use_robust = TRUE, Y, X, S)
  # initial grouping
  g <- initialise_clustering(use_robust = TRUE, Y, g, beta_est, S, k, kg, NA, max_percent_outliers_tkmeans = 0)
  # initial common factorstructure
  temp <- initialise_commonfactorstructure_macropca(use_robust = TRUE, Y, X, beta_est, g, NA, k, kg)
  comfactor <- temp[[1]]
  lambda <- temp[[2]]
  # initial group specific factorstructure
  factor_group <- estimate_factor_group(use_robust = TRUE, Y, X, beta_est, g, NA, NA, NA, S, k, kg, initialise = TRUE)
  lambda_group <- calculate_lambda_group(use_robust = TRUE, Y, X, beta_est, factor_group, g, NA, NA, S, k, kg, initialise = TRUE)


  print("estimate:")
  ######### estimations
  obj_funct_values <- c()
  speed <- 999999 # convergence speed: set to initial high value
  while (iteration < maxit & !check_stopping_rules(iteration, speed, obj_funct_values, verbose = TRUE)) {

    #new errors occuring when using parallel system -> find out why
    temp<- tryCatch(
      iterate(use_robust = TRUE, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor, S, k, kg, verbose = FALSE),
      error = function(e) {
        message(e)
        return(e)
      }
    )
    #temp <- iterate(use_robust = TRUE, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor, S, k, kg, verbose = FALSE)
    beta_est <- temp[[1]]
    g <- temp[[2]]
    comfactor <- temp[[3]]
    lambda <- temp[[4]]
    factor_group <- temp[[5]]
    lambda_group <- temp[[6]]
    value <- temp[[7]]
    print(paste(c("   subset:", subset, "/", max(indices_subset), "config:", i, "/", nrow(configs), "(", S, "-", k, "-", kg[1:S], ")", "iteration", iteration, ":", round(value)), collapse = " "))
    obj_funct_values <- c(obj_funct_values, value)
    iteration <- iteration + 1
    speed <- get_convergence_speed(iteration, obj_funct_values / nrow(Y) / ncol(Y))
  }


  # calculate the estimation errors
  print("e2:")
  pic_e2 <- calculate_error_term(
    Y, X, beta_est, g,
    factor_group, lambda_group,
    comfactor, lambda,
    S, k, kg
  )

  print("pic:")
  pic <- add_pic_parallel(
    use_robust = TRUE, Y, beta_est, g, S, k, kg,
    pic_e2, C_candidates
  )


  # add results of this configuration to df_results
  print("pic_sigma2:")
  pic_sigma2 <- calculate_sigma2(pic_e2, nrow(Y), ncol(Y))

  print("done")


  return(list(S, k, kg, pic, pic_sigma2, g))
}

#' Calculates the PIC for the current configuration.
#'
#' @inheritParams add_pic
add_pic_parallel <- function(use_robust, Y, beta_est, g,
                             S, k, kg, pic_e2, C_candidates, method_estimate_beta = "individual",
                             vars_est = ncol(beta_est), choice_pic = "pic2022") {
  pic_sigma2 <- calculate_sigma2(pic_e2)
  pic <- sapply(C_candidates, function(x) {
    calculate_PIC(x, use_robust, S, k, kg, pic_e2, pic_sigma2,
      NN = nrow(Y), TT = ncol(Y), method_estimate_beta,
      beta_est, g, vars_est, choice_pic
    )
  })
  return(pic)
}

#' Makes a dataframe with information on each configuration.
#'
#' @param x output of the parallel version of the algorithm
#' @inheritParams kg_candidates_expand
#' @export
make_df_results_parallel <- function(x, limit_est_groups = 20) {
  df <- data.frame(
    S = unlist(x %>% purrr::map(1)),
    k_common = unlist(x %>% purrr::map(2)),
    sigma2 = unlist(x %>% purrr::map(5))
  ) %>% cbind(t(matrix(unlist(x %>% purrr::map(3)), nrow = limit_est_groups)))
  names(df)[4:(4 + limit_est_groups - 1)] <- paste0("k", 1:limit_est_groups)
  df$g <- sapply(1:nrow(df), function(y) paste((x %>% purrr::map(6))[[y]], collapse = "-"))
  return(df)
}

#' Makes a dataframe with the PIC for each configuration and each candidate C.
#'
#' @param x output of the parallel version of the algorithm
#' @export
make_df_pic_parallel <- function(x) {
  df <- t(matrix(unlist(x %>% purrr::map(4)), nrow = 2002))
  return(df)
}
