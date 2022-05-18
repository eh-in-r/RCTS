
globalVariables(c("i")) #required to pass R CMD check of a function which uses foreach

#' Wrapper around the non-parallel algorithm, to estimate beta, group membership and the factorstructures.
#'
#' The function estimates beta, group membership and the common and group specific factorstructures for one configuration.
# @param configs dataframe where each row contains one configuration of groups and factors
#' @param config contains one configuration of groups and factors
# @param n_configs total number of configurations
# @param i index of the current configuration
#' @inheritParams calculate_VCsquared
#' @inheritParams estimate_beta
#' @param maxit maximum limit for the number of iterations
run_config <- function(config, C_candidates, Y, X, maxit = 30) {
  print("start:")
  print(config)
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
    #print(paste(c("   subset:", subset, "/", max(indices_subset), "config:", i, "/", n_configs, "(", S, "-", k, "-", kg[1:S], ")", "iteration", iteration, ":", round(value)), collapse = " "))
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

#' Wrapper of the loop over the subsets which in turn use the parallelised algorithm.
#'
#' @param original_data list containing the original data (1: Y, 2: X, 3: true group membership)
#' @param indices_subset vector with indices of the subsets; starts with zero
#' @inheritParams get_best_configuration
#' @inheritParams calculate_VCsquared
#' @export
parallel_algorithm <- function(original_data, indices_subset, S_cand, k_cand, kg_cand, C_candidates, USE_DO = FALSE) {
  rc <- initialise_rc(indices_subset, C_candidates) # dataframe that will contain the optimized number of common factors for each C and subset
  rcj <- initialise_rcj(indices_subset, C_candidates) # dataframe that will contain the optimized number of groups and group specific factors for each C and subset

  for (subset in indices_subset) {
    print(paste("subset:", subset, "/", max(indices_subset)))
    temp <- make_subsamples(original_data, subset)
    Y <- temp[[1]]
    X <- temp[[2]]
    g_true <- temp[[3]]

    # define set of possible combinations of number of group factors
    configs <- define_configurations(S_cand, k_cand, kg_cand)

    # to add a progressbar:
    progress <- function(n) utils::setTxtProgressBar(utils::txtProgressBar(max = nrow(configs), style = 3), n)
    opts <- list(progress = progress)

    message("--maxit is set to 2 for test--")
    `%dopar%` <- foreach::`%dopar%` #to make dopar work within a function/package
    `%do%` <- foreach::`%do%` #for testing purposes
    if(USE_DO) {
      output <- foreach::foreach(
        i = 1:nrow(configs),
        .packages = c("RCTS", "tidyverse"),
        .options.snow = opts,
        .errorhandling = "pass"
      ) %do% {
        tryCatch(
          run_config(configs[i,], C_candidates, Y, X, maxit = 2), #might still be needing an extra trycatch around this function here?
          error = function(e) {
            return("moeilijk ding")
          }
        )
      }
    } else {
      output <- foreach::foreach(
        i = 1:nrow(configs),
        .packages = c("RCTS", "tidyverse"),
        .options.snow = opts,
        .errorhandling = "pass"
      ) %dopar% {
        run_config(configs[i,], C_candidates, Y, X, maxit = 2) #might still be needing an extra trycatch around this function here?
      }
    }

    df_results <- make_df_results_parallel(output)
    df_pic <- make_df_pic_parallel(output)


    pic_sigma2 <- df_results$sigma2[nrow(df_results)]
    df_pic <- adapt_pic_with_sigma2maxmodel(df_pic, df_results, pic_sigma2)


    # calculate for each candidate value for C the best S, k and kg
    all_best_values <- calculate_best_config(df_results, df_pic, C_candidates)
    rc <- fill_rc(rc, all_best_values, subset) # best number of common factors
    rcj <- fill_rcj(rcj, all_best_values, subset, S_cand, kg_cand) # best number of group specific factors and groups
    rm(all_best_values)

    # keep df_results of the full sample (this will contain the final clustering):
    if (subset == 0) df_results_full <- df_results
  } # closes loop over subsets
  return(list(rc, rcj, df_results_full))
}
