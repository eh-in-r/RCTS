
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
#' @inheritParams initialise_beta
#' @inheritParams calculate_PIC
#' @param maxit maximum limit for the number of iterations
#' @return list with the estimators and metrics for this configuration
run_config <- function(robust, config, C_candidates, Y, X, choice_pic, maxit = 30) {
  #print("-----------------------------------------start run_config:-------------------------------------------")
  #print(config)
  # print("-remove sleep again-")
  # Sys.sleep(2)
  S <- config %>%
    dplyr::select(S) %>%
    dplyr::pull()
  k <- config %>%
    dplyr::select(k) %>%
    dplyr::pull()
  kg <- unlist(config %>% dplyr::select(.data$k1:.data$k20)) # must be a vector (class "integer")


  iteration <- 0 # number of the iteration; 0 indicates being in the initialisation phase
  #print("initialise:")
  ########## initialisation
  beta_est <- initialise_beta(robust, Y, X, S)
  # initial grouping
  g <- initialise_clustering(robust, Y, S, k, kg, NA, max_percent_outliers_tkmeans = 0, verbose = FALSE)
  # initial common factorstructure
  temp <- initialise_commonfactorstructure_macropca(robust, Y, X, beta_est, g, NA, k, kg)
  comfactor <- temp[[1]]
  lambda <- temp[[2]]
  # initial group specific factorstructure
  factor_group <- estimate_factor_group(robust, Y, X, beta_est, g, NA, NA, NA, S, k, kg, initialise = TRUE)
  lambda_group <- calculate_lambda_group(robust, Y, X, beta_est, factor_group, g, NA, NA, S, k, kg, initialise = TRUE)


  #print("estimate:")
  ######### estimations
  obj_funct_values <- c()
  speed <- 999999 # convergence speed: set to initial high value
  while (iteration < maxit & !check_stopping_rules(iteration, speed, obj_funct_values, verbose = FALSE)) {

    temp <- iterate(robust, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor, S, k, kg, verbose = FALSE)
    #print("********************************************************(end of iterate()")
    #print(class(temp))
    beta_est <- temp[[1]]
    g <- temp[[2]]
    comfactor <- temp[[3]]
    lambda <- temp[[4]]
    factor_group <- temp[[5]]
    lambda_group <- temp[[6]]
    value <- temp[[7]]
    print(iteration)
    if(iteration == 0) print(paste(S, "-", k, "-", kg[1:S]),collapse = " ")
    #print(paste(c("   subset:", subset, "/", max(indices_subset), "config:", i, "/", "...", "(", S, "-", k, "-", kg[1:S], ")", "iteration", iteration, ":", round(value)), collapse = " "))
    obj_funct_values <- c(obj_funct_values, value)
    iteration <- iteration + 1
    speed <- get_convergence_speed(iteration, obj_funct_values / nrow(Y) / ncol(Y))
    print("--")
  }
  #print("estimation is done")

  # calculate the estimation errors
  est_errors <- calculate_error_term(
    Y, X, beta_est, g,
    factor_group, lambda_group,
    comfactor, lambda,
    S, k, kg
  )

  if(length(choice_pic) == 1) {
    pic <- add_pic_parallel(
      robust, Y, beta_est, g, S, k, kg,
      est_errors, C_candidates, choice_pic
    )
  }
  if(length(choice_pic) > 1) {
    pic = list()
    for( i in 1:length(choice_pic)) {
      pic[[i]] <- add_pic_parallel(
        robust, Y, beta_est, g, S, k, kg,
        est_errors, C_candidates, choice_pic[i]
      )
    }
  }

  pic_sigma2 <- calculate_sigma2(est_errors, nrow(Y), ncol(Y))
  #print("run_config is done")
  return(list(S, k, kg, pic, pic_sigma2, g, iteration, data.frame(beta_est), data.frame(comfactor), data.frame(lambda), factor_group, lambda_group))
}

#' Calculates the PIC for the current configuration.
#'
#' @inheritParams add_pic
#' @return numeric vector with a value for each candidate C
add_pic_parallel <- function(robust, Y, beta_est, g,
                             S, k, kg, est_errors, C_candidates, choice_pic,
                             method_estimate_beta = "individual") {
  if(!is.na(beta_est[1])) {
    vars_est <- ncol(beta_est)
  } else {
    vars_est <- 0
  }
  pic_sigma2 <- calculate_sigma2(est_errors)

  pic <- sapply(C_candidates, function(x) {
    calculate_PIC(x, robust, S, k, kg, est_errors, pic_sigma2,
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
#' @return data.frame
make_df_results_parallel <- function(x, limit_est_groups = 20) {
  df <- data.frame(
    S = unlist(x %>% purrr::map(1)),
    k_common = unlist(x %>% purrr::map(2))
  ) %>% cbind(t(matrix(unlist(x %>% purrr::map(3)), nrow = limit_est_groups))) %>%
    mutate(sigma2 = unlist(x %>% purrr::map(5)))

  names(df)[3:(3 + limit_est_groups - 1)] <- paste0("k", 1:limit_est_groups)
  temp <- x %>% purrr::map(6) #contains the group memberships
  df$g <- sapply(1:nrow(df), function(y) paste(temp[[y]], collapse = "-"))
  df$table_g <- sapply(1:nrow(df), function(y) paste(table(temp[[y]]), collapse = "_"))
  df$number_of_iterations <- unlist(x %>% purrr::map(7))
  df$beta_est <- x %>% purrr::map(8)
  df$comfactor <- x %>% purrr::map(9)
  df$lambda <- x %>% purrr::map(10)
  df$factor_group <- x %>% purrr::map(11)
  df$lambda_group <- x %>% purrr::map(12)
  return(df)
}

#' Makes a dataframe with the PIC for each configuration and each candidate C.
#'
#' @param x output of the parallel version of the algorithm
#' @param C_candidates candidates for C
#' @return data.frame
make_df_pic_parallel <- function(x, C_candidates) {
  df <- t(matrix(unlist(x %>% purrr::map(4)), nrow = length(C_candidates)))
  return(df)
}

#' Wrapper of the loop over the subsets which in turn use the parallelised algorithm.
#'
#' @param original_data list containing the original data (1: Y, 2: X)
#' @param indices_subset vector with indices of the subsets; starts with zero
#' @inheritParams get_best_configuration
#' @inheritParams calculate_VCsquared
#' @inheritParams initialise_beta
#' @inheritParams define_configurations
#' @param choice_pic indicates which PIC to use to estimate the number of groups and factors.
#' Options are "pic2017" (PIC of \insertCite{Ando2017;textual}{RCTS}; works better for large N),
#' "pic2016" (\insertCite{Ando2016;textual}{RCTS}; works better for large T) weighs the fourth term with an extra factor relative to the size of the groups,
#' and "pic2022" which shrinks the NT-space where the number of groups and factors would be over- or underestimated compared to pic2016 and pic2017. This is the default.
#' This parameter can also be a vector with multiple pic's.
#' @param maxit maximum limit for the number of iterations for each configuration; defaults to 30
#' @param USE_DO (for testing purposes) if TRUE, then a serialized version is performed ("do" instead of "dopar")
#' @return Returns a list with three elements.
#' 1. Data.frame with the optimal number of common factors for each candidate C in the rows.
#'    Each column contains the results of one subset of the input data (the first row corresponds to the full dataset).
#' 2. Data.frame with the optimal number of groups and group specific factors for each candidate C in the rows. The structure is the same as in the above.
#'    Each entry is of the form "1_2_3_NA". This is to be interpreted as 3 groups (three non NA values) where group 1 contains 1 group specific factor,
#'    group 2 contains 2 and group 3 contains 3.
#' 3. Data.frame with information about each configuration in the rows.
#' @examples
#' \donttest{
#' #Using a small dataset as an example; this will generate several warnings due to its small size.
#' #Note that this example is run sequentially instead of parallel,
#' #  and consequently will print some intermediate information in the console.
#' #This example uses the classical algorithm instead of the robust algorithm
#' #  to limit its running time.
#' set.seed(1)
#' original_data <- create_data_dgp2(30, 10)
#' #define the number of subsets used to estimate the optimal number of groups and factors
#' indices_subset <- define_number_subsets(2)
#' #define the candidate values for C (this is a parameter in the information criterium
#' #  used to estimate the optimal number of groups and factors)
#' C_candidates <- define_C_candidates()
#'
#' S_cand <- 3:3 # vector with candidate number of groups
#' k_cand <- 0:0 # vector with candidate number of common factors
#' kg_cand <- 1:2 # vector with candidate number of group specific factors
#'
#' #excluding parallel part from this example
#' #cl <- makeCluster(detectCores() - 1)
#' #registerDoSNOW(cl)
#' output <- parallel_algorithm(original_data, indices_subset, S_cand, k_cand, kg_cand,
#'   C_candidates, robust = FALSE, USE_DO = TRUE, maxit = 3)
#' #stopCluster(cl)
#' }
#' @export
parallel_algorithm <- function(original_data, indices_subset, S_cand, k_cand, kg_cand, C_candidates, robust = TRUE, USE_DO = FALSE, choice_pic = "pic2022", maxit = 30) {
  stopifnot(max(kg_cand) > 0)
  df_results_full <- NULL
  rc <- initialise_rc(indices_subset, C_candidates) # dataframe that will contain the optimized number of common factors for each C and subset
  rcj <- initialise_rcj(indices_subset, C_candidates) # dataframe that will contain the optimized number of groups and group specific factors for each C and subset
  if(length(choice_pic) > 1) {
    rc <- list()
    rcj <- list()
    for(i in 1:length(choice_pic)) {
      rc[[i]] <- initialise_rc(indices_subset, C_candidates)
      rcj[[i]] <- initialise_rcj(indices_subset, C_candidates)
    }
  }
  for (subset in indices_subset) {
    print(paste("subset:", subset, "/", max(indices_subset)))
    temp <- make_subsamples(original_data, subset)
    Y <- temp[[1]]
    X <- temp[[2]]
    g_true <- temp[[3]]

    # define set of possible combinations of number of group factors
    configs <- define_configurations(S_cand, k_cand, kg_cand)
    print(paste("There are", nrow(configs), "possible configurations."))

    # to add a progressbar:
    progress <- function(n) utils::setTxtProgressBar(utils::txtProgressBar(max = nrow(configs), style = 3), n)
    opts <- list(progress = progress)

    #message("--maxit is set to 2 for test--")
    `%dopar%` <- foreach::`%dopar%` #to make dopar work within a function/package
    `%do%` <- foreach::`%do%` #for testing purposes
    if(USE_DO) {
      output <- foreach::foreach(
        i = 1:nrow(configs),
        #.packages = c("RCTS", "tidyverse"), #apparently not needed anymore
        .options.snow = opts,
        .errorhandling = "pass"
      ) %do% {
          run_config(robust, configs[i,], C_candidates, Y, X, choice_pic, maxit)

      }
    } else {
      output <- foreach::foreach(
        i = 1:nrow(configs),
        #.packages = c("RCTS", "tidyverse"), #apparently not needed anymore
        .options.snow = opts,
        .errorhandling = "pass"
      ) %dopar% {
        run_config(robust, configs[i,], C_candidates, Y, X, choice_pic, maxit)
      }
    }
    #print("foreach has finished")
    config_groups_plus_errormessages <- purrr::map(output, 1)
    #print(config_groups_plus_errormessages)
    has_error <- unlist(lapply(purrr::map(output, class), function(x) "error" %in% x))
    if(sum(has_error) > 0) {
      message(paste("There are", sum(has_error), "configurations that returned an error."))
    }
    if(sum(has_error) == nrow(configs)) {
      message("1. All possible configurations have produced an error for this subset. Expanding the configurationspace is an option.")
    }
    #Note that these errors is due to trycatch statements not or not properly working within a foreach loop. (this should be solved now)
    #They are in fact solved in the serialized algorithm!
    #Note that these errors are often linked to one estimated group being small, or to many factors fitted to a group.
    #This occurs in configurations that have more groups and factors then the true values.
    #Omitting them from analysis should then have little impact on the overall estimation of the optimal configuration.
    for(i in sort(which(has_error == 1), decreasing = TRUE)) { #first the error with highest index is deleted, then the rest
      #message(i)
      to_print <- paste("Subset", subset, "has a problem with configuration", paste(configs[i,], collapse = " "),
                        "(", config_groups_plus_errormessages[[i]], ") -> no results available -> this configuration is omitted as a candidate.")
      #message(to_print)
      warning(to_print)
      output[[i]] <- NULL
    }
    #new vector showing which configurations have an error; should be all FALSE; can be NULL when all configurations failed
    has_error_new <- unlist(lapply(purrr::map(output, class), function(x) "error" %in% x))
    if (sum(has_error_new) == 0) {
      if (!is.null(has_error_new)) {
        df_results <- make_df_results_parallel(output)
        if(length(choice_pic) > 1) {
          pic_sigma2 <- df_results$sigma2[nrow(df_results)]
          df_pic <- list()
          for( i in 1:length(choice_pic) ) {
            #remaining issue here:
            #print(output %>% purrr::map(4))
            df_pic[[i]] <- t(matrix(unlist(output %>% purrr::map(4) %>% purrr::map(i)), nrow = length(C_candidates)))

            df_pic[[i]] <- adapt_pic_with_sigma2maxmodel(df_pic[[i]], df_results, pic_sigma2)

            # calculate for each candidate value for C the best S, k and kg
            all_best_values <- calculate_best_config(df_results, df_pic[[i]], C_candidates)
            rc[[i]] <- fill_rc(rc[[i]], all_best_values, subset) # best number of common factors
            rcj[[i]] <- fill_rcj(rcj[[i]], all_best_values, subset, S_cand, kg_cand) # best number of group specific factors and groups
            rm(all_best_values)
          }


        } else {
          df_pic <- make_df_pic_parallel(output, C_candidates)
          pic_sigma2 <- df_results$sigma2[nrow(df_results)]
          df_pic <- adapt_pic_with_sigma2maxmodel(df_pic, df_results, pic_sigma2)

          # calculate for each candidate value for C the best S, k and kg
          all_best_values <- calculate_best_config(df_results, df_pic, C_candidates)
          rc <- fill_rc(rc, all_best_values, subset) # best number of common factors
          rcj <- fill_rcj(rcj, all_best_values, subset, S_cand, kg_cand) # best number of group specific factors and groups
          rm(all_best_values)
        }



      } else {
        print(output)
        message("2. All possible configurations have produced an error for this subset. Expanding the configurationspace is an option.")
      }
    } else {
      message("Errors left! This should not happen!")
    }

    # keep df_results of the full sample (this will contain the final clustering):
    if (subset == 0) df_results_full <- df_results
  } # closes loop over subsets
  return(list(rc, rcj, df_results_full))
}
