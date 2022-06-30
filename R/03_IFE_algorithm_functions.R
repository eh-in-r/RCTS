

#' Function used in generating simulated data with non normal errors.
#'
#' Used to include cross-sectional dependence or serial dependence into the simulated panel data.
#' @param parameter amount of cross-sectional dependence
#' @param NN number of time series
#' @return NxN covariance matrix
create_covMat_crosssectional_dependence <- function(parameter, NN) {
  covMat <- matrix(NA, nrow = NN, ncol = NN)
  for (i in 1:NN) {
    for (j in 1:NN) {
      covMat[i, j] <- parameter^abs(i - j)
    }
  }
  return(covMat)
}

#' Helpfunction in create_true_beta() for the option beta_true_heterogeneous_groups. (This is the default option.)
#'
#' @param vars number of observable variables
#' @param S_true true number of groups: should be at least 1 and maximum limit_true_groups
#' @param extra_beta_factor option to multiply the coefficients in true beta; default = 1
#' @param limit_true_groups  Maximum number of true groups in a simulation-DGP for which the code in this package is implemented. Currently equals 12. For application on realworld data this parameter is not relevant.
#' @importFrom stats runif
#' @return matrix where the number of rows equals S_true, and the number of columns equals max(1, vars)
beta_true_heterogroups <- function(vars, S_true, extra_beta_factor = 1, limit_true_groups = 12) {
  stopifnot(S_true <= limit_true_groups) # Code allows up to 12 true groups at this point.
  stopifnot(S_true > 0)

  #######################################################################################
  # Define true values for beta, for each group, when there are 3 or less observable variables
  # These are the values for DGP 3 & 4 (For DGP 1 & 2 beta_true is defined in 08_IFE_create_data_dgp1.R)

  # 1st element is the intercept.
  beta_part1 <- c(0, 4, 3.5, 3, 2.5) * extra_beta_factor
  beta_part2 <- c(0, -2.5, -2, -2.5, -2) * extra_beta_factor
  beta_part3 <- c(0, 1, 0.5, 1.5, 1) * extra_beta_factor

  beta_define_further_true_values <- function(n) {
    # starting from group 4 we randomly generate values for the true values of beta
    return(c(0, (round(runif(n), 1) - 0.5) * 4) * extra_beta_factor) # between -2 and 2
  }
  beta_part4 <- beta_define_further_true_values(4)
  beta_part5 <- beta_define_further_true_values(4)
  beta_part6 <- beta_define_further_true_values(4)
  beta_part7 <- beta_define_further_true_values(4)
  beta_part8 <- beta_define_further_true_values(4)
  beta_part9 <- beta_define_further_true_values(4)
  beta_part10 <- beta_define_further_true_values(4)
  beta_part11 <- beta_define_further_true_values(4)
  beta_part12 <- beta_define_further_true_values(4)

  #######################################################################################
  # Define true values for beta_est, for each group, when there are more than 3 variables
  if (vars > 3) { # add when there are more than 3 observable variables:
    beta_part1 <- c(beta_part1, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part2 <- c(beta_part2, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part3 <- c(beta_part3, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part4 <- c(beta_part4, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part5 <- c(beta_part5, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part6 <- c(beta_part6, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part7 <- c(beta_part7, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part8 <- c(beta_part8, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part9 <- c(beta_part9, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part10 <- c(beta_part10, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part11 <- c(beta_part11, (round(runif(vars - 3), 1) - 0.5) * 4)
    beta_part12 <- c(beta_part12, (round(runif(vars - 3), 1) - 0.5) * 4)
  }

  # Add the values of the different groups together in one object.
  if (S_true >= 1) {
    beta_true <- matrix(beta_part1[1:(vars + 1)])
  }
  if (S_true >= 2) beta_true <- beta_true %>% cbind(beta_part2[1:(vars + 1)])
  if (S_true >= 3) beta_true <- beta_true %>% cbind(beta_part3[1:(vars + 1)])
  if (S_true >= 4) beta_true <- beta_true %>% cbind(beta_part4[1:(vars + 1)])
  if (S_true >= 5) beta_true <- beta_true %>% cbind(beta_part5[1:(vars + 1)])
  if (S_true >= 6) beta_true <- beta_true %>% cbind(beta_part6[1:(vars + 1)])
  if (S_true >= 7) beta_true <- beta_true %>% cbind(beta_part7[1:(vars + 1)])
  if (S_true >= 8) beta_true <- beta_true %>% cbind(beta_part8[1:(vars + 1)])
  if (S_true >= 9) beta_true <- beta_true %>% cbind(beta_part9[1:(vars + 1)])
  if (S_true >= 10) beta_true <- beta_true %>% cbind(beta_part10[1:(vars + 1)])
  if (S_true >= 11) beta_true <- beta_true %>% cbind(beta_part11[1:(vars + 1)])
  if (S_true >= 12) beta_true <- beta_true %>% cbind(beta_part12[1:(vars + 1)])


  return(beta_true)
}

#' Creates beta_true, which contains the true values of beta (= the coefficients of X)
#'
#' @param vars number of observable variables
#' @param NN number of time series
#' @param S_true number of groups
# @param use_real_world_data indicates using realworld data; defaults to FALSE
#' @param extra_beta_factor multiplies coefficients in beta_est; default = 1
#' @param method_true_beta how the true values of beta are defined: "homogeneous" (equal for all individuals),
#' "heterogeneous_groups" (equal within groups, and different between groups) or heterogeneous_individuals (different for all individuals)
#' @param limit_true_groups Maximum number of true groups in a simulation-DGP for which the code in this package is implemented. Currently equals 12. For application on realworld data this parameter is not relevant.
#' @return matrix with number of rows equal to number of observable variables + 1 (the first row contains the intercept) and number of culumns
#' equal to the true number of groups.
#' @importFrom stats rnorm
create_true_beta <- function(vars,
                             NN,
                             S_true,
                             method_true_beta = "heterogeneous_groups",
                             limit_true_groups = 12,
                             extra_beta_factor = 1) {
  stopifnot(method_true_beta %in% c("homogeneous", "heterogeneous_groups", "heterogeneous_individuals"))
  # real world data: beta_true does not exist -> return NA
  # if(use_real_world_data) return(NA)


  if (vars > 0) {
    # common beta_true:
    if (method_true_beta == "homogeneous") {
      beta_true <- c(0, c(1, 2, 3, seq(28, 83, 5))[1:vars]) # intercept 0, and after that values for the true beta's
      beta_true <- matrix(rep(beta_true, S_true), nrow = (vars + 1))
    }
    # group specific beta_true: -> default case
    if (method_true_beta == "heterogeneous_groups") {
      beta_true <- beta_true_heterogroups(
        vars, S_true,
        extra_beta_factor,
        limit_true_groups
      )
    }
    # individualspecific beta_true:
    if (method_true_beta == "heterogeneous_individuals") {
      beta_true <- matrix(rnorm(NN * (vars + 1)), nrow = (vars + 1), ncol = NN)
    }
  } else {
    beta_true <- rep(NA, S_true)
  }

  return(beta_true)
}


#' Creates X (the observable variables) to use in simulations.
#'
#' X is an array with dimensions N, T and number of observable variables.
#' The variables are randomly generated with mean 0 and sd 1.
#' @param NN number of time series
#' @param TT length of time series
#' @param vars number of available observable variables
#' @param scale_robust logical, defines if X will be scaled with robust metrics instead of with non-robust metrics
#' @return array with dimensions N x T x number of observable variables
initialise_X <- function(NN, TT, vars, scale_robust = TRUE) {
  if (vars > 0) {
    X <- array(0, dim = c(NN, TT, vars))
    for (i in 1:NN) {
      for (t in 1:TT) {
        for (k in 1:vars) {
          X[i, t, k] <- rnorm(1, mean = 0, sd = 1)
        }
      }
    }
    X <- scaling_X(X, TRUE, scale_robust, vars)


    return(X)
  } else {
    return(NA)
  }
}


#' Scaling of X.
#'
#' @param X input
#' @param firsttime Scaling before generating Y and before adding outliers: this is always with mean and sd. If this is FALSE, it indicates that
#' we are using the function for a second time, after adding the outliers. In the robust case it uses median and MAD, otherwise again mean and sd.
#' @inheritParams create_true_beta
#' @param robust logical, scaling with robust metrics instead of with non-robust measures
#' @inheritParams initialise_X
#' @return 3D-array with the same dimensions as X
#' @importFrom stats sd
scaling_X <- function(X, firsttime, robust, vars) {
  #
  # replaced this codeblock by a less concise but more clear codeblock
  #
  ##################
  # if(robust & !firsttime) {
  #   message("Scale with median and mad")
  #   med = median(X[,,k])
  #   mad = mad(X[,,k])
  #
  #   X[,,k] = (X[,,k] - med) / mad #Note that this is NOT column-based (timeindex) scaling!
  # } else {
  #   if(use_real_world_data) {
  #     message("Scale with mean and sd of NxT-matrix")
  #     X[,,k] = (X[,,k] - mean(X[,,k])) / sd(X[,,k]) #for eclipzdata, we cannot use scale(),
  #     #  since there are variables (for example age) that have constant columns, so scale() (which is columnbased) would produce errors
  #   } else {
  #     message("Scale with mean and sd (for each t separate)")
  #     X[,,k] = scale(X[,,k])      #Note that this is column-based (timeindex) scaling!
  #   }
  # }
  #################
  if (vars > 0) {
    for (k in 1:vars) {
      # print(paste("sd of variable",k,": Before:",sd(X[,,k])))
      if (mad(X[, , k]) != 0) {
        if (firsttime) {
          # if(use_real_world_data) {
          # #for eclipzdata, we cannot use scale(),
          #   #  since there are variables (for example age) that have constant columns, so scale() (which is columnbased) would produce errors
          #   #message("Scale with mean and sd of NxT-matrix")
          #   X[,,k] = (X[,,k] - mean(X[,,k])) / sd(X[,,k])
          # } else {
          # message("Scale with mean and sd (for each t separate)")
          X[, , k] <- scale(X[, , k]) # Note that this is column-based (timeindex) scaling!
          # }
        } else {
          warning("This part should be obsolete!")
          if (robust) {
            # message("Scale with median and mad")
            med <- median(X[, , k])
            mad <- mad(X[, , k])

            X[, , k] <- (X[, , k] - med) / mad # Note that this is NOT column-based (timeindex) scaling!
          } else {
            # if(use_real_world_data) {
            #   #for eclipzdata, we cannot use scale(),
            #   #  since there are variables (for example age) that have constant columns, so scale() (which is columnbased) would produce errors
            #   #message("Scale with mean and sd of NxT-matrix")
            #   X[,,k] = (X[,,k] - mean(X[,,k])) / sd(X[,,k])
            # } else {
            # message("Scale with mean and sd (for each t separate)")
            X[, , k] <- scale(X[, , k]) # Note that this is column-based (timeindex) scaling!
            # }
          }
        }

        # print(paste("sd of variable",k,": After:",sd(X[,,k])))
      }
    }
    stopifnot(!is.nan(X[1, 1, 1]))
  }
  return(X)
}

#' Restructures X (which is an 3D-array of dimensions (N,T,p) to a 2D-matrix of dimension (NxT,p).
#'
#' @param X input
#' @inheritParams create_true_beta
#' @inheritParams estimate_beta
#' @return The function returns a 2D-array, unless the input X is NA, in which case the output will be NA as well.
restructure_X_to_order_slowN_fastT <- function(X, vars_est) {
  if (!(is.na(X)[1])) {
    if (length(dim(X)) == 2) {
      # occurs when only one element in group
      X <- array(X, dim = c(1, nrow(X), ncol(X)))
    }

    # OLD
    # if(use_real_world_data) {
    #   number_of_vars = dim(X)[3]
    # } else {
    number_of_vars <- vars_est
    # }

    X_local <- matrix(NA, nrow = nrow(X) * ncol(X), ncol = number_of_vars)
    for (k in 1:number_of_vars) {
      X_local[, k] <- c(t(X[, , k]))
    }

    return(X_local)
  } else {
    # when there are no observable variables, X was NA and stays NA
    return(NA)
  }
}



#' Generates the true groupfactorstructure, to use in simulations.
#'
#' Loadings and factors are generated by:
#' factors ~ N(j * fgr_factor_mean, fgr_factor_sd) -> default case will be N(j, 1)
#' loadings ~ N(lgr_factor_mean, j * lgr_factor_sd) -> default case will be N(0, j)
#'
#' @param S true number of groups
#' @param kg_true vector with as length the number of groups, where each element is the true number of groupfactors of that group.
#' @inheritParams estimate_beta
#' @inheritParams generate_Y
#' @param fgr_factor_mean mean of the normal distribution from which the group specific factors are generated (multiplied by a coefficient for each different group)
#' @param fgr_factor_sd sd of the normal distribution from which the group specific factors are generated
#' @param lgr_factor_mean mean of the normal distribution from which the loadings are generated
#' @param lgr_factor_sd sd of the normal distribution from which the loadings are generated (multiplied by a coefficient for each different group)
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @return list: first element contains the true group specific factors and the second element contains the corresponding loadings
generate_grouped_factorstructure <- function(S, kg_true,
                                             TT, g_true,
                                             lgr_factor_mean = 0,
                                             lgr_factor_sd = 1,
                                             fgr_factor_mean = 1,
                                             fgr_factor_sd = 1) {
  if (mean(kg_true) > 0) {
    LGR <- list()
    factor_group_true <- list()
    for (i in 1:S) {
      LGR[[i]] <- matrix(nrow = sum(g_true == i), ncol = kg_true[i])
      factor_group_true[[i]] <- matrix(NA, nrow = kg_true[i], ncol = TT)
    }

    for (j in 1:S) {
      for (i in 1:kg_true[j]) {
        for (k in 1:(sum(g_true == j))) {
          LGR[[j]][k, i] <- rnorm(1, mean = lgr_factor_mean, sd = j * lgr_factor_sd)
        }
        for (k in 1:TT) {
          factor_group_true[[j]][i, k] <- rnorm(1, mean = j * fgr_factor_mean, sd = fgr_factor_sd)
        }
      }
    }
    # lambda_group_true to dataframe (easier to work with)
    lambda_group_true <- data.frame(LGR[[1]]) %>% mutate(group = 1, id = which(g_true == 1))
    if (kg_true[1] == 1) { # change name when only 1 factor (when there are more names are automatically X1,X2,...)
      names(lambda_group_true)[1] <- "X1"
    }
    for (i in 2:S) {
      temporary <- data.frame(LGR[[i]]) %>% mutate(group = i, id = which(g_true == i))
      if (kg_true[i] == 1) {
        names(temporary)[1] <- "X1"
      }
      if (ncol(lambda_group_true) >= ncol(temporary)) {
        lambda_group_true <- lambda_group_true %>% bind_rows(temporary)
      } else {
        lambda_group_true <- temporary %>% bind_rows(lambda_group_true)
      }
    }
    lambda_group_true <- lambda_group_true %>% arrange(group)
    lambda_group_true[is.na(lambda_group_true)] <- 0
  } else {
    factor_group_true <- NA
    lambda_group_true <- NA
  }

  # case of only 1 factor in group: format to matrix:
  for (group in 1:S) {
    if (kg_true[group] == 1) {
      factor_group_true[[group]] <- matrix(factor_group_true[[group]], nrow = kg_true[group])
    }
  }
  return(list(factor_group_true, lambda_group_true))
}



#' Generate panel data Y for simulations.
#'
#' @param NN number of time series
#' @param TT length of time series
#' @param k_true true number of common factors
#' @param kg_true Vector of length the number of groups. Each element contains the true number of group factors for that group.
#' @param g_true vector of length NN with true group memberships
#' @param beta_true true coefficients of the observable variables
#' @param factor_group_true true group specific factors
#' @param lambda_group_true loadings of the true group specific factors
#' @param comfactor_true true common factors
#' @param lambda_true loadings of the true common factors
#' @param eps NN x TT-matrix containing the error term
#' @param X dataframe with the observed variables
#' @return NN x TT matrix
generate_Y <- function(NN, TT, k_true, kg_true,
                       g_true, beta_true, lambda_group_true, factor_group_true,
                       lambda_true, comfactor_true, eps, X) {

  # Define the size of the panel data:
  Y <- matrix(NA, nrow = NN, ncol = TT) # initialisation, later on this gets filled in

  for (i in 1:NN) {

    # if(!use_real_world_data) {
    if (mean(kg_true) > 0) {
      dropvars <- names(lambda_group_true) %in% c("group", "id")
      LAMBDAGROUP <- as.matrix(subset(lambda_group_true, lambda_group_true$id == i)[!dropvars])
      LAMBDAGROUP <- LAMBDAGROUP[1:kg_true[g_true[i]]]
    } else {
      LAMBDAGROUP <- NA
    }
    # }

    if (is.na(X)[1]) {
      vars <- 0
    } else {
      vars <- dim(X)[3]
    }
    for (t in 1:TT) {
      if (vars > 0) {
        XT <- c(1, X[i, t, ]) %*% beta_true[, g_true[i]] # add 1 to X[i,t,] to make room for the intercept (which is in beta_est)
      } else {
        XT <- 0
      }


      # randomly generated data
      if (k_true > 0) {
        LF <- t(lambda_true[, i]) %*% comfactor_true[, t]
      } else {
        LF <- 0
      }
      if (mean(kg_true) > 0) {
        LF_GROUP <- LAMBDAGROUP %*% factor_group_true[[g_true[i]]][, t]
      } else {
        LF_GROUP <- 0
      }

      Y[i, t] <- XT +
        LF +
        LF_GROUP +
        eps[i, t]
      stopifnot(!is.nan(Y[i, t]))
      # }
    }
  }

  return(Y)
}


#' Initialisation of estimation of beta (the coefficients with the observable variables)
#'
#' Note: this needs to be called before the definition of grid.
#' @inheritParams generate_Y
#' @param robust robust or classical estimation
#' @param S estimated number of groups
#' @param g vector with estimated group membership for all individuals
#' @inheritParams define_object_for_initial_clustering_macropca
#' @param nosetting_lmrob option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
#' @param special_case_dgp1 special case for data generated according to dgp 1: it changes the 1st variable in X to 1 (-> intercept). Consequently the estimation of beta needs to be restructured slightly.
#' @return Matrix with number of rows equal to the number of estimated variables plus one. If method_estimate_beta is set to the default ("individual"),
#' the number of columns is equal to the number of time series in Y. If method_estimate_beta is set to "group" or to "homogeneous" the number of columns
#' is equal to the number of groups.
#' @examples
#' \donttest{
#' X <- X_dgp3
#' Y <- Y_dgp3
#' # Set estimations for group factors and its loadings, and group membership
#' #  to the true value for this example.
#' lambda_group <- lambda_group_true_dgp3
#' factor_group <- factor_group_true_dgp3
#' g <- g_true_dgp3
#' beta_init <- initialise_beta(TRUE, Y, X,
#'   S = 3, g
#' )
#' }
#' @importFrom stats lm
#' @export
initialise_beta <- function(robust, Y, X,
                            S, g,
                            method_estimate_beta = "individual",
                            nosetting_lmrob = FALSE,
                            special_case_dgp1 = FALSE) {
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  NN <- nrow(Y)
  TT <- ncol(Y)
  # if(use_real_world_data) {
  #   number_of_vars = vars
  # } else {
  number_of_vars <- vars_est
  # }

  if (number_of_vars > 0) {
    beta_est <- matrix(NA, nrow = (number_of_vars + 1), ncol = S)
    if (method_estimate_beta == "individual") {
      beta_est <- matrix(NA, nrow = (number_of_vars + 1), ncol = NN)
    }


    # X needs to be in the form of (NN*TT x p matrix)

    if (method_estimate_beta == "homogeneous") { # also used in DGP05: BramatiCroux
      X_special <- restructure_X_to_order_slowN_fastT(
        X,
        vars_est
      )
      Y_special <- Y
      # this includes robust estimation of beta_est:
      beta_est <- determine_beta("homogeneous", X_special, Y_special, robust, NN, TT, S, method_estimate_beta,
        initialisation = TRUE,
        indices = 1:NN, vars_est, NA, nosetting_local = nosetting_lmrob,
        special_case_dgp1
      )
    }
    if (method_estimate_beta == "group") {
      for (group in 1:S) { # for each group there must be a column in beta_est.
        # select parts of X and Y of this group
        indices_group <- which(g == group)
        if (length(indices_group) == 0) {
          message("There is an empty group!")
        }

        # X needs to be in the form of (NN*TT x p matrix)
        X_special <- restructure_X_to_order_slowN_fastT(
          array(X[indices_group, , ],
            dim = c(length(indices_group), TT, number_of_vars)
          ),
          vars_est
        )

        Y_special <- as.vector(t(Y[indices_group, ])) # order: N1T1, N1T2,N1T3,...N2T1,...N_endT_end

        beta_est[, group] <- determine_beta("heterogeneous", X_special, Y_special, robust, NN, TT, S, method_estimate_beta,
          initialisation = TRUE,
          indices = indices_group, vars_est, NA, nosetting_local = nosetting_lmrob,
          special_case_dgp1
        )
      }
    }
    if (method_estimate_beta == "individual") {


      # Initialisation in classical case: use of lm instead of ncvreg (as in AndoBai-code (website))
      for (i in 1:NN) {
        X_special <- restructure_X_to_order_slowN_fastT(
          matrix(X[i, , ], ncol = dim(X)[3]),
          vars_est
        )
        Y_special <- as.vector(t(Y[i, ]))

        #####################
        # define model
        #####################
        if (robust) { # robust -> lmrob
          # if(exists("use_bramaticroux")) {
          #   message(paste("bramati croux init",i))
          #   model = RpanFE(Y_special, X_special, TT, 0.20, 20, vars, 1)[[1]]
          #
          # } else {
          if (special_case_dgp1) {
            # for special case of dgp 1 (= ando/bai dgp 2)
            # no intercept, because dgp1_spread_group_centers defines the first variable in X as an intercept
            model <- LMROB(Y_special, X_special, nointercept = TRUE, nosetting = nosetting_lmrob) #-> lmrob(Y_special ~ X_special + 0, setting="KS2014)
          } else {
            model <- LMROB(Y_special, X_special, nosetting = nosetting_lmrob)
          }

          # }
        } else { # classical -> lm
          if (special_case_dgp1) {
            model <- lm(Y_special ~ X_special + 0) # for special case of dgp 1 (= ando/bai dgp 2)
          } else {
            model <- lm(Y_special ~ X_special) # for the rest
          }
        }
        #####################
        # get the coefficients
        #####################

        if (special_case_dgp1) {
          beta_est[, i] <- c(0, model$coefficients) # for dgp 1
        } else {
          beta_est[, i] <- model$coefficients # for dgp 2
        }
      }
    }
    rm(X_special, Y_special)


    return(beta_est)
  } else {
    # when vars == 0, we do not need beta_est
    return(NA)
  }
}


#' Helpfunction used in update_g()
#'
#' This function calculates FgLg (the groupfactorstructure) for all possible groups where individual i can be placed. For each group were
#' the groupfactors (Fg) estimated earlier. Now the grouploadings are needed for each group as well. In the classical case these are calculated by Fg*Y/T. In the robust case
#' these are robust.
#' @param group number of groups
#' @param solve_FG_FG_times_FG This is the same as groupfactor / T. It is only used in the Classical approach
#' @param method_estimate_factors_local specifies the robust algorithm to estimate factors: default is "macro"
#' @param g vector with estimated group membership for all individuals
#' @param NN_local number of time series
#' @param vars_est number of variables that are included in the algorithm and have their coefficient estimated. This is usually equal to vars.
#' @param number_of_group_factors_local number of group factors to be estimated
#' @param number_of_common_factors_local number of common factors to be estimated
#' @inheritParams update_g
#' @return NxT matrix containing the product of virtual groupfactors and virtual loadings
calculate_virtual_factor_and_lambda_group <- function(group, solve_FG_FG_times_FG, robust,
                                                      NN_local,
                                                      method_estimate_factors_local,
                                                      g,
                                                      vars_est,
                                                      number_of_group_factors_local,
                                                      number_of_common_factors_local,
                                                      method_estimate_beta,
                                                      factor_group, lambda, comfactor, Y, X, beta_est,
                                                      # use_real_world_data_local = use_real_world_data,
                                                      verbose = FALSE) {
  FG <- factor_group[[group]]
  indices <- 1:NN_local
  LF <- t(lambda) %*% comfactor
  xbeta <- calculate_XB_estimated(
    X, beta_est, g,
    vars_est, method_estimate_beta, ncol(Y)
  )

  if (!do_we_estimate_common_factors(number_of_common_factors_local)) {
    Y_ster <- Y[indices, ] - xbeta
  } else {
    Y_ster <- Y[indices, ] - xbeta - LF
  }

  # robust grouplambda:
  if (robust) {
    # CHANGE LOCATION 5/7
    if (method_estimate_factors_local %in% c("macro", "pertmm")) {
      # we need a robust version of the virtual factorstructure:
      LG_local <- return_robust_lambdaobject(Y_ster, group,
        type = 1, g, NN_local,
        number_of_common_factors_local,
        number_of_group_factors_local,
        comfactor_rrn = comfactor,
        factor_group_rrn = factor_group,
        verbose
        # use_real_world_data_rrn = use_real_world_data_local,
      )
    } else {
      LG_local <- t(solve_FG_FG_times_FG[[group]] %*% t(Y_ster)) # This equalS to Fg*Y/T
    }
  } else {
    LG_local <- t(solve_FG_FG_times_FG[[group]] %*% t(Y_ster)) # This equalS to Fg*Y/T
  }


  LG_local <- handleNA_LG(LG_local)
  return(LG_local %*% FG)
}





#' Determines parameters of rho-function.
#'
#' Robust updating of group membership is based on a rho function (instead of the non-robust quadratic function) on the norm of the errors.
#' This requires parameters of location and scale.
#' They are defined here (currently as median and madn).
#' This function is applied on the estimated errors: Y - XB - FL - FgLg.
#' This function is used in update_g().
#' @param object input
#' @importFrom stats mad
#' @importFrom stats median
#' @return list
define_rho_parameters <- function(object = NULL) {
  if (is.null(object)) {
    stop("This part should be obsolete (define_rho_parameters()).")
  } else {
    rho_loc <- apply(data.frame(object), 1, function(x) median(x)) # =median over T elements
    rho_scale <- apply(data.frame(object), 1, function(x) mad(x)) # =normalized mad over T elements


    # print(median(unlist(lapply(list(rho_loc,rho_scale),function(x) x[[1]][1]))))
  }

  return(list(rho_loc, rho_scale))
}

#' Helpfunction for update_g(). Calculates the errors for one of the possible groups time series can be placed in.
#'
#' During the updating of group membership, the errorterm is used as the objective function to estimate the group.
#' @param group group
#' @param LF NxT-matrix of the common factorstructure
#' @param virtual_grouped_factor_structure list with length the number of groups; every element of the list contains NxT-matrix
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @return NxT matrix with the errorterms (=Y minus the estimated factorstructure(s) and minus X*beta)
calculate_errors_virtual_groups <- function(group, LF, virtual_grouped_factor_structure,
                                            NN, TT,
                                            k, kg,
                                            vars_est,
                                            method_estimate_beta,
                                            Y, X, beta_est, g) {
  E_prep <- matrix(NA, nrow = NN, ncol = TT)

  a <- do_we_estimate_common_factors(k)
  b <- do_we_estimate_group_factors(kg)
  if( !is.na(X[1]) & !is.null(X[1]) ) {
    X <- adapt_X_estimating_less_variables(X, vars_est)
  }

  for (i in 1:NN) {
    # calculate lambda_group for one individual, based on an hypothetical group membership
    if (b == 1) { #-> we estimate group factors
      virtual_structure <- virtual_grouped_factor_structure[[group]][i, ]
    } else {
      message("Option to not estimate any group factor in any group -> this option is not implemented")
    }
    if (vars_est > 0) {
      if ((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) {
        XT <- cbind(1, X[i, , ]) %*% beta_est[, g[i]] # matrix with TT rows and 1 column
        #-> this should not be beta_est[,group] as only the grouped factorstructure should vary with group (similar to calculate_virtual_factor_and_lambda_group)
      }
      if (method_estimate_beta == "individual") {
        XT <- cbind(1, X[i, , ]) %*% beta_est[, i] # matrix with TT rows and 1 column
      }
    } else {
      XT <- matrix(0, nrow = TT, ncol = 1)
    }

    # calculate the objective function, while making sure that NA's do not have any effect in the sum
    E_prep[i, ] <- Y[i, ] - a * LF[i, ]
    if (b != 0) E_prep[i, ] <- (E_prep[i, ] - b * virtual_structure)
    # print(dim(XT))
    # print(length(E_prep[i,]))
    if (vars_est > 0) E_prep[i, ] <- E_prep[i, ] - XT # note: XT is an TTx1 matrix
  }
  return(E_prep)
}

#' Calculates objective function for individual i and group k in order to estimate group membership.
#'
#' Helpfunction in update_g().
#' Depends on an not yet established group k ( cannot use lgfg_list)
#' @param i individual
#' @param k group
#' @param errors_virtual list with errors for each possible group
#' @param rho_parameters median and madn of the calculated error term
#' @inheritParams initialise_beta
#' @param TT length of time series
#' @return numeric value
calculate_obj_for_g <- function(i, k, errors_virtual, rho_parameters, robust, TT) {
  totalsum <- 0

  if (robust) {
    # define a scaling:
    # It must be individualspecific, and also be equal over all virtual groups;
    # We take here the median over the values of all groups.
    #     (because: use of #rho_parameters[[k]][[2]][i] leads to random grouping)
    # "location" is then defined the same way
    # rho_parameters is a list of 'number_of_groups' elements. Every element has 2 elements with 'NN' values.
    location <- (unlist(lapply(rho_parameters, function(x) x[[1]][i]))) # map over virtual groups; take 1st element (=median) and take individual i
    scaling <- (unlist(lapply(rho_parameters, function(x) x[[2]][i]))) # map over virtual groups; take 2nd element (=mad) and take individual i


    location <- median(location) # median over groups
    scaling <- median(scaling) # median over groups
    if (scaling == 0) { # This should be rare.
      warning("scaling equals 0 and is changed to 0.000001, to evade division by zero")
      scaling <- 0.000001 # make sure no 0/0 will occur
    }
  }

  for (t in 1:TT) {

    # calculate the objective function, while making sure that NA's do not have any effect in the sum

    # calculate the estimation error:
    E_prep <- errors_virtual[[k]][i, t]
    if (robust) { # define the rho-function (bisquare)

      E_prep_loc_scale <- (E_prep - location) / scaling # this is a scalar
      E <- Mpsi(E_prep_loc_scale, cc = 4.685, psi = "bisquare", deriv = -1) # rho-functie on scaled errors
    } else { # classical version:
      E <- E_prep^2
    }


    E <- evade_floating_point_errors(E) # evading floating point errors -> set to zero when values are really small

    totalsum <- totalsum + as.numeric(E)
  }
  return(totalsum)
}

#' Helpfunction in update_g(), to calculate solve(FG x t(FG)) x FG
#'
#' @param TT length of time series
#' @param S number of groups
#' @param kg vector with the estimated number of group specific factors for each group
#' @param testing variable that determines if we are in 'testing phase'; defaults to FALSE (requires Matrix-package if set to TRUE)
#' @inheritParams estimate_beta
#' @return list: the number of elements in this list is equal to S (the number of groups).
#' Each of the elements in this list has a number rows equal to the number of group specific factors, and TT columns.
solveFG <- function(TT, S, kg, factor_group, testing = FALSE) {
  solve_FG_FG_times_FG <- list()
  for (group in 1:S) {
    if (kg[group] > 0) {
      FG <- factor_group[[group]]
      if(testing) {
        if (as.numeric(Matrix::rankMatrix(FG)) != kg[group]) {
          warning("There is an issue in solveFG(): rank of FG should be the same as the number of group factors")
          message(paste("rank of matrix FG: ", as.numeric(Matrix::rankMatrix(FG))))
          message(paste("number_of_group_factors[group]: ", kg[group]))
          print(FG[, 1:3])
        }
      }
      solve_FG_FG_times_FG[[group]] <- solve(FG %*% t(FG)) %*% FG # this is actually the same as FG/T
      rm(FG)
    } else {
      # make 0-matrix with 1 row
      solve_FG_FG_times_FG[[group]] <- matrix(0, nrow = 1, ncol = TT)
    }
  }
  return(solve_FG_FG_times_FG)
}

#' Function that estimates group membership.
#'
#' @inheritParams estimate_beta
# @param use_class_zero if set to TRUE, then individuals with high distance to all possible groups are put in a separate class zero
# @param use_real_world_data_inupdateg Parameter to indicate using real world dataset. Defaults to FALSE.
#' @inheritParams estimate_factor
#' @param robust robust or classical estimation of group membership
#' @param verbose when TRUE, it prints messages
#' @return Returns a list. The first element contains a vector with the estimated group membership for all time series. The second element contains the values
#' which were used to determine the group membership. The third element is only relevant if method_estimate_factors is set to "cz" (non-default) and
#' contains the group membership before moving some of the time series to class zero.
#' @examples
#' \donttest{
#' X <- X_dgp3
#' Y <- Y_dgp3
#' # Set estimations for group factors and its loadings, and group membership to the true value
#' lambda_group <- lambda_group_true_dgp3
#' factor_group <- factor_group_true_dgp3
#' g_true <- g_true_dgp3 # true values of group membership
#' g <- g_true # estimated values of group membership; set in this example to be equal to true values
#' # There are no common factors to be estimated  ->  use placeholder with values set to zero
#' lambda <- matrix(0, nrow = 1, ncol = 300)
#' comfactor <- matrix(0, nrow = 1, ncol = 30)
#' # Choose how coefficients of the observable are estimated
#' beta_est <- estimate_beta(
#'   robust = TRUE, Y, X, NA, g, lambda_group, factor_group,
#'   lambda, comfactor,
#'   S = 3, k = 0, kg = c(3, 3, 3),
#'   vars_est = 3
#' )[[1]]
#' g_new <- update_g(
#'   robust = TRUE, Y, X, beta_est, g,
#'   factor_group, lambda, comfactor,
#'   S = 3,
#'   k = 0,
#'   kg = c(3, 3, 3),
#'   vars_est = 3,
#'   "macro", "individual"
#' )[[1]]
#' }
#' @export
update_g <- function(robust, Y, X, beta_est, g,
                     factor_group,
                     lambda, comfactor,
                     S, k, kg,
                     vars_est,
                     method_estimate_factors, method_estimate_beta,
                     verbose = FALSE) {
  NN <- nrow(Y)
  TT <- ncol(Y)
  if (do_we_estimate_group_factors(kg)) { # if there are groupfactors estimated
    solve_FG_FG_times_FG <- solveFG(TT, S, kg, factor_group)

    # Calculate FgLg (groupfactors times grouploadings) for all the possible groups in which individual i could end up:

    virtual_grouped_factor_structure <- lapply(1:S, function(y) {
      calculate_virtual_factor_and_lambda_group(
        y, solve_FG_FG_times_FG, robust, NN, method_estimate_factors, g,
        vars_est,
        kg,
        k,
        method_estimate_beta,
        factor_group, lambda, comfactor, Y, X, beta_est
      )
    })
    # if (verbose) message("virtual_grouped_factor_structure is created")
  } else {
    virtual_grouped_factor_structure <- NA
  }


  LF <- (t(lambda) %*% comfactor)

  # calculate errors for each possible group
  errors_virtual <- lapply(1:S, function(x) {
    calculate_errors_virtual_groups(
      x, LF, virtual_grouped_factor_structure, NN, TT,
      k,
      kg,
      vars_est,
      method_estimate_beta, Y, X, beta_est, g
    )
  })
  # if (verbose) message("errors_virtual is created")

  if (robust) {
    rho_parameters <- lapply(1:S, function(x) define_rho_parameters(errors_virtual[[x]])) # (parameter object = NA -> returns median and madn of the calculated error term)
  } else {
    rho_parameters <- NA
  }
  # if (verbose) message("rho_parameters is created")

  # init matrix with objectivefunctionvalues for all groups
  matrix_obj_values <- matrix(NA, nrow = NN, ncol = S)
  for (i in 1:NN) {
    obj_values <- sapply(1:S, function(x) calculate_obj_for_g(i, x, errors_virtual, rho_parameters, robust, TT))
    g[i] <- which.min(obj_values)
    matrix_obj_values[i, ] <- obj_values
  } # Note: vectorizing is not faster: g = sapply(1:NN, function(z) which.min(sapply(1:S, function(x) calculate_obj_for_g(z, x, errors_virtual, rho_parameters, robust, TT))))



  g_before_class_zero <- g
  if (method_estimate_factors == "cz") {
    g <- clustering_with_robust_distances(g, S, Y)
  }

  g <- reassign_if_empty_groups(g, S, NN)
  return(list(g, matrix_obj_values, g_before_class_zero))
}

#' Randomly reassign individual(s) if there are empty groups. This can happen if the total number of time series is low compared to the number of desired groups.
#'
#' @param g Vector with group membership for all individuals
#' @param S_true true number of groups
#' @param TT number of time series
#' @return numeric vector with the estimated group membership for all time series
reassign_if_empty_groups <- function(g, S_true, NN) {
  empty_groups <- c()
  while (length(table(g[g != 0])) != S_true) { # class zero (individuals too far from all groups) should not be counted here
    # reassign one individual to empty groups
    for (i in 1:S_true) {
      if (i %in% tibble::as_tibble(table(g))$g) { # if(i %in% broom::tidy(table(g))$g) {
      } else {
        empty_groups <- c(empty_groups, i)
      }
    }

    for (i in 1:length(empty_groups)) {
      randomindex <- ceiling(runif(1) * NN)
      g[randomindex] <- empty_groups[i]
    }
  }
  return(g)
}


#' Function that puts individuals in a separate "class zero", when their distance to all possible groups is bigger then a certain threshold.
#'
#' It starts with defining a robust location and scatter (based on Ma & Genton (2000): Highly robust estimation of the autocovariance function).
#'
#' @param g Vector with group membership for all individuals
#' @param number_of_groups number of groups
#' @param Y Y: the panel data of interest
# @param percent_outliers (only used in plot if it is > 0)
#' @return numeric vector with the new clustering, now including class zero adjustments
#' @importFrom stats qchisq
#' @importFrom stats quantile
# @importFrom tsqn corQn
# @importFrom Matrix rankMatrix
clustering_with_robust_distances <- function(g, number_of_groups, Y) {
  RD <- matrix(NA, nrow = number_of_groups, ncol = nrow(Y)) # every row contains the distances to a certain group, and every column is 1 individual
  for (group in 1:number_of_groups) { # loop over all groups
    mu <- apply(Y[g == group, ], 2, "median")
    s <- apply(Y[g == group, ], 2, "Qn")
    x <- Y[g == group, 2:ncol(Y)]
    laggedx <- Y[g == group, 1:(ncol(Y) - 1)]
    r <- tsqn::corQn(as.vector(x), as.vector(laggedx)) # Computes the robust correlation of x and y

    ar1_cor <- function(n, rho) {
      exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
        (1:n - 1))
      rho^exponent
    }

    if (Matrix::rankMatrix(diag(s))[[1]] < ncol(Y)) {
      message("rank of diag(s) is too small (reason: because Qn(.) equals zero for 1 year)") #-> solve(sig) would generate error
      s[s == 0] <- 1e-6
    }

    R <- ar1_cor(n = ncol(Y), rho = r)
    sig <- diag(s) %*% R %*% diag(s)
    sig_inv <- solve(sig)
    # calculate robust distances for all individuals
    for (i in 1:nrow(Y)) {
      RD[group, i] <- sqrt(t(Y[i, ] - mu) %*% sig_inv %*% (Y[i, ] - mu))
    }
  }

  # define limit as 99%-quantile of chi-squared distribution:
  # if for an individual the robust distance for every group is larger then the limit, the individual is put to class zero
  limit <- sqrt(qchisq(0.99, ncol(Y), ncp = 0))

  minimum_distance <- apply(RD, 2, min)
  limit_empirical <- quantile(minimum_distance, 0.75) #
  cases <- which(minimum_distance > max(limit, limit_empirical)) # those cases go to class zero


  # if(percent_outliers > 0) {
  #   #plot with colored outlierindividuals;
  #   indices_of_outliers = which(apply(Y, 1, function(x) abs(max(x))) > 900) #indices of individuals with at least 1 generated outlier in it
  #   plot(log(apply(RD, 2, min)), col = ((1:nrow(Y)) %in% indices_of_outliers) + 1, main = "minimal log(distance) of i to any group")
  #   abline(h = log(limit), col = "red")
  #   abline(h = log(limit_empirical),col = "orange")
  # } else {
  #   #plot without coloring
  #   plot(apply(RD, 2, min), main = "minimal distance of i to any group", ylim = c(0, limit * 1.05))
  #   abline(h = limit, col = "red")
  #   abline(h = limit_empirical, col = "orange")
  # }
  g[cases] <- 0
  #print(table(g))
  return(g)
}






#' Helpfunction in OF_vectorized3()
#'
#' @param i index of individual
#' @param t index of time
#' @param XBETA matrixproduct of X and beta_est
#' @param LF matrixproduct of common factors and its loadings
#' @param group_memberships Vector with group membership for all individuals
#' @param lgfg_list product of groupfactors and their loadings; list with length the number of groups
#' @inheritParams estimate_beta
#' @param kg vector containing the number of group factors to be estimated for all groups
#' @return numeric: contains the contribution to the objective function of one timepoint for one time series
OF_vectorized_helpfunction3 <- function(i, t, XBETA, LF,
                                        group_memberships,
                                        lgfg_list, Y,
                                        kg) {
  if (do_we_estimate_group_factors(kg) != 0 & group_memberships[i] != 0) {
    if (t > ncol(lgfg_list[[group_memberships[i]]])) { # this is the case when macropca() dropped columns
      warning("If this is seen: there is a unsolved issue: macropca has dropped columns (OF_vectorized_helpfunction3()) for unknown reasons")
      print(t)
      print(dim(lgfg_list[[group_memberships[i]]]))
    } else {
      result <- as.numeric(Y[i, t] - XBETA - LF -
        lgfg_list[[group_memberships[i]]][i, t])^2
    }
  } else {
    result <- as.numeric(Y[i, t] - XBETA -
      LF)^2
  }

  return(result)
}

#' Calculates objective function for the classical algorithm: used in iterate() and in local_search.
#'
#' @param NN number of time series
#' @param TT length of time series
#' @param g Vector with group membership for all individuals
#' @param grid dataframe containing the matrix multiplications XB, FgLg and FL
#' @param beta_est estimated values of beta
#' @param fc estimated common factors
#' @param lc loadings of estimated common factors
#' @param fg estimated groupfactors
#' @param lg estimated grouploadings
#' @inheritParams estimate_beta
#' @importFrom magrittr %>%
#' @return numeric value of the objective function
OF_vectorized3 <- function(NN,
                           TT, g, grid,
                           Y,
                           beta_est,
                           lc, fc,
                           lg, fg,
                           S, k, kg,
                           method_estimate_beta,
                           num_factors_may_vary = TRUE) {
  # this is a list (with as length the number of groups) of the product FgLg (which is the groupfactorstructure)
  lgfg_list <- calculate_lgfg(lg, fg, S, k, kg, num_factors_may_vary, NN, TT)

  if ((method_estimate_beta == "homogeneous") | (method_estimate_beta == "individual")) {
    return(sum(apply(grid, 1, function(x) OF_vectorized_helpfunction3(x[1], x[2], x[3], x[4], g, lgfg_list, Y, kg))))
  }
  if (method_estimate_beta == "group") {
    # construct a vector with XBETA-values depending on the group of the individuals:

    #limit package dependency: note that base R returns a vector if there is only 1 factor, and the colnames also need to be assigned again
    #temp <- grid %>% dplyr::select(starts_with("XBETA"))
    temp <- as.data.frame(grid[, grep("XBETA", names(grid))])
    names(temp) <- paste0("XBETA", 1:ncol(temp))


    XBETA_parameter <- sapply(1:NN, function(x) temp[x, g[x]])

    #prep <- grid %>% dplyr::select(-starts_with("XBETA"))
    prep <- as.data.frame(grid[, -grep("XBETA", names(grid))])

    # used to put LF into the 3rd column -> x[3] in next line
    return(sum(apply(prep, 1, function(x) OF_vectorized_helpfunction3(x[1], x[2], XBETA_parameter, x[3], g, lgfg_list, Y, kg))))
  }
}






#' Calculates the group factor structure: the matrix product of the group factors and their loadings.
#'
#' Returns list (with as length the number of groups) with lgfg (product of grouploadings and groupfactors).
#' Each element of the list with the assumption that all individuals are in the same group k.
#' This function is used to speed up code.
#' @inheritParams estimate_beta
#' @param NN number of time series
#' @param TT length of time series
#' @param S number of groups
#' @param k number of common factors
#' @param kg vector with the number of group specific factors for each group
#' @return list with S elements: each element contains a matrix with NN rows and TT columns with the estimated group factor structure of this particular group
#' @importFrom dplyr arrange
#' @importFrom rlang .data
calculate_lgfg <- function(lambda_group, factor_group, S, k, kg, num_factors_may_vary,
                           NN, TT) {
  lgfg_list <- list()
  # define LgFg for each group (and later on select the correct element (correct group of individual i) of lgfg_list)
  for (gr in 1:S) {
    if (kg[gr] > 0) {
      LG_clean <- (as.matrix(lambda_group %>% arrange(.data$id) %>% dplyr::select(-.data$group, -.data$id)))[, 1:kg[gr]]

      # When using a varying number of groupfactors per group,  and also estimating a positive number of common factors, lgfg_list will contain NA's and crash the algorithm.
      # Reason is that there are NA's in lambda_group (only for the groups that have not the maximum amount of groupfactors), which is fine in itself.
      #-> Since those NA's mean that there is no effect on LGFG (since there is no factorpart and lambdapart for those), we can replace those NA's by zero's.
      if (num_factors_may_vary & k > 0) {
        if (anyNA(LG_clean)) {
          temp <- LG_clean
          temp[is.na(temp)] <- 0
          LG_clean <- temp
        }
      }


      lgfg_list[[gr]] <- LG_clean %*% factor_group[[gr]]
    } else {
      lgfg_list[[gr]] <- NA
    }
  }


  # replace parts of lgfg_list which are NA (because of no groupfactors in that particular group) by 0-matrices
  emptyFL <- sapply(lgfg_list, function(x) is.null(dim(x)))
  for (q in 1:length(lgfg_list)) {
    if (emptyFL[q]) {
      lgfg_list[[q]] <- matrix(0, nrow = NN, ncol = TT)
    }
  }


  return(lgfg_list)
}


#' Helpfunction in estimate_beta() for estimating beta_est.
#'
#' @param string can have values: "homogeneous" (when one beta_est is estimated for all individuals together) or "heterogeneous" (when beta_est is estimated either groupwise or elementwise)
#' @param X_special preprocessed X (2-dimensional matrix with 'var_est' observable variables)
#' @param Y_special preprocessed Y
#' @param initialisation indicator of being in the initialisation phase
#' @param indices individuals for which beta_est is being estimated
# @param optimize_kappa indicates if kappa has to be optimized or not (only relevant for the classical algorithm)
#' @param nosetting_local option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
#' @param kappa_candidates Defines the size of the SCAD-penalty used in the classical algorithm. This vector should contain more than 1 element.
#' @inheritParams initialise_beta
#' @inheritParams generate_Y
#' @inheritParams LMROB
#' @param vars_est number of available observed variables for which a coefficient will be estimated. As default it is equal to the number of available observed variables.
#' @param sigma2 sum of squared error terms, scaled by NT
#' @return The function returns a numeric vector (for the default setting: string == "heterogeneous") or a matrix with the estimated beta (if string == "homogeneous").
#' @importFrom ncvreg ncvreg
determine_beta <- function(string, X_special, Y_special, robust,
                           NN, TT,
                           S,
                           method_estimate_beta,
                           initialisation = FALSE, indices = NA,
                           vars_est, sigma2,
                           nosetting_local = FALSE, kappa_candidates = c(2^(-0:-20), 0),
                           special_case_dgp1 = FALSE) {
  stopifnot(string == "homogeneous" | string == "heterogeneous") # these are the two only options

  optimize_kappa <- FALSE
  stopifnot(optimize_kappa == FALSE) # not implemented in RCTS (1. Needs loop over C_Candidates; 2. Hard to create a robust equivalent.)

  if (!(method_estimate_beta == "group" & initialisation == TRUE)) {
    Y_special <- matrix(Y_special, nrow = length(indices), ncol = TT)
    Y_special <- as.vector(t(Y_special)) # order: N1T1, N2T1,...N1T2,...N_endT_end


    # remove NA's from regression:
    vectorNA <- which(is.na(Y_special))

    if (length(vectorNA) > 0) {
      X_special <- X_special[-vectorNA]
      Y_special <- Y_special[-vectorNA]
    }
  }

  # Regression:
  if (robust) {
    # if(exists("use_bramaticroux")) { #-> keep outside of RCTS
    #   model = RpanFE(Y_special, X_special, TT, 0.20, 20, vars_est, length(Y_special)/TT)[[1]] #function is translation of Bramati/Croux code
    # } else {
    model <- LMROB(Y_special, X_special, nosetting = nosetting_local) #-> lmrob(Y_special ~ X_special, setting = "KS2014")
    # }
  } else {
    # ncvreg, without weights
    model <- ncvreg(X_special, Y_special, family = "gaussian", penalty = "SCAD",
                    lambda = kappa_candidates
                    )
    # print(model$beta)
    # plot(log(abs(model$beta[2,])))
    # print("--")
    if (optimize_kappa) {
    stop("Classical case with optimizing kappa is not implemented in RCTS. This would require a fixed C.")
      # C <- .....
      # beta_temp = estimate_beta_classical_optimizekappa(model, kappa_candidates, Y_special, X_special, C, sigma2, TT)
      # print("--")
      # print(dim(beta_temp))
      # return(beta_temp)
    } else {
      #take kappa = 0 (note that kappa_candidates needs > 1 element according to documentation)
      return(model$beta[, length(kappa_candidates)])
    }
  }


  if (string == "homogeneous") {
    if ("lm" %in% class(model) | "lmrob" %in% class(model)) {
      return(matrix(rep(model$coefficients, S), nrow = (vars_est + 1)))
    } else {
      return(matrix(rep(model$beta[, 1], S), nrow = (vars_est + 1)))
    }
  }
  if (string == "heterogeneous") {
    if ("lm" %in% class(model) | "lmrob" %in% class(model)) {
      if (special_case_dgp1) {
        return(c(0, as.numeric(model$coefficients[-2])))
      } else {
        return(as.numeric(model$coefficients))
      }
    } else { # ncvreg or rpanfe

      if (special_case_dgp1) {
        return(c(0, model$beta[-2, 1]))
      } else {
        return(model$beta[, 1]) # for old dgp 2
      }
    }
  }
}



#' Estimates beta.
#'
#' Update step of algorithm to obtain new estimation for beta. Note that we call it beta_est because beta() exists in base R.
#' @inheritParams determine_beta
# @param use_real_world_data Parameter to indicate using real world dataset. Defaults to FALSE.#' @param lambda_group_true loadings of the group factors
#' @param robust TRUE or FALSE: defines using the classical or robust algorithm to estimate beta
#' @param Y Y: NxT dataframe with the panel data of interest
#' @param X X: NxTxp array containing the observable variables
#' @param beta_est estimated values of beta
#' @param g Vector with estimated group membership for all individuals
#' @param factor_group estimated group specific factors
#' @param lambda_group loadings of the estimated group specific factors
#' @param comfactor estimated common factors
#' @param lambda loadings of the estimated common factors
#' @param S number of estimated groups
#' @param k number of common factors to be estimated
#' @param kg number of group specific factors to be estimated
#' @param vars_est number of variables that will be included in the algorithm and have their coefficient estimated. This is usually equal to the number of observable variables.
#' @param num_factors_may_vary whether or not the number of groupfactors is constant over all groups or not
#' @param nosetting option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
#' @param optimize_kappa indicates if kappa has to be optimized or not (only relevant for the classical algorithm)
#' @inheritParams define_object_for_initial_clustering_macropca
#' @inheritParams initialise_beta
#' @return list: 1st element contains matrix (N columns: 1 for each time series of the panel data) with estimated beta_est's. If vars_est is set to 0, the list contains NA.
#' @examples
#' \donttest{
#' X <- X_dgp3
#' Y <- Y_dgp3
#' # Set estimations for group factors and its loadings, and group membership to the true value
#' lambda_group <- lambda_group_true_dgp3
#' factor_group <- factor_group_true_dgp3
#' g <- g_true_dgp3
#' # There are no common factors to be estimated  -> but needs placeholder
#' lambda <- matrix(0, nrow = 1, ncol = 300)
#' comfactor <- matrix(0, nrow = 1, ncol = 30)
#' #
#' # Choose how coefficients of the observable variables are estimated
#' method_estimate_beta <- "individual"
#' method_estimate_factors <- "macro"
#' beta_est <- estimate_beta(
#'   robust = TRUE, Y, X, NA, g, lambda_group, factor_group,
#'   lambda, comfactor,
#'   S = 3, k = 0, kg = c(3, 3, 3),
#'   vars_est = 3
#' )[[1]]
#' }
#' @importFrom stats filter
#' @importFrom purrr pmap
#' @importFrom purrr map2
#' @importFrom rlang .data
#' @export
estimate_beta <- function(robust, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor,
                          method_estimate_beta = "individual",
                          S, k, kg,
                          vars_est,
                          num_factors_may_vary = TRUE,
                          optimize_kappa = FALSE, nosetting = FALSE,
                          special_case_dgp1 = FALSE) {
  #need sigma2 for classical case, for optimizing kappa
  #(note: original implementation of Ando/Bai calculates this only once: after the initialization)
  #here, this is calculated before each iteration
  if(!robust) {
    sigma2 <- calculate_sigma2(calculate_error_term(Y, X, beta_est, g, factor_group, lambda_group, comfactor, lambda,
                                   S, k, kg))
  } else {
    sigma2 <- NA
  }
  NN <- nrow(Y)
  TT <- ncol(Y)
  vars <- dim(X)[3]
  if (vars_est > 0) {
    if (method_estimate_beta == "homogeneous") {
      # X needs to be in the form of (NN*TT x p matrix)
      X_special <- restructure_X_to_order_slowN_fastT(X, vars_est)
      # define Y* as Y - FcLc - FgLg:
      # and get LgFg out of the loop to speed up:
      lgfg_list <- calculate_lgfg(lambda_group, factor_group, S, k, kg, num_factors_may_vary, NN, TT)

      Y_special <- Y
      for (i in 1:NN) {
        # choose the correct lgfg_list, based on group of i
        for (k in 1:S) {
          if (g[i] == k) temp <- lgfg_list[[k]]
        }

        # subtract common factorstructure
        Y_special[i, ] <- Y_special[i, ] - t(lambda[, i]) %*% comfactor[, ]
        if (max(kg, na.rm = TRUE) > 0) {
          # subtract group factorstructure
          Y_special[i, ] <- Y_special[i, ] - temp[i, ]
        }
      }

      beta_est <- determine_beta("homogeneous", X_special, Y_special, robust,
        NN, TT, S, method_estimate_beta,
        indices = 1:NN, vars_est, sigma2, nosetting_local = nosetting,
        special_case_dgp1
      )
    }
    if (method_estimate_beta == "group") {
      beta_est <- matrix(NA, nrow = (vars + 1), ncol = S)
      for (group in 1:S) {
        # select parts of X and Y of this group
        indices_group <- which(g == group)
        if (length(indices_group) == 0) {
          warning("Empty group! This message should never be seen.")
        }

        # X needs to be in the form of (NN*TT x p matrix)
        X_special <- restructure_X_to_order_slowN_fastT(
          array(X[indices_group, , ], dim = c(length(indices_group), TT, dim(X)[3])),
          vars_est
        )
        # define Y* as Y - FcLc - FgLg:
        Y_special <- matrix(Y[indices_group, ], nrow = length(indices_group)) %>% unlist()


        for (i in 1:length(indices_group)) {
          index <- indices_group[i]
          LAMBDAGROUP <- as.matrix(lambda_group %>% dplyr::filter(.data$id %in% index) %>% dplyr::select(-.data$group, -.data$id))
          Y_special[i, ] <- Y_special[i, ] -
            t(lambda[, index]) %*% comfactor[, ] -
            LAMBDAGROUP %*% factor_group[[g[index]]]
        }


        beta_est[, group] <- determine_beta("heterogeneous", X_special, Y_special, robust,
          NN, TT, S, method_estimate_beta,
          indices = indices_group, vars_est, sigma2, nosetting_local = nosetting,
          special_case_dgp1
        )
      }
    }
    if (method_estimate_beta == "individual") {
      #######################################################################
      # first, X and Y need to be restructured: made into list with N elements.
      #######################################################################

      X_local <- adapt_X_estimating_less_variables(X, vars_est) # makes X smaller when vars_est < vars
      # if(use_real_world_data) {
      #   number_of_vars = vars
      # } else {
      number_of_vars <- vars_est
      # }
      X_special_list <- lapply(1:NN, function(x) {
        restructure_X_to_order_slowN_fastT(
          matrix(X_local[x, , ], ncol = number_of_vars),
          vars_est
        )
      })
      # helpfunction: calculates Y - factorstructures. This is equal to XB + errorterm.
      # The result is a list with number of elements equal to number of time series in Y.
      make_Y_special <- function(i) {
        Y_special <- unlist(matrix(Y[i, ], nrow = 1))
        if (dim(t(Y_special))[2] == TT) {
          Y_special <- t(Y_special)
        }

        to_remove <- c(
          which(names(lambda_group) == "group"),
          which(names(lambda_group) == "id")
        )
        LAMBDAGROUP <- as.matrix(lambda_group[which(lambda_group$id %in% i), ][, -to_remove])
        if (g[i] != 0) {
          LAMBDAGROUP <- LAMBDAGROUP[1:kg[g[i]]]
          Y_special <- Y_special -
            t(lambda[, i]) %*% comfactor[, ] -
            LAMBDAGROUP %*% factor_group[[g[i]]]
        } else { # class zero -> no groupfactorstructure
          LAMBDAGROUP <- NA
          Y_special <- Y_special -
            t(lambda[, i]) %*% comfactor[, ]
        }

        return(Y_special)
      }

      Y_special_list <- lapply(1:NN, function(x) make_Y_special(x))

      #######################################################################
      # secondly, beta_i is estimated for all i
      #######################################################################
      if (optimize_kappa) {
        beta_est <- pmap(
          list(X_special_list, Y_special_list, 1:NN),
          function(x, y, z) {
            determine_beta("heterogeneous", x, y, robust,
              NN, TT, S, method_estimate_beta,
              indices = z, vars_est, sigma2, nosetting_local = nosetting,
              special_case_dgp1
            )
          }
        )
      } else {
        # note that mapply isnt faster than map2
        # print(length(X_special_list))
        # print(length(Y_special_list))
        beta_est <- map2(
          X_special_list, Y_special_list,
          function(x, y) {
            determine_beta("heterogeneous", x, y, robust,
              NN, TT, S, method_estimate_beta,
              indices = NA, vars_est, sigma2, nosetting_local = nosetting,
              special_case_dgp1
            )
          }
        )




        # micro = microbenchmark::microbenchmark(
        #   beta_est = map2(X_special_list, Y_special_list,  function(x,y) determine_beta(...) ),
        #   beta_new = mapply( function(x,y) { determine_beta(...) }, x = X_special_list, y = Y_special_list ),
        #   times = 15
        # )
        # print(summary(micro))
        # expr      min       lq     mean   median       uq      max neval cld
        # 1     beta_est 3.544195 3.585846 3.643604 3.622231 3.672088 3.829970    15   a
        # 2 beta_new 3.543808 3.553852 3.621918 3.597413 3.680516 3.731678    15   a
      }
      ########################################################
      # Possible use of future::map2 instead of map2:
      # (is 2 to 3 times faster (microbenchmark) for (300, 30)-system)
      # BUT:  ERRORS:
      #   Error: 'map2' is not an exported object from 'namespace:future'
      # beta_est = future::map2(X_special_list, Y_special_list,  function(x,y) determine_beta(...) )
      ########################################################
      if (length(unlist(beta_est)) != nrow(Y) * (vars_est + 1)) {
        warning("There are missing beta's for some individuals.")
      }
      # print("***")
      # print(length(unlist(beta_est)))
      # if lmrob did not converge, there are no new coefficients -> this was due to forgotten return in warning()-part of LMROB -> solved now
      # missing_coefficients = which((purrr::map(beta_est, length) %>% unlist) != (vars_est + 1))
      # print(missing_coefficients)
      # print(beta_est[[missing_coefficients-1]])
      # if(length(missing_coefficients) > 0) {
      #   for(q in 1:missing_coefficients) {
      #     beta_est[[q]] = previous_beta_est[,q]
      #   }
      # }
      # Now this should not reach a warning anymore
      beta_est <- tryCatch(
        matrix(unlist(beta_est), ncol = NN),
        warning = function(w) {
          print("------")
          message(w)
          print(unlist(beta_est))
          print(NN)
          # print("sleep")
          # Sys.sleep(36000)
        }
      )
      # beta_est <- matrix(unlist(beta_est), ncol = NN)
    }

    return(list(beta_est))
  } else {
    return(list(NA))
  }
}


#' Calculates W = Y - X*beta_est. It is used in the initialization step of the algorithm, to initialise the factorstructures.
#'
#' @param beta_est estimated values of beta
#' @param g Vector with group membership for all individuals
#' @inheritParams estimate_beta
#' @return NxT matrix
calculate_W <- function(Y, X, beta_est, g,
                        vars_est,
                        method_estimate_beta) {
  NN <- nrow(Y)
  TT <- ncol(Y)
  W <- matrix(0, nrow = NN, ncol = TT) # T x N-matrix in original paper , but I rather define as NxT

  # if vars_est < vars the obsoleterows in beta_est were already erased -> do the same in X
  X <- adapt_X_estimating_less_variables(X, vars_est)

  if (vars_est > 0) {
    if (method_estimate_beta == "homogeneous") {
      # non-dependence on g -> take first column
      for (i in 1:NN) W[i, ] <- Y[i, ] - t(cbind(1, X[i, , ]) %*% as.matrix(beta_est[, 1]))
    }
    if (method_estimate_beta == "group") {
      for (i in 1:NN) W[i, ] <- Y[i, ] - t(cbind(1, X[i, , ]) %*% as.matrix(beta_est[, g[i]]))
    }
    if (method_estimate_beta == "individual") {
      for (i in 1:NN) {
        W[i, ] <- Y[i, ] - t(cbind(1, X[i, , ]) %*% as.matrix(beta_est[, i]))
      }
    }
  } else {
    for (i in 1:NN) W <- as.matrix(Y)
  }

  if (anyNA(W)) {

    # handle NA's: which currently means: remove the rows with NA's
    W <- handleNA(W)[[1]]
  }
  return(W)
}

#' Calculates Z = Y - X*beta_est - LgFg. It is used in the estimate of the common factorstructure.
#'
#' @inheritParams calculate_W
#' @inheritParams estimate_factor
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param initialise indicator of being in the initialisation phase
#' @inheritParams estimate_beta
#' @return NxT matrix
calculate_Z_common <- function(Y, X, beta_est, g, lgfg_list,
                               vars_est,
                               kg,
                               method_estimate_beta, method_estimate_factors,
                               initialise = FALSE) {
  NN <- nrow(Y)
  TT <- ncol(Y)
  Z <- matrix(0, nrow = NN, ncol = TT)

  # if vars_est < vars the obsolete rows in beta_est were already erased -> do the same in X
  X <- adapt_X_estimating_less_variables(X, vars_est)

  for (i in 1:NN) {
    y <- Y[i, ] %>% as.numeric()
    if (do_we_estimate_group_factors(kg)) { # when there are group factors to be estimated

      if (method_estimate_factors == "pertmm" & !("list" %in% class(lgfg_list))) {
        # then lgfg_list is not a list, but a data frame of dimension NxT
        LF_GROUP <- lgfg_list[i, ]
      } else {
        if (g[i] == 0) { # class zero
          LF_GROUP <- 0
        } else {
          LF_GROUP <- lgfg_list[[g[i]]][i, ]
        }
      }
    } else {
      LF_GROUP <- 0
    }

    if (vars_est > 0) {
      if ((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) BETA <- as.matrix(beta_est[, g[i]])
      if (method_estimate_beta == "individual") BETA <- as.matrix(beta_est[, i])
      Z[i, ] <- y - t(cbind(1, X[i, , ]) %*% BETA) - LF_GROUP
    } else {
      Z[i, ] <- y - LF_GROUP
    }
  }


  if (anyNA(Z)) {
    # handle NA's: which currently means: remove the rows with NA's
    Z <- handleNA(Z)[[1]]
  }
  return(Z)
}

#' Calculates Z = Y - X*beta_est - LF. It is used to estimate the groupfactorstructure.
#'
#' @inheritParams calculate_W
#' @inheritParams calculate_Z_common
#' @inheritParams estimate_beta
#' @param group indexnumber of the group
#' @return NxT matrix
calculate_Z_group <- function(Y, X, beta_est, g, lambda, comfactor, group,
                              k,
                              method_estimate_beta,
                              initialise,
                              vars_est) {
  TT <- ncol(Y)
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
    X <- adapt_X_estimating_less_variables(X, vars_est)
  } else {
    vars_est <- 0
  }
  indices_group <- which(g == group)
  if (length(indices_group) == 0) {
    message("empty group (calculate_Z_group())")
  }

  Z <- matrix(0, nrow = length(indices_group), ncol = TT) # Nj x T matrix
  # if vars_est < vars the obsoleterows in beta_est were already erased -> do the same in X

  for (i in 1:length(indices_group)) { # loop over the number of elements in the group
    index <- indices_group[i]
    if(!between(index, 1, nrow(Y))) {
      print(index)
      print(table(g))
      warning("index has an incorrect value")
    }

    # define XT
    if (vars_est > 0) {
      if ((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) BETA <- as.matrix(beta_est[, g[index]])
      if (method_estimate_beta == "individual") BETA <- as.matrix(beta_est[, index])

      XT <- t(cbind(1, X[index, , ]) %*% BETA)
    } else {
      XT <- rep(0, TT)
    }

    # calculate Z = Y - X*beta_est - LF)
    y <- as.numeric(Y[index, ])
    if (initialise) {
      Z[i, ] <- y - XT
    } else {
      a <- do_we_estimate_common_factors(k)
      Z[i, ] <- y - XT -
        a * t(lambda[, index]) %*% comfactor
    }
  }


  if (anyNA(Z)) {
    # handle NA's: which currently means: remove the rows with NA's
    Z <- handleNA(Z)[[1]]
  }

  return(Z)
}

#' Solves a very specific issue with MacroPCA.
#'
#' MacroPCA crashes Rstudio with certain dimensions of the input. Solve this by doubling every row.
#' No information is added by this, so there is no influence on the end result,
#' but crashes of Rstudio are evaded.
#' @param object input
#' @param verbose prints messages
#' @return matrix
evade_crashes_macropca <- function(object, verbose = FALSE) {
  #--------------------MacroPCA seems to make Rstudio crash when the dimension of object = (27, 193) ------------------
  size_with_crashes1 <- c(27, 27, 25, 43, 66, 24)
  size_with_crashes2 <- c(193, 310, 600, 560, 387, 374)
  changed_objectsize <- FALSE

  if (verbose) message("Test size of object (for MacroPCA)")

  for (i in 1:length(size_with_crashes1)) {
    if (changed_objectsize == FALSE & dim(object)[1] == size_with_crashes1[i] & dim(object)[2] == size_with_crashes2[i]) {
      if (verbose) {
        message(paste("Rstudio would crash when applying MacroPCA on an object of size", dim(object)[1], " ", dim(object)[2], " -> double amount of rows of input (no information added, and result does not change)."))
        print(i)
      }
      object <- rbind(object, object)
      if (verbose) message(dim(object))
      changed_objectsize <- TRUE
    }
  }
  if (verbose) message("Test size: done")
  return(object)
}

#' Helpfunction in robustpca().
#'
#' It handles possible thrown errors in MacroPCA.
#' @param object input
#' @param temp this is the result of the trycatch block of using macropca on object
#' @param KMAX parameter kmax in MacroPCA
#' @param number_eigenvectors number of principal components that are needed
#' @inheritParams update_g
#' @return matrix of which the columns contain the chosen amount of eigenvectors of object
handle_macropca_errors <- function(object, temp, KMAX, number_eigenvectors, verbose = FALSE) {
  if ("error" %in% class(temp)) {
    if (verbose) {
      message("*******************************")
      message(paste("MacroPCA with", number_eigenvectors, " eigenvectors fails, -> calculate different (higher) amount of eigenvectors and select afterwards the necessary ones.
                    Start with a couple more and decrease 1 by 1."))
    }
    k_higher_value <- 0
    counter <- 0
    #NOTE: changing k does help now and then (also with lower k, but in that case there is a factor short of course, and that one should then be imputed)
    #-> depreciate this while-loop, and immediately go to the classical estimation as solution
    while ("error" %in% class(temp)) {
      counter <- counter + 1
      # if (verbose) print(paste("counter:", counter))
      # print(paste("set to", max(number_eigenvectors + k_higher_value, number_eigenvectors) - counter))
      # temp <- tryCatch(
      #   cellWise::MacroPCA(object, k = max(number_eigenvectors + k_higher_value, number_eigenvectors) - counter, MacroPCApars = list(kmax = KMAX, silent = TRUE)),
      #   error = function(e) {
      #     print(e)
      #     return(e)
      #   },
      #   finally = {
      #     #print("new macropca is done")
      #   }
      # )

      # sometimes, when there occurred too often errors in macropca, we end up with only a limited amount of possible eigenvectors that are calculated.
      # This happens when a group has very little elements.
      # if (verbose) {
      #   print("----")
      #   print(class(temp))
      # }
      # number_columns <- 999
      # if(!("error" %in% class(temp))) { #otherwise temp$loadings would return NULL (and no error)
      #   number_columns <- tryCatch(
      #     ncol(temp$loadings),
      #     error = function(e) {
      #       print(e)
      #       return(999)
      #     }
      #   )
      #
      #   temp <- tryCatch(
      #     temp$loadings[, 1:min(number_eigenvectors, number_columns)],
      #     error = function(e) {
      #       print(e)
      #       return(e)
      #     }
      #   )
      # }

      # if (number_columns < number_eigenvectors) {
      #   message(paste("add", number_eigenvectors - number_columns, "columns to temp"))
      #   print(temp)
      # }

      if (counter >= k_higher_value) { # counter >= number_eigenvectors
        #-> then infinite loop (MacroPCA does not work with any k)

        # cause unknown:
        # unlikely to be caused by "small" amount of units in this group : (occurs up to 51 units have I found)
        # can happen during any iteration at (150,30), but happens only in iteration 0 or 1 at (300,200)
        # it does seem to occur more with small N / small T
        # happens with higher amount of groupfactors to be estimated

        # note: warning() does not print in console here -> use message()
        if (verbose) message("\nMacroPCA failed to estimate factors for a group with ", nrow(object), " units -> use in this iteration for this group the classical estimation of the factors.")
        # eigen() needs always square matrix -> take covmatrix
        temp <- eigen(t(object) %*% object)$vectors[, 1:number_eigenvectors]
      }
    }
    # print(paste("endcounter:", counter))
  }
  if(is.null(temp)) {
    message("--temp should not be null!--")
  }
  return(temp)
}

#' Function that uses robust PCA and estimates robust factors and loadings.
#'
#' Contains call to MacroPCA()
#'
#' Notes:
#'
#' Different values for kmax give different factors, but the product lambda*factor stays constant. Note that
#'  this number needs to be big enough, otherwise eigen() will be used. Variation in k does give different results for lambda*factor
#'
#' MacroPCA() crashes with specific values of dim(object). For example when dim(object) = c(193,27).
#'   This is solved with evade_crashes_macropca(), for those problematic dimensions that are already encountered during tests.
#' @param object input
#' @param number_eigenvectors number of eigenvectors to extract
#' @param KMAX The maximal number of principal components to compute. This is a parameter in cellWise::MacroPCA()
#' @param verbose_robustpca when TRUE, it prints messages: used for testing (requires Matrix-package when set to TRUE)
#' @return list with as the first element the robust factors and as the second element the robust factor loadings
robustpca <- function(object, number_eigenvectors, KMAX = 20, verbose_robustpca = FALSE) {
  if (verbose_robustpca) {
    print(paste("*************************************************robust PCA with:", number_eigenvectors))
    print("dimension of input:")
    print(dim(object))
  }


  ######################
  # MacroPCA
  ######################
  error_macropca <- FALSE
  object <- evade_crashes_macropca(object, verbose = verbose_robustpca)

  # print(number_eigenvectors)
  if (number_eigenvectors > KMAX) {
    # Note that when k > kmax, k gets the value of kmax.
    message("MacroPCA is (through KMAX) limited to 20 factors.")
  }

  macropca_kmax <- number_eigenvectors

  # note that in the documentation it says the default for scale is FALSE
  # But:
  #
  # MacroPCA(M, 1, MacroPCApars=(list(scale = TRUE)))$loadings
  # MacroPCA(M, 1, MacroPCApars=(list(scale = FALSE)))$loadings
  # MacroPCA(M, 1)$loadings
  # -> 1st and 3rd give the same result
  # -> The actual default for scale seems to be TRUE


  if (verbose_robustpca) {
    message("--start of macropca")
    print(paste("rank of input:", Matrix::rankMatrix(object)))
    print(paste("required number of eigenvectors:", number_eigenvectors))
  }

  temp <- tryCatch(
    {
      # print("test1")
      # print(dim(object))
      cellWise::MacroPCA(object, k = max(macropca_kmax, number_eigenvectors), MacroPCApars = list(kmax = KMAX, silent = TRUE))
    },
    error = function(e) {
      #message(e) #THIS IS FORBIDDEN, AS IT FAILS THE PARALLEL SYSTEM SOMEHOW!!!!!!!!! (using message() itself is ok though)
      #print(e) #using print is ok, but the text gets too long, and looks too unhelpful (and I know the cause anyway)
      return(e)
    },
    finally = {
      # print("test2")
    }
  )
#--HERE USED TO BE AN ISSUE IN THE PARALLEL SYSTEM (both with "do" and with "dopar"): SOMETIMES IT STOPS HERE!!!! The serialized system works normal.--
#--This seems also to be linked with small amount of units in one of the groups.
#--Increasing ROBUST_THRESHOLD (e.g. to 25) works as a hack around it, (works as a hybrid classical/robust estimation), but it really should not come to this...
#-----> SOLUTION: removing 'message(e)' seems to work...
  # print("-") #this should always print, after "test2".
  # print(paste("is.null: ",is.null(temp)))
  # print("test3")
  if (verbose_robustpca) print(class(temp))
  if ("error" %in% class(temp)) {
    error_macropca <- TRUE
  }
  ###############################
  # sometimes the output of MacroPCA() has too little columns (should be equal to "number_eigenvectors").
  # Possible related to performing MacroPCA on a very small group, and/or the elements of the group are very similar.
  # Solution: add one-column(s) to the factors (-> temp$loadings), and zero-column(s) to factor loadings (-> temp$scores).
  # This ensures that the product of factor and factor loadings does not get altered.
  #  (Note: adding zero-column to temp$scores would have the result that the rank of factor_group[...] is too low, therefore use one-column(s).)
  if (verbose_robustpca) {
    message("----dimension of output of macropca (factorloadings and factors):----")
    print(dim(temp$scores))
    print(dim(temp$loadings))
    print(paste("required number of eigenvectors:", number_eigenvectors))
  }
  if (!is.null(dim(temp$scores))) {
    if (dim(temp$loadings)[2] != number_eigenvectors) {
      # note: this only works if macropca returns maximum one factor too little

      # the warning message that is returned:
      # In svd::propack.svd(Y, neig = ncomp) :
      #   Only ... singular triplets converged within ... iterations.
      message(paste(
        "There is an issue with the dimensions of the factors: MacroPCA returns only", dim(temp$loadings)[2],
        "eigenvectors, instead of", number_eigenvectors, "."
      ))

      # for (i in 1:(number_eigenvectors - dim(temp$loadings)[2])) {
      #   temp$scores <- cbind(temp$scores, 0) # add 0 to factor loadings
      #   temp$loadings <- cbind(temp$loadings, 1) # add 1 to factors
      #   if (verbose_robustpca) {
      #     message("resulting factors:")
      #     print(temp$loadings[1:5, ])
      #   }
      # }

      temp$scores <- NULL
      temp$loadings <- NULL
    }
  }
  ##############################
  scores <- temp$scores[, 1:number_eigenvectors] # these are the factor loadings
  factors_macropca <- temp$loadings[, 1:number_eigenvectors] # these are the factors

  #note: new function is() requires as dependency package "methods" -> use %in% instead
  # suppressWarnings(
  #   if (class(factors_macropca) == "numeric") { # case of estimating 1 factor -> numeric -> make matrix
  #     factors_macropca <- matrix(factors_macropca)
  #   }
  # )
  # if (is(factors_macropca, "numeric")) { # case of estimating 1 factor -> numeric -> make matrix
  #        factors_macropca <- matrix(factors_macropca)
  # }
  if ("numeric" %in% class(factors_macropca)) { # case of estimating 1 factor -> numeric -> make matrix
    factors_macropca <- matrix(factors_macropca)
  }

  if (!error_macropca & !is.null(factors_macropca)) {
    # rare issue with macropca
    if (nrow(factors_macropca) != ncol(object)) {
      print(dim(object))
      print(dim(factors_macropca))
      warning("--MacroPCA has dropped a column for unknown reasons---") # This leads to wrong dimensions in the factors, and gives error in rstudio.
    }
  }
  if (error_macropca) {
    factors_macropca <- handle_macropca_errors(object, temp, KMAX, number_eigenvectors, verbose_robustpca)
  }
  if (verbose_robustpca) print("end of macropca")
  return(list(factors_macropca, scores))
}




#' Estimates common factor(s) F.
#'
#' The estimator for F, see Anderson (1984), is equal to the first k eigenvectors (multiplied by sqrt(T) due to the restriction F'F/T = I)
#' associated with first r largest eigenvalues of the matrix WW' (which is of size TxT).
#' @inheritParams calculate_W
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param method_estimate_factors defines method of robust estimaton of the factors: "macro", "pertmm" or "cz"
#' @inheritParams estimate_beta
#' @inheritParams calculate_virtual_factor_and_lambda_group
#' @inheritParams calculate_Z_common
#' @inheritParams update_g
#' @return Return a list. The first element contains the k x T matrix with the k estimated common factors. The second element contains either the
#' robust MacroPCA-based loadings or NA.
estimate_factor <- function(robust, Y, X, beta_est, g, lgfg_list,
                            k, kg,
                            method_estimate_beta,
                            method_estimate_factors,
                            initialise = FALSE,
                            verbose = FALSE) {

  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  NN <- nrow(Y)
  TT <- ncol(Y)
  if (!do_we_estimate_common_factors(k)) {
    # return matrix with only zero's
    return(list(t(matrix(rep(0, TT))), NA))
  }
  # initialisation: has no grouped factorstructure yet
  if (initialise) {
    W <- calculate_W(
      Y, X, beta_est, g,
      vars_est,
      method_estimate_beta
    ) # Y - XT
  } else {
    W <- calculate_Z_common(
      Y, X, beta_est, g, lgfg_list,
      vars_est,
      kg,
      method_estimate_beta, method_estimate_factors
    ) # this is a "rows_without_NA x T - matrix"
  }
  if (verbose) message("W is constructed")
  # if there are individuals in "class zero", then they are considered outliers
  #  and need to be excluded in the estimation of the common factors as well
  if (verbose) message("removing class zero individuals from estimation of common factors")
  W <- W[g != 0, ]



  # Define the object on which (robust or classical) PCA will be performed
  if (robust) {
    # CHANGE LOCATION 4/7
    if (method_estimate_factors %in% c("macro", "pertmm")) {
      temp <- prepare_for_robpca(W, NN, TT)
    } else {
      temp <- t(W) %*% W / (NN * TT)
    }
  } else {
    # Classical case
    # calculate 1/NT * W'W
    temp <- t(W) %*% W / (NN * TT) # TxT matrix (this division by NT does not actually matter for the calculation of the eigenvectors)
  }


  # first k eigenvectors * sqrt(T)

  # take k eigenvectors
  message(paste0("Estimate ", k, " common factors"))
  if (robust) {
    # CHANGE LOCATION 1/7
    if (method_estimate_factors %in% c("macro", "pertmm")) {
      temp2 <- robustpca(temp, k, verbose_robustpca = verbose)
      estimatorF <- t(sqrt(TT) * temp2[[1]])
      scores <- temp2[[2]] # =factor loadings coming out of macropca
      rm(temp2)
    } else {
      estimatorF <- t(sqrt(TT) * eigen(temp)$vectors[, 1:k])
      scores <- NA # because loadings will be calculated later
    }
  } else {
    estimatorF <- t(sqrt(TT) * eigen(temp)$vectors[, 1:k])
    scores <- NA # because loadings will be calculated later
  }
  # OLD:
  # if(!do_we_estimate_common_factors(k) ) {
  #   #then estimatorF has length TT, which is still an appropriate size to use in this case -> just set values to zero
  #   estimatorF = estimatorF - estimatorF
  # }


  return(list(estimatorF, scores))
}

#' Helpfunction: prepares object to perform robust PCA on.
#'
#' It contains options to use the classical or robust covmatrix or no covariance matrix at all.
#' @param object this is the object of which we may take the covariance matrix and then to perform robust PCA on
#' @param NN N
#' @param TT T
#' @param option 1 (robust covmatrix), 2 (classical covmatrix), 3 (no covmatrix)
#' @return matrix
prepare_for_robpca <- function(object, NN, TT, option = 3) {
  if (option == 1) {
    message("not used anymore")
    # temp = get_robust_covmatrix(object)
  }

  if (option == 2) temp <- t(object) %*% object / (NN * TT)

  if (option == 3) temp <- object

  return(temp)
}



#' Estimates group factors Fg.
#'
#' @inheritParams calculate_W
#' @inheritParams calculate_Z_common
#' @inheritParams estimate_beta
#' @inheritParams calculate_virtual_factor_and_lambda_group
#' @inheritParams estimate_factor
#' @inheritParams update_g
#' @return Return a list with an element for each estimated group. Each element of the list is a matrix with the group specific factors as rows.
#' @export
#' @examples
#' #example with data generated with DGP 2
#' data <- create_data_dgp2(30, 10)
#' Y <- data[[1]]
#' X <- data[[2]]
#' g <- data[[3]] #true group membership
#' set.seed(1)
#' beta_est <- matrix(rnorm(4 * nrow(Y)), nrow = 4)
#' factor_group <- data[[5]] #true values of group specific factors
#' comfactor <- matrix(0, nrow = 1, ncol = ncol(Y))
#' lambda <- matrix(0, nrow = 1, ncol = nrow(Y))
#' estimate_factor_group(TRUE, Y, X, beta_est, g, lambda, comfactor, factor_group,
#' 3, 0, c(3, 3, 3))
estimate_factor_group <- function(robust, Y, X, beta_est, g, lambda, comfactor, factor_group,
                                  S, k, kg,
                                  method_estimate_beta = "individual", method_estimate_factors = "macro",
                                  initialise = FALSE,
                                  verbose = FALSE
                                  # returnscores = FALSE
) {
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  NN <- nrow(Y)
  TT <- ncol(Y)
  factor_group_previous <- factor_group # save for if macropca fails

  estimatorF <- list()
  scores <- list()
  for (group in 1:S) {
    if (verbose) print(paste("Estimate factors for group", group))
    if (kg[group] > 0) {
      if (TT < kg[group]) {
        warning("There are too many factors to be estimated, compared to TT.")
      }
      Wj <- calculate_Z_group(
        Y, X, beta_est, g, lambda, comfactor, group,
        k,
        method_estimate_beta,
        initialise,
        vars_est
      )


      # use limit on nrow(Wj) due to error ("The input data must have at least 3 rows (cases)")
      # When the number of rows within the group is bigger than 3, macropca runs, but in some cases can drop columns (when values within this column are equal).
      # So, in practice it is better to increase the threshold a bit.
      ROBUST_THRESHOLD <- 5
      if (robust & nrow(Wj) > ROBUST_THRESHOLD) {


        # CHANGE LOCATION 2/7
        if (method_estimate_factors %in% c("macro", "pertmm")) {
          temp <- prepare_for_robpca(Wj, NN, TT)
          temp2 <- robustpca(temp, kg[group], verbose_robustpca = verbose)
          if (is.null(temp2[[1]])) {
            message(paste("MacroPCA returned a wrong amount of factors, so skip the update of the group factors for group", group, "(has ", as.numeric(table(g)[group]), " units) in this iteration."))
            if (!initialise) {
              estimatorF[[group]] <- factor_group_previous[[group]] #-1 #placeholder
            } else {
              # if MacroPCA() fails in the initial step, there is no previous value available.
              message("MacroPCA() fails in the initial step for this group, there is thus no previous value available. Initialize factor at random for this group.")
              estimatorF[[group]] <- matrix(rnorm(kg[group] * TT), nrow = kg[group], ncol = TT)
            }
          } else {
            estimatorF[[group]] <- t(sqrt(TT) * temp2[[1]]) # robust pca
            scores[[group]] <- temp2[[2]]
            rm(temp2)
          }
        } else {
          ## CZ
          temp <- t(Wj) %*% Wj # dividing by NT does not make any difference for the eigenvectors
          scores[[group]] <- NA # because loadings will be calculated later
          # print(kg)
          # print(group)
          estimatorF[[group]] <- t(sqrt(TT) * eigen(temp)$vectors[, 1:kg[group]])
        }
      } else {
        if (nrow(Wj) < ROBUST_THRESHOLD & robust & verbose) message("One or more groups contain a very small number of elements. The robust estimation of factors with macropca cannot work -> use non-robust estimation of factors for these groups.")

        temp <- t(Wj) %*% Wj # dividing by NT does not make any difference for the eigenvectors
        scores[[group]] <- NA # because loadings will be calculated later
        estimatorF[[group]] <- t(sqrt(TT) * eigen(temp)$vectors[, 1:kg[group]])
      }
    } else {
      estimatorF[[group]] <- matrix(0, 1, TT)
      scores[[group]] <- NA # because loadings will be calculated later
    }
  }



  # as test:
  # if(returnscores) {
  #   return(scores)
  # }
  return(estimatorF)
}

#' calculates factor loadings of common factors
#' @inheritParams calculate_W
#' @param comfactor common factors
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @inheritParams estimate_beta
#' @inheritParams calculate_Z_common
#' @inheritParams estimate_factor
#' @return Returns a matrix where each row contains a common factor. If the number of estimated common factors equals zero, it returns a matrix with 1 row, containing zero's.
calculate_lambda <- function(robust, Y, X, beta_est, comfactor, factor_group, g, lgfg_list,
                             k,
                             kg,
                             method_estimate_beta,
                             method_estimate_factors,
                             verbose = FALSE,
                             initialise = FALSE) {
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  NN <- nrow(Y)
  TT <- ncol(Y)
  if (!do_we_estimate_common_factors(k)) {
    return(t(matrix(rep(0, NN))))
  }
  if (initialise) {
    W <- calculate_W(
      Y, X, beta_est, g,
      vars_est,
      method_estimate_beta
    )
  } else {
    W <- calculate_Z_common(
      Y, X, beta_est, g, lgfg_list,
      vars_est,
      kg,
      method_estimate_beta, method_estimate_factors
    )
  }

  if (robust) {
    # CHANGE LOCATION 6/7
    if (method_estimate_factors %in% c("macro", "pertmm")) {
      lambda <- return_robust_lambdaobject(W, NA,
        type = 2, g, NN,
        nrow(comfactor),
        kg,
        comfactor_rrn = comfactor,
        factor_group_rrn = factor_group,
        verbose
        # use_real_world_data_rrn = use_real_world_data,
      )
    } else {
      lambda <- t(W %*% t(comfactor) / TT)
    }
  } else {
    lambda <- t(W %*% t(comfactor) / TT)
  }




  ############################
  # for time-series with NA's,
  # the lambda's could not be determined (dimension of lambda is smaller) -> set NA's to zero
  ############################
  if (anyNA(Y)) {
    rows_with_NA <- which(apply(Y, 1, anyNA))
    rows_without_NA <- which(!apply(Y, 1, anyNA))
    temp <- data.frame(matrix(NA, nrow = max(k, 1), ncol = NN)) # new dataframe
    temp[, rows_without_NA] <- lambda # put calculated lambda's in df
    if (sum(rows_with_NA, na.rm = TRUE) > 0) temp[, rows_with_NA] <- 0 # set the rest to zero
    lambda <- as.matrix(temp)
  }


  return(lambda)
}

#' calculates factor loadings of groupfactors
#'
#' returns object which includes group and id of the individuals
#' @param beta_est estimated values of beta
# @param UPDATE1 option to indicate the number of groupfactors is updated during the algorithm; defaults to FALSE
# @param UPDATE2 option to indicate the number of common factors is updated during the algorithm; defaults to FALSE
#' @inheritParams estimate_beta
#' @inheritParams calculate_W
#' @inheritParams calculate_Z_common
#' @inheritParams estimate_factor
#' @importFrom rlang .data
#' @return Returns a data.frame with a row for each time series. The first number of columns contain the individual loadings to the group specific factors.
#' Furthermore "group" (group membership) and id (the order in which the time series appear in Y) are added.
#' @export
#' @examples
#' #' #example with data generated with DGP 2
#' data <- create_data_dgp2(30, 10)
#' Y <- data[[1]]
#' X <- data[[2]]
#' g <- data[[3]] #true group membership
#' set.seed(1)
#' beta_est <- matrix(rnorm(4 * nrow(Y)), nrow = 4)
#' factor_group <- data[[5]] #true values of group specific factors
#' comfactor <- matrix(0, nrow = 1, ncol = ncol(Y))
#' lambda <- matrix(0, nrow = 1, ncol = nrow(Y))
#' calculate_lambda_group(TRUE, Y, X, beta_est, factor_group, g, lambda, comfactor,
#' 3, 0, c(3, 3, 3) )
calculate_lambda_group <- function(robust, Y, X, beta_est, factor_group, g, lambda, comfactor,
                                   S,
                                   k,
                                   kg,
                                   method_estimate_beta = "individual", method_estimate_factors = "macro",
                                   verbose = FALSE,
                                   initialise = FALSE
                                   # UPDATE1 = FALSE, UPDATE2 = FALSE
) {
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  NN <- nrow(Y)
  TT <- ncol(Y)

  lambda_local <- list()

  for (group in 1:S) {
    if (kg[group] > 0) {
      Wj <- calculate_Z_group(
        Y, X, beta_est, g, lambda, comfactor, group,
        k,
        method_estimate_beta,
        initialise,
        vars_est
      )

      # robust things:
      if (robust) {

        # CHANGE LOCATION 7/7
        if (method_estimate_factors %in% c("macro", "pertmm")) {
          # In the classical approach each lambda is the mean of a set of products of F and Y (Z),
          #   (lambda_N1 = (F_T1 * Y_N1T1 + F_T2 * Y_N1T2 + ...) / TT)
          #   We replace this mean by an M-estimator in the robust approach.

          lambda_local[[group]] <- return_robust_lambdaobject(Wj, group,
            type = 3, g, NN,
            k,
            unlist(lapply(factor_group, nrow)),
            comfactor,
            factor_group,
            verbose
            # use_real_world_data_rrn = use_real_world_data
          )
        } else {
          FGG <- factor_group[[group]]
          if (dim(t(FGG))[1] == 1) {
            FGG <- matrix(FGG, nrow = 1)
          }
          lambda_local[[group]] <- Wj %*% t(FGG) / TT
        }
      } else {
        FGG <- factor_group[[group]]
        if (dim(t(FGG))[1] == 1) {
          FGG <- matrix(FGG, nrow = 1)
        }
        lambda_local[[group]] <- Wj %*% t(FGG) / TT
      }
    } else { # estimate zero group factors for this group -> matrix with 1 column
      # lambda_local[[group]] = matrix(0,
      #                                length(rows_without_NA[(rows_without_NA %in% which(g == group))]),
      #                                1)
      lambda_local[[group]] <- matrix(0, table(g)[group], 1)
    }
  }


  # for income-series with NA's, the lambda's cannot be determined -> set to zero with this function
  add_lambdas_from_NArows <- function(DF, group) {
    rows_with_NA <- which(apply(Y, 1, anyNA))
    rows_without_NA <- which(!apply(Y, 1, anyNA))
    ids <- rows_with_NA[(rows_with_NA %in% which(g == group))]
    if (length(ids) > 0) {
      temp_init <- rep(0, length(ids)) # zero's
      temp <- temp_init
      if (kg[group] > 1) {
        for (i in 2:kg[group]) {
          temp <- temp %>% cbind(temp_init)
        }
      }
      temp <- data.frame(temp %>% cbind(group) %>% cbind(ids))
      names(temp) <- names(DF)
      return(temp)
    } else {
      return(NULL)
    }
  }

  # change list to dataframe:
  lambda_local2 <- data.frame(lambda_local[[1]])
  names(lambda_local2) <- paste0("X", 1:ncol(lambda_local2))  #str_c("X", 1:ncol(lambda_local2))
  rows_without_NA <- which(!apply(Y, 1, anyNA))
  lambda_local2 <- lambda_local2 %>% mutate(group = 1, id = rows_without_NA[(rows_without_NA %in% which(g == 1))])


  # fill up with "rows_with_NA" of this group:
  # lambda_group can not be defined for these rows, so put to 0
  lambda_local2 <- lambda_local2 %>%
    rbind(add_lambdas_from_NArows(lambda_local2, 1)) %>%
    arrange(.data$id)

  if (S > 1) {
    for (i in 2:S) {
      df_to_add <- data.frame(lambda_local[[i]])
      names(df_to_add) <- paste0("X", 1:ncol(df_to_add))

      df_to_add <- df_to_add %>% mutate(group = i, id = rows_without_NA[(rows_without_NA %in% which(g == i))])

      df_to_add <- df_to_add %>%
        rbind(add_lambdas_from_NArows(df_to_add, i)) %>%
        arrange(.data$id)

      lambda_local2 <- lambda_local2 %>%
        bind_rows(df_to_add)
    }
  }



  # get rid of NA's when using update1 or update2
  # if (UPDATE1 | UPDATE2) {
  #   if (anyNA(lambda_local2)) {
  #     lambda_local2[is.na(lambda_local2)] <- 0
  #   }
  # }

  # add rows to lambda_local2 for individuals of class zero (=individuals that are too far from all groups)
  # put those values to zero
  if (sum(g == 0) > 0) {
    indices_class_zero <- which(g == 0)
    for (ind in indices_class_zero) {
      extrarow <- c(
        rep(0, max(kg, na.rm = TRUE)),
        0, # this is group membership
        ind # this is id
      )
      lambda_local2 <- rbind(lambda_local2, extrarow)
    }
  }
  stopifnot(nrow(lambda_local2) == NN)

  return(lambda_local2)
}


#' Function which is used to have a dataframe (called "grid") with data (individualindex, timeindex, XT and LF) available.
#'
#' It is used in iterate().
#' @param grid dataframe containing values for X*beta_est and LF (product of common factor and its loadings)
#' @param beta_est estimated values of beta
#' @param limit_est_groups_heterogroups maximum amount of groups that can be estimated when method_estimate_beta is set to "group"
#' @inheritParams estimate_beta
#' @return data.frame
grid_add_variables <- function(grid, Y, X, beta_est, g, lambda, comfactor, method_estimate_beta,
                               vars_est,
                               S,
                               limit_est_groups_heterogroups = 15) {
  NN <- nrow(Y)
  TT <- ncol(Y)
  if( !is.na(X[1]) & !is.null(X[1]) ) {
    vars <- dim(X)[3]
  } else {
    vars <- 0
  }
  if (vars > 0) {
    # for homogeneous beta_est (1 -> 4 at this moment), we only need 1 column as all columns are the same
    if (method_estimate_beta == "homogeneous") {
      beta_used <- as.matrix(beta_est[, 1])
      # calculate matrix multiplications outside the OF-functions to speed up:
      grid$XBETA <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_used))
    } else {
      if (method_estimate_beta == "group") {
        stopifnot((S >= 0 & S < limit_est_groups_heterogroups)) # code exists up to 14 groups

        # for each group: define XT
        if (S > 0) grid$XBETA1 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 1]))
        if (S > 1) grid$XBETA2 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 2]))
        if (S > 2) grid$XBETA3 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 3]))
        if (S > 3) grid$XBETA4 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 4]))
        if (S > 4) grid$XBETA5 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 5]))
        if (S > 5) grid$XBETA6 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 6]))
        if (S > 6) grid$XBETA7 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 7]))
        if (S > 7) grid$XBETA8 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 8]))
        if (S > 8) grid$XBETA9 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 9]))
        if (S > 9) grid$XBETA10 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 10]))
        if (S > 10) grid$XBETA11 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 11]))
        if (S > 11) grid$XBETA12 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 12]))
        if (S > 12) grid$XBETA13 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 13]))
        if (S > 13) grid$XBETA14 <- (apply(grid, 1, function(x) c(1, X[x[1], x[2], ]) %*% beta_est[, 14]))
      }
      if (method_estimate_beta == "individual") {
        grid$XBETA <- c(calculate_XB_estimated(X, beta_est, g, vars_est, method_estimate_beta, ncol(Y)))
      }
    }
  } else {
    if ((method_estimate_beta == "homogeneous") | (method_estimate_beta == "individual")) {
      grid$XBETA <- 0
    }
    if (method_estimate_beta == "group") {
      if (S > 0) grid$XBETA1 <- 0
      if (S > 1) grid$XBETA2 <- 0
      if (S > 2) grid$XBETA3 <- 0
      if (S > 3) grid$XBETA4 <- 0
      if (S > 4) grid$XBETA5 <- 0
      if (S > 5) grid$XBETA6 <- 0
      if (S > 6) grid$XBETA7 <- 0
      if (S > 7) grid$XBETA8 <- 0
      if (S > 8) grid$XBETA9 <- 0
      if (S > 9) grid$XBETA10 <- 0
      if (S > 10) grid$XBETA11 <- 0
      if (S > 11) grid$XBETA12 <- 0
      if (S > 12) grid$XBETA13 <- 0
      if (S > 13) grid$XBETA14 <- 0
    }
  }
  LF <- t(lambda) %*% comfactor
  grid$LF <- c(LF)

  return(grid)
}
















#' Helpfunction to shorten code: are common factors being estimated.
#'
#' @param k number of common factors
#' @return numeric: 0 or 1
do_we_estimate_common_factors <- function(k) {
  if (k == 0) {
    a <- 0
  } else {
    a <- 1
  }
  return(a)
}
#' Helpfunction to shorten code: are group factors being estimated.
#'
#' @param kg number of group factors to be estimated
#' @return numeric: 0 or 1
do_we_estimate_group_factors <- function(kg) {
  if (mean(kg, na.rm = TRUE) == 0) {
    b <- 0
  } else {
    b <- 1
  }
  return(b)
}


#' Calculates the error term Y - X*beta_est - LF - LgFg.
#' @param g Vector with estimated group membership for all individuals
#' @param no_common_factorstructure if there is a common factorstructure being estimated
#' @param no_group_factorstructure if there is a group factorstructure being estimated
#' @inheritParams grid_add_variables
#' @inheritParams estimate_beta
#' @inheritParams initialise_beta
#' @importFrom stringr str_detect
#' @importFrom rlang .data
#' @return NxT matrix
#' @examples
#' X <- X_dgp3
#' Y <- Y_dgp3
#' # Set estimations for group factors and its loadings, and group membership
#' #  to the true value for this example.
#' lambda_group <- lambda_group_true_dgp3
#' factor_group <- factor_group_true_dgp3
#' g <- g_true_dgp3
#' set.seed(1)
#' beta_est <- matrix(rnorm(nrow(Y) * 4), ncol = nrow(Y)) #random values for beta
#' comfactor <- matrix(0, ncol = ncol(Y))
#' lambda  <- matrix(0, ncol = nrow(Y))
#' calculate_error_term(Y, X, beta_est, g, factor_group, lambda_group, comfactor, lambda,
#'   3, 0, c(3, 3, 3))
#' @export
calculate_error_term <- function(Y, X, beta_est, g, factor_group, lambda_group, comfactor, lambda,
                                 S,
                                 k,
                                 kg,
                                 method_estimate_beta = "individual",
                                 no_common_factorstructure = FALSE, no_group_factorstructure = FALSE) {
  stopifnot(table(lambda_group$group) == table(g)) # in this case lambda_group would need to be updated
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  NN <- nrow(Y)
  TT <- ncol(Y)

  u <- matrix(NA, nrow = NN, ncol = TT)
  e <- matrix(NA, nrow = NN, ncol = TT)
  lf <- t(lambda) %*% comfactor

  if (vars_est > 0) {
    if ((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) {
      xt <- sapply(
        1:NN,
        function(y) sapply(1:TT, function(x) c(1, X[y, x, ]) %*% beta_est[, g[y]])
      )
    }
    if (method_estimate_beta == "individual") {
      xt <- t(calculate_XB_estimated(X, beta_est, g, vars_est, method_estimate_beta, ncol(Y))) # TxN matrix
    }
  } else {
    xt <- matrix(0, nrow = TT, ncol = NN) # TxN
  }

  lf_group <- list()
  group_membership <- list()
  for (gr in 1:S) {
    LGclean <- lambda_group %>%
      arrange(.data$id) %>%
      dplyr::filter(.data$group == gr)

    # replaced dplyr::select because of using ".data" in combination with "starts_with"
    # %>%  dplyr::select(starts_with("X")))

    selectcolumns <- which(str_detect(colnames(LGclean), "X"))
    LGclean <- as.matrix(LGclean[, selectcolumns])

    lf_group[[gr]] <- LGclean[, 1:kg[gr]] %*% factor_group[[gr]]
    group_membership[[gr]] <- data.frame(as.matrix(lambda_group %>% arrange(.data$id) %>%
      dplyr::filter(.data$group == gr) %>%
      dplyr::select(.data$group, .data$id)))
  }

  for (i in 1:NN) {
    if (g[i] != 0) {
      lf_group_i <- lf_group[[g[i]]]
      index <- which(group_membership[[g[i]]]$id == i)
    } else {
      lf_group_i <- NA
      index <- NA
    }
    for (t in 1:TT) {
      u[i, t] <- Y[i, t] - xt[t, i]

      a <- do_we_estimate_common_factors(k)
      b <- do_we_estimate_group_factors(kg)

      part2 <- a * lf[i, t]
      part3 <- ifelse(g[i] != 0,
        b * lf_group_i[index, t],
        0
      )
      if (no_common_factorstructure) part2 <- 0
      if (no_group_factorstructure) part3 <- 0
      e[i, t] <- u[i, t] - part2 - part3
    }
  }
  return(e)
}





#' Calculates sum of squared errors, divided by NT
#'
#' @param e matrix with error terms
#' @param NN N
#' @param TT T
#' @return numeric
#' @examples
#' Y <- Y_dgp3
#' set.seed(1)
#' e <- matrix(rnorm(nrow(Y) * ncol(Y)), nrow = nrow(Y))
#' calculate_sigma2(e)
#' @export
calculate_sigma2 <- function(e, NN = nrow(e), TT = ncol(e)) {
  if (anyNA(e)) {
    warning("There are NA's in e (calculate_sigma2(). This should not happen. ")
  }
  sigma2 <- sum(e * e) / (NN * TT)
  return(sigma2)
}



#' Function to calculate the norm of a matrix.
#'
#' @param mat input matrix
#' @return numeric
matrixnorm <- function(mat) {
  return(sqrt(sum(diag(t(mat) %*% mat))))
}



#' Function to calculate the first term of PIC (panel information criterium)
#'
#' This is used in calculate_PIC()
#' @param e NxT matrix with the error terms
#' @inheritParams initialise_beta
#' @importFrom robustbase Mpsi
#' @return numeric
calculate_PIC_term1 <- function(e, robust) {
  if (robust) {
    # This replaces the classical sum(z^2) by sum(rho(z)) with rho the bisquare function.

    # NOTE: dividing by mad(e) scales the errors of configurations with bad estimations to the results of good estimations
    #  -> this leads to wrong estimation of S & k
    rho_e_scaled <- Mpsi((e - median(e)), cc = 4.685, psi = "bisquare", deriv = -1)
    term1_new <- sum(rho_e_scaled) / (nrow(e) * ncol(e))
  } else {
    term1_new <- sum(e^2) / (nrow(e) * ncol(e))
  }

  return(term1_new)
}



#' Helpfunction. Calculates part of the 4th term of the PIC.
#'
#' @param TT length of time series
#' @param Nj number of time series in group j
#' @return numeric
calculate_TN_factor <- function(TT, Nj) {
  return((TT + Nj) / (TT * Nj) * log(TT * Nj))
}

#' Function to determine PIC (panel information criterium)
#'
#' This depends on kappa1 -> kappaN, through p (=number of nonzero elements of beta_est).
#' The parameter 'sigma2' is the non-robust sigma2. As it is only used in term 2 to 4, it does not actually matter what its value is (needs to be > 0).
#' It could be set to 1 as well.
#' @param C determines relative contribution of the penalty terms compared to the estimation error term
#' @inheritParams estimate_beta
#' @inheritParams OF_vectorized3
#' @param NN number of time series
#' @param TT length of time series
#' @param e2 NxT matrix with error terms
#' @param sigma2 scalar: sum of squared error terms, scaled by NT
#' @param choice_pic indicates which PIC to use to estimate the number of groups and factors:
#' options are "pic2017" (uses the PIC of \insertCite{Ando2017;textual}{RCTS}; works better for large N),
#' "pic2016" (\insertCite{Ando2016;textual}{RCTS}; works better for large T) weighs the fourth term with an extra factor relative to the size of the groups,
#' and "pic2022" which shrinks the NT-space where the number of groups and factors would be over- or underestimated compared to pic2016 and pic2017.
#' @return numeric
#' @examples
#' set.seed(1)
#' NN <- 30
#' TT <- 10
#' e <- matrix(rnorm(NN * TT), nrow = NN)
#' beta_est <- matrix(rnorm(NN * 4), ncol = NN) #random values for beta
#' g <- round(runif(NN, 1, 3))
#' calculate_PIC(0.51, TRUE, 3, 0, c(3, 3, 3), e, e^2/(NN*TT), NN, TT, "individual", beta_est, g, 3)
#' @export
calculate_PIC <- function(C, robust, S, k, kg, e2, sigma2,
                          NN,
                          TT,
                          method_estimate_beta,
                          beta_est, g,
                          vars_est,
                          choice_pic = "pic2022") {
  term1 <- calculate_PIC_term1(e2, robust)

  if (vars_est > 0) {
    if ((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) {
      # p is number of nonzero elements of beta_est; we need the sum of all p's (including the intercept)
      p_sum <- sum(sapply(1:NN, function(x) sum(beta_est[, g[x]] != 0)))
    } else {
      p_sum <- sum(sapply(1:NN, function(x) sum(beta_est[, x] != 0)))
    }
  } else {
    p_sum <- 0
  }

  term2 <- (C / NN) * sigma2 * log(TT) * p_sum

  # term3: penalty on number of common factors
  term3 <- C * k * sigma2 * (TT + NN) / (TT * NN) * log(TT * NN)

  # term4: penalty on number of groupfactors
  term4 <- 0
  for (j in 1:S) {
    Nj <- sum(g == j)
    TN_factor <- calculate_TN_factor(TT, Nj)
    temp <- C * kg[j] * sigma2 * TN_factor

    if (choice_pic == "pic2017") {
      # (from Ando/Bai2017-paper)
      term4 <- term4 + C * kg[j] * sigma2 * TN_factor
    }
    if (choice_pic == "pic2016") {
      # (from Ando/Bai2016-paper")
      term4 <- term4 + (C * kg[j] * sigma2 * (Nj / NN) * TN_factor) # =weigh with relative size of groups
    }
    if (choice_pic == "pic2022") {
      # note: any finite param_pic2022 > 0 will do. The larger it is, the smaller the problematic NT-region becomes.
      # 1e308 is the max for floating point stuff and should never be reached in a multiplication
      # -> param_pic2022 lower: 1e60 is quite conservative
      param_pic2022 <- 1e60
      term4 <- term4 + (C * kg[j] * sigma2 * (Nj / NN) * (TT + Nj) / (TT * Nj) * log(log((TT * Nj) * param_pic2022)))
    }
  }

  return(term1 + term2 + term3 + term4)
}







#' Calculates the product of X*beta_true .
#'
#' @inheritParams estimate_beta
#' @inheritParams generate_Y
#' @param g Vector with estimated group membership for all individuals
#' @param g_true true group membership
#' @return Returns a NxT matrix (if method_estimate_beta == "individual"), and otherwise NA.
calculate_XB_true <- function(X, beta_true, g, g_true, method_estimate_beta) {
  NN <- dim(X)[1]
  TT <- dim(X)[2]
  vars <- dim(X)[3]
  if (method_estimate_beta == "homogeneous") { # This is only relevant for DGP5 (BramatiCroux)
    XB_true <- beta_true[1, ] + X[, , 1] * beta_true[1, ]
    stopifnot(vars == 1)
  }
  if (method_estimate_beta == "group") {
    if (vars > 0) {
      if (!is.na(g[1])) {
        XB_true <- t(sapply(
          1:NN,
          function(x) matrix(cbind(1, X[x, , ]) %*% beta_true[, g[x]], nrow = 1)
        ))
      } else {
        XB_true <- NA
      }
    } else {
      XB_true <- NA
    }
  }
  if (method_estimate_beta == "individual") {
    if (vars > 0) {
      XB_true <- (t(sapply(
        1:NN,
        function(y) sapply(1:TT, function(x) c(1, X[y, x, 1:vars]) %*% beta_true[, g_true][, y])
      )))
    } else {
      XB_true <- NA
    }
  }

  return(XB_true)
}

#' When running the algorithm with a different number of observable variables then the number that is available, reformat X. (Mainly used for testing)
#'
#' @inheritParams create_true_beta
#' @inheritParams generate_Y
#' @inheritParams initialise_beta
#' @param vars_est number of available observed variables for which a coefficient will be estimated
#' @return Returns a 3D-array. If vars_est is set to 0, it returns NA.
adapt_X_estimating_less_variables <- function(X, vars_est) {
  vars <- dim(X)[3]
  # if vars_est < vars, then the obsolete rows in beta_est are already erased -> do the same in X
  if (vars_est < vars) {
    if (vars_est > 0) {
      X <- X[, , 1:vars_est]
    } else {
      X <- NA
    }
  }
  return(X)
}

#' Calculates (the estimated value of) the matrix X*beta_est.
#'
#' @inheritParams estimate_beta
#' @inheritParams OF_vectorized3
#' @return Returns a NxT matrix. If vars_est is set to 0, it returns NA.
calculate_XB_estimated <- function(X, beta_est, g,
                                   vars_est,
                                   method_estimate_beta, TT) {
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars <- dim(X)[3] # do before adapt_X_estimating_less_variables()
    # if vars_est < vars, then the obsolete rows in beta_est are already erased -> now do the same in X
    X <- adapt_X_estimating_less_variables(X, vars_est)
    NN <- dim(X)[1]
  } else {
    NN <- length(g)
    return(matrix(0, nrow = NN, ncol = TT))
  }

  if (vars > 0 & !is.na(X[1])) {
    if (method_estimate_beta == "homogeneous") { # currently only designed for DGP 5
      XT_estimated <- beta_est[1] + X[, , 1] * beta_est[1]
      stopifnot(vars == 1)
    }
    if (method_estimate_beta == "group") {
      if (vars_est > 0) {
        XT_estimated <- t(sapply(
          1:NN,
          function(x) matrix(cbind(1, X[x, , ]) %*% beta_est[, g[x]], nrow = 1)
        ))
      } else {
        XT_estimated <- NA
      }
    }
    if (method_estimate_beta == "individual") {
      if (vars_est > 0) {
        XT_estimated <- t(sapply(
          1:NN,
          function(x) matrix(cbind(1, X[x, , ]) %*% beta_est[, x], nrow = 1)
        ))
      } else {
        XT_estimated <- NA
      }
    }
  } else {
    XT_estimated <- matrix(0, nrow = NN, ncol = TT)
  }
  return(XT_estimated)
}


#' Calculate the true groupfactorstructure.
#'
#' @param lgt true group factor loadings
#' @param fgt true group factors
#' @param S_true true number of groups
#' @param kg_true true number of group factors for each group
#' @param k_true true number of common factors
# @param using_bramaticroux parameter to indicate that we are using data generated with dgp 5
#' @inheritParams estimate_beta
#' @inheritParams generate_Y
#' @param dgp1_AB_local gives information about which DGP we use; TRUE of FALSE
#' @return list with NjxT matrices
calculate_FL_group_true <- function(lgt, fgt, g_true, NN, TT,
                                    S_true,
                                    k_true,
                                    kg_true,
                                    num_factors_may_vary = TRUE,
                                    dgp1_AB_local = FALSE) {
  if (!dgp1_AB_local) { # not using for DGP 1, because true values are not available given current coding
    temp <- calculate_lgfg(
      lgt, fgt,
      S_true,
      k_true,
      kg_true,
      num_factors_may_vary,
      NN, TT
    )
    FL_group_true <- t(lapply(1:NN, function(x) temp[[g_true[x]]][x, ]) %>% unlist() %>% matrix(nrow = TT))
  } else {
    FL_group_true <- NA
  }
  return(FL_group_true)
}


#' Returns the estimated groupfactorstructure.
#'
#' @param fg estimated group factors
#' @param lg loadings of estimated group factors
#' @inheritParams calculate_FL_group_true
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @return list with NjxT matrices
calculate_FL_group_estimated <- function(lg, fg, g,
                                         NN,
                                         TT,
                                         S, k, kg,
                                         num_factors_may_vary = TRUE) {
  temp <- calculate_lgfg(lg, fg, S, k, kg, num_factors_may_vary, NN, TT)
  # helpfunction to return the groupfactorstructure of x. For elements in class zero, for which no factors are estimated, this is set to 0.
  FL_helpf <- function(x, g) {
    if (g[x] != 0) {
      return(temp[[g[x]]][x, ])
    } else {
      return(rep(0, TT))
    }
  }
  suppressWarnings( # the condition has length > 1 and only the first element will be used
    if (!is.na(temp)) {
      FL_group_est <- lapply(1:NN, function(x) FL_helpf(x, g))
      FL_group_est <- t(matrix(unlist(FL_group_est), nrow = TT))
    } else {
      FL_group_est <- NA
    }
  )
  return(FL_group_est)
}

#' Function to calculate the mean squared error of beta_est.
#'
#' For DGP 1 & 2: When the true number of variables in X is not equal to the standard of 3 it currently returns NA.
#' @param beta_est estimated values of beta
#' @inheritParams generate_Y
#' @param without_intercept TRUE of FALSE: to remove the intercept in the calculation of the MSE
#' @inheritParams estimate_beta
#' @inheritParams calculate_FL_group_true
#' @return numeric, or NA if the true number of variables is not equal to the standard of 3
#' @examples
#' set.seed(1)
#' beta_est <- matrix(rnorm(30 * 4), ncol = 30) #random values for beta
#' beta_true <- matrix(rnorm(4 * 3), nrow = 4)
#' g_true <- round(runif(30, 1,3)) #random values for true group membership
#' calculate_mse_beta(beta_est, beta_true, 30, 10, g_true, "individual")
#' @export
calculate_mse_beta <- function(beta_est, beta_true, NN, TT, g_true, method_estimate_beta,
                               # number_of_variables,
                               without_intercept = FALSE,
                               special_case_dgp1 = FALSE) {
  if (method_estimate_beta == "homogeneous") { # relevant for DGP05 (Bramati-Croux)
    mse <- mean((beta_est - beta_true)^2)
    if (without_intercept) mse <- mean((beta_est[-1, ] - beta_true[-1, ])^2)
    return(mse)
  }

  if (method_estimate_beta == "group") {
    pre_mse <- mse_heterogeneous_groups(without_intercept, special_case_dgp1, TT)
  }

  if (method_estimate_beta == "individual") { # Default case; In case of DGP 1 & 2 it returns NA when number of variables is not equal to 3.
    if (without_intercept) {
      if (special_case_dgp1) { # calculate MSE without the normal, and without the de facto intercept
        beta_est <- matrix(beta_est[c(-1, -2), ], ncol = NN)
        beta_true <- beta_true[c(-1, -2), ]
      } else { # calculate MSE without the normal intercept
        beta_est <- beta_est[-1, ]
        beta_true <- beta_true[-1, ]
      }
    }


    pre_mse <- rep(NA, NN)
    for (i in 1:NN) {
      # (note: length of this vector depends on without_intercept and special_case_dgp1)
      afw <- beta_est[, i] - beta_true[, g_true[i]]

      pre_mse[i] <- mean(afw^2) # this is (beta_est - beta_true)^2 (mean() goes over the number of variables)
    }
  }

  return(mean(pre_mse)) # this is E[(beta_est - beta_true)^2]
}

#' Helpfunction in calculate_mse_beta(), when method_estimate_beta == "group" (beta is then estimated for each group separately).
#'
#' @param without_intercept boolean to remove the intercept in the calculation
#' @inheritParams calculate_mse_beta
#' @inheritParams initialise_beta
#' @param TT length of time series
#' @param g vector with estimated group membership for all individuals
#' @return numeric vector
mse_heterogeneous_groups <- function(beta_est, beta_true, TT, g, g_true, without_intercept, special_case_dgp1) {

  # 1. The order of the estimated groups is not necessarily the same as the order of the true groups. -> possible necessary to permutate g to get the MSE?
  #-> this is countered by the order in the object beta_est
  #-> we do not need permutation here

  # if(sd(v1) != 0) { #when extra noise is added to the DGP (this is a non-standard case)
  #   message("mse_heterogeneous_groups(): not implemented when v1 contains values (=extra noise in dgp)")
  #   return(NA)
  # }

  if (without_intercept) {
    if (special_case_dgp1) { #->calculate MSE without true and without de facto intercept of DGP1
      beta_est <- beta_est[c(-1, -2), ]
      beta_true <- beta_true[c(-1, -2), ]
    } else { # calculate MSE without true intercept
      beta_est <- beta_est[-1, ]
      beta_true <- beta_true[-1, ]
    }
  }

  pre_mse <- rep(NA, TT)
  for (i in 1:TT) {
    afw <- beta_est[, g[i]] - (beta_true[, g_true[i]])
    pre_mse[i] <- mean(afw^2)
  }
  return(pre_mse)
}



#' Calculates VC², to determine the stability of the found number of groups and factors over the subsamples.
#'
#' VC² depends on C (this is the scale parameter in PIC). When VC² is equal to zero, the found number of groups and factors
#' are the same over the subsamples.
#' @param rcj dataframe containg the numer of groupfactors for all candidate C's and all subsamples
#' @param rc dataframe containg the numer of common factors for all candidate C's and all subsamples
#' @param C_candidates candidates for C (parameter in PIC)
#' @param indices_subset all indices of the subsets
#' @param Smax maximum allowed number of estimated groups
#' @inheritParams calculate_lambda_group
#' @inheritParams kg_candidates_expand
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @return numeric vector with the VC2-value for each candidate C
#' @examples
#' rcj <- data.frame(X1 = rep("3_3_3", 5), X2 = rep("3_2_1", 5))
#' rc <- data.frame(X1 = rep(1, 5), X2 = rep(0, 5))
#' calculate_VCsquared(rcj, rc, 1:5, 0:1, 3)
#' @export
calculate_VCsquared <- function(rcj, rc, C_candidates, indices_subset,
                                Smax,
                                # UPDATE1 = FALSE, UPDATE2 = FALSE,
                                limit_est_groups = 20) {
  VC_squared <- rep(NA, length(C_candidates))

  # vector with max number of groups
  Smax_local <- rep(Smax, length(C_candidates))

  number_subsets <- length(indices_subset) # number of subsamples used

  for (C_local in C_candidates) {
    C_index <- which(C_candidates == C_local)
    part1 <- 0
    for (i in 1:number_subsets) {
      part1 <- part1 + (rc[C_index, i] - sum(rc[C_index, ]) / number_subsets)^2
    }
    part1 <- part1 / number_subsets

    part2 <- 0
    for (j in 1:Smax_local[C_index]) {
      if (j <= limit_est_groups) {
        # if (j <= limit_est_groups | UPDATE1 | UPDATE2) {
        part2_part <- 0
        for (aA in 1:number_subsets) {
          temp <- 0
          for (bA in 1:number_subsets) {
            if (j <= length(unlist(str_split(rcj[C_index, bA], "_")))) {
              temp <- sum(temp, as.numeric(
                gsub("NA", NA, purrr::map(str_split(rcj[C_index, bA], "_"), j)) # add gsub() to evade warnings "NAs introduced by coercion"
              ), na.rm = TRUE)
            }
          }
          if (j <= length(unlist(str_split(rcj[C_index, aA], "_")))) {
            rcjNT <- as.numeric(
              gsub("NA", NA, purrr::map(str_split(rcj[C_index, aA], "_"), j))
            )
          } else {
            rcjNT <- 0
          }
          if (is.na(rcjNT)) rcjNT <- 0
          part2_part <- part2_part + (rcjNT - temp / number_subsets)^2
        }
        part2_part <- part2_part / number_subsets
        part2 <- part2 + part2_part
      } else {
        stop("error in calculating VC2")
      }
    }
    VC_squared[C_index] <- part1 + part2
  }
  return(unlist(VC_squared))
}









#' Wrapper around lmrob.
#'
#' Desgined to make sure the following error does not happen anymore:
#' Error in if (init$scale == 0)  : missing value where TRUE/FALSE needed.
#' KS2014 is the recommended setting (use "nosetting = FALSE").
#' @param parameter_y dependent variable in regression
#' @param parameter_x independent variables in regression
#' @param nointercept if TRUE it performs regression without an intercept
#' @param nosetting option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
# @param newkmax sets parameter k.max in lmrob.control
#' @importFrom robustbase lmrob
#' @importFrom robustbase lmrob.control
#' @return An object of class lmrob. If something went wrong it returns an object of class error.
LMROB <- function(parameter_y, parameter_x, nointercept = FALSE, nosetting = FALSE) {
  if (is.na(parameter_x)[1]) { # when there are no independent variables
    result2 <- tryCatch(
      {
        if (nosetting) {
          result <- lmrob(parameter_y ~ 1)
          return(result)
        } else {
          result <- lmrob(parameter_y ~ 1, setting = "KS2014")
          return(result)
        }
      },
      error = function(e) {
        print(e)
        print("error, therefore use lmrob without 'setting'")
        result <- lmrob(parameter_y ~ 1)
        return(result)
      },
      finally = {}
    )
  } else {
    if (nosetting) {
      if (nointercept) {
        result2 <- lmrob(parameter_y ~ parameter_x + 0)
      } else {
        result2 <- lmrob(parameter_y ~ parameter_x)
      }
    } else {
      if (nointercept) {
        result2 <- lmrob(parameter_y ~ parameter_x + 0, setting = "KS2014")
      } else {
        # sometimes Error in if (init$scale == 0) happens. In that case: run without setting = KS2014.

        result2 <- tryCatch(
          {
            # note: I cannot use both "setting" and "control" (to set k.max higher),
            # therefore use the parameters of 'setting="KS2014"' listed in documentation of lmrob.control
            # The option setting="KS2011" alters the default arguments. They are changed to method = "SMDM", psi = "lqq", max.it = 500, k.max = 2000, cov = ".vcov.w".
            # The defaults of all the remaining arguments are not changed.
            # The option setting="KS2014" builds upon setting="KS2011". More arguments are changed to best.r.s = 20, k.fast.s = 2, nResample = 1000.

            #I also increased k.max and maxit.scale
            lmrobcontrol_ks2014 <- lmrob.control(method = "SMDM", psi = "lqq", max.it = 500, k.max = 3000, cov = ".vcov.w", best.r.s = 20, k.fast.s = 2, nResample = 1000, maxit.scale = 600)
            result <- lmrob(parameter_y ~ parameter_x, control = lmrobcontrol_ks2014)
            return(result)
          },
          # possible warnings are:
          # "S refinements did not converge (to refine.tol=1e-07) in 2000 (= k.max) steps"
          #-> solution: i increased k.max to 3000

          # find_scale() did not converge in 'maxit.scale' (= 200) iterations with tol=1e-10, last rel.diff=0
          #-> increased to 600

          # but note that if code reaches the warning, there is no access to "result", and i won't calculate it again -> comment out the warning part
          # (this thus only lists the warnings after the algorithm has run, and not anymore per iteration)

          # warning = function(w) {
          #   message(w)
          #   message(paste("\nlength of input:",length(parameter_y)))
          #   return(result)
          # },
          error = function(e) {
            print(e)
            print("error, therefore use lmrob without 'setting'") # does this still occur?
            result <- lmrob(parameter_y ~ parameter_x)
            return(result)
          },
          finally = {}
        )
      }
    }
  }
  return(result2)
}






#' Calculates sigma2maxmodel
#'
#' Sigma2 is the sum of the squared errors, divided by NT. We need the sigma2 of the maxmodel to use (in term 2,3,4 of the PIC) instead of the configuration-dependent sigma2. (See paper AndoBai 2016).
#' sigma2_max_model could actually be set to 1 as well, as it can be absorbed in parameter C of the PIC.
#' @param e NxT-matrix containing the estimated error term
#' @param kg vector with the estimated number of group specific factors for each group
#' @param kg_max scalar: maximum allowed number of estimated factors for any group
#' @param S estimated number of groups
#' @param S_cand vector with candidate values for the number of groups
#' @param k estimated number of common factors
#' @param k_cand vector with candidate value for the number of common factors
#' @return numeric
calculate_sigma2maxmodel <- function(e, kg_max,
                                     S, S_cand,
                                     kg,
                                     k, k_cand
                                     # UPDATE1 = FALSE, UPDATE2 = FALSE
) {


  if ( # !UPDATE1 & !UPDATE2 &
    mean(kg, na.rm = TRUE) == kg_max &
    S == max(S_cand) &
    k == max(k_cand)) {
    sigma2_max_model <- calculate_sigma2(e)
  }
  # if (UPDATE1 & UPDATE2) {
  #   # update number of group- and common factors during algorithm
  #   sigma2_max_model <- NA
  # }

  return(sigma2_max_model)
}

#' Adapts the object that contains PIC for all candidate C's and all subsamples with sigma2_max_model.
#'
#' The PIC is calculated with a sigma2 specific to the configuration (= number of groups and factors).
#' Because the method to estimate the number of groups and factors requires sigma2 to be equal over all configurations
#' (see proofs of different papers of Ando/Bai) we replace sigma2
#' by the sigma2 of the configuration with maximum number of groups and factors (this is the last one that was executed).
#' @param df contains PIC for all candidate C's and all subsamples
#' @param sigma2_max_model sigma2 of model with maximum number of groups and factors
#' @inheritParams calculate_lambda_group
#' @inheritParams add_configuration
#' @return data.frame of same size as df
#' @examples
#' set.seed(1)
#' df_pic <- data.frame(matrix(rnorm(4 * 50), nrow = 4)) #4 configuration / 50 candidate values for C
#' df_results <- data.frame(sigma2 = rnorm(4))
#' pic_sigma2 <- 3.945505
#' adapt_pic_with_sigma2maxmodel(df_pic, df_results, pic_sigma2)
#' @export
adapt_pic_with_sigma2maxmodel <- function(df, df_results, sigma2_max_model
                                          # UPDATE1 = FALSE, UPDATE2 = FALSE
) {
  if (is.null(df)) warning("Warning in adapt_allpic_with_sigma2maxmodel(): df is empty")
  # if (!UPDATE1 & !UPDATE2) {
  for (i in 1:nrow(df)) {
    sig2 <- as.numeric(df_results$sigma2[i])
    df[i, ] <- (unlist(df[i, ]) - sig2) / sig2 * sigma2_max_model + sig2
  }
  # }
  return(df)
}

#' Creates an instance of DGP 2, as defined in \insertCite{BoudtHeyndels2021;textual}{RCTS}.
#'
#' The default has 3 groups with each 3 group specific factors. Further it contains 0 common factors and 3 observed variables.
#' The output is a list where the first element is the simulated panel dataset (a dataframe with N (amount of time series) rows and T (length of time series) columns).
#' The second element contains the NxTxp array with the p observed variables. The third element contains the true group membership.
#' The fourth element contains the true beta's (this has p+1 rows and one column for each group).
#' The fifth element contains a list with the true group specific factors.
#' The sixth element contains a dataframe with N rows where each row contains the group specific factor loadings that corresponds to the group specific factors.
#' Further it contains the true group membership and an index (this corresponds to the rownumber in Y and X).
#' The seventh and eighth elements contain the true common factor(s) and its loadings respectively.
#' @importFrom Rdpack reprompt
#' @param N number of time series
#' @param TT length of time series
#' @param S_true true number of groups
#' @param vars number of available observed variables
#' @param k_true true number of common_factors
#' @param kg_true vector with the true number of group factors for each group
#' @return list
#' @examples
#' create_data_dgp2(30, 10)
#' @export
create_data_dgp2 <- function(N, TT, S_true = 3, vars = 3, k_true = 0, kg_true = c(3, 3, 3)) {
  # true group membership
  g_true <- ceiling(runif(N) * S_true)
  # true beta
  beta_true <- create_true_beta(vars, N, S_true, "heterogeneous_groups")

  # true group-specific factors and loadings
  temp <- generate_grouped_factorstructure(S_true, kg_true, TT, g_true)
  factor_group_true <- temp[[1]]
  lambda_group_true <- temp[[2]]
  rm(temp)

  # true common factor structure
  comfactor_true <- matrix(t(rnorm(TT * k_true)), nrow = k_true)
  lambda_true <- matrix(t(rnorm(N * k_true)), nrow = k_true)

  # make object X: gives random values to the variables
  X <- initialise_X(N, TT, vars)

  # errorterm
  eps <- matrix(rnorm(N * TT, sd = 1), nrow = N)

  Y <- generate_Y(
    N, TT,
    k_true,
    kg_true,
    g_true, beta_true, lambda_group_true, factor_group_true,
    lambda_true, comfactor_true, eps, X
  )
  return(list(Y, X, g_true, beta_true, factor_group_true, lambda_group_true, comfactor_true, lambda_true))
}

#' Selects a subsample of the time series, and of the length of the time series.
#' Based on this it returns a list with a subsample of Y, the corresponding subsample of X and of the true group membership and factorstructures if applicable.
#'
#' @param original_data list containing the true data: Y, X, g_true, beta_true, factor_group_true, lambda_group_true, comfactor_true, lambda_true
#' @param subset index of the subsample: this defines how many times stepsize_N is subtracted from the original N time series. Similar for stepsize_T.
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @return
#' Y, X, g_true, comfactor_true, lambda_true, factor_group_true, lambda_group_true, sampleN, sampleT
#' The output is a list where the first element is a subset of the panel dataset.
#' The second element contains a subsetted 3D-array with the p observed variables.
#' The third element contains the subsetted true group membership.
#' The fourth and fifth elements contain the subsetted true common factor(s) and its loadings respectively.
#' The sixth element contains a list with the subsetted true group specific factors.
#' The seventh element contains a dataframe where each row contains the group specific factor loadings that corresponds to the group specific factors.
#' The eighth and ninth element contain the indices of N and T respectively, which were used to create the subsets.
#' @examples
#' set.seed(1)
#' original_data <- create_data_dgp2(30, 10)
#' make_subsamples(original_data, 1)
#' @export
make_subsamples <- function(original_data,
                            subset,
                            verbose = TRUE) {
  #determine whether input data is simulated or real world data, based on size of input
  if(length(original_data) == 2) {
    use_real_world_data = TRUE
  } else {
    use_real_world_data = FALSE
  }
  Y <- original_data[[1]]
  N_fulldata <- nrow(Y)
  T_fulldata <- ncol(Y)
  stepsize_N <- round(N_fulldata / 10)
  stepsize_T <- round(T_fulldata / 30)

  X <- original_data[[2]]
  if(!use_real_world_data) {
    g_true <- original_data[[3]]
    factor_group_true <- original_data[[5]]
    lambda_group_true <- original_data[[6]]
    comfactor_true <- original_data[[7]]
    lambda_true <- original_data[[8]]
  }


  # define size of the subsample
  subN <- N_fulldata - subset * stepsize_N
  subT <- T_fulldata - subset * stepsize_T
  if (!(subN > 0 & subT > 0)) {
    stop("subN or subT < 0 -> stop")
  }

  # define which N and T are in the subsample
  sampleN <- sort(sample(1:N_fulldata, subN))
  sampleT <- sort(sample(1:T_fulldata, subT))

  # subsample of Y
  Y <- Y[sampleN, sampleT]

  # subsample of X
  if(!is.na(X[1]) & !is.null(X[1])) {
    if (dim(X)[3] > 0) {
      X_temp <- array(NA, dim = c(subN, subT, dim(X)[3])) # subT*T_FACTOR
      X_temp[, , ] <- X[sampleN, sampleT, ]
      X <- X_temp
      rm(X_temp)
    }
  }

  if (!use_real_world_data) {
    # subsample of true group membership
    g_true <- g_true[sampleN]

    S <- length(factor_group_true)
    k <- nrow(comfactor_true)
    kg <- unlist(lapply(factor_group_true, nrow))

    # subsample of the factors and their loadings
    if (k > 0) {
      comfactor_true <- matrix(comfactor_true[, sampleT], nrow = k)
      lambda_true <- matrix(lambda_true[, sampleN], nrow = k)
    }
    if (max(kg, na.rm = TRUE) > 0) {
      for (ii in 1:S) {
        factor_group_true[[ii]] <- matrix(factor_group_true[[ii]][, sampleT], nrow = kg[ii])
        lambda_group_true <- lambda_group_true %>% dplyr::filter(.data$id %in% sampleN)
      }
    }
  } else {
    g_true <- NA
    comfactor_true <- NA
    lambda_true <- NA
    factor_group_true <- NA
    lambda_group_true <- NA
  }
  return(list(Y, X, g_true, comfactor_true, lambda_true, factor_group_true, lambda_group_true, sampleN, sampleT))
}

#' Defines the object that will be used to define a initial clustering.
#'
#' This is a short version of define_object_for_initial_clustering() which only contains implementations for robust macropca case and classical case.
#' @inheritParams estimate_beta
# @param g Vector with group membership for all individuals
#' @param method_estimate_beta defines how beta is estimated. Default case is an estimated beta for each individual. Default value is "individual." Possible values are "homogeneous", "group" or "individual".
#' @param method_estimate_factors specifies the robust algorithm to estimate factors: default is "macro". The value is not used when robust is set to FALSE.
#' @param verbose when TRUE, it prints messages
#' @return matrix with N rows and 10 columns
define_object_for_initial_clustering_macropca <- function(robust, Y, k, kg, comfactor, method_estimate_beta = "individual", method_estimate_factors = "macro", verbose = FALSE) {
  TT <- ncol(Y)
  number_of_initial_factors <- 10
  stopifnot(method_estimate_beta == "individual")
  stopifnot(method_estimate_factors == "macro")

  if (robust) {
    if (method_estimate_factors == "macro") {
      ev <- robustpca(Y, number_of_initial_factors, verbose_robustpca = verbose) # these are the eigenvectors
      if (is.null(ev)) {
        # then macropca returned a wrong amount of factors
        message("define_object_for_initial_clustering_macropca() fails due to macoPCA returning wrong amount of factors")
      }
      factor_for_grouping <- sqrt(TT) * (ev[[1]])[, 1:number_of_initial_factors]
      lambda_for_grouping <- t(return_robust_lambdaobject(Y, NA, type = 4, NA, nrow(Y), k, kg, comfactor, factor_for_grouping, verbose))
    }
  } else {
    # classical
    ev <- eigen(t(Y) %*% Y)$vectors
    factor_for_grouping <- sqrt(TT) * (ev)[, 1:number_of_initial_factors]
    lambda_for_grouping <- t(factor_for_grouping) %*% t(Y) / TT
  }

  to_divide <- t(lambda_for_grouping)
  rm(factor_for_grouping, lambda_for_grouping)

  return(to_divide)
}

#' Function that clusters time series in a dataframe with kmeans (classical algorithm) or trimmed kmeans(robust algorithms).
#'
#' If a time series contains NA's a random cluster will be assigned to that time series.
#' @inheritParams define_object_for_initial_clustering_macropca
#' @inheritParams initialise_beta
# @param g vector with estimated group membership for all individuals
#' @param S the desired number of groups
#' @param max_percent_outliers_tkmeans the proportion of observations to be trimmed
#' @importFrom stats kmeans
#' @return numeric vector
#' @examples
#' Y <- Y_dgp3
#' comfactor <- matrix(0, nrow = ncol(Y))
#' initialise_clustering(TRUE, Y, 3, 0, c(3, 3, 3), comfactor)
#' @export
initialise_clustering <- function(robust, Y, S, k, kg, comfactor, max_percent_outliers_tkmeans = 0, verbose = FALSE) {
  df <- define_object_for_initial_clustering_macropca(robust, Y, k, kg, comfactor, verbose = verbose)
  # in the kmeans-call NA's are not allowed:
  # -> function handleNA() will drop time series with NA's.
  NN <- nrow(df)
  if (anyNA(df)) {
    temp <- handleNA(df)
    df <- temp[[1]]
    rows_with_NA <- temp[[2]]
    rows_without_NA <- which(!(1:NN %in% rows_with_NA))
  } else {
    rows_with_NA <- NA
    rows_without_NA <- 1:NN
  }

#Note: to limit dependencies on other packages -> use also kmeans for robust algorithm
  #(note: tkmeans also uses 50 initializations, but since the initial clustering of the RCTS-algorithm is not really influential, kmeans is fine enough)
  # if (robust) { # only use trimmed kmeans in the robust version
  #   counter <- 0
  #   Km <- tkmeans(df, S, max_percent_outliers_tkmeans)
  #
  #   # solve possible issue with empty clusters
  #   while (max_percent_outliers_tkmeans > 0 & length(table(Km$cluster)) != (S + 1) & counter < 5) {
  #     message("try again until there are no empty clusters (which means that there is no warning message")
  #     counter <- counter + 1
  #     Km <- tkmeans(df, S, max_percent_outliers_tkmeans)
  #   }
  #   # trimmed kmeans identifies outliers and puts them in zero (if parameter > 0) -> they would thus need to be replaced
  # } else {
    # use kmeans
    Km <- tryCatch(kmeans(df, S), error = function(e) e)
    if ("error" %in% class(Km)) {
      print(Km)
      # more cluster centers than distinct data points.
      warning("Km does not exist now: take as initial kmeans one group less + change 1 element to the empty group")
      counter <- 0
      while ("error" %in% class(Km)) {
        counter <- counter + 1
        Km <- tryCatch(kmeans(df, S - counter), error = function(e) e)
        if (counter > 50) {
          stop("infinite loop -> stop")
        }
      }


      indices_to_change_randomly <- ceiling(runif(counter, max = nrow(df)))
      for (tel in 0:(counter - 1)) {
        Km$cluster[indices_to_change_randomly[tel + 1]] <- S - tel
      }
    # }
  }

  g <- rep(NA, NN)
  g[rows_without_NA] <- Km$cluster
  # time series with NA do not get any result with kmeans (since they were deleted up front) -> assign a random initial group for those
  g[rows_with_NA] <- ceiling(runif(length(rows_with_NA)) * S)
  return(g)
}

#' Initialises the estimation of the common factors and their loadings.
#'
#' This is a short version of initialise_commonfactorstructure() which only contains implementations for the robust macropca case and the classical case.
#' @inheritParams define_object_for_initial_clustering_macropca
#' @inheritParams generate_Y
#' @inheritParams initialise_beta
#' @inheritParams estimate_beta
#' @param beta_est estimated values of beta
#' @param g Vector with estimated group membership for all individuals
#' @param k number of estimated common factors
#' @param kg vector with the number of estimated group specific factors
# @param vars_est number of available observed variables for which a coefficient will be estimated. As default it is equal to the number of available observed variables.
#' @return list: 1st element contains the common factor(s) and the second element contains the factor loadings
#' @examples
#' set.seed(1)
#' original_data <- create_data_dgp2(30, 20)
#' Y <- original_data[[1]]
#' X <- original_data[[2]]
#' g <- original_data[[3]]
#' beta_est <- matrix(rnorm(4 * ncol(Y)), nrow = 4)
#' initialise_commonfactorstructure_macropca(TRUE, Y, X, beta_est, g, NA, 0, c(3, 3, 3))
#' @export
initialise_commonfactorstructure_macropca <- function(robust, Y, X, beta_est, g, factor_group,
                                                      k, kg,
                                                      method_estimate_beta = "individual", method_estimate_factors = "macro",
                                                      verbose = FALSE) {
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  NN <- nrow(Y)
  TT <- ncol(Y)
  if ((method_estimate_beta == "individual")) {
    if (k == 0) {
      comfactor <- t(matrix(rep(0, TT)))
      lambda <- t(matrix(rep(0, NN)))
    } else {
      comfactor <- estimate_factor(robust, Y, X, beta_est, g, NA,
        k, kg,
        method_estimate_beta, method_estimate_factors,
        initialise = TRUE
      )[[1]][1:k, , drop = FALSE]
      lambda <- calculate_lambda(robust, Y, X, beta_est, comfactor, factor_group, g, NA,
        k, kg,
        method_estimate_beta, method_estimate_factors,
        verbose,
        initialise = TRUE
      )[1:k, , drop = FALSE]
    }
  }
  return(list(comfactor, lambda))
}

#' Wrapper around estimate_beta(), update_g(), and estimating the factorstructures.
#'
#' @inheritParams estimate_beta
#' @param S number of groups to estimate
#' @param k number of common factors to estimate
#' @param kg vector with length S. Each element contains the number of group specific factors to estimate.
#' @inheritParams define_object_for_initial_clustering_macropca
#' @param special_case_dgp1 TRUE or FALSE: whether data is generated from dgp1 and has the extra spread in group centers. Default is FALSE.
# @param vars_est number of variables for which a beta is estimated. Usually equal to the number of variables available in X.
#' @inheritParams update_g
#' @return list with
#' 1. estimated beta
#' 2. vector with group membership
#' 3. matrix with the common factor(s) (contains zero's if there are none estimated)
#' 4. loadings to the common factor(s)
#' 5. list with the group specific factors for each of the groups
#' 6. data.frame with loadings to the group specific factors augmented with group membership and id (to have the order of the time series)
#' 7. the value of the objective function
#' @examples
#' set.seed(1)
#' original_data <- create_data_dgp2(30, 10)
#' Y <- original_data[[1]]
#' X <- original_data[[2]]
#' g <- original_data[[3]]
#' beta_est <- matrix(rnorm(4 * ncol(Y)), nrow = 4)
#' factor_group <- original_data[[5]]
#' lambda_group <- original_data[[6]]
#' comfactor <- matrix(0, nrow = 1, ncol = ncol(Y))
#' lambda <- matrix(0, nrow = 1, ncol = nrow(Y))
#' iterate(TRUE, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor, 3, 0, c(3, 3, 3),
#'   verbose = FALSE)
#' @export
iterate <- function(robust, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor, S, k, kg,
                    method_estimate_beta = "individual", method_estimate_factors = "macro",
                    special_case_dgp1 = FALSE, verbose = FALSE) {
  if(!is.na(X[1]) & !is.null(X[1])) {
    vars_est <- dim(X)[3]
  } else {
    vars_est <- 0
  }
  if (verbose) message("update beta")
  beta_est <- estimate_beta(robust, Y, X, beta_est, g, lambda_group, factor_group, lambda, comfactor,
    method_estimate_beta,
    S, k, kg,
    vars_est,
    num_factors_may_vary = TRUE, special_case_dgp1 = special_case_dgp1
  )[[1]]
  if (verbose) message("update group membership")
  g <- update_g(
    robust, Y, X, beta_est, g,
    factor_group,
    lambda, comfactor,
    S, k, kg,
    vars_est,
    method_estimate_factors, method_estimate_beta
  )[[1]]

  if (verbose) message("update common factorstructure")
  lgfg_list <- calculate_lgfg(lambda_group, factor_group,
    S, k, kg,
    num_factors_may_vary = TRUE, NN = nrow(Y), TT = ncol(Y)
  )
  comfactor <- estimate_factor(robust, Y, X, beta_est, g, lgfg_list,
    k, kg,
    method_estimate_beta, method_estimate_factors,
    verbose = FALSE
  )[[1]]
  lambda <- calculate_lambda(
    robust, Y, X, beta_est, comfactor, factor_group, g, lgfg_list,
    k, kg,
    method_estimate_beta, method_estimate_factors
  )

  if (verbose) message("update group specific factorstructure")
  factor_group <- estimate_factor_group(
    robust, Y, X, beta_est, g, lambda, comfactor, factor_group,
    S, k, kg,
    method_estimate_beta, method_estimate_factors,
    verbose = verbose
  )
  if (verbose) message("factor_group is estimated")
  lambda_group <- calculate_lambda_group(
    robust, Y, X, beta_est, factor_group, g, lambda, comfactor,
    S, k, kg,
    method_estimate_beta, method_estimate_factors,
    verbose = verbose # , UPDATE1 = FALSE, UPDATE2 = FALSE
  )
  if (verbose) message("lambda_group is estimated")

  if (verbose) message("calculate objective function")
  grid <- expand.grid(1:nrow(Y), 1:ncol(Y))
  grid <- grid_add_variables(
    grid, Y, X, beta_est, g, lambda, comfactor, method_estimate_beta,
    vars_est,
    S
  )

  # value to minimize:
  if (robust) {
    # robust: sum of "rho-functioned" (in context of M-estimators) errors / NT
    value <- calculate_PIC_term1(
      calculate_error_term(Y, X, beta_est,
        g,
        factor_group,
        lambda_group,
        comfactor,
        lambda,
        S,
        k,
        kg,
        method_estimate_beta,
      ),
      robust
    ) * nrow(Y) * ncol(Y)
  } else {
    # classical: sum of squared errors / NT
    value <- OF_vectorized3(nrow(Y), ncol(Y), g, grid, Y,
      beta_est,
      lc = lambda, fc = comfactor,
      lg = lambda_group, fg = factor_group,
      S, k, kg,
      method_estimate_beta,
      num_factors_may_vary = TRUE
    )
  }
  if (verbose) print(value)
  return(list(beta_est, g, comfactor, lambda, factor_group, lambda_group, value))
}



#' Defines the convergence speed.
#'
#' @param iteration number of iteration
#' @param of objective function
# @param verbose when TRUE, it prints more information
#' @return numeric if iteration > 3, otherwise NA
#' @examples
#' get_convergence_speed(5, 10:1)
#' @export
get_convergence_speed <- function(iteration, of) { # , verbose = FALSE) {
  if (iteration > 3) {
    evo_min <- min(of[(iteration - 3):iteration], na.rm = T)
    evo_max <- max(of[(iteration - 3):iteration], na.rm = T)
    speed <- round((evo_max - evo_min), 4)
    # if(robust) {
    #   if(verbose) print(paste("Range:",round(evo_min, 2), " -> ", round(evo_max, 2)))
    #   if(verbose) message(paste("Speed (robust): ", speed))
    # } else {
    #   if(verbose) print(paste("Range:", round(evo_min), " -> ", round(evo_max)))
    #   if(verbose) message(paste("Speed (non-robust): ", speed))
    # }
  } else {
    speed <- NA
  }

  return(speed)
}

#' Checks the rules for stopping the algorithm, based on its convergence speed.
#'
#' @param iteration number of iteration
#' @param speed convergence speed
#' @param speedlimit if the convergence speed falls under this limit the algorithm stops
#' @param all_OF_values vector containing the values of the objective function from previous iterations
#' @param verbose if TRUE, more information is printed
#' @return logical
#' @examples
#' check_stopping_rules(4, 1.7, 5:1)
#' @export
check_stopping_rules <- function(iteration, speed, all_OF_values, speedlimit = 0.01, verbose = FALSE) {
  if (is.na(speed)) {
    return(FALSE)
  } else {
    #########################################################################
    # stop when range of OF is small enough during last 3 + current iteration
    if (speed < speedlimit) {
      if (verbose) print("stop (rule 1)")
      return(TRUE)
    }
    #################################################
    # stop when there are recurring periods in the OF
    ###########################################################
    # stop when the last XX steps do not give a better mean(OF)
    if (iteration > 14) {
      stoplimit <- 7
      OF_9_15 <- mean(all_OF_values[(iteration - stoplimit + 1):iteration], na.rm = T) # mean of objective function in iteration 9->15, 10->16,...
      OF_1_8 <- mean(all_OF_values[(iteration - (2 * stoplimit)):(iteration - stoplimit)], na.rm = T) # mean of objective function in iteration 1->8, 2->9,...
      if (OF_9_15 >= OF_1_8) {
        if (verbose) print("stop (rule 3)")
        return(TRUE)
      }
    }
    return(FALSE)
  }
}


#' Function that returns the set of combinations of groupfactors for which the algorithm needs to run.
#'
#' @param S number of groups
#' @param kg_min minimum value for number of group specific factors
#' @param kg_max minimum value for number of group specific factors
#' @param limit_est_groups maximum allowed number of groups that can be estimated
#' @return data.frame where each row contains a possible combination of group specific factors for each of the groups
kg_candidates_expand <- function(S, kg_min, kg_max, limit_est_groups = 20) {
  stopifnot(S <= limit_est_groups)
  if (S == 1) kg_cand <- expand.grid(kg_min:kg_max)
  if (S == 2) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max)
  if (S == 3) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 4) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 5) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 6) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 7) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 8) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 9) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 10) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 11) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 12) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 13) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 14) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 15) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 16) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 17) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 18) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 19) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S == 20) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  # if (S == 21) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  # if (S == 22) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  # if (S == 23) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  # if (S == 24) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  # if (S == 25) kg_cand <- expand.grid(kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max, kg_min:kg_max)
  if (S > limit_est_groups) {
    stop(paste("The number of estimated groups can be maximum", limit_est_groups, "."))
  }

  # When the number of groupfactors is allowed to vary among the groups the number of possible combinations increases very fast.
  # Since the order of the groupfactors does NOT matter for the result, (1,2,3) factors and (3,2,1) factors will give the same results.
  ####### (see "testen_volgorde_aantal_geschatte_factoren.R"),
  #    Therefore we can filter out some of the combinations.
  if (S >= 2) kg_cand <- kg_cand %>% dplyr::filter(.data$Var2 >= .data$Var1)
  if (S >= 3) kg_cand <- kg_cand %>% dplyr::filter(.data$Var3 >= .data$Var2)
  if (S >= 4) kg_cand <- kg_cand %>% dplyr::filter(.data$Var4 >= .data$Var3)
  if (S >= 5) kg_cand <- kg_cand %>% dplyr::filter(.data$Var5 >= .data$Var4)
  if (S >= 6) kg_cand <- kg_cand %>% dplyr::filter(.data$Var6 >= .data$Var5)
  if (S >= 7) kg_cand <- kg_cand %>% dplyr::filter(.data$Var7 >= .data$Var6)
  if (S >= 8) kg_cand <- kg_cand %>% dplyr::filter(.data$Var8 >= .data$Var7)
  if (S >= 9) kg_cand <- kg_cand %>% dplyr::filter(.data$Var9 >= .data$Var8)
  if (S >= 10) kg_cand <- kg_cand %>% dplyr::filter(.data$Var10 >= .data$Var9)
  if (S >= 11) kg_cand <- kg_cand %>% dplyr::filter(.data$Var11 >= .data$Var10)
  if (S >= 12) kg_cand <- kg_cand %>% dplyr::filter(.data$Var12 >= .data$Var11)
  if (S >= 13) kg_cand <- kg_cand %>% dplyr::filter(.data$Var13 >= .data$Var12)
  if (S >= 14) kg_cand <- kg_cand %>% dplyr::filter(.data$Var14 >= .data$Var13)
  if (S >= 15) kg_cand <- kg_cand %>% dplyr::filter(.data$Var15 >= .data$Var14)
  if (S >= 16) kg_cand <- kg_cand %>% dplyr::filter(.data$Var16 >= .data$Var15)
  if (S >= 17) kg_cand <- kg_cand %>% dplyr::filter(.data$Var17 >= .data$Var16)
  if (S >= 18) kg_cand <- kg_cand %>% dplyr::filter(.data$Var18 >= .data$Var17)
  if (S >= 19) kg_cand <- kg_cand %>% dplyr::filter(.data$Var19 >= .data$Var18)
  if (S >= 20) kg_cand <- kg_cand %>% dplyr::filter(.data$Var20 >= .data$Var19)

  if (S > limit_est_groups) {
    stop(paste("The number of estimated groups can be maximum", limit_est_groups, "."))
  }


  return(kg_cand)
}

#' Defines the set of combinations of group specific factors.
#'
#' @param S number of estimated groups
#' @param kg_min minimum value for number of group specific factors
#' @param kg_max minimum value for number of group specific factors
#' @param nfv logical; whether the number of group specific factors is allowed to change among the groups
#' @inheritParams kg_candidates_expand
#' @return Returns a data frame where each row contains the number of group specific factors for all the estimated groups.
#' The number of columns is set to 20 (the current maximum amount of group that can be estimated)
#' @examples
#' define_kg_candidates(3, 2, 4)
#' @export
define_kg_candidates <- function(S, kg_min, kg_max, nfv = TRUE, limit_est_groups = 20) {
  if (nfv) {
    kg_candidates <- kg_candidates_expand(S, kg_min, kg_max)
  } else {
    kg_candidates <- matrix(rep(kg_min:kg_max, S), ncol = S)
  }
  while (ncol(kg_candidates) < limit_est_groups) kg_candidates <- cbind(kg_candidates, NA)
  colnames(kg_candidates) <- paste0("k", 1:limit_est_groups)
  return(kg_candidates)
}

#' Initialises a dataframe that will contain an overview of metrics for each estimated configuration (for example adjusted randindex).
#'
#' @inheritParams initialise_beta
#' @param limit_est_groups maximum allowed number of groups that can be estimated
#' @return Returns an empty data.frame.
#' @examples
#' initialise_df_results(TRUE)
#' @export
initialise_df_results <- function(robust, limit_est_groups = 20) {
  df_results <- data.frame(
    kappa = NA, S = NA, k_common = NA,
    k1 = NA, k2 = NA, k3 = NA, k4 = NA, k5 = NA, k6 = NA, k7 = NA, k8 = NA, k9 = NA,
    k10 = NA, k11 = NA, k12 = NA, k13 = NA, k14 = NA, k15 = NA, k16 = NA, k17 = NA, k18 = NA, k19 = NA, k20 = NA,
    g = NA, g_true = NA, adjustedrandindex = NA, sigma2 = NA
  )
  df_results <- df_results[-1, ] # delete empty row
  if (robust) df_results <- df_results %>% dplyr::select(-kappa)
  return(df_results)
}

#' Defines the candidate values for C.
#'
#' @return numeric vector
#' @examples
#' define_C_candidates()
#' @export
define_C_candidates <- function() {
  C_candidates <- c(0, 10^(seq(-10, 10, by = 0.01)))
  return(C_candidates)
}

#' Initialises a dataframe which will contain the PIC for each configuration and for each value of C.
#'
#' @inheritParams calculate_VCsquared
#' @return Returns an empty data.frame.
#' @examples
#' initialise_df_pic(1:10)
#' @export
initialise_df_pic <- function(C_candidates) {
  df_pic <- data.frame(t(matrix(rep(NA, length(C_candidates)))))
  return(df_pic)
}

#' Adds the current configuration (number of groups and factors) to df_results.
#'
#' @param df_results dataframe with results for each estimated configuration
#' @param S estimated number of groups in current configuration
#' @param k estimated number of common factors in current configuration
#' @param kg vector with the estimated number of group specific factors in current configuration (augmented with NA's to reach a length of 20)
#' @return data.frame
#' @examples
#' add_configuration(initialise_df_results(TRUE), 3, 0, c(3, 3, 3, rep(NA, 17)))
#' @export
add_configuration <- function(df_results, S, k, kg) {
  colnames <- colnames(df_results) # using rbind removes/changes colnames, so save these here
  if (!("kappa" %in% colnames)) {
    # then robust = TRUE
    configurationdata <- c(S, k, kg, rep(NA, 4))
  } else {
    configurationdata <- c(NA, S, k, kg, rep(NA, 4))
  }
  suppressWarnings( #configurationdata should be added in the first .. columns; rest stays empty -> warning not important
    df_results <- df_results %>% rbind(configurationdata)
  )
  colnames(df_results) <- colnames

  return(df_results)
}

#' Fills in df_pic: adds a row with the calculated PIC for the current configuration.
#'
#' @param df input data frame
#' @param index_configuration index of the configuration of groups and factors
#' @param pic_e2 NxT matrix with the error terms
#' @inheritParams initialise_beta
#' @inheritParams estimate_beta
#' @param S number of estimated groups
#' @param k estimated number of common factors
#' @param kg vector with the estimated number of group specific factors for each group
#' @inheritParams calculate_VCsquared
#' @inheritParams OF_vectorized3
#' @param choice_pic parameter that defines which PIC is used to select the best configuration of groups and factors.
#' Options are "pic2017" (uses the PIC of \insertCite{Ando2017;textual}{RCTS}),
#' "pic2016" (\insertCite{Ando2016;textual}{RCTS}) weighs the fourth term with an extra factor relative to the size of the groups, and "pic2022".
#' They differ in the penalty they perform on the number of group specific factors (and implicitly on the number of groups). They also differ in the sense that they have
#' different NT-regions (where N is the number of time series and T is the length of the time series) where the estimated number of groups, and thus group specific factors will be wrong.
#' Pic2022 is the default (this PIC shrinks the problematic NT-region to very large N / very small T).
#' @return data.frame
#' @examples
#' set.seed(1)
#' original_data <- create_data_dgp2(30, 10)
#' Y <- original_data[[1]]
#' g <- original_data[[3]]
#' beta_est <- matrix(rnorm(4 * nrow(Y)), nrow = 4)
#' df_pic <- initialise_df_pic(1:5)
#' e <- matrix(rnorm(nrow(Y) * ncol(Y)), nrow(Y))
#' add_pic(df_pic, 1, TRUE, Y, beta_est, g, 3, 0, c(3, 3, 3), e, 1:5)
#' @export
add_pic <- function(df, index_configuration, robust, Y, beta_est, g, S, k, kg,
                    pic_e2, C_candidates,
                    method_estimate_beta = "individual",
                    choice_pic = "pic2022") {
  if(!is.na(beta_est[1])) {
    vars_est <- ncol(beta_est)
  } else {
    vars_est <- 0
  }
  pic_sigma2 <- calculate_sigma2(pic_e2)

  df[index_configuration, ] <- sapply(C_candidates, function(x) {
    calculate_PIC(x, robust, S, k, kg, pic_e2, pic_sigma2,
      NN = nrow(Y), TT = ncol(Y),
      method_estimate_beta, beta_est, g,
      vars_est,
      choice_pic
    )
  })
  if (anyNA(df[index_configuration, ])) {
    stop("There are NA's in the PIC for this configuration.")
  }
  if (min(df[index_configuration, ]) < 0) {
    stop("PIC Values should not be < 0.")
  }
  return(df)
}

#' Returns a vector with the indices of the subsets. Must start with zero.
#'
#' @param n number of subsets
#' @return numeric
#' @examples
#' define_number_subsets(3)
#' @export
define_number_subsets <- function(n) {
  return(0:(n - 1))
}

#' Function that returns for each candidate C the best number of groups and factors, based on the PIC.
#'
#' @inheritParams add_configuration
#' @inheritParams calculate_VCsquared
#' @param df_pic dataframe with the PIC for each configuration and for each candidate C
#' @return Returns a matrix with a row for each candidate value for C.
#' The first column contains the optimized number of groups (for each candidate C).
#' The second columns does the same for the number of common factors.
#' Column 3 until 22 do the same for the number of group specific factors. This is set to NA if the configuration has less than 20 groups estimated.
#' @examples
#' df_results <- add_configuration(initialise_df_results(TRUE),
#'   3, 0, c(3, 3, 3, rep(NA, 17))) #data.frame with one configuration
#' calculate_best_config(df_results, data.frame(t(1:5)), 1:5)
#' @export
calculate_best_config <- function(df_results, df_pic, C_candidates, limit_est_groups = 20) {
  indices_minimal_pic <- apply(df_pic, 2, which.min)


  result <- cbind(
    df_results$S[indices_minimal_pic],
    df_results$k_common[indices_minimal_pic]
  )

  temp <- df_results %>% dplyr::select(.data$k1:.data$k20)
  if (limit_est_groups > 20) {
    warning("change k1:k20 in line above")
  }
  for (i in 1:limit_est_groups) {
    result <- cbind(result, temp[, i][indices_minimal_pic])
  }
  return(as.matrix(result, nrow = (length(C_candidates))))
}

#' Fills in the optimized number of common factors for each C.
#'
#' @param df input
#' @param all_best_values data frame with the optimal number of groups, common factors and group specific factors
#' @param subset index of the subsample
#' @return data.frame
#' @examples
#' df_results <- add_configuration(initialise_df_results(TRUE),
#'   3, 0, c(3, 3, 3, rep(NA, 17))) #data.frame with one configuration
#' all_best_values <- calculate_best_config(df_results, data.frame(t(1:5)), 1:5)
#' rc <- fill_rc(initialise_rc(0:2, 1:5), all_best_values, 1)
#' @export
fill_rc <- function(df, all_best_values, subset) {
  df[, subset + 1] <- all_best_values[, 2] # fill in the best number of common factors
  return(df)
}

#' Fills in the optimized number of groups and group specific factors for each C.
#' @param df input
#' @param all_best_values data frame with the optimal number of groups, common factors and group specific factors
#' @param subset index of the subsample
#' @param S_cand vector with candidate values for the number of estimated groups
#' @param kg_cand vector with candidate values for the number of estimated group specific factors
#' @return data.frame
#' @examples
#' df_results <- add_configuration(initialise_df_results(TRUE),
#'   3, 0, c(3, 3, 3, rep(NA, 17))) #data.frame with one configuration
#' all_best_values <- calculate_best_config(df_results, data.frame(t(1:5)), 1:5)
#' rcj <- fill_rcj(initialise_rcj(0:2, 1:5) , all_best_values, 1, 2:4, 2:4)
#' @export
fill_rcj <- function(df, all_best_values, subset, S_cand, kg_cand) {
  if (max(kg_cand, na.rm = TRUE) > 0) {
    if (max(S_cand) > 1) {
      df[, subset + 1] <- apply(all_best_values[, 3:(3 + max(S_cand) - 1)], 1, function(x) paste(x, collapse = "_"))
    }
  } else {
    df[, subset + 1] <- 0
  }
  return(df)
}

#' Initialises rc.
#'
#' This function initialises a data frame which will eventually be filled with the optimized number of common factors for each C and for each subset of the original dataset.
#' @param indices_subset all indices of the subsets
#' @inheritParams calculate_VCsquared
#' @return data.frame
#' @examples
#' initialise_rc(0:2, 1:5)
#' @export
initialise_rc <- function(indices_subset, C_candidates) {
  return(data.frame(matrix(NA, length(C_candidates), length(indices_subset))))
}

#' Initialises rcj.
#'
#' This function initialises a data frame which will eventually be filled with the optimized number of groups and group specific factors for each C and for each subset of the original dataset.
#' @param indices_subset all indices of the subsets
#' @inheritParams calculate_VCsquared
#' @return data.frame
#' @examples
#' initialise_rcj(0:2, 1:5)
#' @export
initialise_rcj <- function(indices_subset, C_candidates) {
  return(data.frame(matrix(NA, length(C_candidates), length(indices_subset))))
}

#' Adds several metrics to df_results.
#'
#' @inheritParams add_configuration
#' @inheritParams add_pic
#' @inheritParams generate_Y
#' @inheritParams estimate_beta
#' @param pic_sigma2 sum of squared errors divided by NT
#' @param iteration number of iteration
# @importFrom mclust adjustedRandIndex
#' @importFrom tidyr nest
#' @importFrom dplyr everything
#' @importFrom dplyr tibble
#' @param add_rand adds the adjusted randindex to the df (requires the mclust package); used for simulations
#' @return data.frame with final estimations of each configuration
#' @examples
#' df_results <- add_configuration(initialise_df_results(TRUE),
#'   3, 0, c(3, 3, 3, rep(NA, 17))) #data.frame with one configuration
#' add_metrics(df_results, 1, 3.94, NA, round(runif(30, 1, 3)), NA, NA, NA, NA, 9)
#' @export
add_metrics <- function(df_results, index_configuration, pic_sigma2, beta_est, g, comfactor, lambda, factor_group, lambda_group,
                        iteration, g_true = NA, add_rand = FALSE) {
  df_results$sigma2[index_configuration] <- pic_sigma2
  if (add_rand) {
    if (!is.na(g_true[1])) df_results$adjustedrandindex[index_configuration] <- mclust::adjustedRandIndex(g, g_true)
  }
  df_results$table_g[index_configuration] <- paste(table(g), collapse = "_")
  if (!is.na(g_true[1])) df_results$table_g_true[index_configuration] <- paste(table(g_true), collapse = "_")
  df_results$g[index_configuration] <- paste(g, collapse = "-")
  df_results$number_of_iterations[index_configuration] <- iteration
  #note: use data.frame instead of tibble -> size is more clear then
  df_results$beta_est[index_configuration] <- nest(data.frame(beta_est), data = everything())
  df_results$comfactor[index_configuration] <- nest(data.frame(comfactor), data = everything())
  df_results$lambda[index_configuration] <- nest(data.frame(lambda), data = everything())
  df_results$factor_group[index_configuration] <- nest(tibble(factor_group), data = everything()) #do use tibble here (factor_group is a list)
  df_results$lambda_group[index_configuration] <- nest(data.frame(lambda_group), data = everything())
  return(df_results)
}



#' Finds the first stable interval after the first unstable point. It then defines the value for C for the begin, middle and end of this interval.
#'
#' @param list_vc list with resulting expression(VC^2) for each run
#' @param list_rc list with resulting rc for each run
#' @param list_rcj list with resulting rcj for each run
#' @param C_candidates candidates for C
#' @param S_cand candidates for S (number of groups)
# @param k_cand candidates for k (number of common factors)
# @param kg_cand candidates for kg (number of group specific factors)
#' @param return_short if TRUE, the function returns the dataframe filtered for several specified potential candidates for C
#' @param verbose when TRUE, it prints messages
#' @importFrom dplyr lead
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_count
#' @importFrom stringr str_sub
#' @importFrom stringr str_locate
#' @return data.frame with the optimized configuration for each candidate C (if return_short is FALSE) and
#' for each of the selected C's in the chosen stable interval (if return_short is TRUE).
#' @examples
#' set.seed(1)
#' all_best_values <- calculate_best_config(add_configuration(initialise_df_results(TRUE),
#'   3, 0, c(3, 3, 3, rep(NA, 17))),
#'   data.frame(t(1:5)), 1:5)
#' rc <- fill_rc(initialise_rc(0:1, 1:5), all_best_values, 0)
#' rc <- fill_rc(rc, all_best_values, 1)
#' rcj <- fill_rcj(initialise_rcj(0:1, 1:5) , all_best_values, 0, 2:4, 2:4)
#' rcj <- fill_rcj(rcj, all_best_values, 1, 2:4, 2:4)
#' get_best_configuration(sort(runif(5)), rc, rcj, 1:5, 2:4, return_short = FALSE)
#' @export
get_best_configuration <- function(list_vc, list_rc, list_rcj, C_candidates, S_cand, return_short = FALSE, verbose = FALSE) {

  # take the results of the full samples:
  #REPLACED if (class(list_rcj) == "list") { # if there were multiple runs
  #note: new function is() requires as dependency package "methods" -> use %in% instead
  if ("list" %in% class(list_rcj)) { # if there were multiple runs

    runs <- length(list_rcj)
    list_rc <- list_rc %>% map(1)
    list_rcj <- list_rcj %>% map(1)
  } else { #-> if there was one run -> matrices
    runs <- 1
    # note these are matrices, but use the same name as in the list-case above
    list_rc <- list_rc[, 1] # take first column which contains the rsults from the full dataset
    list_rcj <- list_rcj[, 1] # take first column which contains the rsults from the full dataset
    # list_rc <- matrix(unlist(list_rc), nrow = length(C_candidates))[, 1] #take first column which contains the rsults from the full dataset
    # list_rcj <- matrix(unlist(list_rcj), nrow = length(C_candidates))[, 1] #take first column which contains the rsults from the full dataset
  }

  # lists to matrices
  matrix_VC <- matrix(unlist(list_vc), nrow = length(C_candidates)) # 1 column for each run
  matrix_RC <- matrix(unlist(list_rc), nrow = length(C_candidates)) # 1 column for each run
  matrix_RCJ_allgroups <- matrix(unlist(list_rcj), nrow = length(C_candidates)) # 1 column for each run


  # get the number of factors of the first group
  matrix_und <- apply(matrix_RCJ_allgroups, 1:2, function(x) str_locate(x, "_")[1] - 1) # Matrix with index of first occurence of "_"
  matrix_RCJ <- matrix(as.numeric(str_sub(matrix_RCJ_allgroups, 1, matrix_und)), nrow = length(C_candidates))
  matrix_groups <- matrix(unlist(list_rcj), nrow = length(C_candidates))
  matrix_groups <- matrix(as.numeric(max(S_cand) - str_count(matrix_groups, "NA")), nrow = length(C_candidates)) # 1 column for each run

  # combine into df
  df_plot_median <- data.frame(
    C = C_candidates,
    median_VC = apply(data.frame(t(matrix_VC)), 2, median),
    median_groups = apply(data.frame(t(matrix_groups)), 2, median),
    median_RC = apply(data.frame(t(matrix_RC)), 2, median),
    median_RCJ = apply(data.frame(t(matrix_RCJ)), 2, median)
    # matrix_RCJ_allgroups = matrix_RCJ_allgroups,
  )
  # add the number of group factors for the other groups as well
  for (i in 2:max(S_cand)) {
    temp <- (matrix_RCJ_allgroups %>% apply(1:2, function(x) str_split(x, "_")))
    # get number of group factors for group i for all runs
    groupfactors <- matrix(as.numeric(temp %>% purrr::map(1) %>% purrr::map(i) %>% str_replace_all("NA", "")), nrow = length(C_candidates))
    df_plot_median <- cbind(df_plot_median, apply(groupfactors, 1, median))
  }

  # get the first unstable point of VC2 and the configurations corresponding to that point
  first_unstable_point <- which(df_plot_median$median_VC > 0)[1] # first point where (median of) VC2 > 0 .
  if (max(df_plot_median$median_VC) == 0) {
    print("VC2 is always 0, set first unstable point to 1.")
    first_unstable_point <- 1 # set to one in this case
  }
  # find the first stable point (beginning of first interval) after the first unstable point
  if (verbose) print(paste("first_unstable_point = ", first_unstable_point))
  beginpoint <- which(df_plot_median$median_VC == 0) # ...min(df_plot_median$median_VC))
  beginpoint <- beginpoint[beginpoint >= first_unstable_point][1]
  chosen_C_beginpoint <- C_candidates[beginpoint]
  if (verbose) {
    print(paste("beginpoint:", beginpoint))
    print(paste("C (beginpoint) = ", chosen_C_beginpoint))
  }

  endpoint <- which(df_plot_median$median_VC == 0 & dplyr::lead(df_plot_median$median_VC) > 0 &
    dplyr::lead(df_plot_median$median_VC, 2) > 0) # make sure the next 2 points have VC2 > 0: (just in case the first point would be an abberation)
  endpoint <- endpoint[endpoint > beginpoint][1]
  # Check if endpoint is defined (not the case when VC2 stays at zero for all C's starting from 'chosen_C_beginpoint' of if there is only one halfopen stable interval)
  if (length(endpoint) > 0 & !is.na(endpoint)) {
    chosen_C_endpoint <- C_candidates[endpoint]
    middlepoint <- which.min(abs(C_candidates - (chosen_C_beginpoint + chosen_C_endpoint) / 2))
    if (length(middlepoint) > 0) {
      chosen_C_middle <- C_candidates[middlepoint]
    } else {
      chosen_C_middle <- NA
    }
    middlepoint_log <- round(mean(c(endpoint, beginpoint))) # logarithmic middlepoint

    if (verbose) {
      print(paste("middlepoint:", middlepoint))
      print(paste("C (middle) =", round(C_candidates[middlepoint], 2)))
      print(paste("middlepoint_log:", middlepoint_log))
      print(paste("C (middle_log) =", round(C_candidates[middlepoint_log], 2)))
      print(paste("endpoint:", endpoint))
      print(paste("C (endpoint) =", round(C_candidates[endpoint], 2)))
    }
  } else {
    if (verbose) message("No endpoint could be defined")
    middlepoint <- NA
    chosen_C_middle <- -100
    middlepoint_log <- NA
    endpoint <- NA
    chosen_C_endpoint <- -100
  }

  colnames(df_plot_median)[5:ncol(df_plot_median)] <- c(paste0("k", 1:max(S_cand)))
  if (return_short) {
    return(tabulate_potential_C(df_plot_median, runs, beginpoint, middlepoint_log, middlepoint, endpoint, S_cand))
  } else {
    return(df_plot_median)
  }
}

#' Shows the configurations for potential C's of the first stable interval (beginpoint, middlepoint and endpoint)
#'
#' @param df input dataframe
#' @param runs number of panel data sets for which the algorithm has run. If larger than one, the median VC2 is used to determine C.
#' @param beginpoint first C of the chosen stable interval
#' @param middlepoint_log middle C (on a logscale) of the chosen stable interval
#' @param middlepoint middle C of the chosen stable interval
#' @param endpoint last C of the chosen stable interval
#' @param S_cand candidate number for the number of groups
#' @return data.frame
tabulate_potential_C <- function(df, runs, beginpoint, middlepoint_log, middlepoint, endpoint, S_cand) {
  if (runs == 1) {
    colnames(df)[1:5] <- c("C", "VC2", "S", "k", "k1")
    for (i in 2:max(S_cand)) colnames(df)[4 + i] <- paste0("k", i)
  } else {
    colnames(df)[1:5] <- c("C", "median VC2", "median S", "median k", "median k1")
    for (i in 2:max(S_cand)) colnames(df)[4 + i] <- paste0("median k", i)
  }
  rownames(df) <- NULL
  if (is.na(middlepoint_log) & is.na(middlepoint) & is.na(endpoint)) { # when there is only one stable interval after the first unstable point and this interval is halfopen -> only beginpoint can be defined
    return(df[c(beginpoint), ] %>% mutate(location = c("beginpoint")))
  } else {
    return(df[c(beginpoint, middlepoint_log, middlepoint, endpoint), ] %>% mutate(location = c("beginpoint", "middlepoint (logscale)", "middlepoint", "endpoint")))
  }
}

#' Plots expression(VC^2) along with the corresponding number of groups (orange), common factors (darkblue) and group factors of the first group (lightblue).
#'
#' @inheritParams calculate_VCsquared
#' @param VC_squared measure of variability in the optimal configuration between the subsets
#' @param S_cand candidate numbers for the number of groups
# @param k_cand candidate numbers for the number of common factors
# @param kg_cand candidate numbers for the number of group specific factors
#' @param xlim_min starting point of the plot
#' @param xlim_max end point of the plot
#' @param add_true_lines if set to TRUE, for each C the true number of groups, common factors, and group specific factors of group 1 will be added to the plot
#' @param verbose if TRUE, more details are printed
#' @return A ggplot object.
#' @examples
#' \donttest{
#' set.seed(1)
#' #requires filled in dataframes rc and rcj
#' all_best_values <- calculate_best_config(add_configuration(initialise_df_results(TRUE),
#'   3, 0, c(3, 3, 3, rep(NA, 17))),
#'   data.frame(t(1:20)), 1:20)
#' rc <- fill_rc(initialise_rc(0:1, 1:20), all_best_values, 0)
#' rc <- fill_rc(rc, all_best_values, 1)
#' rcj <- fill_rcj(initialise_rcj(0:1, 1:20) , all_best_values, 0, 2:4, 2:4)
#' rcj <- fill_rcj(rcj, all_best_values, 1, 2:4, 2:4)
#' plot_VCsquared(c(runif(9), 0, 0, runif(9)), rc, rcj, 1:20, 2:4)
#' }
#' @export
plot_VCsquared <- function(VC_squared, rc, rcj, C_candidates, S_cand,
                           xlim_min = 1e-3,
                           xlim_max = 1e2, add_true_lines = FALSE, verbose = FALSE) {
  # get data frame collecting for each C the VC2 and the best configuration
  df <- get_best_configuration(VC_squared, rc, rcj, C_candidates, S_cand, verbose = verbose)
  df <- df %>%
    dplyr::filter(.data$C != 0)
  colnames(df) <- c("C", "VC2", "S", "k", paste0("k", 1:max(S_cand)))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$C, y = .data$VC2)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    # geom_vline(aes(xintercept = C_candidates[beginpoint]), col = "green", linetype = "dashed") +
    ggplot2::ylab(expression("VC"^2)) +
    ggplot2::scale_x_log10() +
    ggplot2::ggtitle("") +
    ggplot2::coord_cartesian(xlim = c(xlim_min, xlim_max))

  if (add_true_lines) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(y = .data$k), col = "blue") +
      ggplot2::geom_point(ggplot2::aes(y = .data$k1), col = "lightblue", size = 0.7) +
      ggplot2::geom_line(ggplot2::aes(y = .data$S), col = "orange")
  }

  return(p)
}



#' Function that returns the final clustering, based on the estimated number of groups and common and group specific factors.
#'
#' @param df input dataframe (this will be df_results_full)
#' @param opt_groups the optimal number of groups
#' @param k the optimal number of common factors
#' @param kg vector with the optimal number of group specific factors
#' @param type defines which estimation to return: options are "clustering", "beta", "fg" (group specific factors), "lg" (loadings corresponding to fg), "f" (common factors), "l" (loadings corresponding to f),
#' @param limit_est_groups maximum allowed number of groups that can be estimated
#' @return This function returns the estimations of the chosen configuration.
#' If type is "clustering" it returns a numeric vector with the estimated group membership for all time series.
#' If type is "beta", "lg" the function returns a data.frame.
#' If type is "f" or "l" the function also returns a data.frame. If no common factors were estimated in the optimized configuration, then NA is returned.
#' If type is "fg" the function returns a list.
#' @examples
#' get_final_estimation(df_results_example, 3, 0, c(3, 3, 3), "clustering")
#' get_final_estimation(df_results_example, 3, 0, c(3, 3, 3), "beta")
#' get_final_estimation(df_results_example, 3, 0, c(3, 3, 3), "fg")
#' get_final_estimation(df_results_example, 3, 0, c(3, 3, 3), "lg")
#' @export
get_final_estimation <- function(df, opt_groups, k, kg, type, limit_est_groups = 20) {
  stopifnot(length(kg) <= limit_est_groups) # code is implemented up to 20 estimated groups
  df <- df %>%
    dplyr::filter(.data$S == opt_groups, .data$k_common == k, .data$k1 == kg[1]) %>%
    final_estimations_filter_kg(kg)
  if(type == "clustering") {
    final_g <- as.numeric(unlist(df$g %>% str_split("-")))
    return(final_g)
  }
  if(type == "fg") {
    object <- df$factor_group[[1]]
    #note: use object[[1]] in while instead of object, since output must be a list of groupfactors for each group
    while("list" %in% class(object[[1]])) object <- object[[1]] #for the serialized algorithm i need 3 times [[1]]. parallel algorithm only 1 time
    return(object)
  }
  if(type == "f") {
    if(k == 0) return(NA)
    object <- df$comfactor[[1]]
    while("list" %in% class(object)) object <- object[[1]] #for the serialized algorithm i need 3 times [[1]]. parallel algorithm only 1 time
    return(object)
  }
  if(type == "lg") {
    object <- df$lambda_group[[1]]
    while("list" %in% class(object)) object <- object[[1]] #for the serialized algorithm i need 3 times [[1]]. parallel algorithm only 1 time
    return(object %>% arrange(.data$id))
  }
  if(type == "l") {
    if(k == 0) return(NA)
    object <- df$lambda[[1]]
    while("list" %in% class(object)) object <- object[[1]] #for the serialized algorithm i need 3 times [[1]]. parallel algorithm only 1 time
    return(object)
  }
  if(type == "beta") {
    object <- df$beta_est[[1]]
    while("list" %in% class(object)) object <- object[[1]] #for the serialized algorithm i need 3 times [[1]]. parallel algorithm only 1 time
    return(object)
  }
}

#' Filters dataframe on the requested group specific factors configuration.
#'
#' @param df input dataframe
#' @param kg vector with number of group specific factors for each group, on which should be filtered
#' @return data.frame
final_estimations_filter_kg <- function(df, kg) {
  if (length(kg) > 1) df <- df %>% dplyr::filter(.data$k2 == kg[2])
  if (length(kg) > 2) df <- df %>% dplyr::filter(.data$k3 == kg[3])
  if (length(kg) > 3) df <- df %>% dplyr::filter(.data$k4 == kg[4])
  if (length(kg) > 4) df <- df %>% dplyr::filter(.data$k5 == kg[5])
  if (length(kg) > 5) df <- df %>% dplyr::filter(.data$k6 == kg[6])
  if (length(kg) > 6) df <- df %>% dplyr::filter(.data$k7 == kg[7])
  if (length(kg) > 7) df <- df %>% dplyr::filter(.data$k8 == kg[8])
  if (length(kg) > 8) df <- df %>% dplyr::filter(.data$k9 == kg[9])
  if (length(kg) > 9) df <- df %>% dplyr::filter(.data$k10 == kg[10])
  if (length(kg) > 10) df <- df %>% dplyr::filter(.data$k11 == kg[11])
  if (length(kg) > 11) df <- df %>% dplyr::filter(.data$k12 == kg[12])
  if (length(kg) > 12) df <- df %>% dplyr::filter(.data$k13 == kg[13])
  if (length(kg) > 13) df <- df %>% dplyr::filter(.data$k14 == kg[14])
  if (length(kg) > 14) df <- df %>% dplyr::filter(.data$k15 == kg[15])
  if (length(kg) > 15) df <- df %>% dplyr::filter(.data$k16 == kg[16])
  if (length(kg) > 16) df <- df %>% dplyr::filter(.data$k17 == kg[17])
  if (length(kg) > 17) df <- df %>% dplyr::filter(.data$k18 == kg[18])
  if (length(kg) > 18) df <- df %>% dplyr::filter(.data$k19 == kg[19])
  if (length(kg) > 19) df <- df %>% dplyr::filter(.data$k20 == kg[20])
  if (nrow(df) != 1) warning("only 1 configuration should be found")
  return(df)
}
#' Constructs dataframe where the rows contains all configurations that are included and for which the estimators will be estimated.
#'
#' @param S_cand candidates for S (number of groups)
#' @param k_cand candidates for k (number of common factors)
#' @param kg_cand candidates for kg (number of group specific factors)
#' @return data.frame
#' @examples
#' define_configurations(2:4, 0, 2:3)
#' @export
define_configurations <- function(S_cand, k_cand, kg_cand) {
  config_cand <- define_kg_candidates(min(S_cand), min(kg_cand), max(kg_cand), nfv = TRUE) %>% mutate(S = min(S_cand), k = min(k_cand))

  for (i in S_cand) {
    for (j in k_cand) {
      if (i != min(S_cand) | j != min(k_cand)) {
        config_cand <- rbind(config_cand, define_kg_candidates(i, min(kg_cand), max(kg_cand), nfv = TRUE) %>% mutate(S = i, k = j))
      }
    }
  }

  return(config_cand %>% dplyr::select(.data$S, .data$k, dplyr::everything()))
}
