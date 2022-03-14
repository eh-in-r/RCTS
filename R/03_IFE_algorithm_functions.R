
utils::globalVariables(c("use_robust",
                         "g", "g_true",
                         "beta_est", "beta_true",
                         "number_of_time_series", "length_of_time_series",
                         "length_of_time_series_fulldata",
                          "number_of_groups_fixedvalue", "number_of_groups_candidates_fixedvalue",
                         "number_of_common_factors_fixedvalue", "k_common_candidates",
                         "number_of_group_factors_fixedvalue", "k_max",
                         "comfactor", "lambda",
                         "factor_group", "factor_group_true",
                         "lambda_group", "lambda_group_true",
                         "num_factors_may_vary_fixedvalue",
                         "number_of_variables_fixedvalue", "number_vars_estimated_fixedvalue",
                         "ALTERNATIVE_PIC", "C",
                         #"FGR_FACTOR","FGR_FACTOR_sd","LGR_FACTOR","LGR_FACTOR_mean",
                         #"iteration",
                         #"beta_true_heterogeen_groups", "beta_true_heterogeen_individueel", "beta_true_homogeen",
                         "method_estimate_factors",
                         "Smax", "X", "Y",
                         "percent_outliers",
                         "df_results", "sigma2_max_model",
                         "indices_subset",
                         #"kappa_candidates",
                         "grid",
                         #"v1", "v2", "v3",
                         #"usecoviddata_cases","usecoviddata_deaths",
                         "method_estimate_beta",
                         "limit_true_groups"


))



#' Function used in generating simulated data with non normal errors.
#'
#' Used to include cross-sectional dependence or serial dependence into the simulated panel data.
#' @param parameter amount of cross-sectional dependence
#' @inheritParams estimate_beta
#' @return NxN covariance matrix
#' @examples
#' create_covMat_crosssectional_dependence(0.3, 300)
#' @export
create_covMat_crosssectional_dependence <- function(parameter,NN) {
  covMat = matrix(NA, nrow = NN, ncol = NN)
  for(i in 1:NN) {
    for(j in 1:NN) {
      covMat[i,j] = parameter^abs(i-j)
    }
  }
  return(covMat)
}

#' Helpfunction in create_true_beta() for the option beta_true_heterogeen_groups. (This is the default option.)
#'
#' @param number_of_variables number of observable variables
#' @param true_number_of_groups true number of groups
#' @param extra_beta_factor_bth option to multiply the coefficients in true beta; default = 1
#' @param limit_true_groups_bth  Maximum number of true groups in a simulation-DGP for which the code in this package is implemented. Currently equals 12. For application on realworld data this parameter is not relevant.
#' @importFrom stats runif
#' @importFrom magrittr %>%
beta_true_heterogroups <- function(number_of_variables, true_number_of_groups, extra_beta_factor_bth = 1, limit_true_groups_bth = limit_true_groups) {
  stopifnot(true_number_of_groups < limit_true_groups_bth) #Code allows up to 12 true groups at this point.

  #######################################################################################
  # Define true values for beta, for each group, when there are 3 or less observable variables
  #These are the values for DGP 3 & 4 (For DGP 1 & 2 beta_true is defined in 08_IFE_create_data_dgp1.R)

  #1st element is the intercept.
  beta_part1 = c(0, 4,  3.5,  3,  2.5)  * extra_beta_factor_bth
  beta_part2 = c(0,-2.5, -2, -2.5, -2)  * extra_beta_factor_bth
  beta_part3 = c(0, 1,  0.5,  1.5,  1)  * extra_beta_factor_bth

  beta_define_further_true_values <- function(number_of_values) {
    #starting from group 4 we randomly generate values for the true values of beta
    return(c(0,(round(runif(number_of_values),1) - 0.5) * 4) * extra_beta_factor_bth) #between -2 and 2
  }
  beta_part4 = beta_define_further_true_values(4)
  beta_part5 = beta_define_further_true_values(4)
  beta_part6 = beta_define_further_true_values(4)
  beta_part7 = beta_define_further_true_values(4)
  beta_part8 = beta_define_further_true_values(4)
  beta_part9 = beta_define_further_true_values(4)
  beta_part10 = beta_define_further_true_values(4)
  beta_part11 = beta_define_further_true_values(4)
  beta_part12 = beta_define_further_true_values(4)

  #######################################################################################
  # Define true values for beta_est, for each group, when there are more than 3 variables
  if(number_of_variables > 3) { #add when there are more than 3 observable variables:
    beta_part1 = c(beta_part1, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part2 = c(beta_part2, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part3 = c(beta_part3, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part4 = c(beta_part4, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part5 = c(beta_part5, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part6 = c(beta_part6, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part7 = c(beta_part7, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part8 = c(beta_part8, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part9 = c(beta_part9, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part10 = c(beta_part10, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part11 = c(beta_part11, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    beta_part12 = c(beta_part12, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
  }

  #Add the values of the different groups together in one object.
  if(true_number_of_groups >= 1) {
    beta_true = matrix(beta_part1[1:(number_of_variables + 1)])
  }
  if(true_number_of_groups >= 2) beta_true = beta_true %>% cbind(beta_part2[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 3) beta_true = beta_true %>% cbind(beta_part3[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 4) beta_true = beta_true %>% cbind(beta_part4[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 5) beta_true = beta_true %>% cbind(beta_part5[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 6) beta_true = beta_true %>% cbind(beta_part6[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 7) beta_true = beta_true %>% cbind(beta_part7[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 8) beta_true = beta_true %>% cbind(beta_part8[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 9) beta_true = beta_true %>% cbind(beta_part9[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 10) beta_true = beta_true %>% cbind(beta_part10[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 11) beta_true = beta_true %>% cbind(beta_part11[1:(number_of_variables + 1)])
  if(true_number_of_groups >= 12) beta_true = beta_true %>% cbind(beta_part12[1:(number_of_variables + 1)])


  return(beta_true)
}

#' Creates beta_true, which contains the true values of beta (= the coefficients of X)
#'
#' beta_true_heterogeen_groups is the default case
#' @inheritParams estimate_beta
#' @param true_number_of_groups number of groups
# @param use_real_world_data indicates using realworld data; defaults to FALSE
#' @param extra_beta_factor multiplies coefficients in beta_est; default = 1
#' @param beta_true_homogeneous whether true beta is equal for all individuals
#' @param beta_true_heterogeneous_groups whether true beta is equal within groups, and different between groups
#' @param beta_true_heterogeneous_individuals whether true beta is different for all individuals
#' @param limit_true_groups Maximum number of true groups in a simulation-DGP for which the code in this package is implemented. Currently equals 12. For application on realworld data this parameter is not relevant.
#' @return matrix with number of rows equal to number of observable variables + 1 (the first row contains the intercept) and number of culumns
#' equal to the true number of groups.
#' @examples
#' library(tidyverse)
#' #Decide if beta_est is common, or specific to groups or individuals: Choose 1 of the following 3.
#' create_true_beta(3, NN = 300, true_number_of_groups = 3,
#'   beta_true_homogeneous = FALSE, beta_true_heterogeneous_groups = TRUE,
#'   beta_true_heterogeneous_individuals = FALSE)
#' @importFrom stats rnorm
#' @export
create_true_beta <- function(number_of_variables,
                             NN,
                             true_number_of_groups,
                             beta_true_homogeneous,
                             beta_true_heterogeneous_groups,
                             beta_true_heterogeneous_individuals,
                             limit_true_groups = 12,
                             extra_beta_factor = 1) {
  stopifnot((beta_true_homogeneous + beta_true_heterogeneous_groups + beta_true_heterogeneous_individuals) == 1)
  #real world data: beta_true does not exist -> return NA
  #if(use_real_world_data) return(NA)


  if(number_of_variables > 0) {
    #common beta_true:
    if(beta_true_homogeneous) {
      beta_true = c(0, c(1,2,3,28,33,38,43,48,53,58,63,68,73,78,83)[1:number_of_variables]) #intercept 0, and after that values for the true beta's
      beta_true = matrix(rep(beta_true, true_number_of_groups), nrow = (number_of_variables + 1))
    }
    #groupsspecific beta_true: -> default case
    if(beta_true_heterogeneous_groups) {
      beta_true = beta_true_heterogroups(number_of_variables, true_number_of_groups, extra_beta_factor_bth = extra_beta_factor, limit_true_groups_bth = limit_true_groups)
    }
    #individualspecific beta_true:
    if(beta_true_heterogeneous_individuals) {
      beta_true = matrix(rnorm(NN * (number_of_variables + 1)), nrow = (number_of_variables + 1), ncol = NN)
    }
  } else {
    beta_true = rep(NA,true_number_of_groups)
  }

  return(beta_true)
}


#' Creates X (the observable variables) to use in simulations.
#'
#' X is an array with dimensions N, T and number of observable variables.
#' The variables are randomly generated with mean 0 and sd 1.
#' @inheritParams create_covMat_crosssectional_dependence
#' @inheritParams estimate_beta
#' @return array with dimensions N x T x number of observable variables
#' @examples
#' initialise_X(300,30, number_of_variables = 3)
#' @export
initialise_X <- function(NN,TT, number_of_variables = number_of_variables_fixedvalue) {
  if(number_of_variables > 0) {
    X = array(0,dim=c(NN, TT, number_of_variables))
    for(i in 1:NN) {
      for(t in 1:TT) {
        for(k in 1:number_of_variables) {
          X[i,t,k] = rnorm(1, mean = 0, sd = 1)
        }
      }
    }
    X = scaling_X(X, firsttime = TRUE, number_of_variables = number_of_variables)


    return(X)
  } else {
    return(NA)
  }
}


#'Scaling of X.
#'
#' @param X input
#' @param firsttime Scaling before generating Y and before adding outliers: this is always with mean and sd. If this is FALSE, it indicates that
#' we are using the function for a second time, after adding the outliers. In the robust case it uses median and MAD, otherwise again mean and sd.
# @param use_real_world_data Parameter to indicate using real world data. Defaults to FALSE.
#' @inheritParams create_true_beta
#' @examples
#' X = initialise_X(300,30, number_of_variables = 3)
#' use_robust = TRUE
#' scaling_X(X,TRUE, number_of_variables = 3)
#' @importFrom stats sd
#' @export
scaling_X <- function(X, firsttime, number_of_variables = number_of_variables_fixedvalue) {
  #
  # replaced this codeblock by a less concise but more clear codeblock
  #
  ##################
  # if(use_robust & !firsttime) {
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
  if(number_of_variables > 0) {
    for(k in 1:number_of_variables) {
      #print(paste("sd of variable",k,": Before:",sd(X[,,k])))
      if(mad(X[,,k]) != 0) {
        if(firsttime) {
          # if(use_real_world_data) {
          # #for eclipzdata, we cannot use scale(),
          #   #  since there are variables (for example age) that have constant columns, so scale() (which is columnbased) would produce errors
          #   #message("Scale with mean and sd of NxT-matrix")
          #   X[,,k] = (X[,,k] - mean(X[,,k])) / sd(X[,,k])
          # } else {
            #message("Scale with mean and sd (for each t separate)")
            X[,,k] = scale(X[,,k])      #Note that this is column-based (timeindex) scaling!
          # }
        } else {
          if(use_robust) {
            #message("Scale with median and mad")
            med = median(X[,,k])
            mad = mad(X[,,k])

            X[,,k] = (X[,,k] - med) / mad #Note that this is NOT column-based (timeindex) scaling!
          } else {
            # if(use_real_world_data) {
            #   #for eclipzdata, we cannot use scale(),
            #   #  since there are variables (for example age) that have constant columns, so scale() (which is columnbased) would produce errors
            #   #message("Scale with mean and sd of NxT-matrix")
            #   X[,,k] = (X[,,k] - mean(X[,,k])) / sd(X[,,k])
            # } else {
              #message("Scale with mean and sd (for each t separate)")
              X[,,k] = scale(X[,,k])      #Note that this is column-based (timeindex) scaling!
            # }
          }
        }

        #print(paste("sd of variable",k,": After:",sd(X[,,k])))
      }
    }
    stopifnot(!is.nan(X[1,1,1]))
  }
  return(X)
}

#' This function restructures X (which is an 3D-array of dimensions (N,T,p) to a 2D-matrix of dimension (NxT,p).
#'
#' This 2D-matrix is what Ando/Bai worked with.
#' @param X input
#' @inheritParams create_true_beta
#' @inheritParams estimate_beta
#' @examples
#' X = initialise_X(300, 30, number_of_variables = 3)
#' X_restructured = restructure_X_to_order_slowN_fastT(X,
#'   number_of_variables = 3, number_vars_estimated = 3)
#' @export
restructure_X_to_order_slowN_fastT <- function(X,
                                               number_of_variables = number_of_variables_fixedvalue,
                                               number_vars_estimated = number_vars_estimated_fixedvalue) {

  if(number_of_variables > 0) {
    if( length(dim(X)) == 2) {
      #occurs when only one element in group
      X = array(X, dim=c(1,nrow(X),ncol(X)))
    }
    # if(use_real_world_data) {
    #   number_of_vars = number_of_variables
    # } else {
      number_of_vars = number_vars_estimated
    # }
    X_local = matrix(NA, nrow = nrow(X)*ncol(X), ncol = number_of_vars)
    for(k in 1:number_of_vars) {
      X_local[,k] = c(t(X[,,k]))
    }

    return(X_local)
  } else {
    #when there are no observable variables, X was NA and stays NA
    return(NA)
  }
}



#' Generates the true groupfactorstructure, to use in simulations.
#'
#' Loadings and factors are generated by:
#' loadings ~ N(LGR_FACTOR_mean, j * LGR_FACTOR) -> default case will be N(0,j)
#' factors ~ N(j * FGR_FACTOR, FGR_FACTOR_sd) -> default case will be N(j,1)
#'
#' @param S true number of groups
#' @param true_number_of_group_factors vector with as length the number of groups, where each element is the true number of groupfactors of that group.
#' @inheritParams estimate_beta
#' @param g_true true group membership
#' @return list: first element contains true groupfactors and second element contains true groupfactor loadings
#' @importFrom dplyr mutate
#' @examples
#' library(tidyverse)
#' #For 3 groups, each with 3 groupfactors:
#' g_true = ceiling(runif(300) * 3)
#' generate_grouped_factorstructure(3, c(3, 3, 3), TT = 30, g_true)
#' @importFrom dplyr bind_rows
#' @export
generate_grouped_factorstructure <- function(S, true_number_of_group_factors, TT = length_of_time_series_fulldata, g_true) {
  if(!exists("LGR_FACTOR") | !exists("FGR_FACTOR")) {
    #set these to their default value
    LGR_FACTOR_mean = 0
    LGR_FACTOR = 1
    FGR_FACTOR_sd = 1
    FGR_FACTOR = 1
  }
  if(mean(true_number_of_group_factors) > 0) {


    LGR = list()
    factor_group_true = list()
    for(i in 1:S) {
      LGR[[i]] = matrix(nrow = sum(g_true == i), ncol = true_number_of_group_factors[i])
      factor_group_true[[i]] = matrix(NA, nrow = true_number_of_group_factors[i], ncol = TT)
    }

    for(j in 1:S) {
      for(i in 1:true_number_of_group_factors[j]) {
        for(k in 1:(sum(g_true == j))) {
          LGR[[j]][k,i] = rnorm(1, mean = LGR_FACTOR_mean, sd = j * LGR_FACTOR)
        }
        for(k in 1:TT) {
          factor_group_true[[j]][i,k] = rnorm(1, mean = j * FGR_FACTOR, sd = FGR_FACTOR_sd)
        }
      }
    }
    #lambda_group_true to dataframe (easier to work with)
    test = data.frame(LGR[[1]])
    test %>% mutate(groep = 1)
    lambda_group_true = data.frame(LGR[[1]]) %>% mutate(groep = 1, id = which(g_true == 1))
    if(true_number_of_group_factors[1] == 1) { #change name when only 1 factor (when there are more names are automatically X1,X2,...)
      names(lambda_group_true)[1] = "X1"
    }
    for(i in 2:S) {
      temporary = data.frame(LGR[[i]]) %>% mutate(groep = i, id =  which(g_true == i))
      if(true_number_of_group_factors[i] == 1) {
        names(temporary)[1] = "X1"
      }
      if(ncol(lambda_group_true) >= ncol(temporary)) {
        lambda_group_true = lambda_group_true %>% bind_rows(temporary)
      } else {
        lambda_group_true = temporary %>% bind_rows(lambda_group_true)
      }
    }
    lambda_group_true = lambda_group_true %>% arrange(groep)
    lambda_group_true[is.na(lambda_group_true)] <- 0
  } else {
    factor_group_true = NA
    lambda_group_true = NA
  }

  #geval van 1 factor in group: format to matrix:
  for(groep in 1:S) {
    if(true_number_of_group_factors[groep] == 1) {
      factor_group_true[[groep]] = matrix(factor_group_true[[groep]], nrow = true_number_of_group_factors[groep])
    }
  }
  return(list(factor_group_true, lambda_group_true))
}


#' Generate panel data Y for simulations.
#' @inheritParams estimate_beta
#' @param true_number_of_common_factors true number of common factors
#' @param true_number_of_group_factors Vector of length the number of groups. Each element contains the true number of group factors for that group.
#' @param g_true vector of length N with true group memberships
#' @param beta_true true coefficients with the observable variables
#' @param lambda_group_true loadings of the group factors
#' @param factor_group_true groupfactors
#' @param lambda_true loadings of the common factors
#' @param comfactor_true common factors
#' @param epsilon NN x TT-matrix containing the error term
#' @param X dataframe with the observed variables
# @param ando_bai loads Ando/Bai-dataset; only for testing purposes. Defaults to FALSE.
# @param ando_bai_2017 loads Ando/Bai-dataset; only for testing purposes. Defaults to FALSE.
#' @inheritParams create_true_beta
#' @return N x T matrix
#' @export
generate_Y <- function(NN, TT, true_number_of_common_factors, true_number_of_group_factors,
                       g_true, beta_true, lambda_group_true, factor_group_true,
                       lambda_true, comfactor_true, epsilon, X,
                       # ando_bai = FALSE,
                       # ando_bai_2017 = FALSE,
                       # use_real_world_data = FALSE,
                       number_of_variables = number_of_variables_fixedvalue) {

  #Define the size of the panel data:
  Y = matrix(NA, nrow = NN, ncol = TT) #initialisation, later on this gets filled in


    for(i in 1:NN) {

      #if(!ando_bai & !ando_bai_2017 & !use_real_world_data) {
        if(mean(true_number_of_group_factors) > 0) {
          dropvars <- names(lambda_group_true) %in% c("groep","id")
          LAMBDAGROUP = as.matrix(subset(lambda_group_true, lambda_group_true$id == i)[!dropvars])
          LAMBDAGROUP = LAMBDAGROUP[1:true_number_of_group_factors[g_true[i]]]
        } else {
          LAMBDAGROUP = NA
        }
      #}


      for(t in 1:TT) {

        if(number_of_variables > 0) {
          XT = c(1, X[i,t,]) %*% beta_true[,g_true[i]] #add 1 to X[i,t,] to make room for the intercept (which is in beta_est)

        } else {
          XT = 0
        }


        #interactive fixed effects: both common and group-specific factors
        # if(ando_bai) {
        #   Y[i,t] = XT +
        #     XL[t,i] + ERR[t,i]
        # } else if(ando_bai_2017) {
        #   Y[i,t] = XT +
        #     XG[t,i] +
        #     XL[t,i] + ERR[t,i]
        # } else {
          #randomly generated data
          if(true_number_of_common_factors > 0) {
            LF = t(lambda_true[,i]) %*% comfactor_true[,t]
          } else {
            LF = 0
          }
          if(mean(true_number_of_group_factors) > 0) {
            LF_GROUP = LAMBDAGROUP %*% factor_group_true[[g_true[i]]][,t]
          } else {
            LF_GROUP = 0
          }

          Y[i,t] = XT +
            LF +
            LF_GROUP +
            epsilon[i,t]
          stopifnot(!is.nan(Y[i,t]))
        # }
      }
    }

  return(Y)
}


#' Initialisation of estimation of beta (the coefficients with the observable variables)
#'
#' Note: this needs to be called before the definition of grid.
#' @inheritParams generate_Y
#' @param use_robust robust or classical estimation
#' @param number_vars_estimated number of variables from which the coefficients are estimated
#' @param number_of_groups number of groups
#' @inheritParams define_object_for_initial_clustering_macropca
#' @param nosetting_lmrob option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
#' @param special_case_dgp1 special case for data generated according to dgp 1: it changes the 1st variable in X to 1 (-> intercept). Consequently the estimation of beta needs to be restructured slightly.
#' @examples
#' library(RCTS)
#' library(tidyverse)
#' library(robustbase)
#'
#' X = RCTS::X_dgp3
#' Y = RCTS::Y_dgp3
#' #Set estimations for group factors and its loadings, and group membership to the true value
#' lambda_group = RCTS::lambda_group_true_dgp3
#' factor_group = RCTS::factor_group_true_dgp3
#' g = RCTS::g_true_dgp3
#' #There are no common factors to be estimated  -> but needs placeholder
#' lambda = matrix(0,nrow=1,ncol=300)
#' comfactor = matrix(0,nrow=1,ncol=30)
#'
#' #Choose how coefficients of the observable are estimated
#' method_estimate_beta = "individual" #estimating beta_i for every individual
#' beta_init = initialise_beta(use_robust = TRUE, NN = 300, TT = 30,
#'   number_of_variables = 3, number_vars_estimated = 3, number_of_groups = 3)
#' @importFrom stats lm
#' @export
initialise_beta <- function(use_robust, NN,
                            TT,
                            number_of_variables,
                            number_vars_estimated,
                            number_of_groups,
                            method_estimate_beta = "individual",
                            nosetting_lmrob = FALSE,
                            special_case_dgp1 = FALSE) {


  # if(use_real_world_data) {
  #   number_of_vars = number_of_variables
  # } else {
    number_of_vars = number_vars_estimated
  # }

  if(number_of_vars > 0) {
    beta_est = matrix(NA, nrow = (number_of_vars + 1), ncol = number_of_groups)
    if(method_estimate_beta == "individual") {
      beta_est = matrix(NA, nrow = (number_of_vars + 1), ncol = NN)
    }


    #X needs to be in the form of (NN*TT x p matrix)

    if(method_estimate_beta == "homogeneous") { #also used in DGP05: BramatiCroux
      X_special = restructure_X_to_order_slowN_fastT(X)
      Y_special = Y
      #this includes robust estimation of beta_est:
      beta_est = determine_beta("homogeneous", X_special, Y_special, use_robust, initialisation = TRUE,
                                indices = 1:NN,  TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting_lmrob,
                                special_case_dgp1)

    }
    if(method_estimate_beta == "group") {
      for(group in 1:number_of_groups) { #for each group there must be a column in beta_est.
        #select parts of X and Y of this group
        indices_group = which(g == group)
        if(length(indices_group) == 0) {
          message("There is an empty group!")
        }

        #X needs to be in the form of (NN*TT x p matrix)
        X_special = restructure_X_to_order_slowN_fastT(array(X[indices_group,,], dim = c(length(indices_group), TT, number_of_vars)),
                                                       number_of_variables = number_of_variables,
                                                       number_vars_estimated = number_vars_estimated)

        Y_special = as.vector(t(Y[indices_group,])) #order: N1T1, N1T2,N1T3,...N2T1,...N_endT_end

        beta_est[,group] = determine_beta("heterogeneous", X_special, Y_special, use_robust, initialisation = TRUE,
                                          indices = indices_group,  TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting_lmrob,
                                          special_case_dgp1)

      }

    }
    if(method_estimate_beta == "individual") {


      #Initialisation in classical case: use of lm instead of ncvreg (as in AndoBai-code)
      for(i in 1:NN) {
        X_special = restructure_X_to_order_slowN_fastT(matrix(X[i,,], ncol = dim(X)[3]),
                                                       number_of_variables = number_of_variables,
                                                       number_vars_estimated = number_vars_estimated)
        Y_special = as.vector(t(Y[i,]))

        #####################
        # define model
        #####################
        if(use_robust) { #robust -> lmrob
          # if(exists("use_bramaticroux")) {
          #   message(paste("bramati croux init",i))
          #   model = RpanFE(Y_special, X_special, TT, 0.20, 20, number_of_variables, 1)[[1]]
          #
          # } else {
          if(special_case_dgp1)  {
            #for special case of dgp 1 (= ando/bai dgp 2)
            #no intercept, because dgp1_spread_group_centers defines the first variable in X as an intercept
            model <- LMROB(Y_special, X_special, nointercept = TRUE, nosetting = nosetting_lmrob)  #-> lmrob(Y_special ~ X_special + 0, setting="KS2014)
          } else {
            model <- LMROB(Y_special, X_special, nosetting = nosetting_lmrob)
          }

          # }
        } else {  #classical -> lm
          if(special_case_dgp1) {
            model <- lm(Y_special ~ X_special + 0) #for special case of dgp 1 (= ando/bai dgp 2)
          } else {
            model <- lm(Y_special ~ X_special) #for the rest
          }

        }
        #####################
        #get the coefficients
        #####################

        if(special_case_dgp1) {
          beta_est[,i] = c(0, model$coefficients) #for dgp 1
        } else {
          beta_est[,i] = model$coefficients #for dgp 2
        }

      }

      #Note: in case of using USE_MEDIAN_OF_INDIVIDUAL_THETA:
      #this is not possible in this initialisation step, because g has not been initialized yet.
      #-> estimation of beta_est starts with the individual values

    }
    rm(X_special, Y_special)


    return(beta_est)
  } else {
    #when number_of_vars == 0, we do not need beta_est
    return(rep(NA, number_of_groups))
  }
}


#' Helpfunction used in update_g()
#'
#' We calculate FgLg (the groupfactorstructure) for all possible groups where individual i can be placed. For each group we have before estimated
#' the groupfactors (Fg). Now we need the grouploadings for each group as well. In the classical case these are calculated by Fg*Y/T. In the robust case
#' these are robust.
#' @param group number of groups
#' @param solve_FG_FG_times_FG This is the same as groupfactor / T. It is only used in the Classical approach
#' @param method_estimate_factors_local specifies the robust algorithm to estimate factors: default is "macro".
#' @param g vector with group memberships
#' @param NN_local N
#' @param TT_local T
#' @param number_of_variables_local number of observable variables
#' @param number_vars_estimated_local number of variables that are included in the algorithm and have their coefficient estimated. This is usually equal to number_of_variables.
#' @param number_of_group_factors_local number of group factors to be estimated
#' @param number_of_common_factors_local number of common factors to be estimated
#' @inheritParams update_g
#' @return NxT matrix containing the product of virtual groupfactors and virtual loadings

calculate_virtual_factor_and_lambda_group <- function(group, solve_FG_FG_times_FG, use_robust,
                                                      NN_local, TT_local,
                                                      method_estimate_factors_local,
                                                      g,
                                                      number_of_variables_local = number_of_variables_fixedvalue,
                                                      number_vars_estimated_local = number_vars_estimated_fixedvalue,
                                                      number_of_group_factors_local = number_of_group_factors_fixedvalue,
                                                      number_of_common_factors_local = number_of_common_factors_fixedvalue,
                                                      #use_real_world_data_local = use_real_world_data,
                                                      verbose = FALSE
                                                      ) {

  FG = factor_group[[group]]
  indices = 1:NN_local
  LF = t(lambda) %*% comfactor
  xbeta = calculate_XB_estimated(NN = NN_local, TT = TT_local, number_of_variables = number_of_variables_local, number_vars_estimated = number_vars_estimated_local)


  if(!do_we_estimate_common_factors(number_of_common_factors_local)) {
    Y_ster = Y[indices,] - xbeta
  } else {
    Y_ster = Y[indices,] - xbeta - LF
  }

  #robust grouplambda:
  if(use_robust) {
    #CHANGE LOCATION 5/7
    if(method_estimate_factors_local %in% c("macro", "pertmm")) {
      #we need a robust version of the virtual factorstructure:
      LG_local = return_robust_lambdaobject(Y_ster, group, type = 1, g,
                                            NN_rrn = NN_local,
                                             number_of_group_factors_rrn = number_of_group_factors_local,
                                             number_of_common_factors_rrn = number_of_common_factors_local,
                                            #use_real_world_data_rrn = use_real_world_data_local,
                                            #application_covid = exists("usecoviddata_cases") | exists("usecoviddata_deaths")
                                            )

    } else {
      LG_local = t(solve_FG_FG_times_FG[[group]] %*% t(Y_ster)) #This equalS to Fg*Y/T
    }
  } else {
    LG_local = t(solve_FG_FG_times_FG[[group]] %*% t(Y_ster)) #This equalS to Fg*Y/T
  }


  LG_local = handleNA_LG(LG_local)
  return(LG_local %*% FG)
}



#' Function to evade floating point errors.
#'
#' Sets values that should be zero but are >0 (e.g. 1e-13) on zero.
#' @param LIMIT limit under which value is set to 0
#' @param A input
#' @export
evade_floating_point_errors <- function(A, LIMIT = 1e-13) {
  for(i in 1:length(A)) {
    if(abs(A[i]) < LIMIT) A[i] = 0
  }
  return(A)
}

#' Determines parameters of rho-function.
#'
#' In robust updating of group membership we use a rho function (instead of the non-robust quadratic function) on the norm of the errors.
#' For this we need parameters of location and scale.
#' They are defined here (currently as median and madn).
#' This function is applied on the estimated errors: Y - XB - FL - FgLg.
#' This function is used in update_g().
#' @param object input
#' @importFrom stats mad
#' @importFrom stats median
define_rho_parameters <- function(object = NULL) {
  if(is.null(object)) {
    stop("This part should be obsolete (define_rho_parameters()). Sleep if encountered.")

  } else {
    rho_loc =  apply(data.frame(object), 1, function(x) median(x)) #=median over T elements
    rho_scale = apply(data.frame(object), 1, function(x) mad(x)) #=normalized mad over T elements


    #print(median(unlist(lapply(list(rho_loc,rho_scale),function(x) x[[1]][1]))))

  }

  return(list(rho_loc, rho_scale))
}

#' Helpfunction for update_g(). Calculates the errors for one of the possible groups time series can be placed in.
#'
#' As we are updating group membership, we use the errorterm as objective function to estimate the group. We assume group membership equals 1,...,NG
#' (with NG the total number of groups) and calculate the error term.
#' @param k group
#' @param LF NxT-matrix of the commonfactorstructure
#' @param virtual_grouped_factor_structure list with length the number of groups; every element of the list contains NxT-matrix
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @return NxT matrix with the errorterms (=dataset minus the estimated factorstructure and minus X*beta)
calculate_errors_virtual_groups <- function(k, LF, virtual_grouped_factor_structure,
                                            NN,
                                            TT,
                                            number_of_variables,
                                            number_of_common_factors,
                                            number_of_group_factors,
                                            number_vars_estimated = number_vars_estimated_fixedvalue
                                            ) {
  E_prep = matrix(NA, nrow = NN, ncol = TT)

  a = do_we_estimate_common_factors(number_of_common_factors)
  b = do_we_estimate_group_factors(number_of_group_factors)
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated)

  for(i in 1:NN) {
    #calculate lambda_group for one individual, based on an hypothetical groupmembership
    if(b == 1) { #-> we estimate group factors
      virtual_structure = virtual_grouped_factor_structure[[k]][i,]
    } else {
      message("Option to not estimate any group factor in any group -> this option is not implemented")
    }
    if(number_vars_estimated > 0) {
      if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) {
        XT = cbind(1, X[i,,]) %*% beta_est[,g[i]] #matrix with TT rows and 1 column
        #-> this should not be beta_est[,k] as only the grouped factorstructure should vary with k (similar to calculate_virtual_factor_and_lambda_group)
      }
      if(method_estimate_beta == "individual") {
        XT = cbind(1, X[i,,]) %*% beta_est[,i] #matrix with TT rows and 1 column
      }
    } else {
      XT = matrix(0, nrow = TT, ncol = 1)
    }


    for(t in 1:TT) {
      #calculate the objective function, while making sure that NA's do not have any effect in the sum
      #define the estimationerror:
      #print(paste(i,t))
      E_prep[i,t] = Y[i,t] - a * LF[i,t]
      # print("--")
      # print(E_prep[i,t])
      # print( virtual_structure[t])
      if(b != 0) E_prep[i,t] = (E_prep[i,t] - b * virtual_structure[t])
      # print(E_prep[i,t])
      #print(XT[t,])
      if(number_vars_estimated > 0) E_prep[i,t] = E_prep[i,t] - XT[t,] #note: XT is an TTx1 matrix
    }
  }
  return(E_prep)
}

#' Calculates objective function for individual i and group k in order to estimate group membership.
#'
#' Helpfunction in update_g().
#' Depends on an not yet established group k ( cannot use lgfg_list)
#' @param i individual
#' @param k group
#' @param ERRORS_VIRTUAL list with errors for each possible group
#' @param rho_parameters median and madn of the calculated error term
#' @param use_robust robust or classical estimation
#' @param TT T
calculate_obj_for_g <- function(i, k, ERRORS_VIRTUAL, rho_parameters, use_robust, TT = length_of_time_series) {

  totalsum = 0

  if(use_robust) {
    #define a scaling:
    #It must be individualspecific, and also be equal over all virtual groups;
    #We take here the median over the values of all groups.
    #     (because: use of #rho_parameters[[k]][[2]][i] leads tot random grouping)
    #"location" is then defined the same way
    #rho_parameters is a list of 'number_of_groups' elements. Every element has 2 elements with 'NN' values.
    #Note: do not use FUTURE_MAP() here: this is way slower.
    location = (unlist(lapply(rho_parameters,function(x) x[[1]][i]))) #map over virtual groups; take 1st element (=median) and take individual i
    scaling = (unlist(lapply(rho_parameters,function(x) x[[2]][i]))) #map over virtual groups; take 2nd element (=mad) and take individual i


    location = median(location) #median over groups
    scaling = median(scaling) #median over groups
    if(scaling == 0) { #This should be rare.
      warning("scaling equals 0 and is changed to 0.000001, to evade division by zero")
      scaling = 0.000001 #make sure no 0/0 will occur
    }
  }

  for(t in 1:TT) {

    #calculate the objective function, while making sure that NA's do not have any effect in the sum

    #calculate the estimationerror:
    E_prep = ERRORS_VIRTUAL[[k]][i,t]
    if(use_robust) { #define the rho-function (bisquare)

      E_prep_loc_scale = (E_prep - location) / scaling #this is a scalar
      E = Mpsi(E_prep_loc_scale, cc = 4.685, psi = "bisquare", deriv = -1) #rho-functie on scaled errors

    } else { #classical version:
      E = E_prep^2
    }


    E = evade_floating_point_errors(E) #evading floating point errors -> set to zero when values are really small

    totalsum = totalsum + as.numeric(E)
  }
  return(totalsum)
}

#' Helpfunction in update_g(), to calculate solve(FG x t(FG)) x FG
#'
#' @param TT length of time series
#' @param number_of_groups number of groups
#' @param number_of_group_factors number of group factors to be estimated
solveFG <- function(TT, number_of_groups, number_of_group_factors){
  solve_FG_FG_times_FG = list()
  for(group in 1:number_of_groups) {
    if(number_of_group_factors[group] > 0) {
      FG = factor_group[[group]]
      if(as.numeric(rankMatrix(FG)) != number_of_group_factors[group]) {
        message("___There is an issue in solveFG(): rank of FG should be the same as the number of group factors")
        message(paste("rank of matrix FG: ",as.numeric(rankMatrix(FG))))
        message(paste("number_of_group_factors[group]: ",number_of_group_factors[group]))
        print(FG[,1:3])
      }
      solve_FG_FG_times_FG[[group]] = solve(FG%*%t(FG))%*%FG #this is actually the same as FG/T
      rm(FG)
    } else {
      #make 0-matrix with 1 row
      solve_FG_FG_times_FG[[group]] = matrix(0, nrow = 1, ncol = TT)
    }
  }
  return(solve_FG_FG_times_FG)
}

#' Function that estimates group membership.
#'
#' @return list: 1st element contains group membership and second element contains the values which are used to determine group membership
#' @inheritParams estimate_beta
# @param use_class_zero if set to TRUE, then individuals with high distance to all possible groups are put in a separate class zero
# @param use_real_world_data_inupdateg Parameter to indicate using real world dataset. Defaults to FALSE.
#' @inheritParams estimate_factor
#' @param use_robust robust or classical estimation of group membership
#' @param verbose when TRUE, it prints messages
#' @examples
#' #This function needs several initial parameters to be initialized in order to work on itself
#' library(RCTS)
#' library(tidyverse)
#' library(robustbase)
#'
#' X = RCTS::X_dgp3
#' Y = RCTS::Y_dgp3
#' #Set estimations for group factors and its loadings, and group membership to the true value
#' lambda_group = RCTS::lambda_group_true_dgp3
#' factor_group = RCTS::factor_group_true_dgp3
#' g_true = RCTS::g_true_dgp3 #true values of group membership
#' g = g_true #estimated values of group membership; set in this example to be equal to true values
#' #There are no common factors to be estimated  ->  use placeholder with values set to zero
#' lambda = matrix(0,nrow=1,ncol=300)
#' comfactor = matrix(0,nrow=1,ncol=30)
#
#' grid = expand.grid(1:300,1:30)
#' #Choose how coefficients of the observable are estimated
#' method_estimate_beta = "individual"
#' method_estimate_factors = "macro"
#' number_vars_estimated_fixedvalue = 3
#' beta_est = estimate_beta(use_robust = TRUE, Y, X, lambda_group, factor_group,
#' lambda, comfactor, NN = 300, TT = 30,
#'   number_of_groups = 3, number_of_group_factors = c(3, 3, 3), number_of_common_factors = 0,
#'   number_of_variables = 3, number_vars_estimated = number_vars_estimated_fixedvalue,
#'   num_factors_may_vary = TRUE)[[1]]
#' grid = grid_add_variables(grid, beta_est, lambda, comfactor, method_estimate_beta,
#'   NN = 300, TT = 30,
#'   number_of_variables = 3, number_vars_estimated = number_vars_estimated_fixedvalue,
#'   number_of_groups = 3)
#' g_new = update_g(use_robust = TRUE, NN = 300, TT = 30, method_estimate_factors,
#'   number_of_groups = 3, number_of_variables = 3,
#'   number_vars_estimated = number_vars_estimated_fixedvalue,
#'   number_of_group_factors = c(3, 3, 3),
#'   number_of_common_factors = 0)[[1]]
#' @importFrom purrr map_dbl
#' @export
update_g <- function(use_robust, NN, TT, method_estimate_factors,
                     number_of_groups = number_of_groups_fixedvalue,
                     number_of_variables = number_of_variables_fixedvalue,
                     number_vars_estimated = number_vars_estimated_fixedvalue,
                     number_of_group_factors = number_of_group_factors_fixedvalue,
                     number_of_common_factors = number_of_common_factors_fixedvalue,
                     #use_real_world_data_inupdateg = use_real_world_data,
                     verbose = FALSE) {



  if(do_we_estimate_group_factors(number_of_group_factors)) { #if there are groupfactors estimated
    solve_FG_FG_times_FG = solveFG(TT, number_of_groups, number_of_group_factors)

    #we calculate FgLg (groupfactors times grouploadings) for all the possible groups in which individual i could end up:

    virtual_grouped_factor_structure = lapply(1:number_of_groups, function(y) calculate_virtual_factor_and_lambda_group(y, solve_FG_FG_times_FG, use_robust, NN, TT, method_estimate_factors, g,
                                                                                                                        number_of_variables_local = number_of_variables,
                                                                                                                        number_vars_estimated_local = number_vars_estimated,
                                                                                                                        number_of_group_factors_local = number_of_group_factors,
                                                                                                                        number_of_common_factors_local = number_of_common_factors #,
                                                                                                                        #use_real_world_data_local = use_real_world_data_inupdateg
                                                                                                                        ))
    if(verbose) message("virtual_grouped_factor_structure is created")


  } else {
    virtual_grouped_factor_structure = NA
  }


  LF = (t(lambda) %*% comfactor)


  #calculate errors for each possible group
  ERRORS_VIRTUAL = lapply(1:number_of_groups, function(x) calculate_errors_virtual_groups(x, LF, virtual_grouped_factor_structure, NN, TT,
                                                                                          number_of_variables,
                                                                                          number_of_common_factors,
                                                                                          number_of_group_factors,
                                                                                          number_vars_estimated = number_vars_estimated
                                                                                          ))
  if(verbose) message("ERRORS_VIRTUAL is created")


  if(use_robust) {
    rho_parameters = lapply(1:number_of_groups,function(x) define_rho_parameters(ERRORS_VIRTUAL[[x]])) #(parameter object = NA -> returns median and madn of the calculated error term)
  } else {
    rho_parameters = NA
  }
  if(verbose) message("rho_parameters is created")

  #init matrix with objectivefunctionvalues for all groups
  matrix_obj_values = matrix(NA, nrow = NN, ncol = number_of_groups)
  for(i in 1:NN) {
    obj_values = map_dbl(1:number_of_groups, function(x) calculate_obj_for_g(i, x, ERRORS_VIRTUAL, rho_parameters, use_robust, TT = TT))
    g[i] = which.min(obj_values)
    matrix_obj_values[i,] = obj_values
  } #Note: vectorizing is not faster: g = sapply(1:NN, function(z) which.min(map_dbl(1:number_of_groups, function(x) calculate_obj_for_g(z, x, ERRORS_VIRTUAL, rho_parameters, use_robust))))

  if(verbose) {
    print("current g-table is:")
    if(!is.na(g_true[1])) {
      print(table(g,g_true))
    } else {
      print(table(g))
    }
  }

  g_before_class_zero = g
  if(method_estimate_factors == "cz") {
    g = clustering_with_robust_distances(g, number_of_groups)
  }

  g = reassign_if_empty_groups(g, number_of_groups, TT)

  return(list(g, matrix_obj_values, g_before_class_zero))

}

#' Randomly reassign individual(s) if there are empty groups. This can happen if the total number of time series is low compared to the number of desired groups.
#'
#' @param g Vector of group membership for all individuals.
#' @param S_true true number of groups
#' @param TT length of time series
reassign_if_empty_groups <- function(g, S_true, TT) {
  empty_groups = c()
  while(length(table(g[g != 0])) != S_true) { #class zero (individuals too far from all groups) should not be counted here
    #reassign one individual to empty groups
    for(i in 1:S_true) {
      if(i %in% tibble::as_tibble(table(g))$g) { #if(i %in% broom::tidy(table(g))$g) {

      } else {
        empty_groups = c(empty_groups, i)
      }
    }

    for(i in 1:length(empty_groups)) {
      randomindex = ceiling(runif(1) * TT)
      g[randomindex] = empty_groups[i]
    }
  }
  return(g)
}


#' Function that puts individuals in a separate "class zero", when their distance to all possible groups is bigger then a threshold.
#'
#' @param g vector with group membership
#' @param number_of_groups number of groups
#' It starts with defining a robust location and scatter (based on Ma & Genton (2000): Highly robust estimation of the autocovariance function)
#' @return new clustering, including class zero
#' @importFrom graphics abline
#' @importFrom graphics plot
#' @importFrom stats qchisq
#' @importFrom stats quantile
#' @importFrom tsqn corQn
#' @importFrom Matrix rankMatrix
clustering_with_robust_distances <- function(g, number_of_groups) {
  RD = matrix(NA, nrow = number_of_groups, ncol = nrow(Y)) #every row contains the distances to a certain group, and every column is 1 individual
  for(group in 1:number_of_groups) { #loop over all groups
    mu = apply(Y[g==group,], 2, "median")
    s = apply(Y[g==group,], 2, "Qn")
    x = Y[g==group, 2:ncol(Y)]
    laggedx = Y[g==group,1:(ncol(Y)-1)]
    r = corQn( as.vector(x), as.vector(laggedx)) #Computes the robust correlation of x and y

    ar1_cor <- function(n, rho) {
      exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                        (1:n - 1))
      rho^exponent
    }

    if(rankMatrix(diag(s))[[1]] < ncol(Y)) {
      message("rank of diag(s) is too small (reason: because Qn(.) equals zero for 1 year)")  #-> solve(sig) would generate error
      s[s==0] <- 1e-6
    }

    R = ar1_cor(n=ncol(Y),rho=r)
    sig = diag(s)%*%R%*%diag(s)
    sig_inv = solve(sig)
    #calculate robust distances for all individuals
    for(i in 1:nrow(Y)) {
      RD[group,i] = sqrt(t(Y[i,] - mu) %*% sig_inv %*% (Y[i,] - mu))
    }
  }

  #define limit as 99%-quantile of chi-squared distribution:
  #if for an individual the robust distance for every group is larger then the limit, the individual is put to class zero
  limit = sqrt(qchisq(0.99, ncol(Y), ncp = 0))

  minimum_distance = apply(RD, 2, min)
  limit_empirical = quantile(minimum_distance, 0.75) #
  cases = which(minimum_distance > max(limit,limit_empirical)) #those cases go to class zero

  # message("limits are:")
  # print(limit)
  # print(limit_empirical)
  if(percent_outliers > 0) {
    #plot with colored outlierindividuals;
    indices_of_outliers = which(apply(Y, 1, function(x) abs(max(x))) > 900) #indices of individuals with at least 1 generated outlier in it
    plot(log(apply(RD, 2, min)), col = ((1:nrow(Y)) %in% indices_of_outliers) + 1, main = "minimal log(distance) of i to any group")
    abline(h = log(limit), col="red")
    abline(h = log(limit_empirical),col="orange")
  } else {
    #plot without coloring
    plot(apply(RD,2,min),main = "minimal distance of i to any group", ylim = c(0, limit * 1.05))
    abline(h = limit, col = "red")
    abline(h = limit_empirical, col="orange")
  }
  g[cases] = 0
  print(table(g))
  return(g)
}





#' Returns the estimated groupfactorstructure.
#'
#' This is the same function as calculate_FL_group_estimated(), but with adjustable parameters
#' @param f groupfactors
#' @param l grouploadings
#' @param j Number of groupfactors that are being used in the calculation.
#' @inheritParams calculate_FL_group_true
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @export
calculate_FL_group_estimated2 <- function(f,l,j,
                                          number_of_groups = number_of_groups_fixedvalue,
                                          number_of_common_factors = number_of_common_factors_fixedvalue,
                                          num_factors_may_vary = num_factors_may_vary_fixedvalue,
                                          NN = number_of_time_series,
                                          TT = length_of_time_series) {

  temp = calculate_lgfg(l,f,number_of_groups, rep(j, number_of_groups), number_of_common_factors, num_factors_may_vary)
  FL_GROUP_geschat = lapply(1:NN, function(x) temp[[g[x]]][x,]) %>% unlist %>% matrix(nrow = TT) %>% t

  return(FL_GROUP_geschat)
}










#' Helpfunction in OF_vectorized3()
#'
#' @param i index of individual
#' @param t index of time
#' @param XBETA matrixproduct of X and beta_est
#' @param LF matrixproduct of common factors and its loadings
#' @param group_memberships vector with group memberships
#' @param lgfg_list product of groupfactors and their loadings; list with length the number of groups
#' @param number_of_group_factors vector containing the number of group factors to be estimated for all groups
OF_vectorized_helpfunction3 <- function(i, t, XBETA, LF,
                                        group_memberships,
                                        lgfg_list,
                                        number_of_group_factors) {


  if(do_we_estimate_group_factors(number_of_group_factors) != 0 & group_memberships[i] != 0) {

    if(t > ncol(lgfg_list[[group_memberships[i]]])) { #this is the case when macropca() dropped columns
      warning("Unsolved issue: macropca has dropped columns (OF_vectorized_helpfunction3())")
      print(t)
      print(dim(lgfg_list[[group_memberships[i]]]))
    } else {
      result = as.numeric(Y[i,t] - XBETA - LF -
                            lgfg_list[[group_memberships[i]]][i,t]
      )^2
    }


  } else {
    result = as.numeric(Y[i,t] - XBETA -
                          LF
    )^2
  }

  return(result)
}

#' Calculates objective function: used in local_search + to determine "best_result".
#'
#' @param group_memberships Vector containing the group membership for all individuals.
#' @param grid dataframe containing XB, FgLg and FL
#' @param beta_est estimated values of beta
#' @param fc estimated common factors
#' @param lc loadings of estimated common factors
#' @param fg estimated groupfactors
#' @param lg estimated grouploadings
#' @inheritParams estimate_beta
#' @importFrom tidyselect starts_with
#' @export
OF_vectorized3 <- function(group_memberships, grid, beta_est = beta_est,
                           lc = lambda, fc = comfactor,
                           lg = lambda_group, fg = factor_group,
                           NN = number_of_time_series,
                           number_of_groups = number_of_groups_fixedvalue,
                           number_of_common_factors = number_of_common_factors_fixedvalue,
                           number_of_group_factors = number_of_group_factors_fixedvalue,
                           num_factors_may_vary = num_factors_may_vary_fixedvalue) {
  #this is a list (length number of groups) of the product FgLg (which is the groupfactorstructure)
  lgfg_list = calculate_lgfg(lg, fg, number_of_groups, number_of_group_factors, number_of_common_factors, num_factors_may_vary)

  if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "individual")) {
    return(sum(apply(grid, 1, function(x) OF_vectorized_helpfunction3(x[1], x[2], x[3], x[4], group_memberships, lgfg_list, number_of_group_factors))))
  }
  if(method_estimate_beta == "group") {
    #construct a vector with XBETA-values depending on the group of the individuals:
    temp = grid %>% dplyr::select(starts_with("XBETA"))
    XBETA_parameter = sapply(1:NN, function(x) temp[x,g[x]])

    prep = grid %>%
      dplyr::select(-starts_with("XBETA"))
    #used to put LF into the 3rd column -> x[3] in next line
    return(sum(apply(prep, 1, function(x) OF_vectorized_helpfunction3(x[1],x[2],XBETA_parameter,x[3],group_memberships, lgfg_list, number_of_group_factors))))
  }


}






#' Returns list (with as length the number of groups) with lgfg (product of grouploadings a&nd groupfactors).
#' Each element of the list with the assumption that all individuals are in the same group k.
#'
#' This function is used to speed up code.
#' @param lambda_group grouploadings
#' @param factor_group groupfactors
#' @inheritParams estimate_beta
#' @importFrom dplyr arrange
#' @examples
#' library(tidyverse)
#' lambda_group = RCTS::lambda_group_true_dgp3
#' factor_group = RCTS::factor_group_true_dgp3
#' calculate_lgfg(lambda_group,factor_group, 3, c(3, 3, 3), 0, FALSE)
#' @importFrom rlang .data
#' @export
calculate_lgfg <- function(lambda_group, factor_group, number_of_groups, number_of_group_factors, number_of_common_factors, num_factors_may_vary,
                           NN = number_of_time_series, TT = length_of_time_series) {
  lgfg_list = list()
  #define LgFg for each group (and later on select the correct element (correct group of individual i) of lgfg_list)
  for(k in 1:number_of_groups) {
    if(number_of_group_factors[k] > 0) {
      LG_clean = (as.matrix(lambda_group %>% arrange(.data$id) %>% dplyr::select(-.data$groep, -.data$id)))[, 1:number_of_group_factors[k]]

      #When using a varying number of groupfactors per group,  and also estimating a positive number of common factors, lgfg_list will contain NA's and crash the algorithm.
      #Reason is that there are NA's in lambda_group (only for the groups that have not the maximum amount of groupfactors), which is fine in itself.
      #-> Since those NA's mean that there is no effect on LGFG (since there is no factorpart and lambdapart for those), we can replace those NA's by zero's.
      if(num_factors_may_vary & number_of_common_factors > 0) {
          if(anyNA(LG_clean)) {
            temp = LG_clean
            temp[is.na(temp)] <- 0
            LG_clean = temp
          }
      }


      lgfg_list[[k]] = LG_clean %*% factor_group[[k]]
    } else {
      lgfg_list[[k]] = NA
    }
  }


  #replace parts of lgfg_list which are NA (because of no groupfactors in that particular group) by 0-matrices
  emptyFL = sapply(lgfg_list,function(x) is.null(dim(x)))
  for(q in 1:length(lgfg_list)) {
    if(emptyFL[q]) {
      lgfg_list[[q]] = matrix(0, nrow=NN, ncol=TT)
    }
  }


  return(lgfg_list)
}


#' Helpfunction in estimate_beta() for estimating beta_est.
#'
#' @param string can have values: "homogeneous" (when one beta_est is estimated for all individuals together) or "heterogeneous" (when beta_est is estimated either groupwise or elementwise)
#' @param X_special preprocessed X (observable variables)
#' @param Y_special preprocessed Y
#' @param use_robust robust or classical estimation
#' @param initialisation indicator of being in the initialisation phase
#' @param indices individuals for which beta_est is being estimated
#' @param optimize_kappa indicates if kappa has to be optimized or not (only relevant for the classical algorithm); defaults to FALSE
#' @param nosetting_local option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
#' @param kappa_candidates defines the size of the SCAD-penalty used in the classical algorithm
#' @inheritParams initialise_beta
#' @importFrom ncvreg ncvreg
#' @inheritParams generate_Y
#' @inheritParams LMROB
determine_beta <- function(string, X_special, Y_special, use_robust, initialisation = FALSE, indices = NA, optimize_kappa = FALSE,
                           TT = length_of_time_series, NN = number_of_time_series,
                           number_of_variables = number_of_variables_fixedvalue, nosetting_local = FALSE, kappa_candidates = c(0.1),
                           special_case_dgp1 = FALSE) {
  stopifnot(string == "homogeneous" | string == "heterogeneous")

  if(!(method_estimate_beta == "group" & initialisation == TRUE)) {

    Y_special = matrix(Y_special, nrow = length(indices), ncol = TT)

    if(initialisation) {
      #initialisation of beta_est, so no factorstructure in Y_special
      Y_special = as.vector(t(Y)) #order: N1T1, N2T1,...N1T2,...N_endT_end
    } else {
      Y_special = as.vector(t(Y_special)) #order: N1T1, N2T1,...N1T2,...N_endT_end
    }

    #remove NA's from regression:
    vectorNA = which(is.na(Y_special))

    if(length(vectorNA) > 0) {
      X_special = X_special[-vectorNA]
      Y_special = Y_special[-vectorNA]
    }
  }

  #Regression:
  if(use_robust) {
    # if(exists("use_bramaticroux")) {
    #   model = RpanFE(Y_special, X_special, TT, 0.20, 20, number_of_variables, length(Y_special)/TT)[[1]]
    # } else {
      model <- LMROB(Y_special, X_special, nosetting = nosetting_local) #-> lmrob(Y_special ~ X_special, setting = "KS2014")
    # }
  } else {


    #ncvreg, without weights
    model <- ncvreg(X_special, Y_special, family = "gaussian", penalty = "SCAD", lambda = (kappa_candidates))
    if(optimize_kappa) {
      minBIC <- 10^10
      best_lami = 0
      #print("indices")
      #print(indices)
      for(LAMi in 1:length(kappa_candidates)){
        beta_temp <- model$beta[, LAMi] #c(0,model$beta)[-2,LAMi]
        BIC <- sum( (Y[indices,] - cbind(1, X[indices,,]) %*% beta_temp )^2 )/(TT) + C * sigma2_max_model * log(TT) * sum(beta_temp != 0)/(TT)

        if(BIC <= minBIC){
          best_lami = LAMi
          minBIC <- BIC
          beta_temp <- model$beta[, LAMi] #c(0,model$beta)[-2,LAMi]
        }
      }
      #print(kappa_candidates[best_lami]) #-> this is kappa[i]
      #Sys.sleep(0.1)
      return(beta_temp)
    }

  }


  if(string == "homogeneous") {
    if(class(model) == "lm" | class(model) == "lmrob") return(matrix(rep(model$coefficients, number_of_groups_fixedvalue), nrow = (number_of_variables + 1)))
    else return(matrix(rep(model$beta[,1], number_of_groups_fixedvalue), nrow = (number_of_variables + 1)))
  }
  if(string == "heterogeneous") {
    if(class(model) == "lm" | class(model) == "lmrob") {
      if(special_case_dgp1) {
        return(c(0,as.numeric(model$coefficients[-2])))
      } else {
        return(as.numeric(model$coefficients))
      }

    } else { #ncvreg or rpanfe

      # if(exists("use_bramaticroux")) {
      #   return( c(0,model) )
      # }
      if(special_case_dgp1) {
        return(c(0,model$beta[-2,1]))
      } else {
        return(model$beta[,1])  #for old dgp 2
      }

    }
  }

}



#' Estimates beta.
#'
#' Update step of algorithm to obtain new estimation for beta. Note that we call it beta_est because beta() exists in base R.
#' @inheritParams determine_beta
# @param use_real_world_data Parameter to indicate using real world dataset. Defaults to FALSE.#' @param lambda_group_true loadings of the group factors
#' @param use_robust TRUE or FALSE: defines using the classical, robust algorithm to estimate beta
#' @param Y Y
#' @param X X
#' @param factor_group estimated group specific factors
#' @param lambda_group loadings of the group specific factors
#' @param comfactor estimated common factors
#' @param lambda loadings of the common factors
#' @param NN number of individuals
#' @param TT length of time series
#' @param number_of_groups number of groups estimated
#' @param number_of_group_factors number of groupfactors to be estimated
#' @param number_of_common_factors number of common factors to be estimated
#' @param number_of_variables number of observable variables
#' @param number_vars_estimated number of variables that are included in the algorithm and have their coefficient estimated. This is usually equal to number_of_variables.
#' @param num_factors_may_vary whether or not the number of groupfactors is constant over all groups or not
#' @param nosetting option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
#' @inheritParams define_object_for_initial_clustering_macropca
#' @inheritParams initialise_beta
#' @return list: 1st element contains matrix (N columns: 1 for each element of the panel data) with estimated beta_est's.
#' @examples
#' #This function needs several initial parameters to be initialized in order to work on itself.
#' library(RCTS)
#' library(tidyverse)
#' library(robustbase)
#'
#' X = RCTS::X_dgp3
#' Y = RCTS::Y_dgp3
#' #Set estimations for group factors and its loadings, and group membership to the true value
#' lambda_group = RCTS::lambda_group_true_dgp3
#' factor_group = RCTS::factor_group_true_dgp3
#' g = RCTS::g_true_dgp3
#' #There are no common factors to be estimated  -> but needs placeholder
#' lambda = matrix(0,nrow=1,ncol=300)
#' comfactor = matrix(0,nrow=1,ncol=30)
#
#' #Choose how coefficients of the observable variables are estimated
#' method_estimate_beta = "individual"
#' method_estimate_factors = "macro"
#' beta_est = estimate_beta(use_robust = TRUE, Y, X, lambda_group, factor_group,
#'   lambda, comfactor, NN = 300, TT = 30,
#'   number_of_groups = 3, number_of_group_factors = c(3, 3, 3), number_of_common_factors = 0,
#'   number_of_variables = 3,number_vars_estimated = 3, num_factors_may_vary = FALSE)[[1]]
#' @importFrom stats filter
#' @importFrom purrr pmap
#' @importFrom purrr map2
#' @importFrom rlang .data
#' @export
estimate_beta <- function(use_robust, Y, X, lambda_group, factor_group, lambda, comfactor,
                          method_estimate_beta = "individual",
                          NN = nrow(Y),
                           TT = ncol(Y),
                           number_of_groups = number_of_groups_fixedvalue,
                           number_of_group_factors = number_of_group_factors_fixedvalue,
                           number_of_common_factors = number_of_common_factors_fixedvalue,
                           number_of_variables = number_of_variables_fixedvalue,
                           number_vars_estimated = number_vars_estimated_fixedvalue,
                           num_factors_may_vary = num_factors_may_vary_fixedvalue,
                           optimize_kappa = FALSE, nosetting = FALSE,
                          special_case_dgp1 = FALSE) {
  if(number_vars_estimated > 0) {


    if(method_estimate_beta == "homogeneous") {
      #X needs to be in the form of (NN*TT x p matrix)
      X_special = restructure_X_to_order_slowN_fastT(X)
      #define Y* as Y - FcLc - FgLg:
      #and get LgFg out of the loop to speed up:
      lgfg_list = calculate_lgfg(lambda_group, factor_group, number_of_groups, number_of_group_factors, number_of_common_factors, num_factors_may_vary)

      Y_special = Y
      for(i in 1:NN) {
        #choose the correct lgfg_list, based on group of i
        for(k in 1:number_of_groups) {
          if(g[i] == k) temp = lgfg_list[[k]]
        }

        #subtract common factorstructure
        Y_special[i,] = Y_special[i,] - t(lambda[,i]) %*% comfactor[,]
        if(max(number_of_group_factors, na.rm=T) > 0) {
          #subtract group factorstructure
          Y_special[i,] = Y_special[i,] - temp[i,]
        }
      }

      beta_est = determine_beta("homogeneous", X_special, Y_special, use_robust,
                                indices = 1:NN, TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting,
                                special_case_dgp1)

    }
    if(method_estimate_beta == "group") {
      beta_est = matrix(NA, nrow = (number_of_variables + 1), ncol = number_of_groups)
      for(group in 1:number_of_groups) {
        #select parts of X and Y of this group
        indices_group = which(g == group)
        if(length(indices_group) == 0) {
          warning("Empty group! This message should never be seen.")
        }

        #X needs to be in the form of (NN*TT x p matrix)
        X_special = restructure_X_to_order_slowN_fastT(array(X[indices_group,,], dim = c(length(indices_group), TT, number_of_variables)),
                                                       number_of_variables = number_of_variables,
                                                       number_vars_estimated = number_vars_estimated)
        #define Y* as Y - FcLc - FgLg:
        Y_special = matrix(Y[indices_group,], nrow = length(indices_group)) %>% unlist


        for(i in 1:length(indices_group)) {
          index = indices_group[i]
          LAMBDAGROUP = as.matrix(lambda_group %>% dplyr::filter(.data$id %in% index) %>% dplyr::select(-.data$groep, -.data$id))
          Y_special[i,] = Y_special[i,] -
            t(lambda[,index]) %*% comfactor[,] -
            LAMBDAGROUP %*% factor_group[[g[index]]]
        }


        beta_est[,group] = determine_beta("heterogeneous", X_special, Y_special, use_robust,
                                          indices = indices_group, TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting,
                                          special_case_dgp1)

      }
    }
    if(method_estimate_beta == "individual") {
      #######################################################################
      #first, X and Y need to be restructured: made into list with N elements.
      #######################################################################

      X_local = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated) #makes X smaller when number_vars_estimated < number_of_variables
      # if(use_real_world_data) {
      #   number_of_vars = number_of_variables
      # } else {
        number_of_vars = number_vars_estimated
      # }
      X_special_list = lapply(1:NN, function(x) restructure_X_to_order_slowN_fastT(matrix(X_local[x,,], ncol = number_of_vars),
                                                                                   number_of_variables = number_of_variables,
                                                                                   number_vars_estimated = number_vars_estimated))
      #helpfunction: calculates Y - factorstructures. This is equal to XB + errorterm.
      #The result is a list with number of elements equal to number of time series in Y.
      make_Y_special <- function(i) {
        Y_special = unlist(matrix(Y[i,], nrow = 1))
        if(dim(t(Y_special))[2] == TT) {
          Y_special = t(Y_special)
        }

        to_remove = c(
          which(names(lambda_group) == "groep"),
          which(names(lambda_group) == "id")
        )
        LAMBDAGROUP = as.matrix(lambda_group[which(lambda_group$id %in% i),][,-to_remove])
        if(g[i] != 0) {
          LAMBDAGROUP = LAMBDAGROUP[1:number_of_group_factors[g[i]]]
          Y_special = Y_special -
            t(lambda[,i]) %*% comfactor[,] -
            LAMBDAGROUP %*% factor_group[[g[i]]]
        } else { #class zero -> no groupfactorstructure
          LAMBDAGROUP = NA
          Y_special = Y_special -
            t(lambda[,i]) %*% comfactor[,]
        }

        return(Y_special)
      }

      Y_special_list = lapply(1:NN, function(x) make_Y_special(x) )

      #######################################################################
      #secondly, beta_i is estimated for all i
      #######################################################################
      if(optimize_kappa) {
        beta_est = pmap(list(X_special_list, Y_special_list, 1:NN),
                        function(x, y, z) determine_beta("heterogeneous", x, y, use_robust,
                                                         indices = z,  TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting,
                                                         special_case_dgp1) )
      } else {
        #note that mapply isnt faster than map2
        beta_est = map2(X_special_list, Y_special_list,
                        function(x, y) determine_beta("heterogeneous", x, y, use_robust,
                                                      indices = NA,  TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting,
                                                      special_case_dgp1) )
        # print(summary(c( matrix(unlist(beta_est), ncol = NN)  - beta_new)))
        # print(sum(c( matrix(unlist(beta_est), ncol = NN)  - beta_new)))
        #
        # micro = microbenchmark::microbenchmark(
        #   beta_est = map2(X_special_list, Y_special_list,  function(x,y) determine_beta("heterogeneous", x, y, use_robust, indices = NA,  TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting, special_case_dgp1) ),
        #   beta_new = mapply( function(x,y) { determine_beta("heterogeneous", x, y, use_robust, indices = NA,  TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting, special_case_dgp1) }, x = X_special_list, y = Y_special_list ),
        #   times = 15
        # )
        # print(summary(micro))
        # expr      min       lq     mean   median       uq      max neval cld
        # 1     beta_est 3.544195 3.585846 3.643604 3.622231 3.672088 3.829970    15   a
        # 2 beta_new 3.543808 3.553852 3.621918 3.597413 3.680516 3.731678    15   a
      }
      ########################################################
      #Possible use of future::map2 instead of map2:
      #(is 2 to 3 times faster (microbenchmark) for (300, 30)-system)
      #BUT:  ERRORS:
      #   Error: 'map2' is not an exported object from 'namespace:future'
      #beta_est = future::map2(X_special_list, Y_special_list,  function(x,y) determine_beta("heterogeneous", x, y, use_robust, indices = NA,  TT = TT, number_of_variables = number_of_variables, nosetting_local = nosetting, special_case_dgp1) )
      ########################################################

      beta_est = matrix(unlist(beta_est), ncol = NN)

    }


    return(list(beta_est))
  } else {
    return(list(NA))
  }
}


#' Calculates W = Y - X*beta_est. It is used in the initialization step of the algorithm, to initialise the factorstructures.
#'
#' @param beta_est estimated values of beta
#' @param g vector with group membership
#' @inheritParams estimate_beta
#' @return NxT matrix
#' @export
calculate_W <- function(beta_est, g ,
                        NN = number_of_time_series,
                        TT = length_of_time_series,
                        number_of_variables = number_of_variables_fixedvalue,
                        number_vars_estimated = number_vars_estimated_fixedvalue) {
  W = matrix(0, nrow = NN, ncol = TT) #T x N-matrix in original paper , but I rather define as NxT

  #if number_vars_estimated < number_of_variables the obsoleterows in beta_est were already erased -> do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated)

  if(number_vars_estimated > 0) {
    if(method_estimate_beta == "homogeneous") {
      #non-dependence on g -> take first column
      for(i in 1:NN) W[i,] = Y[i,] - t(cbind(1,X[i,,]) %*% as.matrix(beta_est[,1]))
    }
    if(method_estimate_beta == "group") {
      for(i in 1:NN) W[i,] = Y[i,] - t(cbind(1,X[i,,]) %*% as.matrix(beta_est[,g[i]]))
    }
    if(method_estimate_beta == "individual") {
      for(i in 1:NN) {
          W[i,] = Y[i,] - t(cbind(1,X[i,,]) %*% as.matrix(beta_est[,i]))
      }
    }
  } else {
    for(i in 1:NN) W = as.matrix(Y)
  }

  if(anyNA(W)) {

    #handle NA's: which currently means: remove the rows with NA's
    W = handleNA(W)[[1]]
  }
  return(W)
}

#' Calculates Z = Y - X*beta_est - LgFg. It is used in the estimate of the common factorstructure.
#' @inheritParams calculate_W
#' @inheritParams estimate_factor
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param initialise boolean
#' @inheritParams estimate_beta
#' @export
calculate_Z_common <- function(beta_est, g, lgfg_list,
                               method_estimate_factors,
                               initialise = FALSE,
                               NN = number_of_time_series,
                               TT = length_of_time_series,
                               number_of_variables = number_of_variables_fixedvalue,
                               number_vars_estimated = number_vars_estimated_fixedvalue,
                               number_of_group_factors = number_of_group_factors_fixedvalue
                               ) {
  Z = matrix(0, nrow = NN, ncol = TT) #T x N-matrix in paper , maar ik definieer liever als NxT

  #if number_vars_estimated < number_of_variables the obsolete rows in beta_est were already erased -> do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated)

  for(i in 1:NN) {
    y = Y[i,] %>% as.numeric
    if(do_we_estimate_group_factors(number_of_group_factors)) { #when there are group factors to be estimated

      if(method_estimate_factors == "pertmm" & class(lgfg_list) != "list") {
        #then lgfg_list is not a list, but a data frame of dimension NxT
        LF_GROUP = lgfg_list[i,]
      } else {
        if(g[i] == 0) { #class zero
          LF_GROUP = 0
        } else {
          LF_GROUP = lgfg_list[[g[i]]][i,]
        }
      }
    } else {
      LF_GROUP = 0
    }

    if(number_vars_estimated > 0) {
      if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) BETA = as.matrix(beta_est[,g[i]])
      if(method_estimate_beta == "individual") BETA = as.matrix(beta_est[,i])
      Z[i,] = y - t(cbind(1,X[i,,]) %*% BETA) - LF_GROUP
    } else {
      Z[i,] = y - LF_GROUP
    }
  }


  if(anyNA(Z)) {
    #handle NA's: which currently means: remove the rows with NA's
    Z = handleNA(Z)[[1]]
  }
  return(Z)
}

#' Calculates Z = Y - X*beta_est - LF. It is used to estimate the groupfactorstructure.
#' @inheritParams calculate_W
#' @inheritParams estimate_beta
#' @param lambda common factor loadings
#' @param comfactor common factors
#' @param group the number of the group
#' @param initialise boolean
#' @export
calculate_Z_group <- function(beta_est, g, lambda, comfactor, group, initialise,
                              TT = length_of_time_series,
                              number_of_variables = number_of_variables_fixedvalue,
                              number_vars_estimated = number_vars_estimated_fixedvalue,
                              number_of_common_factors = number_of_common_factors_fixedvalue) {

  indices_group = which(g == group)
  if(length(indices_group) == 0) {
    message("empty group (calculate_Z_group())")
    Sys.sleep(1)
  }

  Z = matrix(0, nrow = length(indices_group), ncol = TT) #Nj x T matrix
  #if number_vars_estimated < number_of_variables the obsoleterows in beta_est were already erased -> do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated)

  for(i in 1:length(indices_group)) { #loop over the number of elements in the group
    index = indices_group[i]

    #define XT
    if(number_vars_estimated > 0) {
      if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group"))  BETA = as.matrix(beta_est[,g[index]])
      if(method_estimate_beta == "individual") BETA = as.matrix(beta_est[,index])

      XT = t(cbind(1,X[index,,]) %*% BETA)
    } else {
      XT = rep(0,TT)
    }

    #calculate Z = Y - X*beta_est - LF)
    y = Y[index,] %>% as.numeric
    if(initialise) {
      Z[i,] = y - XT
    } else {
      a = do_we_estimate_common_factors(number_of_common_factors)
      Z[i,] = y - XT -
        a * t(lambda[,index]) %*% comfactor
    }
  }


  if(anyNA(Z)) {
    #handle NA's: which currently means: remove the rows with NA's
    Z = handleNA(Z)[[1]]
  }

  return(Z)
}

#' Solves a very specific issue with MacroPCA.
#'
#' MacroPCA crashes Rstudio with certain dimensions of the input. Solve this by doubling every row.
#' No information is added by this, so no influence on end result,
#' but crashes of Rstudio are evaded.
#' @param object input
#' @param verbose prints messages
evade_crashes_macropca <- function(object, verbose = FALSE) {
  #--------------------MacroPCA seems to make Rstudio crash when the dimension of object = (27,193) ------------------
  size_with_crashes1 = c( 27,  27,  25,  43,  66,  24)#,  44,  41)
  size_with_crashes2 = c(193, 310, 600, 560, 387, 374)#, 201, 201)
  changed_objectsize = FALSE

  if(verbose) message("Test size of object (for MacroPCA)")

  for(i in 1:length(size_with_crashes1)) {
    if(changed_objectsize == FALSE & dim(object)[1] == size_with_crashes1[i] & dim(object)[2] == size_with_crashes2[i]) {
      message("Rstudio would crash -> double amount of rows")
      print(i)
      object = rbind(object, object)
      message(dim(object))
      changed_objectsize = TRUE
    }
  }
  if(verbose) message("Test size: done")
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
handle_macropca_errors <- function(object, temp, KMAX, number_eigenvectors, verbose = FALSE) {

  if("error" %in% class(temp)) {
    message("*******************************")
    message(paste("MacroPCA with .... fails, -> use different amount of eigenvectors. Start with a couple more and decrease 1 by 1."))

    teller = 0
    while("error" %in% class(temp)) {
      teller = teller + 1
      if(verbose) print(paste("teller:",teller))


      temp =  tryCatch(
          cellWise::MacroPCA(object, k = max(12, number_eigenvectors) - teller, MacroPCApars = list(kmax=KMAX)),
          error = function(e) { message(e); return(e) }
      )

      #sometimes, when there occurred too often errors in macropca, we end up with only a limited amount of possible eigenvectors that are calculated.
      #This happens when a group has very little elements.
      if(verbose) {
        print("----")
        print(class(temp))
      }
      number_columns = 999
      temp =  tryCatch(
        number_columns = ncol(temp$loadings),
        error = function(e) { message(e); return(e) }
      )
      temp =  tryCatch(
          temp = temp$loadings[,1:min(number_eigenvectors, number_columns)],
          error = function(e) { message(e); return(e) }
      )
      if(number_columns < number_eigenvectors) {
        message("---- This should be solved. Does this still occur? ----")
        message(paste("add",number_eigenvectors - number_columns,"columns to temp"))
        print(temp)
        message("sleep...")
        Sys.sleep(9000)
      }

      if(teller >= 10) {
        warning("-------------infinite loop (MacroPCA does not work with any k) -> use (classical) eigen(): has to be squared matrix -> take covariance matrix of object-------------")
        #(reason: Error in svd::propack.svd(Y, neig = min(n, d)) : BLAS/LAPACK routine 'DLASCL' gave error code -4)
        #eigen() needs always square matrix -> take covmatrix
        temp = eigen(t(object)%*%object)$vectors[,1:number_eigenvectors]

      }
    }
  }
  return(temp)
}

#' Function that uses robust PCA and estimates robust factors and loadings.
#'
#' Contains MacroPCA
#'
#' Notes for MACROPCA:
#'
#' Different values for kmax give different factors, but the product lambda*factor stays constant. Note that
#'  this number needs to be big enough, otherwise eigen() will be used. Variation in k does give different results for lambda*factor
#'
#' NOTE: this function may (not certain) crash when ncol(object) >> nrow(object)
#'   Actually it crashes with specific values of dim(object). For example when dim(object) = c(193,27).
#'   This is solved with evade_crashes_macropca(), for those problematic dimensions that are already encountered.
#' @param object input
#' @param number_eigenvectors number of eigenvectors to extract
#' @param KMAX The maximal number of principal components to compute. This is a paramater in cellWise::MacroPCA()
#' @param verbose_robustpca when TRUE, it prints messages
#' @export
robustpca <- function(object, number_eigenvectors, KMAX = 20, verbose_robustpca = FALSE) {

  if(verbose_robustpca) {
    print(paste("*************************************************robust PCA with:",number_eigenvectors))
    print("dimension of input:")
    print(dim(object))
  }


  ######################
  # MacroPCA
  ######################
  error_macropca = FALSE
  object = evade_crashes_macropca(object, verbose = verbose_robustpca)

    #print(number_eigenvectors)
    if(number_eigenvectors > KMAX) {
      #Note that when k > kmax, k gets the value of kmax.
      message("MacroPCA is (through KMAX) limited to 20 factors.")
    }

    macropca_kmax = number_eigenvectors

    #note that in the documentation it says the default for scale is FALSE
    #But:
    #
    # MacroPCA(M,1,MacroPCApars=(list(scale=T)))$loadings
    # MacroPCA(M,1,MacroPCApars=(list(scale=F)))$loadings
    # MacroPCA(M,1)$loadings
    # -> 1st and 3rd give the same result
    # -> The actual default for scale seems to be TRUE

    # if(exists("macropca_screeplot")) { #plots the screeplot
    #   cumul_var = cellWise::MacroPCA(object)$cumulativeVar
    #   Sys.sleep(5)
    # }

    if(verbose_robustpca) {
      message("--start of macropca")
      print(paste("rank of input:", rankMatrix(object)))
      print(paste("required number of eigenvectors:",number_eigenvectors))
    }
    temp =  tryCatch(
      cellWise::MacroPCA(object, k = max(macropca_kmax, number_eigenvectors), MacroPCApars = list(kmax=KMAX)),
      error = function(e) { message(e); return(e) }
    )
    if("error" %in% class(temp)) {
      error_macropca = TRUE
    }
    ###############################
    #sometimes the output of MacroPCA() has too little columns (should be equal to "number_eigenvectors").
    #Possible related to performing MacroPCA on a very small group, and/or the elements of the group are very similar.
    #Solution: add one-column(s) to the factors (-> temp$loadings), and zero-column(s) to factor loadings (-> temp$scores).
    #This ensures that the product of factor and factor loadings does not get altered.
    #  (Note: adding zero-column to temp$scores would have the result that the rank of factor_group[...] is too low, therefore use one-column(s).)
    if(verbose_robustpca) {
      message("----dimension of output of macropca (factorloadings and factors):----")
      print(dim(temp$scores))
      print(dim(temp$loadings))
      print(paste("number of eigenvectors:", number_eigenvectors))
    }
    if(!is.null(dim(temp$scores))) {
      if(dim(temp$loadings)[2] != number_eigenvectors) {
        message(paste("There is an issue with the dimensions of the factors: MacroPCA returns only",dim(temp$loadings)[2],
        "eigenvectors, instead of", number_eigenvectors,  "-> add one-column(s)"))
        for(i in 1:(number_eigenvectors - dim(temp$loadings)[2])) {
          temp$scores = cbind(temp$scores, 0) #add 0 to factor loadings
          temp$loadings = cbind(temp$loadings, 1) #add 1 to factors
          if(verbose_robustpca) {
            message("resulting factors:")
            print(temp$loadings[1:5,])
          }

        }
      }
    }
    ##############################
    scores = temp$scores[,1:number_eigenvectors] #these are the factor loadings
    factors_macropca = temp$loadings[,1:number_eigenvectors] #these are the factors
    suppressWarnings(
      if(class(factors_macropca) == "numeric") { #case of estimating 1 factor -> numeric -> make matrix
        factors_macropca = matrix(factors_macropca)
      }
    )
    if(!error_macropca) {
        #rare issue with macropca
        if(nrow(factors_macropca) != ncol(object)) {
          print(dim(object))
          print(dim(factors_macropca))
          message("--MacroPCA has dropped a column for unknown reasons---") #This leads to wrong dimensions in the factors, and gives error in rstudio.
        }
    }
    if(error_macropca) {
      factors_macropca = handle_macropca_errors(object, temp, KMAX,number_eigenvectors)

    }
    if(verbose_robustpca) print("end of macropca")
    return(list(factors_macropca,scores))


}




#'Estimates common factor F
#'
#'The estimator for F, see Anderson (1984), is equal to the first r eigenvectors (multiplied by sqrt(T) due to the restriction F'F/T = I)
#'associated with first r largest eigenvalues of the matrix WW' (which is TxT)
#'W is called Z in Ando and Bai (2017)
#' @inheritParams calculate_W
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param initialise boolean
#' @param method_estimate_factors defines method of robust estimaton of the factors: "macro", "pertmm" or "cz"
#' @inheritParams estimate_beta
#' @inheritParams calculate_virtual_factor_and_lambda_group
#' @inheritParams calculate_Z_common
#' @inheritParams update_g
#' @return r x T matrix
#' @importFrom stringr str_c
#' @export
estimate_factor <- function(use_robust, beta_est, g, lgfg_list,
                            method_estimate_factors,
                            initialise = FALSE,
                            NN = number_of_time_series,
                            TT = length_of_time_series,
                            number_of_common_factors = number_of_common_factors_fixedvalue,
                            number_of_variables = number_of_variables_fixedvalue,
                            number_vars_estimated = number_vars_estimated_fixedvalue,
                            verbose = FALSE) {


  if(!do_we_estimate_common_factors(number_of_common_factors) ) {
    #return matrix with only zero's
    return(list(t(matrix(rep(0, TT))), NA))
  }
  #initialisation: has no grouped factorstructure yet
  if(initialise) {
    W = calculate_W(beta_est, g,
                    NN = NN,
                    TT = TT,
                    number_of_variables = number_of_variables,
                    number_vars_estimated = number_vars_estimated) #Y - XT
  } else {
    W = calculate_Z_common(beta_est, g, lgfg_list, method_estimate_factors,
                           NN = NN,
                           TT = TT,
                           number_of_variables = number_of_variables,
                           number_vars_estimated = number_vars_estimated) #this is a "rows_without_NA x T - matrix"
  }
  if(verbose) message("W is constructed")
  #if there are individuals in "class zero", then they are considered outliers
  #  and need to be excluded in the estimation of the common factors as well
  if(verbose) message("removing class zero individuals from estimation of common factors")
  W = W[g != 0, ]



  #Define the object on which (robust or classical) PCA will be performed
  if(use_robust) {
    #CHANGE LOCATION 4/7
    if(method_estimate_factors %in% c("macro","pertmm")) {
      temp = prepare_for_robpca(W)
    } else {
      temp = t(W)%*%W / (NN * TT)
    }
  } else {
    #Classical case
    #calculate 1/NT * W'W
    temp = t(W)%*%W / (NN * TT) #TxT matrix (this division by NT does not actually matter for the calculation of the eigenvectors)
  }


  #first r eigenvectors * sqrt(T)

  #take number_of_common_factors eigenvectors
  message(str_c("Estimate ",number_of_common_factors," common factors"))
  if(use_robust) {
    #CHANGE LOCATION 1/7
    if(method_estimate_factors %in% c("macro","pertmm")) {
      temp2 = robustpca(temp, number_of_common_factors, verbose_robustpca = verbose)
      estimatorF = t(sqrt(TT) * temp2[[1]])
      scores = temp2[[2]] #=factor loadings coming out of macropca
      rm(temp2)
    } else {
      estimatorF = t(sqrt(TT) * eigen(temp)$vectors[,1:number_of_common_factors])
      scores = NA #because loadings will be calculated later
    }
  } else {
    estimatorF = t(sqrt(TT) * eigen(temp)$vectors[,1:number_of_common_factors])
    scores = NA #because loadings will be calculated later
  }
  #OLD:
  # if(!do_we_estimate_common_factors(number_of_common_factors) ) {
  #   #then estimatorF has length TT, which is still an apropriate size to use in this case -> just set values to zero
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
prepare_for_robpca <- function(object, NN = number_of_time_series, TT = length_of_time_series, option = 3) {
  if(option == 1) {
    message("not used anymore")
    #temp = get_robust_covmatrix(object)
  }

  if(option == 2) temp = temp = t(object)%*%object / (NN * TT)

  if(option == 3) temp = object

  return(temp)
}



#' Estimates group factors
#'
#' @inheritParams calculate_W
#' @param lambda common factor loadings
#' @param comfactor common factors
#' @param initialise boolean
#' @inheritParams estimate_beta
#' @inheritParams calculate_virtual_factor_and_lambda_group
#' @inheritParams estimate_factor
#' @inheritParams update_g
#' @export
estimate_factor_group <- function(use_robust, beta_est, g, lambda, comfactor,
                                  method_estimate_factors,
                                  initialise = FALSE,
                                  NN = number_of_time_series,
                                  TT = length_of_time_series,
                                  number_of_groups = number_of_groups_fixedvalue,
                                  number_of_group_factors = number_of_group_factors_fixedvalue,
                                  number_of_common_factors = number_of_common_factors_fixedvalue,
                                  number_of_variables = number_of_variables_fixedvalue,
                                  number_vars_estimated = number_vars_estimated_fixedvalue,
                                  verbose = FALSE
                                  #returnscores = FALSE
                                  ) {

  estimatorF = list()
  scores = list()

  for(group in 1:number_of_groups) {
    if(verbose) print(paste("Estimate factors for group", group))
    if(number_of_group_factors[group] > 0) {
      if(TT < number_of_group_factors[group]) {
        warning("There are too many factors to be estimated, compared to TT.")
      }

      Wj = calculate_Z_group(beta_est, g, lambda, comfactor, group, initialise,
                             TT = TT,
                             number_of_variables = number_of_variables,
                             number_vars_estimated = number_vars_estimated,
                             number_of_common_factors = number_of_common_factors)


      #use limit on nrow(Wj) due to error ("The input data must have at least 3 rows (cases)")
      #When the number of rows within the group is bigger than 3, macropca runs, but in some cases can drop columns (when values within this column are equal).
      #So, in practice it is better to increase the threshold a bit.
      ROBUST_THRESHOLD = 5
      if(use_robust & nrow(Wj) > ROBUST_THRESHOLD) {


        #CHANGE LOCATION 2/7
        if(method_estimate_factors %in% c("macro","pertmm")) {
          temp = prepare_for_robpca(Wj)
          temp2 = robustpca(temp, number_of_group_factors[group], verbose_robustpca = verbose)
          estimatorF[[group]] = t(sqrt(TT) * temp2[[1]]) #robust pca
          scores[[group]] = temp2[[2]]
          rm(temp2)
        } else {
          ##CZ
          temp = t(Wj)%*%Wj #dividing by NT does not make any difference for the eigenvectors
          scores[[group]] = NA #because loadings will be calculated later
          estimatorF[[group]] = t(sqrt(TT) * eigen(temp)$vectors[,1:number_of_group_factors[group]])
        }

      } else {
        if(nrow(Wj) < ROBUST_THRESHOLD & use_robust) message("This group contains a very small number of elements. -> macropca cannot work -> use non-robust estimation of factors")

        temp = t(Wj)%*%Wj #dividing by NT does not make any difference for the eigenvectors
        scores[[group]] = NA #because loadings will be calculated later
        estimatorF[[group]] = t(sqrt(TT) * eigen(temp)$vectors[,1:number_of_group_factors[group]])
      }



    } else {
      estimatorF[[group]] = matrix(0, 1, TT)
      scores[[group]] = NA #because loadings will be calculated later
    }
  }



  #as test:
  # if(returnscores) {
  #   return(scores)
  # }
  return(estimatorF)
}

#' calculates factor loadings of common factors
#' @inheritParams calculate_W
#' @param comfactor common factors
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param initialise boolean
#' @inheritParams estimate_beta
#' @inheritParams calculate_Z_common
#' @inheritParams estimate_factor
#' @export
calculate_lambda <- function(use_robust, beta_est, comfactor, g, lgfg_list,
                             method_estimate_factors,
                             initialise = FALSE,
                             NN = number_of_time_series, TT = length_of_time_series,
                             number_of_common_factors = number_of_common_factors_fixedvalue,
                             number_of_variables = number_of_variables_fixedvalue,
                             number_vars_estimated = number_vars_estimated_fixedvalue) {

  if(!do_we_estimate_common_factors(number_of_common_factors)) {
    return(t(matrix(rep(0, 300))))
  }
  if(initialise) {
    W = calculate_W(beta_est, g,
                    NN = NN,
                    TT = TT,
                    number_of_variables = number_of_variables,
                    number_vars_estimated = number_vars_estimated)
  }  else {
    W = calculate_Z_common(beta_est, g, lgfg_list, method_estimate_factors,
                           NN = NN,
                           TT = TT,
                           number_of_variables = number_of_variables,
                           number_vars_estimated = number_vars_estimated)
  }

  if(use_robust) {
    #CHANGE LOCATION 6/7
    if(method_estimate_factors %in% c("macro","pertmm")) {
      lambda = return_robust_lambdaobject(W, NA, type = 2, g, comfactor_rrn = comfactor, number_of_common_factors_rrn = nrow(comfactor),
                                          NN_rrn = NN #,
                                          #use_real_world_data_rrn = use_real_world_data,
                                          #application_covid = exists("usecoviddata_cases") | exists("usecoviddata_deaths")
                                          )
    } else {
      lambda = t(W %*% t(comfactor) / TT)
    }
  } else {
    lambda = t(W %*% t(comfactor) / TT)
  }

  #OLD:
  # if(!do_we_estimate_common_factors(number_of_common_factors) & !initialise) {
  #   lambda = lambda - lambda
  # }



  ############################
  # for time-series with NA's,
  # the lambda's could not be determined (dimension of lambda is smaller) -> set NA's to zero
  ############################
  if(anyNA(Y)) {
    rows_with_NA = which(apply(Y,1,anyNA))
    rows_without_NA = which(!apply(Y,1,anyNA))
    temp = data.frame(matrix(NA, nrow = max(number_of_common_factors, 1), ncol = NN)) #new dataframe
    temp[,rows_without_NA] = lambda #put calculated lambda's in df
    if(sum(rows_with_NA,na.rm=T) > 0)  temp[,rows_with_NA] = 0 #set the rest to zero
    lambda = as.matrix(temp)
  }


  return(lambda)
}

#'calculates factor loadings of groupfactors
#'
#'returns object which includes group and id of the individuals
#' @param beta_est estimated values of beta
#' @param factor_group group factors
#' @param lambda loadings
#' @param comfactor common factors
#' @param initialise boolean
#' @param UPDATE1 option to indicate the number of groupfactors is updated during the algorithm; defaults to FALSE
#' @param UPDATE2 option to indicate the number of common factors is updated during the algorithm; defaults to FALSE
#' @inheritParams estimate_beta
#' @inheritParams calculate_W
#' @inheritParams estimate_factor
#' @importFrom rlang .data
#' @export
calculate_lambda_group <- function(use_robust, beta_est, factor_group, g, lambda, comfactor,
                                   method_estimate_factors,
                                   initialise = FALSE,
                                   UPDATE1 = FALSE, UPDATE2 = FALSE,
                                   NN = number_of_time_series,
                                   TT = length_of_time_series,
                                   number_of_groups = number_of_groups_fixedvalue,
                                   number_of_group_factors = number_of_group_factors_fixedvalue,
                                   number_of_common_factors = number_of_common_factors_fixedvalue,
                                   number_of_variables = number_of_variables_fixedvalue,
                                   number_vars_estimated = number_vars_estimated_fixedvalue) {


  lambda_local = list()

  for(group in 1:number_of_groups) {
    if(number_of_group_factors[group] > 0) {
      Wj = calculate_Z_group(beta_est, g, lambda, comfactor, group, initialise,
                             TT = TT,
                             number_of_variables = number_of_variables,
                             number_vars_estimated = number_vars_estimated,
                             number_of_common_factors = number_of_common_factors)

      #robust things:
      if(use_robust) {

        #CHANGE LOCATION 7/7
        if(method_estimate_factors %in% c("macro","pertmm")) {
          #In the classical approach each lambda is the mean of a set of products of F and Y (Z),
          #   (lambda_N1 = (F_T1 * Y_N1T1 + F_T2 * Y_N1T2 + ...) / TT)
          #   We replace this mean by an M-estimator in the robust approach.
          if(UPDATE1) {
            lambda_local[[group]] = return_robust_lambdaobject(Wj, group, type = 3, g, factor_group_rrn = factor_group,
                                                               number_of_group_factors_rrn = unlist(lapply(factor_group, nrow))#,
                                                               #use_real_world_data_rrn = use_real_world_data,
                                                               #application_covid = exists("usecoviddata_cases") | exists("usecoviddata_deaths")
                                                               )
          } else {
            lambda_local[[group]] = return_robust_lambdaobject(Wj, group, type = 3, g,
                                                               number_of_group_factors_rrn = number_of_group_factors#,
                                                               #use_real_world_data_rrn = use_real_world_data,
                                                               #application_covid = exists("usecoviddata_cases") | exists("usecoviddata_deaths")
                                                               )
          }
        } else {
          FGG = factor_group[[group]]
          if(dim(t(FGG))[1] == 1) {
            FGG = matrix(FGG, nrow = 1)
          }
          lambda_local[[group]] = Wj %*% t(FGG) / TT
        }


      } else {

        FGG = factor_group[[group]]
        if(dim(t(FGG))[1] == 1) {
          FGG = matrix(FGG, nrow = 1)
        }
        lambda_local[[group]] = Wj %*% t(FGG) / TT
      }

    } else { #estimate zero group factors for this group -> matrix with 1 column
        # lambda_local[[group]] = matrix(0,
        #                                length(rows_without_NA[(rows_without_NA %in% which(g == group))]),
        #                                1)
        lambda_local[[group]] = matrix(0, table(g)[group], 1)
    }

  }


  #for income-series with NA's, the lambda's cannot be determined -> set to zero with this function
  add_lambdas_from_NArows <- function(DF, groep) {
    rows_with_NA = which(apply(Y,1,anyNA))
    rows_without_NA = which(!apply(Y,1,anyNA))
    ids = rows_with_NA[(rows_with_NA %in% which(g == groep))]
    if(length(ids) > 0) {
      temp_init = rep(0, length(ids)) #zero's
      temp = temp_init
      if(number_of_group_factors[groep] > 1) {
        for(i in 2:number_of_group_factors[groep]) {
          temp = temp %>% cbind(temp_init)
        }
      }
      temp = data.frame(temp %>% cbind(groep) %>% cbind(ids))
      names(temp) = names(DF)
      return(temp)
    } else {
      return(NULL)
    }
  }

  #change list to dataframe:
  lambda_local2 = data.frame(lambda_local[[1]])
  names(lambda_local2) = str_c("X",1:ncol(lambda_local2))  #str_c("X",1:max(1,number_of_group_factors[1]))
  rows_without_NA = which(!apply(Y,1,anyNA))
  lambda_local2 = lambda_local2 %>% mutate(groep = 1, id = rows_without_NA[(rows_without_NA %in% which(g == 1))])


  #fill up with "rows_with_NA" of this groep:
  #lambda_group can not be defined for these rows, so put to 0
  lambda_local2 = lambda_local2 %>% rbind(add_lambdas_from_NArows(lambda_local2,1)) %>% arrange(.data$id)

  if(number_of_groups > 1) {
    for(i in 2:number_of_groups) {
      df_to_add = data.frame(lambda_local[[i]])
      names(df_to_add) = str_c("X",1:ncol(df_to_add))

      df_to_add = df_to_add %>% mutate(groep = i, id = rows_without_NA[(rows_without_NA %in% which(g == i))])

      df_to_add = df_to_add %>% rbind(add_lambdas_from_NArows(df_to_add,i)) %>% arrange(.data$id)

      lambda_local2 = lambda_local2 %>%
        bind_rows(df_to_add)
    }
  }



  #get rid of NA's when using update1 or update2
  if(UPDATE1 | UPDATE2) {
    if(anyNA(lambda_local2)) {
      lambda_local2[is.na(lambda_local2)] <- 0
    }
  }

  #add rows to lambda_local2 for individuals of class zero (=individuals that are too far from all groups)
  #put those values to zero
  if(sum(g==0) > 0) {
    indices_class_zero = which(g==0)
    for(ind in indices_class_zero) {
      extrarow = c(rep(0, max(number_of_group_factors, na.rm = TRUE)),
                 0, #this is group membership
                 ind #this is id
                 )
      lambda_local2 = rbind(lambda_local2, extrarow)
    }
  }
  stopifnot(nrow(lambda_local2) == NN)

  return(lambda_local2)
}


#' Function which is used to have a dataframe (called "grid") with data (individualindex, timeindex, XT and LF) available.
#'
#' (The dataframe "grid" used to be required in update_g(), via OF_vectorized3(), but not anymore.)
#' @param grid dataframe containing values for X*beta_est and LF (product of common factor and its loadings)
#' @param beta_est estimated values of beta
#' @param lambda loadings of the common factors
#' @param comfactor common factors
#' @param limit_est_groups_heterogroups maximum amount of groups that can be estimated when method_estimate_beta is set to "group"
#' @inheritParams estimate_beta
#' @export
grid_add_variables <- function(grid, beta_est, lambda, comfactor, method_estimate_beta,
                               NN = number_of_time_series,
                               TT = length_of_time_series,
                               number_of_variables = number_of_variables_fixedvalue,
                               number_vars_estimated = number_vars_estimated_fixedvalue,
                               number_of_groups = number_of_groups_fixedvalue,
                               limit_est_groups_heterogroups = 15) {
  if(number_of_variables > 0) {
    #for homogeneous beta_est (1 -> 4 at this moment), we only need 1 column as all columns are the same
    if(method_estimate_beta == "homogeneous") {
      beta_used = as.matrix(beta_est[,1])
      #calculate matrix multiplications outside the OF-functions to speed up:
      grid$XBETA = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_used))
    } else {
      if(method_estimate_beta == "group") {
        stopifnot((number_of_groups >= 0 & number_of_groups < limit_est_groups_heterogroups)) #code exists up to 14 groups

        #for each group: define XT
        if(number_of_groups > 0) grid$XBETA1 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,1]))
        if(number_of_groups > 1) grid$XBETA2 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,2]))
        if(number_of_groups > 2) grid$XBETA3 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,3]))
        if(number_of_groups > 3) grid$XBETA4 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,4]))
        if(number_of_groups > 4) grid$XBETA5 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,5]))
        if(number_of_groups > 5) grid$XBETA6 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,6]))
        if(number_of_groups > 6) grid$XBETA7 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,7]))
        if(number_of_groups > 7) grid$XBETA8 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,8]))
        if(number_of_groups > 8) grid$XBETA9 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,9]))
        if(number_of_groups > 9) grid$XBETA10 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,10]))
        if(number_of_groups > 10) grid$XBETA11 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,11]))
        if(number_of_groups > 11) grid$XBETA12 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,12]))
        if(number_of_groups > 12) grid$XBETA13 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,13]))
        if(number_of_groups > 13) grid$XBETA14 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% beta_est[,14]))

      }
      if(method_estimate_beta == "individual") {
        grid$XBETA = c(calculate_XB_estimated(NN = NN, TT = TT, number_of_variables = number_of_variables, number_vars_estimated = number_vars_estimated))
      }
    }
  } else {
    if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "individual")) {
      grid$XBETA = 0
    }
    if(method_estimate_beta == "group") {
      if(number_of_groups > 0) grid$XBETA1 = 0
      if(number_of_groups > 1) grid$XBETA2 = 0
      if(number_of_groups > 2) grid$XBETA3 = 0
      if(number_of_groups > 3) grid$XBETA4 = 0
      if(number_of_groups > 4) grid$XBETA5 = 0
      if(number_of_groups > 5) grid$XBETA6 = 0
      if(number_of_groups > 6) grid$XBETA7 = 0
      if(number_of_groups > 7) grid$XBETA8 = 0
      if(number_of_groups > 8) grid$XBETA9 = 0
      if(number_of_groups > 9) grid$XBETA10 = 0
      if(number_of_groups > 10) grid$XBETA11 = 0
      if(number_of_groups > 11) grid$XBETA12 = 0
      if(number_of_groups > 12) grid$XBETA13 = 0
      if(number_of_groups > 13) grid$XBETA14 = 0


    }

  }
  LF = t(lambda) %*% comfactor
  grid$LF = c(LF)

  return(grid)
}













#' Function with as input a dataframe. (this will be "Y" or "to_divide") It filters out rows with NA.
#'
#' @param df input
#' @importFrom stats na.omit
#' @export
handleNA <- function(df) {
  rows_with_NA = which(apply(df,1,function(x) sum(is.na(x)) != 0)) #rownumbers of rows containing NA
  #message(paste("There are ",length(rows_with_NA),"of",nrow(df), "rows with NA's."))

  return(list(na.omit(df), rows_with_NA))
}

#' Removes NA's in LG (in function calculate_virtual_factor_and_lambda_group() )
#' @param df input
#' @importFrom dplyr mutate_if
#' @importFrom tidyr replace_na
#' @export
handleNA_LG <- function(df) {
  result = (as.matrix(data.frame(df) %>% mutate_if(is.numeric , replace_na, replace = 0)))
  return(result)
}



#' Helpfunction to shorten code: use a and b
#'
#' @param number_of_common_factors number of common factors
do_we_estimate_common_factors <- function(number_of_common_factors) {
  if(number_of_common_factors == 0) {
    a = 0
  } else {
    a = 1
  }
  return(a)
}
#' Helpfunction to shorten code: use a and b
#' @param number_of_group_factors number of group factors to be estimated
do_we_estimate_group_factors <- function(number_of_group_factors) {
  if(mean(number_of_group_factors, na.rm = T) == 0) {
    b = 0
  } else {
    b = 1
  }
  return(b)
}


#' Calculates the error term Y - X*beta_est - LF - LgFg.
#' @return NxT matrix#'
#' @param g vector with group memberships
#' @param no_common_factorstructure if there is a common factorstructure being estimated
#' @param no_group_factorstructure if there is a group factorstructure being estimated
#' @inheritParams grid_add_variables
#' @inheritParams estimate_beta
#' @inheritParams initialise_beta
#' @importFrom stringr str_detect
#' @importFrom rlang .data
#' @export
calculate_error_term <- function(Y, X, beta_est, g, factor_group, lambda_group, comfactor, lambda,
                                 method_estimate_beta,
                                 NN = number_of_time_series,
                                 TT = length_of_time_series,
                                 number_of_variables = number_of_variables_fixedvalue,
                                 number_vars_estimated = number_vars_estimated_fixedvalue,
                                 number_of_groups = number_of_groups_fixedvalue,
                                 number_of_group_factors = number_of_group_factors_fixedvalue,
                                 number_of_common_factors = number_of_common_factors_fixedvalue,
                                 no_common_factorstructure = FALSE, no_group_factorstructure = FALSE) {
  stopifnot(table(lambda_group$groep) == table(g)) #in this case lambda_groep would need to be updated

  u = matrix(NA, nrow = NN, ncol = TT)
  e = matrix(NA, nrow = NN, ncol = TT)
  lf = t(lambda) %*% comfactor

  if(number_vars_estimated > 0) {
    if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) {
      xt = sapply(1:NN,
                  function(y) sapply(1:TT, function(x) c(1, X[y, x, ]) %*% beta_est[, g[y]]))
    }
    if(method_estimate_beta == "individual") {
      xt = t(calculate_XB_estimated(NN = NN, TT = TT, number_of_variables = number_of_variables, number_vars_estimated = number_vars_estimated)) #TxN matrix
    }
  } else {
    xt = matrix(0, nrow = TT, ncol = NN) #TxN
  }
  print("2")
  lf_group = list()
  group_membership = list()
  for(k in 1:number_of_groups) {
    #print(k)
    LGclean = lambda_group %>% arrange(.data$id) %>%
                          dplyr::filter(.data$groep == k)

    #replaced dplyr::select because of using ".data" in combination with "starts_with"
    #%>%  dplyr::select(starts_with("X")))

    selectcolumns = which(str_detect(colnames(LGclean), "X"))
    LGclean = as.matrix(LGclean[,selectcolumns])

    lf_group[[k]] = LGclean[,1:number_of_group_factors[k]] %*% factor_group[[k]]
    group_membership[[k]] = data.frame(as.matrix(lambda_group %>% arrange(.data$id) %>%
                                                   dplyr::filter(.data$groep == k) %>%
                                                   dplyr::select(.data$groep, .data$id)))
  }

  for(i in 1:NN) {
    if(g[i] != 0) {
      lf_group_i = lf_group[[g[i]]]
      index = which(group_membership[[g[i]]]$id == i)
    } else {
      lf_group_i = NA
      index = NA
    }
    for(t in 1:TT) {
      u[i,t] = Y[i,t] - xt[t,i]

      a = do_we_estimate_common_factors(number_of_common_factors)
      b = do_we_estimate_group_factors(number_of_group_factors)

      part2 = a * lf[i,t]
      part3 = ifelse(g[i] != 0,
                     b * lf_group_i[index, t],
                     0)
      if(no_common_factorstructure) part2 = 0
      if(no_group_factorstructure) part3 = 0
      e[i,t] = u[i,t] - part2 - part3
    }
  }
  return(e)
}


#' Calculates the error term Y - X*beta_est - LF - LgFg.
#'
#' Gives same value as calculate_error_term().
#' @inheritParams estimate_beta
#' @return NxT matrix
#' @importFrom rlang .data
#' @export
calculate_error_term_individuals <- function(NN = number_of_time_series,
                                             TT = length_of_time_series,
                                             number_of_variables = number_of_variables_fixedvalue,
                                             number_vars_estimated = number_vars_estimated_fixedvalue,
                                             number_of_groups = number_of_groups_fixedvalue,
                                             number_of_group_factors = number_of_group_factors_fixedvalue,
                                             number_of_common_factors = number_of_common_factors_fixedvalue) {
  u = matrix(NA,nrow = NN, ncol = TT)
  e = matrix(NA,nrow = NN, ncol = TT)
  lf = t(lambda) %*% comfactor

  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated)
  if(number_vars_estimated > 0) {
    if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) {
      xt = sapply(1:NN,
                  function(y) sapply(1:TT, function(x) c(1, X[y,x,]) %*% beta_est[,g[y]]))
    }
    if(method_estimate_beta == "individual") {
      xt = sapply(1:NN,
                  function(y) sapply(1:TT, function(x) c(1, X[y,x,]) %*% beta_est[,y]))
    }

  } else {
    xt = matrix(0, nrow = TT, ncol = NN) #TxN
  }

  lf_group = list()
  group_membership = list()
  for(k in 1:number_of_groups) {
    LGclean = lambda_group %>%
      arrange(.data$id) %>%
      dplyr::filter(.data$groep == k)

    #replaced dplyr::select because of using ".data" in combination with "starts_with"
    #%>% dplyr::select(starts_with("X")))

    selectcolumns = which(str_detect(colnames(LGclean), "X"))
    LGclean = as.matrix(LGclean[,selectcolumns])


    lf_group[[k]] = LGclean[,1:number_of_group_factors[k]] %*% factor_group[[k]]
    group_membership[[k]] = data.frame(as.matrix(lambda_group %>% arrange(.data$id) %>%
                                                   dplyr::filter(.data$groep == k) %>%
                                                   dplyr::select(.data$groep, .data$id)))

  }

  #define a and b to shorten the following code
  a = do_we_estimate_common_factors(number_of_common_factors)
  b = do_we_estimate_group_factors(number_of_group_factors)

  for(i in 1:NN) {
    if(g[i] != 0) {
      lf_group_i = lf_group[[g[i]]]
      index = which(group_membership[[g[i]]]$id == i)
    } else {
      lf_group_i = NA
      index = NA
    }
    u[i,] = (Y[i,] - xt[,i]) %>% as.numeric
    if(g[i] != 0) {
      e[i,] = u[i,] - a * lf[i,] - b * lf_group_i[index,]
    } else {
      e[i,] = u[i,] - a * lf[i,]
    }


  }
  return(e)
}


#' Calculates sum of squared errors, divided by NT
#'
#' @param e matrix with error terms
#' @param NN N
#' @param TT T
#' @export
calculate_sigma2 <- function(e, NN = number_of_time_series, TT = length_of_time_series) {
  if(anyNA(e)) {
    message("There are NA's in e (calculate_sigma2(). This should not happen. ")
  }
  sigma2 = sum(e*e)/(NN * TT)
  return(sigma2)
}



#' Function to calculate the norm of a matrix.
#'
#' @param A input matrix
#' @export
matrixnorm <- function(A) {
  return(sqrt(sum(diag(t(A) %*% A))))
}



#' Function to calculate the first term of PIC (panel information criterium)
#'
#' This is used in calculate_PIC()
#' @param e2 NxT matrix with the error terms
#' @importFrom robustbase Mpsi
#' @export
calculate_PIC_term1 <- function(e2) {
  if(use_robust) {
    #This replaces the classical sum(z^2) by sum(rho(z)) with rho the bisquare function.

    #NOTE: dividing by mad(e2) scales the errors of configurations with bad estimations to the results of good estimations
    #  -> this leads to wrong estimation of S & k
    RHO_E_SCALED = Mpsi((e2 - median(e2)) , cc = 4.685, psi = "bisquare", deriv = -1)
    term1_new = sum(RHO_E_SCALED) / (nrow(e2) * ncol(e2))

  } else {
    term1_new = sum(e2^2) / (nrow(e2) * ncol(e2))

  }

  return(term1_new)
}

#' Function to calculate term4 of PIC for a group j; includes options for testing alternative PIC's.
#'
#' Used in calculate_PIC(), in a loop over all groups.
#' @param temp term4 of PIC for group n
#' @param term4 term 4 of PIC, updated until group n-1; this function updates this value with group n
#' @param Nj number of individuals in group j
#' @param NN N
#' @param use_alterpic indicates using an alternative version of the PIC (the one of AB2016 instead of AB2017). It weighs the fourth term with an extra factor relative to the size of the groups.
calculate_PIC_term4 <- function(temp, term4, Nj, NN = number_of_time_series, use_alterpic = ALTERNATIVE_PIC) {


  if(exists("ALTERNATIVE_PIC_2")) { #PIC of AB2017-code (includes 3 other changes (see test_alternative_PIC()))
    term4 = term4 + (temp * Nj/NN)
  } else {

    # if(exists("ALTERNATIVE_PIC")) {
      if(use_alterpic == TRUE) {
        #message("We use an Alternative PIC -> Ando/Bai2016-paper")
        term4 = term4 + (temp * Nj/NN) #=weigh with relative size of groups
      } else {
        #=Standard PIC (from Ando/Bai2017-paper)
        term4 = term4 + temp
      }
    # } else {
    #   #=Standard PIC (from Ando/Bai2017-paper)
    #   term4 = term4 + temp
    # }
  }
  return(term4)
}



#' Function to determine PIC (panel information criterium)
#'
#' This depends on kappa1 -> kappaN, through p (=number of nonzero elements of beta_est).
#' The parameter 'sigma2' is the non-robust sigma2. As it is only used in term 2 to 4, it does not actually matter what its value is (needs to be > 0).
#' Could be set to 1 as well.
#' @param C determines relative contribution of the penamlty terms
#' @inheritParams estimate_beta
#' @param e2 matrix with error terms
#' @param sigma2 scalar: sum of squared error terms, scaled by NT
#' @export
calculate_PIC <- function(C, number_of_common_factors, number_of_group_factors, e2, sigma2,
                          NN = number_of_time_series,
                          TT = length_of_time_series,
                          number_of_variables = number_of_variables_fixedvalue,
                          number_vars_estimated = number_vars_estimated_fixedvalue,
                          number_of_groups = number_of_groups_fixedvalue) {


  term1 = calculate_PIC_term1(e2)

  if(number_vars_estimated > 0) {
    if((method_estimate_beta == "homogeneous") | (method_estimate_beta == "group")) {
      #p is number of nonzero elements of beta_est; we need the sum of all p's (including the intercept)
      p_sum = sum(sapply(1:NN, function(x) sum(beta_est[,g[x]] != 0)))
    } else{
      p_sum = sum(sapply(1:NN, function(x) sum(beta_est[,x] != 0)))
    }
  } else {
    p_sum = 0
  }

  term2 = C/NN * sigma2 * log(TT) * p_sum

  #term3: penalty on number of common factors
  term3 = C * number_of_common_factors * sigma2 * (TT + NN)/(TT*NN) * log(TT*NN)

  #term4: penalty on number of groupfactors
  term4 = 0
  for(j in 1:number_of_groups) {
    Nj = sum(g == j)
    TN_factor = (TT + Nj)/(TT*Nj) * log(TT*Nj)
    AFG = number_of_group_factors[j]
    temp = C * AFG * sigma2 * TN_factor

    term4 = calculate_PIC_term4(temp,term4,Nj)

  }


  newterms = test_alternative_PIC(C, term2,term3,term4) #returns input if there are no alternative PIC's used, which is the default
  term2 = newterms[[1]]
  term3 = newterms[[2]]
  term4 = newterms[[3]]


  return(term1 + term2 + term3 + term4)
}

#' This function tests alternative PIC's.
#'
#' @inheritParams calculate_PIC
#' @param term2 term2 of PIC
#' @param term3 term3 of PIC
#' @param term4 term4 of PIC
test_alternative_PIC <- function(C, term2, term3, term4) {
  #For testing purposes:
  if(exists("ALTERNATIVE_PIC_2")) {
    #message("use PIC of AB-code (4 things are different (part2: 1/C and 2xC)")
    #note that Nj/N was added earlier already
    term2 = term2 / length_of_time_series
    term3 = C * term3
    term4 = C * term4
  }
  return(list(term2,term3,term4))
}















#' Calculates the product of X*beta_true .
#' @inheritParams estimate_beta
#' @param g estimated group membership
#' @export
calculate_XB_true <- function(NN = number_of_time_series, TT = length_of_time_series, number_of_variables = number_of_variables_fixedvalue, g = g) {
  if(method_estimate_beta == "homogeneous") { #This is only relevant for DGP5 (BramatiCroux)
    XB_true = beta_true[1,] + X[,,1] * beta_true[1,]
    stopifnot(number_of_variables == 1)
  }
  if(method_estimate_beta == "group") {
    if(number_of_variables > 0 & exists("g")) {
      if(!is.na(g)) {
        XB_true = t(sapply(1:NN,
                              function(x) matrix(cbind(1, X[x,,]) %*% beta_true[,g[x]], nrow = 1)))
      } else {
        XB_true = NA
      }
    } else {
      XB_true = NA
    }
  }
  if(method_estimate_beta == "individual") {
    if(number_of_variables > 0 ) {
      XB_true = (t(sapply(1:NN,
                          function(y) sapply(1:TT, function(x) c(1, X[y,x,1:number_of_variables_fixedvalue]) %*% beta_true[,g_true][,y]))))
    } else {
      XB_true = NA
    }
  }

  return(XB_true)
}

#' When we run the algorithm with a different number of observable variables then the number we actually have, we need to reformat X.
#'
#' @inheritParams create_true_beta
#' @inheritParams generate_Y
#' @inheritParams initialise_beta
adapt_X_estimating_less_variables <- function(number_of_variables,
                                              number_vars_estimated
                                              ) {
  #if number_vars_estimated < number_of_variables, then the obsolete rows in beta_est are already erased -> do the same in X
  if(number_vars_estimated < number_of_variables) {
    if(number_vars_estimated > 0) {
      X = X[,,1:number_vars_estimated]
    } else {
      X = NA
    }
  }
  return(X)
}

#' Calculates (the estimated value of) the matrix X*beta_est.
#'
#' @inheritParams estimate_beta
#' @export
calculate_XB_estimated <- function(NN = number_of_time_series, TT = length_of_time_series,
                                   number_of_variables = number_of_variables_fixedvalue,
                                   number_vars_estimated = number_vars_estimated_fixedvalue) {
  #if number_vars_estimated < number_of_variables, then the obsolete rows in beta_est are already erased -> now do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated)

  if(number_of_variables > 0 & !is.na(X[1])) {
    if(method_estimate_beta == "homogeneous") { #currently only designed for DGP 5
      XT_geschat = beta_est[1] + X[,,1] * beta_est[1]
      stopifnot(number_of_variables == 1)
    }
    if(method_estimate_beta == "group") {
      if(number_vars_estimated > 0) {
        XT_geschat = t(sapply(1:NN,
                 function(x) matrix(cbind(1, X[x,,]) %*% beta_est[,g[x]], nrow = 1)))
      } else {
        XT_geschat = NA
      }
    }
    if(method_estimate_beta == "individual") {
      if(number_vars_estimated > 0) {
        XT_geschat = t(sapply(1:NN,
                 function(x) matrix(cbind(1, X[x,,]) %*% beta_est[,x], nrow = 1)))

      } else {
        XT_geschat = NA
      }
    }
  } else {
    XT_geschat = matrix(0, nrow = NN, ncol = TT)
  }
  return(XT_geschat)
}


#' Calculate true groupfactorstructure.
#' @return list with NjxT matrices
#' @param lgt true group factor loadings
#' @param fgt true group factors
#' @param true_number_of_groups true number of groups
#' @param true_number_of_group_factors true number of group factors for each group
#' @param true_number_of_common_factors true number of common factors
# @param using_bramaticroux parameter to indicate that we are using data generated with dgp 5
#' @inheritParams estimate_beta
#' @inheritParams generate_Y
#' @param dgp1_AB_local gives information about which DGP we use; TRUE of FALSE
#' @export
calculate_FL_group_true <- function(NN, TT,
                                    lgt, fgt,
                                    true_number_of_groups,
                                    true_number_of_group_factors,
                                    true_number_of_common_factors,
                                    num_factors_may_vary,
                                    dgp1_AB_local = FALSE) {



  if(!dgp1_AB_local) { #not using for DGP 1, because true values are not available given current coding
    temp = calculate_lgfg(lgt, fgt,
                          true_number_of_groups,
                          true_number_of_group_factors,
                          true_number_of_common_factors,
                          num_factors_may_vary)
    FL_GROUP_true = t(lapply(1:NN, function(x) temp[[g_true[x]]][x,]) %>% unlist %>% matrix(nrow = TT))
  } else {
    FL_GROUP_true = NA
  }
  return(FL_GROUP_true)
}


#' Returns the estimated groupfactorstructure.
#' @param fg estimated group factors
#' @param lg loadings of estimatd group factors
#' @inheritParams calculate_FL_group_true
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @return list with NjxT matrices
#' @export
calculate_FL_group_estimated <- function(lg = lambda_group, fg = factor_group,
                                         NN = number_of_time_series,
                                         TT = length_of_time_series,
                                         number_of_groups = number_of_groups_fixedvalue,
                                         number_of_group_factors = number_of_group_factors_fixedvalue,
                                         number_of_common_factors = number_of_common_factors_fixedvalue,
                                         num_factors_may_vary = num_factors_may_vary_fixedvalue) {
  temp = calculate_lgfg(lg, fg, number_of_groups, number_of_group_factors, number_of_common_factors, num_factors_may_vary)
  #helpfunction to return the groupfactorstructure of x. For elements in class zero, for which no factors are estimated, this is set to 0.
  FL_helpf <- function(x) {
    if(g[x] != 0) {
      return(temp[[g[x]]][x,])
    } else {
      return(rep(0, TT))
    }
  }
  suppressWarnings( #the condition has length > 1 and only the first element will be used
    if(!is.na(temp)) {
      FL_group_est = lapply(1:NN, function(x) FL_helpf(x))
      FL_group_est = t(matrix(unlist(FL_group_est), nrow = TT))
    } else {
      FL_group_est = NA
    }
  )
  return(FL_group_est)
}

#' Function to calculate the mean squared error of beta_est.
#'
#' It is calculated with the formula on p19 of the supplement of \insertCite{Ando2017;textual}{RCTS}.
#' For DGP 1 & 2: When the true number of variables in X is not equal to the standard of 3 it currently returns NA.
#' @param beta_est estimated values of beta
#' @param beta_true true values of beta
#' @param without_intercept TRUE of FALSE: to remove the intercept in the calculation of the MSE
#' @inheritParams estimate_beta
#' @inheritParams calculate_FL_group_true
#' @export
calculate_mse_beta <- function(beta_est, beta_true, without_intercept = FALSE,
                                NN = number_of_time_series,
                                number_of_variables = number_of_variables_fixedvalue,
                                special_case_dgp1 = FALSE) {

  if(method_estimate_beta == "homogeneous") { #relevant for DGP05 (Bramati-Croux)
    mse = mean( (beta_est - beta_true)^2 )
    if(without_intercept) mse = mean( (beta_est[-1,] - beta_true[-1,])^2 )
    return(mse)
  }

  if(method_estimate_beta == "group") {
    pre_mse = mse_heterogeneous_groups(without_intercept, special_case_dgp1)
  }

  if(method_estimate_beta == "individual") { #Default case; In case of DGP 1 & 2 it returns NA when number of variables is not equal to 3.
    if(without_intercept) {
      if(special_case_dgp1) { #calculate MSE without the normal, and without the de facto intercept
        beta_est = matrix(beta_est[c(-1, -2),], ncol = NN)
        beta_true = beta_true[c(-1, -2),]
      } else { #calculate MSE without the normal intercept
        beta_est = beta_est[-1,]
        beta_true = beta_true[-1,]
      }
    }


    pre_mse = rep(NA, NN)
    for(i in 1:NN) {
      #Case of DGP01: These dgp's (can) contain added noise: v1, v2, ...
      ##########################################################################################################################
      #NOTE: replaced this part since allowing for extra noise (v1, v2, ...) complicates code and secondly has never been used.
      # if(dgp1_AB_local) {
      #   if(number_of_variables == 3) { #if not TRUE, bvb if > 3, then v4,v5,... should be added
      #     if(without_intercept) {
      #       if(dgp1_spread_group_centers) {
      #         #stopifnot(length(beta_est[,i]) == 2)
      #         afw = beta_est[,i] - (beta_true[,g_true[i]] + c(v2[i], v3[i])) #v1 is not necessary (it belongs to the de facto intercept)
      #       } else {
      #         afw = beta_est[,i] - (beta_true[,g_true[i]] + c(v1[i], v2[i], v3[i])) #beta_true[,...] has number_of_variables elements
      #       }
      #     } else {
      #       afw = beta_est[,i] - (beta_true[,g_true[i]] + c(0, v1[i], v2[i], v3[i])) #beta_true[,...] has (number_of_variables + 1) elements
      #     }
      #   } else { #when number_of_variables != 3
      #
      #     message("number_of_variables != 3 -> return NA")
      #     return(NA)
      #
      #   }
      # } else { #Case of DGP2
      #   afw = beta_est[,i] - beta_true[,g_true[i]]
      # }
      ##########################################################################################################################
      #replacement is shorter:
      #(note: length of this vector depends on without_intercept and special_case_dgp1)
      afw = beta_est[, i] - beta_true[, g_true[i]]
      ##########################################################################################################################

      pre_mse[i] = mean(afw^2) #this is (beta_est - beta_true)^2 (mean() goes over the number of variables)
    }
  }


  return(mean(pre_mse)) #this is E[(beta_est - beta_true)^2]
}

#' Helpfunction in calculate_mse_beta, when method_estimate_beta == "group".
#' (beta is estimated for each group separately).
#'
#' @param without_intercept boolean to remove the intercept in the calculation
#' @inheritParams calculate_mse_beta
#' @inheritParams initialise_beta
mse_heterogeneous_groups <- function(without_intercept, special_case_dgp1) {

  #1. The order of the estimated groups is not necessarily the same as the order of the true groups. -> possible necessary to permutate g to get the MSE?
  #-> this is countered by the order in the object beta_est
  #-> we do not need permutation here

  # if(sd(v1) != 0) { #when extra noise is added to the DGP (this is a non-standard case)
  #   message("mse_heterogeneous_groups(): not implemented when v1 contains values (=extra noise in dgp)")
  #   return(NA)
  # }

  if(without_intercept) {
    if(special_case_dgp1) { #->calculate MSE without true and without de facto intercept of DGP1
      beta_est = beta_est[c(-1, -2),]
      beta_true = beta_true[c(-1, -2),]
    } else { #calculate MSE without true intercept
      beta_est = beta_est[-1,]
      beta_true = beta_true[-1,]
    }
  }

  pre_mse = rep(NA,number_of_time_series)
  for(i in 1:number_of_time_series) {

    afw = beta_est[,g[i]] - (beta_true[,g_true[i]])
    pre_mse[i] = mean(afw^2)
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
#' @inheritParams calculate_lambda_group
#' @importFrom stringr str_split
#' @importFrom purrr map
#' @export
calculate_VCsquared <- function(rcj, rc, C_candidates, UPDATE1 = FALSE, UPDATE2 = FALSE,
                                number_of_group_factors = number_of_group_factors_fixedvalue) {
  VC_squared = rep(NA, length(C_candidates))

  #vector with max number of groups
  Smax_local = rep(Smax, length(C_candidates))

  number_subsets = length(indices_subset) #number of subsamples used

  for(C_local in C_candidates) {
    C_index = which(C_candidates == C_local)
    part1 = 0
    for(i in 1:number_subsets) {
      part1 = part1 + (rc[C_index, i] - sum(rc[C_index,]) / number_subsets)^2
    }
    part1 = part1 / number_subsets

    part2 = 0
    for(j in 1:Smax_local[C_index]) {
      #print(paste("---",j))
      if(j <= length(number_of_group_factors) | UPDATE1 | UPDATE2) {
        part2_part = 0
        for(aA in 1:number_subsets) {
          temp = 0
          for(bA in 1:number_subsets) {
            if(j <= length(unlist(rcj[C_index, bA] %>% str_split("_")))) {
              temp = sum(temp, as.numeric(map(str_split(rcj[C_index, bA],"_"), j)), na.rm=T)
            }
          }
          if(j <= length(unlist(rcj[C_index, aA] %>% str_split("_")))) {
            rcjNT = as.numeric(map(str_split(rcj[C_index, aA], "_"), j))
          } else {
            rcjNT = 0
          }
          if(is.na(rcjNT)) rcjNT = 0
          part2_part = part2_part + (rcjNT - temp/number_subsets)^2
        }
        part2_part = part2_part / number_subsets
        part2 = part2 + part2_part
      } else {
        print("-*-")
      }
    }
    VC_squared[C_index] = part1 + part2
  }
  return(unlist(VC_squared))
}









#' Wrapper around lmrob.
#'
#' Desgined to make sure the following error does not happen anymore:
#' Error in if (init$scale == 0)  : missing value where TRUE/FALSE needed.
#' KS2014 is the recommended setting. It however does lead to a strongly increased computational time.
#' @param parameter_y dependent variable in regression
#' @param parameter_x independent variables in regression
#' @param nointercept if TRUE it performs regression without an intercept
#' @param nosetting option to remove the recommended setting in lmrob(). It is much faster. Defaults to FALSE.
#' @importFrom robustbase lmrob
#' @export
LMROB <- function(parameter_y, parameter_x, nointercept = FALSE, nosetting = FALSE) {


  if(is.na(parameter_x)[1]) { #when there are no independent variables
    result2 = tryCatch(
      {
        if(nosetting) {
          result = lmrob(parameter_y ~ 1)
        } else {
          result = lmrob(parameter_y ~ 1, setting = "KS2014")
        }
      },
      error=function(e) {
        message(e)
        print("error, therefore use lmrob without 'setting'")
        result = lmrob(parameter_y ~ 1)
        return(result)
      },
      finally={}
    )


  } else {
    if(nosetting) {
      if(nointercept) {
        result2 = lmrob(parameter_y ~ parameter_x + 0)
      } else {
        result2 = lmrob(parameter_y ~ parameter_x)
      }
    } else {
      if(nointercept) {

          result2 = lmrob(parameter_y ~ parameter_x + 0, setting = "KS2014")

      } else {
        #sometimes Error in if (init$scale == 0) happens. In that case: run without setting = KS2014.

          result2 = tryCatch(
            {
              lmrob(parameter_y ~ parameter_x, setting = "KS2014")
            },
            error=function(e) {
              message(e)
              print("error, therefore use lmrob without 'setting'")
              result2 = lmrob(parameter_y ~ parameter_x)
              return(result2)
            },
            finally={}
          )

      }
    }
  }
  return(result2)
}






#' Calculates sigma2maxmodel
#'
#' Sigma2 is the sum of the squared errors, divided by NT. We need the sigma2 of the maxmodel to use (in term 2,3,4 of the PIC) instead of the configuration-dependent sigma2. (See text in paper AndoBai 2016).
#' sigma2_max_model could actually be set to 1 as well, as it can be absorbed in parameter C of the PIC.
#' @inheritParams determine_beta
#' @inheritParams calculate_lambda_group
#' @param number_of_groups_candidates vector with candidate values for the number of groups
#' @param number_of_common_factors number of common factors
#' @export
calculate_sigma2maxmodel <- function(optimize_kappa = FALSE, UPDATE1 = FALSE, UPDATE2 = FALSE,
                                     number_of_groups = number_of_groups_fixedvalue,
                                     number_of_groups_candidates = number_of_groups_candidates_fixedvalue,
                                     number_of_group_factors = number_of_group_factors_fixedvalue,
                                     number_of_common_factors = number_of_common_factors_fixedvalue
                                     ) {
  if(!optimize_kappa) {
    if(!UPDATE1 & !UPDATE2 &
       mean(number_of_group_factors, na.rm=T) == k_max &
       number_of_groups == max(number_of_groups_candidates) &
       number_of_common_factors == max(k_common_candidates) ) {
      sigma2_max_model = calculate_sigma2(calculate_error_term())
    }
    if(UPDATE1 & UPDATE2) {
      #update number of group- and common factors during algorithm
      sigma2_max_model = NA
    }
  } else {
    message("Not defined.")
  }
  return(sigma2_max_model)
}

#' Adapts allpic (object that contains PIC for all candidate C's and all subsamples) with sigma2_max_model.
#'
#' The PIC is calculated with a sigma2 specific to the configuration (= number of groups and factors).
#' Because the method requires sigma2 to be equal over all configurations (see proofs of different papers of Ando/Bai) we replace sigma2 by the sigma2 of the configuration
#' with maximum number of groups and factors (this is the last one that is run).
#' @param all_pic contains PIC for all candidate C's and all subsamples
#' @param sigma2_max_model sigma2 of model with maximum number of groups and factors
#' @inheritParams calculate_lambda_group
#' @export
adapt_allpic_with_sigma2maxmodel <- function(all_pic, sigma2_max_model, UPDATE1 = FALSE, UPDATE2 = FALSE) {
  if(is.null(all_pic))  message("Warning in adapt_allpic_with_sigma2maxmodel(): all_pic is empty")
  if(!UPDATE1 & !UPDATE2) {
    #message("Use sigma2_max_model in PIC (in object all_pic).")
    for(i in 1:nrow(all_pic)) {
      sig2 = as.numeric(df_results$sigma2[i])
      all_pic[i,] = (unlist(all_pic[i,]) - sig2) / sig2 * sigma2_max_model + sig2
    }
    #message("...done")
  }
  return(all_pic)
}

#' This function creates an instance of DGP 2, as defined in \insertCite{BoudtHeyndels2021;textual}{RCTS}.
#'
#' The output is a dataframe with N (amount of time series) rows and T (length of time series) columns.
#' @importFrom Rdpack reprompt
#' @param N number of time series
#' @param TT length of time series
#' @param true_number_of_groups true number of groups
#' @param number_external_variables number of available observed variables
#' @param true_number_of_common_factors true number of common_factors
#' @param true_number_of_group_factors true number of group factors
#' @export
create_data_dgp2 <- function(N, TT, true_number_of_groups = 3, number_external_variables = 3, true_number_of_common_factors = 0, true_number_of_group_factors = c(3,3,3)) {
  #true group membership
  g_true <- ceiling(runif(N) * true_number_of_groups)
  #true beta
  beta_true <- create_true_beta(number_external_variables, N, true_number_of_groups,
                                FALSE, TRUE,
                                FALSE, limit_true_groups = 12)

  #true group-specific factors and loadings
  temp <- generate_grouped_factorstructure(true_number_of_groups, true_number_of_group_factors, TT, g_true)
  factor_group_true <- temp[[1]]
  lambda_group_true <- temp[[2]]
  rm(temp)

  #true common factor structure
  comfactor_true <- matrix(t(rnorm(TT * true_number_of_common_factors)), nrow = true_number_of_common_factors)
  lambda_true <- matrix(t(rnorm(N * true_number_of_common_factors)), nrow = true_number_of_common_factors)

  #make object X: gives random values to the number_external_variables variables
  X <- initialise_X(N, TT, number_external_variables)

  #errorterm
  epsilon <- matrix(rnorm(N * TT, sd = 1), nrow = N)

  Y <- generate_Y(N, TT,
                  true_number_of_common_factors,
                  true_number_of_group_factors ,
                  g_true, beta_true, lambda_group_true, factor_group_true,
                  lambda_true, comfactor_true, epsilon, X,
                  number_of_variables = number_external_variables)
  return(list(Y, X, g_true, beta_true, factor_group_true, lambda_group_true, comfactor_true, lambda_true))
}

#' This function selects a subsample of the time series, and of the length of the time series. Based on this it returns a subsample of Y.
#' It also returns the corresponding subsamples X and the group membership.
#' @param Y input Y
#' @param X input X
# @param sampleN subsample of time series
# @param sampleT subsample of time moments
#' @param number_of_time_series_fulldata number of time series of the original dataset
#' @param length_of_time_series_fulldata length of time series of the original dataset
#' @param stepsize_N size of the decrease in N; multiplied with stepsize
#' @param stepsize_T size of the decrease in T; multiplied with stepsize
#' @param stepsize index of the subsample: this defines how many times stepsize_N is subtracted from the original N time series. Similar for stepsize_T.
#' @param k true number of common factors
#' @param kg vector with for each group the true number of group specific factors
#' @inheritParams generate_Y
#' @param subsamples_factors determines whether subsamples of the (unobservable) factors and loadings also need to be constructed. This cannot be done for real world data. Default is TRUE.
#' @export
make_subsamples <- function(Y, X, factor_true, lambda_true, factor_group_true, lambda_group_true,
                            number_of_time_series_fulldata, length_of_time_series_fulldata, stepsize,
                            stepsize_N = round(number_of_time_series_fulldata / 10),
                            stepsize_T = round(length_of_time_series_fulldata / 30), k, kg, subsamples_factors = TRUE,
                            verbose = TRUE) {

  #define size of the subsample
  subN = number_of_time_series_fulldata - stepsize * stepsize_N
  subT = length_of_time_series_fulldata - stepsize * stepsize_T
  if(!(subN > 0 & subT > 0)) {
    stop("subN or subT < 0 -> stop")
  }

  #define which N and T are in the subsample
  sampleN = sort(sample(1:number_of_time_series_fulldata, subN))
  sampleT = sort(sample(1:length_of_time_series_fulldata, subT))

  #subsample of Y
  Y = Y[sampleN, sampleT]

  #subsample of X
  if(dim(X)[3] > 0) {
    X_temp = array(NA, dim = c(subN, subT, dim(X)[3]))   #subT*T_FACTOR
    X_temp[,,] = X[sampleN, sampleT,]
    X = X_temp
    rm(X_temp)
  }

  #subsample of true group membership
  g_true = g_true[sampleN]

  if(subsamples_factors) {
    #subsample of the factors and their loadings
    if(k > 0) {
      factor_true = matrix(factor_true[, sampleT], nrow = k)
      lambda_true = matrix(lambda_true[, sampleN], nrow = k)
    }
    if(max(kg) > 0) {
      for(ii in 1:aantalgroepen_real) {
        factor_group_true[[ii]] = matrix(factor_group_true[[ii]][, sampleT], nrow = kg[ii])
        lambda_group_true = lambda_group_true %>% filter(id %in% sampleN)
      }
    }
  } else {
    factor_true = NA
    lambda_true = NA
    factor_group_true = NA
    lambda_group_true = NA
  }
  return(list(Y, X, g_true, factor_true, lambda_true, factor_group_true, lambda_group_true))
}

#' Defines the object that will be used to define a initial clustering.
#'
#' This is a short version of define_object_for_initial_clustering() which only contains implementations for robust macropca case and classical case.
#' @param Y Y
#' @param beta_est estimated values of beta
#' @param use_robust TRUE or FALSE: defines using the classical, or one of the robust algorithms
#' @param method_estimate_beta defines how beta is estimated. Default case is an estimated beta for each individual. Default value is "individual." Possible values are "homogeneous", "group" or "individual".
#' @param method_estimate_factors specifies the robust algorithm to estimate factors: default is "macro". The value is not used when use_robust is set to FALSE.
#' @export
define_object_for_initial_clustering_macropca <- function(Y, beta_est, use_robust = TRUE, method_estimate_beta = "individual", method_estimate_factors = "macro") {

  length_of_time_series = ncol(Y)
  number_of_initial_factors = 10
  stopifnot(method_estimate_beta == "individual")
  stopifnot(method_estimate_factors == "macro")

    if(use_robust) {
      if(method_estimate_factors == "macro") {
        VEC <- robustpca(Y, number_of_initial_factors) #these are the eigenvectors
        message("-")
        factor_for_grouping <- sqrt(length_of_time_series) * (VEC[[1]])[, 1:number_of_initial_factors]
        message("--")
        print(dim(factor_for_grouping))
        lambda_for_grouping = t(return_robust_lambdaobject(Y, NA, type = 4, g, factor_group_rrn = factor_for_grouping))
        message("---")
      }
    } else {
      #classical
      VEC = eigen(t(Y) %*% Y)$vectors
      factor_for_grouping <- sqrt(length_of_time_series) * (VEC)[, 1:number_of_initial_factors]
      lambda_for_grouping <- t(factor_for_grouping) %*% t(Y)/length_of_time_series
    }

    to_divide = t(lambda_for_grouping)
    rm(factor_for_grouping, lambda_for_grouping)

  return(to_divide)
}

#' Function that clusters time series in a dataframe with kmeans (classical algorithm) or trimmed kmeans(robust algorithms).
#'
#' If a time series contains NA's a random cluster will be assigned to that time series.
#' @inheritParams define_object_for_initial_clustering_macropca
#' @param df dataframe to cluster
#' @param k the desired number of groups
#' @param max_percent_outliers_tkmeans The proportion of observations to be trimmed.
#' @importFrom stats kmeans
#' @importFrom tclust tkmeans
#' @export
initialise_clustering <- function(df, use_robust, k, max_percent_outliers_tkmeans = 0) {
  #in the kmeans-call NA's are not allowed:
  # -> function handleNA() will drop time series with NA's.
  number_of_time_series = nrow(df)
  if(anyNA(df)) {
    temp = handleNA(df)
    df = temp[[1]]
    rows_with_NA = temp[[2]]
    rows_without_NA = which(!(1:number_of_time_series %in% rows_with_NA))
  } else {
    rows_with_NA = NA
    rows_without_NA = 1:number_of_time_series
  }


  if(use_robust) { #only use trimmed kmeans in the robust version
    counter = 0
    Km = tkmeans(df, k, max_percent_outliers_tkmeans)

    #solve possible issue with empty clusters
    while(max_percent_outliers_tkmeans > 0 & length(table(Km$cluster)) != (k + 1) & counter < 5) {
      message("try again until there are no empty clusters (which means that there is no warning message")
      counter = counter + 1
      Km = tkmeans(df, k, max_percent_outliers_tkmeans)
    }
    #trimmed kmeans identifies outliers and puts them in zero (if parameter > 0) -> they would thus need to be replaced
  } else {
    #use kmeans
    AA = tryCatch(Km <- kmeans(df, k), error = function(e) e)
    if("error" %in% class(AA)) {
      print(AA)
      #more cluster centers than distinct data points.
      message("Km does not exist now: take as initial kmeans one group less + change 1 element to the empty group")
      counter = 0
      while("error" %in% class(AA)) {
        counter = counter + 1
        AA = tryCatch( Km <- kmeans(df, k - counter),error = function(e) e)
        if(counter > 50) {
          stop("infinite loop -> stop")
        }
      }


      indices_to_change_randomly = ceiling(runif(counter, max = nrow(df)))
      for(tel in 0:(counter-1)) {
        Km$cluster[indices_to_change_randomly[tel + 1]] = k - tel
      }
    }
  }

  g = rep(NA, number_of_time_series)
  g[rows_without_NA] = Km$cluster
  #time series with NA do not get any result with kmeans (since they were deleted up front) -> assign a random initial group for those
  g[rows_with_NA] = ceiling(runif(length(rows_with_NA)) * k)
  return(g)
}

#' This is a short version of initialise_commonfactorstructure() which only contains implementations for robust macropca case and classical case.
#'
#' @inheritParams define_object_for_initial_clustering_macropca
#' @inheritParams generate_Y
#' @param g estimated group membership
#' @param NCF_est number of estimated common factors
#' @export
initialise_commonfactorstructure_macropca <- function(use_robust, Y, beta_est, g, method_estimate_beta, method_estimate_factors, NCF_est) {

  stopifnot(method_estimate_beta == "individual")
  stopifnot(method_estimate_factors == "macro")
  number_of_time_series = nrow(Y)
  length_of_time_series = ncol(Y)
  if((method_estimate_beta == "individual")) {
    if(NCF_est == 0) {
      comfactor = t(matrix(rep(0, length_of_time_series)))
      lambda = t(matrix(rep(0, number_of_time_series)))
    } else {
      comfactor = estimate_factor(use_robust, beta_est, g, NA, method_estimate_factors, initialise = TRUE)[[1]][1:NCF_est,, drop=FALSE]
      lambda = calculate_lambda(use_robust, beta_est, comfactor, g, NA, method_estimate_factors, initialise = TRUE)[1:NCF_est,, drop=FALSE]

    }
  }
  return(list(comfactor, lambda))
}

#' Wrapper around estimate_beta, update_g, and estimating the factorstructures.
#'
#' @inheritParams estimate_beta
#' @param S number of groups to estimate
#' @param k number of common factors to estimate
#' @param kg vector with length S. Each element contains the number of group specific factors to estimate.
#' @inheritParams define_object_for_initial_clustering_macropca
#' @param special_case_dgp1 TRUE or FALSE: whether data is generated from dgp1 and has the extra spread in group centers. Default is FALSE.
#' @param vars number of variables available (in X)
#' @param vars_est number of variables for which a beta is estimated. Usually equal to vars.
#' @inheritParams update_g
iterate <- function(use_robust, Y, X, lambda_group, factor_group, lambda, comfactor, S, k, kg, method_estimate_beta,
                    vars = 3, vars_est = 3, special_case_dgp1 = FALSE, verbose = FALSE) {
  if(verbose) message("update beta")
  beta <- estimate_beta(use_robust, Y, X, lambda_group, factor_group, lambda, comfactor, NN = nrow(Y), TT = ncol(Y),
                        number_of_groups = S, number_of_group_factors = kg, number_of_common_factors = k,
                        number_of_variables = vars, number_vars_estimated = vars_est, num_factors_may_vary = TRUE, special_case_dgp1)[[1]]
  if(verbose) message("update group membership")
  g <- update_g(use_robust, NN = nrow(Y), TT = ncol(Y), number_of_groups = S, number_of_variables = vars, number_vars_estimated = vars_est,
                number_of_group_factors = kg, number_of_common_factors = k,
                method_estimate_factors)[[1]]
  if(verbose) print(table(g))
  if(verbose) message("update common factorstructure")
  lgfg_list <- calculate_lgfg(lambda_group, factor_group, NN = nrow(Y), TT = ncol(Y),  number_of_groups = S, number_of_group_factors = kg,
                              number_of_common_factors = k, num_factors_may_vary = TRUE)
  comfactor = estimate_factor(use_robust, beta_est, g, lgfg_list, method_estimate_factors,
                              NN = nrow(Y), TT = ncol(Y),
                              number_of_common_factors = k,
                              number_of_variables = vars,
                              number_vars_estimated = vars_est,
                              verbose = FALSE)[[1]]
  lambda = calculate_lambda(use_robust, beta_est, comfactor, g, lgfg_list, method_estimate_factors,
                            NN = nrow(Y), TT = ncol(Y),
                            number_of_common_factors = k,
                            number_of_variables = vars,
                            number_vars_estimated = vars_est)

  if(verbose) message("update group specific factorstructure")
  factor_group <- estimate_factor_group(use_robust, beta, g, lambda, comfactor, method_estimate_factors,
                                        NN = nrow(Y),
                                        TT = ncol(Y),
                                        number_of_groups = S,
                                        number_of_group_factors = kg,
                                        number_of_common_factors = k,
                                        number_of_variables = vars,
                                        number_vars_estimated = vars_est)
  lambda_group <- calculate_lambda_group(use_robust, beta, factor_group, g, lambda, comfactor, method_estimate_factors,
                                         NN = nrow(Y), TT = ncol(Y), number_of_groups = S, number_of_group_factors = kg,
                                         number_of_common_factors = k,
                                         number_of_variables = vars, number_vars_estimated = vars_est, UPDATE1 = FALSE, UPDATE2 = FALSE)

  if(verbose) message("calculate objective function")
  grid <- expand.grid(1:nrow(Y), 1:ncol(Y))
  grid <- grid_add_variables(grid, beta, lambda, comfactor, method_estimate_beta, NN = nrow(Y), TT = ncol(Y),
                             number_of_variables = vars, number_vars_estimated = vars_est)

  #value to minimize:
  value = OF_vectorized3(g, grid,
                         beta_est = beta,
                         lc = lambda, fc = comfactor,
                         lg = lambda_group, fg = factor_group,
                         NN = nrow(Y),
                         number_of_groups = S,
                         number_of_common_factors = k,
                         number_of_group_factors = kg,
                         num_factors_may_vary = TRUE)
  if(verbose) print(value)
  return(list(beta, g, comfactor, lambda, factor_group, lambda_group, value))
}
