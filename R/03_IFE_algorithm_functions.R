

#' Function used in generating simulated data with non normal errors.
#'
#' Used to include cross-sectional dependence or serial dependence into the simulated panel data. It is based on Ando/Bai-code.
#' @param parameter amount of cross-sectional dependence
#' @inheritParams estimate_theta
#' @return NxN covariance matrix
#' @examples
#' create_covMat_crosssectional_dependence(0.3,300)
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

#' Helpfunction in create_theta_real() for the option theta_real_heterogeen_groups. (This is the default option.)
#'
#' @param number_of_variables number of variables
#' @param number_of_groups_real real numbr of groups
#' @param EXTRA_THETA_FACTOR multiplies coefficients in theta; default = 1
theta_real_heterogroups <- function(number_of_variables, number_of_groups_real, EXTRA_THETA_FACTOR = 1) {
  stopifnot(number_of_groups_real < 10) #Code allows up to 9 real groups at this point. IF higher limit needed, add more code in this function.

  #These are the values for DGP 3 & 4 (For DGP 1 & 2 theta_real is defined in 08_IFE_make_AB_DGP1.R)

  #1st element is the intercept.
  theta_part1 = c(0, 4,  3.5,  3,  2.5)  * EXTRA_THETA_FACTOR
  theta_part2 = c(0,-2.5, -2, -2.5, -2)  * EXTRA_THETA_FACTOR
  theta_part3 = c(0, 1,  0.5,  1.5,  1)  * EXTRA_THETA_FACTOR
  theta_part4 = c(0,(round(runif(4),1) - 0.5) * 4) * EXTRA_THETA_FACTOR#between -2 and 2
  theta_part5 = c(0,(round(runif(4),1) - 0.5) * 4) * EXTRA_THETA_FACTOR#between -2 and 2
  theta_part6 = c(0,(round(runif(4),1) - 0.5) * 4) * EXTRA_THETA_FACTOR#between -2 and 2
  theta_part7 = c(0,(round(runif(4),1) - 0.5) * 4) * EXTRA_THETA_FACTOR#between -2 and 2
  theta_part8 = c(0,(round(runif(4),1) - 0.5) * 4) * EXTRA_THETA_FACTOR#between -2 and 2
  theta_part9 = c(0,(round(runif(4),1) - 0.5) * 4) * EXTRA_THETA_FACTOR#between -2 and 2

  if(number_of_variables > 3) { #add when there are more than 3 observable variables:
    theta_part1 = c(theta_part1, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part2 = c(theta_part2, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part3 = c(theta_part3, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part4 = c(theta_part4, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part5 = c(theta_part5, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part6 = c(theta_part6, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part7 = c(theta_part7, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part8 = c(theta_part8, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
    theta_part9 = c(theta_part9, (round(runif(number_of_variables - 3),1) - 0.5) * 4)
  }

  #Add the values of the different groups together in one object.
  if(number_of_groups_real >= 1) {
    theta_real = matrix(theta_part1[1:(number_of_variables + 1)])
  }
  if(number_of_groups_real >= 2) theta_real = theta_real %>% cbind(theta_part2[1:(number_of_variables + 1)])
  if(number_of_groups_real >= 3) theta_real = theta_real %>% cbind(theta_part3[1:(number_of_variables + 1)])
  if(number_of_groups_real >= 4) theta_real = theta_real %>% cbind(theta_part4[1:(number_of_variables + 1)])
  if(number_of_groups_real >= 5) theta_real = theta_real %>% cbind(theta_part5[1:(number_of_variables + 1)])
  if(number_of_groups_real >= 6) theta_real = theta_real %>% cbind(theta_part6[1:(number_of_variables + 1)])
  if(number_of_groups_real >= 7) theta_real = theta_real %>% cbind(theta_part7[1:(number_of_variables + 1)])
  if(number_of_groups_real >= 8) theta_real = theta_real %>% cbind(theta_part8[1:(number_of_variables + 1)])
  if(number_of_groups_real >= 9) theta_real = theta_real %>% cbind(theta_part9[1:(number_of_variables + 1)])


  return(theta_real)
}

#' Creates theta_real, which contains the real values of theta (= the coefficients of X)
#'
#' theta_real_heterogeen_groups is the default case
#' @inheritParams estimate_theta
#' @param number_of_groups_real number of groups
#' @param eclipz boolean to indicate using eClipzdatabase; defaults to FALSE
#' @param extra_theta_factor multiplies coefficients in theta; default = 1
#' @return matrix with number of rows equal to number of observable variables + 1 (the first row contains the intercept) and number of culumns
#' equal to the real number of groups.
#' @examples
#' library(tidyverse)
#' #Decide if theta is common, or specific to groups or individuals: Choose 1 of the following 3.
#' theta_real_homogeen = FALSE
#' theta_real_heterogeen_groups = TRUE
#' theta_real_heterogeen_individueel = FALSE
#' create_theta_real(3,NN=300,number_of_groups_real = 3)
#' @export
create_theta_real <- function(number_of_variables,
                              NN = aantal_N,
                              number_of_groups_real = aantalgroepen_real,
                              eclipz = FALSE,
                              extra_theta_factor = 1) {
  #real world eclipzdata: theta_real does not exist -> return NA
  if(eclipz) return(NA)


  if(number_of_variables > 0) {
    #common theta_real:
    if(theta_real_homogeen) {
      theta_real = c(0, c(1,2,3,28,33,38,43,48,53,58,63,68,73,78,83)[1:number_of_variables]) #intercept 0, and after that values for the real theta's
      theta_real = matrix(rep(theta_real, number_of_groups_real), nrow = (number_of_variables + 1))
    }
    #groupsspecific theta_real: -> default case
    if(theta_real_heterogeen_groups) {
      theta_real = theta_real_heterogroups(number_of_variables, number_of_groups_real, EXTRA_THETA_FACTOR = extra_theta_factor)
    }
    #individualspecifiec theta_real:
    if(theta_real_heterogeen_individueel) {
      theta_real = matrix(rnorm(NN * (number_of_variables + 1)), nrow = (number_of_variables + 1), ncol = NN)
    }
  } else {
    theta_real = rep(NA,number_of_groups_real)
  }

  return(theta_real)
}


#' Creates X (the observable variables) to use in simulations.
#'
#' X is an array with dimensions N, T and number of observable variables.
#' The variables are randomly generated with mean 0 and sd 1.
#' @inheritParams create_covMat_crosssectional_dependence
#' @inheritParams estimate_theta
#' @return array with dimensions N x T x number of observable variables
#' @examples
#' initialise_X(300,30, number_of_variables = 3)
#' @export
initialise_X <- function(NN,TT, number_of_variables = aantalvars) {
  if(number_of_variables > 0) {
    X = array(0,dim=c(NN, TT, number_of_variables))
    for(i in 1:NN) {
      for(t in 1:TT) {
        for(k in 1:number_of_variables) {
          X[i,t,k] = rnorm(1, mean = 0, sd = 1)
        }
      }
    }

    return(X)
  } else {
    return(NA)
  }
}


#'Scaling of X.
#'
#' @param X input
#' @param firsttime Scaling before generating Y and before adding outliers -> this is always with mean and sd. If this is FALSE, it indicates that
#' @param eclipz parameter to indicate using real world Eclipzdataset. Defaults to FALSE.
#' we are using the function for a second time, after adding the outliers. In the robust case it uses median and MAD, otherwise again mean and sd.
#' @inheritParams create_theta_real
#' @examples
#' X = initialise_X(300,30, number_of_variables = 3)
#' use_robust = TRUE
#' scaling_X(X,TRUE, number_of_variables = 3)
#' @export
scaling_X <- function(X, firsttime, eclipz = FALSE, number_of_variables = aantalvars) {

  if(number_of_variables > 0) {
    for(k in 1:number_of_variables) {
      #print(paste("sd of variable",k,": Before:",sd(X[,,k])))
      if(mad(X[,,k]) != 0) {
        if(use_robust & !firsttime) {
          message("Scale with median and mad")
          med = median(X[,,k])
          mad = mad(X[,,k])

          X[,,k] = (X[,,k] - med) / mad #Note that this is NOT column-based (timeindex) scaling!
        } else {
          if(eclipz) {
            message("Scale with mean and sd of NxT-matrix")
            X[,,k] = (X[,,k] - mean(X[,,k])) / sd(X[,,k]) #for eclipzdata, we cannot use scale(),
            #  since there are variables (for example age) that have constant columns, so scale() (which is columnbased) would produce errors
          } else {
            message("Scale with mean and sd (for each t separate)")
            X[,,k] = scale(X[,,k])      #Note that this is column-based (timeindex) scaling!
          }
        }
        print(paste("sd of variable",k,": After:",sd(X[,,k])))
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
#' @inheritParams create_theta_real
#' @inheritParams estimate_theta
#' @examples
#' X = initialise_X(300,30, number_of_variables = 3)
#' X_restructured = restructure_X_to_order_slowN_fastT(X, FALSE,
#'   number_of_variables = 3, number_vars_estimated = 3)
#' @export
restructure_X_to_order_slowN_fastT <- function(X, eclipz,
                                               number_of_variables = aantalvars,
                                               number_vars_estimated = SCHATTEN_MET_AANTALVARS) {

  if(number_of_variables > 0) {
    if( length(dim(X)) == 2) {
      #occurs when only one element in group
      X = array(X, dim=c(1,nrow(X),ncol(X)))
    }
    if(eclipz) {
      number_of_vars = number_of_variables
    } else {
      number_of_vars = number_vars_estimated
    }
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


#' Generates the real groupfactorstructure, to use in simulations.
#'
#' Loadings and factors are generated by:
#' loadings ~ N(LGR_FACTOR_mean, j * LGR_FACTOR) -> default case will be N(0,j)
#' factors ~ N(j * FGR_FACTOR, FGR_FACTOR_sd) -> default case will be N(j,1)
#' @param S real number of groups
#' @param number_of_group_factors_real vector with as length the number of groups, where each element is the number of real groupfactors of that group.
#' @inheritParams estimate_theta
#' @return list: first element contains real groupfactors and second element contains real groupfactor loadings
#' @examples
#' library(tidyverse)
#' #3 groups, each with 3 groupfactors:
#' g_real = ceiling(runif(300) * 3)
#' LGR_FACTOR_mean = 0
#' LGR_FACTOR = 1
#' FGR_FACTOR_sd = 1
#' FGR_FACTOR = 1
#' generate_grouped_factorstructure(3,c(3,3,3), TT = 30)
#' @export
generate_grouped_factorstructure <- function(S, number_of_group_factors_real, TT = aantal_T_fulldata) {
  if(mean(number_of_group_factors_real) > 0) {


    LGR = list()
    factor_group_real = list()
    for(i in 1:S) {
      LGR[[i]] = matrix(nrow = sum(g_real == i), ncol = number_of_group_factors_real[i])
      factor_group_real[[i]] = matrix(NA, nrow = number_of_group_factors_real[i], ncol = TT)
    }

    for(j in 1:S) {
      for(i in 1:number_of_group_factors_real[j]) {
        for(k in 1:(sum(g_real == j))) {
          LGR[[j]][k,i] = rnorm(1, mean = LGR_FACTOR_mean, sd = j * LGR_FACTOR)
        }
        for(k in 1:TT) {
          factor_group_real[[j]][i,k] = rnorm(1, mean = j * FGR_FACTOR, sd = FGR_FACTOR_sd)
        }
      }
    }
    #lambda_group_real to dataframe (easier to work with)
    lambda_group_real = data.frame(LGR[[1]]) %>% mutate(groep = 1, id =  which(g_real == 1))
    if(number_of_group_factors_real[1] == 1) { #change name when only 1 factor (when there are more names are automatically X1,X2,...)
      names(lambda_group_real)[1] = "X1"
    }
    for(i in 2:S) {
      temporary = data.frame(LGR[[i]]) %>% mutate(groep = i, id =  which(g_real == i))
      if(number_of_group_factors_real[i] == 1) {
        names(temporary)[1] = "X1"
      }
      if(ncol(lambda_group_real) >= ncol(temporary)) {
        lambda_group_real = lambda_group_real %>% bind_rows(temporary)
      } else {
        lambda_group_real = temporary %>% bind_rows(lambda_group_real)
      }
    }
    lambda_group_real = lambda_group_real %>% arrange(groep)
    lambda_group_real[is.na(lambda_group_real)] <- 0
  } else {
    factor_group_real = NA
    lambda_group_real = NA
  }

  #geval van 1 factor in group: format to matrix:
  for(groep in 1:S) {
    if(number_of_group_factors_real[groep] == 1) {
      factor_group_real[[groep]] = matrix(factor_group_real[[groep]], nrow = number_of_group_factors_real[groep])
    }
  }
  return(list(factor_group_real, lambda_group_real))
}


#' Generate panel data Y for simulations.
#' @inheritParams estimate_theta
#' @param number_of_common_factors_real real number of common factors
#' @param number_of_group_factors_real Vector of length the number of groups. Each element contains the real number of group factors for that group.
#' @param g_real vector with real group memberships
#' @param theta_real real coefficients with the observable variables
#' @param lambda_group_real loadings of the group factors
#' @param factor_group_real groupfactors
#' @param lambda_real loadings of the common factors
#' @param factor_real common factors
#' @param epsilon NN x TT-matrix containing the error term
#' @param ando_bai loads Ando/Bai-dataset; only for testing purposes. Defaults to FALSE.
#' @param ando_bai_2017 loads Ando/Bai-dataset; only for testing purposes. Defaults to FALSE.
#' @inheritParams create_theta_real
#' @return N x T matrix
#' @export
generate_Y <- function(NN, TT, number_of_common_factors_real,number_of_group_factors_real,
                       g_real, theta_real, lambda_group_real, factor_group_real,
                       lambda_real, factor_real, epsilon,
                       ando_bai = FALSE,
                       ando_bai_2017 = FALSE,
                       eclipz = FALSE,
                       number_of_variables = aantalvars) {

  #Define the size of the panel data:



  Y = matrix(NA, nrow = NN, ncol = TT) #initialisation, later on this gets filled in
  # if(abctypes) { #
  #   Y = create_Y_abctypes(Y, proporties_types)
  # } else {

    for(i in 1:NN) {

      if(!ando_bai & !ando_bai_2017 & !eclipz) {
        if(mean(number_of_group_factors_real) > 0) {
          dropvars <- names(lambda_group_real) %in% c("groep","id")
          LAMBDAGROUP = as.matrix(subset(lambda_group_real, lambda_group_real$id == i)[!dropvars])
          LAMBDAGROUP = LAMBDAGROUP[1:number_of_group_factors_real[g_real[i]]]
        } else {
          LAMBDAGROUP = NA
        }
      }


      for(t in 1:TT) {

        if(number_of_variables > 0) {
          XT = c(1, X[i,t,]) %*% theta_real[,g_real[i]] #add 1 to X[i,t,] to make room for the intercept (which is in theta)

        } else {
          XT = 0
        }


        #interactive fixed effects: both common and group-speficic factors
        if(ando_bai) {
          Y[i,t] = XT +
            XL[t,i] + ERR[t,i]
        } else if(ando_bai_2017) {
          Y[i,t] = XT +
            XG[t,i] +
            XL[t,i] + ERR[t,i]
        } else {
          #randomly generated data
          if(number_of_common_factors_real > 0) {
            LF = t(lambda_real[,i]) %*% factor_real[,t]
          } else {
            LF = 0
          }
          if(mean(number_of_group_factors_real) > 0) {
            LF_GROUP = LAMBDAGROUP %*% factor_group_real[[g_real[i]]][,t]
          } else {
            LF_GROUP = 0
          }

          Y[i,t] = XT +
            LF +
            LF_GROUP +
            epsilon[i,t]
          stopifnot(!is.nan(Y[i,t]))
        }
      }
    }

  return(Y)
}


#' Initialisation of estimation of beta (the coefficients with the observable variables)
#'
#' Note: this needs to be called before the definition of grid.
#' @inheritParams generate_Y
#' @param number_vars_estimated number of variables from which the coefficients are estimated
#' @param number_of_groups number of groups
#' @examples
#' library(RobClustTimeSeries)
#' library(tidyverse)
#' library(robustbase)
#'
#' X = RobClustTimeSeries::X_dgp3
#' Y = RobClustTimeSeries::Y_dgp3
#' #Set estimations for group factors and its loadings, and group membership to the real value
#' lambda_group = RobClustTimeSeries::lambda_group_real_dgp3
#' factor_group = RobClustTimeSeries::factor_group_real_dgp3
#' g = RobClustTimeSeries::g_real_dgp3
#' #There are no common factors to be estimated  -> but needs placeholder
#' lambda = matrix(0,nrow=1,ncol=300)
#' comfactor = matrix(0,nrow=1,ncol=30)
#'
#' use_robust = TRUE
#' #Choose how coefficients of the observable are estimated
#' homogeneous_coefficients = FALSE
#' heterogeneous_coefficients_groups = FALSE
#' heterogeneous_coefficients_individuals = TRUE #estimating theta_i for every individual
#' ABDGP1 = FALSE
#' ABintercept = TRUE
#' theta_init = initialise_theta(NN = 300, TT = 30,
#'   number_of_variables = 3, number_vars_estimated = 3, number_of_groups = 3)
#' @export
initialise_theta <- function(eclipz = FALSE,
                            NN = aantal_N,
                            TT = aantal_T,
                            number_of_variables = aantalvars,
                            number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                            number_of_groups = aantalgroepen) {

  stopifnot((homogeneous_coefficients + heterogeneous_coefficients_groups + heterogeneous_coefficients_individuals) == 1)

  if(eclipz) {
    number_of_vars = number_of_variables
  } else {
    number_of_vars = number_vars_estimated
  }

  if(number_of_vars > 0) {
    theta = matrix(NA, nrow = (number_of_vars + 1), ncol = number_of_groups)
    if(heterogeneous_coefficients_individuals) {
      theta = matrix(NA, nrow = (number_of_vars + 1), ncol = NN)
    }


    #X needs to be in the form of (NN*TT x p matrix)

    if(homogeneous_coefficients) { #also used in DGP05: BramatiCroux
      X_special = X_restructured
      Y_special = Y
      #this includes robuust estimation of theta:
      theta = determine_theta("homogeen", X_special, Y_special, TRUE, initialisatie = TRUE, indices = 1:NN,  TT = TT, number_of_variables = number_of_variables)

    }
    if(heterogeneous_coefficients_groups) {
      for(group in 1:number_of_groups) { #for each  real group there must be a column in theta.
        #select parts of X and Y of this group
        indices_group = which(g == group)
        if(length(indices_group) == 0) {
          message("There is an empty group!")
        }

        #X needs to be in the form of (NN*TT x p matrix)
        X_special = restructure_X_to_order_slowN_fastT(array(X[indices_group,,], dim = c(length(indices_group), TT, number_of_vars)),
                                                       eclipz,
                                                       number_of_variables = number_of_variables,
                                                       number_vars_estimated = number_vars_estimated)

        Y_special = as.vector(t(Y[indices_group,])) #order: N1T1, N1T2,N1T3,...N2T1,...N_endT_end

        theta[,group] = determine_theta("heterogeen", X_special, Y_special, TRUE, initialisatie = TRUE, indices = indices_group,  TT = TT, number_of_variables = number_of_variables)

      }

    }
    if(heterogeneous_coefficients_individuals) {


      #Initialisation: use of lm here instead of ncvreg (cf AndoBai-code)
      for(i in 1:NN) {
        X_special = restructure_X_to_order_slowN_fastT(matrix(X[i,,], ncol = dim(X)[3]),
                                                       eclipz,
                                                       number_of_variables = number_of_variables,
                                                       number_vars_estimated = number_vars_estimated)
        Y_special = as.vector(t(Y[i,]))

        if(use_robust) {
          if(ABDGP1 & ABintercept) { #This is DGP 1
            #no intercept, because ABDGP1 defines the first variable in X as an intercept
            model <- LMROB(Y_special,X_special, nointercept = TRUE)  #-> lmrob(Y_special ~ X_special + 0,setting="KS2014)
          } else {
            model <- LMROB(Y_special,X_special)
          }
        } else {
          if(ABDGP1 & ABintercept) { #This is DGP 1
            model <- lm(Y_special ~ X_special + 0)
          } else {
            model <- lm(Y_special ~ X_special)

          }
        }
        if(ABDGP1 & ABintercept) {
          values = c(0,model$coefficients)
          theta[,i] = values
        } else {
          theta[,i] = model$coefficients
        }
      }
    }
    rm(X_special, Y_special)


    return(theta)
  } else {
    #when number_of_vars == 0, we do not need theta
    return(rep(NA,number_of_groups))
  }
}


#' Helpfunction used in update_g()
#'
#' We calculate FgLg (the groupfactorstructure) for all  possible groups where individual i can be placed. For each group we have before estimated
#' the groupfactors (Fg). Now we need to the grouploadings for each group as well. In the classical case these are calculated by Fg*Y/T. In the robust case
#' these are robust.
#' @param group number of group
#' @param solve_FG_FG_times_FG This is the same as groupfactor / T. It is only used in the Classical approach
#' @inheritParams estimate_theta
#' @return NxT matrix containing the product of virtual groupfactors and virtual loadings

calculate_virtual_factor_and_lambda_group <- function(group, solve_FG_FG_times_FG,
                                                      NN, TT,
                                                      number_of_variables = aantalvars,
                                                      number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                                                      number_of_group_factors = aantalfactoren_groups,
                                                      number_of_common_factors = aantalfactoren_common) {
  FG = factor_group[[group]]
  indices = 1:NN
  LF = t(lambda) %*% comfactor
  xtheta = calculate_XT_estimated(NN = NN, TT = TT, number_of_variables = number_of_variables, number_vars_estimated = number_vars_estimated)


  if(number_of_common_factors == 0) {
    Y_ster = Y[indices,] - xtheta
  } else {
    Y_ster = Y[indices,] - xtheta - LF
  }


  #robust grouplambda:
  if(use_robust) {
    #we need a robust version of the virtual factorstructure:
    LG_local = return_robust_lambdaobject(Y_ster, group, type = 1,
                                          NN = NN,
                                          number_of_group_factors = number_of_group_factors,
                                          number_of_common_factors = number_of_common_factors)
  } else {
    LG_local = t(solve_FG_FG_times_FG[[group]] %*% t(Y_ster)) #This equal to Fg*Y/T
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
define_rho_parameters <- function(object = NULL) {
  if(is.null(object)) {
    print("This part is probably obsolete (define_rho_parameters). Sleep if encountered.")
    Sys.sleep(108000)
    # e = calculate_error_term()
    # all_norms = apply(e, 1, function(x) matrixnorm(x))
    #
    # #use median and madn in the rho-function:
    # rho_loc = c(median(all_norms[g==1]),
    #             median(all_norms[g==2]),
    #             median(all_norms[g==3])
    #             )
    # rho_scale = c(mad(all_norms[g==1]),
    #             mad(all_norms[g==2]),
    #             mad(all_norms[g==3]))  #this is actually 1.4826 * mad (or mad / 0.675), so this is madn (normalized mad)
  } else {
    rho_loc =  apply(data.frame(object),1, function(x) median(x)) #=median over T elements
    rho_scale = apply(data.frame(object),1, function(x) mad(x)) #=normalized mad over T elements


    #print(median(unlist(lapply(list(rho_loc,rho_scale),function(x) x[[1]][1]))))

  }

  return(list(rho_loc,rho_scale))
}

#' Helpfunction for update_g(). Calculates the errors for every virtual group.
#'
#' As we are updating group membership, we use the errorterm as objective function to estimate the group. We assume group membership equals 1,...,NG
#' (with NG the total number of groups) and calculate the error term.
#' @param k group
#' @param LF NxT-matrix of the commonfactorstructure
#' @param virtual_grouped_factor_structure list with length the number of groups; every element of the list contains NxT-matrix
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @return NxT matrix with errors (=Y - XB - FgLg - FL)
calculate_errors_virtual_groups <- function(k,LF,virtual_grouped_factor_structure,
                                            NN,
                                            TT,
                                            number_of_variables,
                                            number_of_common_factors,
                                            number_of_group_factors) {
  E_prep = matrix(NA, nrow = NN, ncol = TT)

  a = do_we_estimate_common_factors(number_of_common_factors)
  b = do_we_estimate_group_factors(number_of_group_factors)

  for(i in 1:NN) {
    #calculate lambda_group for one individual, based on an hypothetical groupmembership
    if(b == 1) { #-> we estimate group factors
      virtual_structure = virtual_grouped_factor_structure[[k]][i,]
    } else {
      message("Option to not estimate any group factor in any group -> this option is not implemented")
    }
    if(number_of_variables > 0) {
      if(homogeneous_coefficients | heterogeneous_coefficients_groups) {
        XT = cbind(1, X[i,,]) %*% theta[,g[i]]
        #-> this should not be theta[,k] as only the grouped factorstructure should vary with k (similar to calculate_virtual_factor_and_lambda_group)
      }
      if(heterogeneous_coefficients_individuals) {
        XT = cbind(1, X[i,,]) %*% theta[,i]
      }
    }


    for(t in 1:TT) {
      #calculate the objective function, while making sure that NA's do not have any effect in the sum
      #define the estimationerror:
      E_prep[i,t] = Y[i,t] - a * LF[i,t]
      if(b != 0) E_prep[i,t] = (E_prep[i,t] - b * virtual_structure[t])
      if(number_of_variables > 0) E_prep[i,t] = E_prep[i,t] - XT[t,]
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
#' @param ERRORS_VIRTUAL ...
#' @param rho_parameters ...
#' @param TT T
calculate_obj_for_g <- function(i, k, ERRORS_VIRTUAL, rho_parameters, TT = aantal_T) {

  totalsum = 0

  if(use_robust) {
    #define a scaling:
    #It must be individualspecific, and also be equal over all virtual groups;
    #We take here the median over the values of all groups.
    #     (because: use of #rho_parameters[[k]][[2]][i] leads tot random grouping)
    #"location" is then defined the same way
    #rho_parameters is a list of 'number_of_groups' elements. Every element has 2 elements with 'NN' values.
    # DO NOT USE FUTURE_MAP() HERE: this is way slower for some reason!
    location = (unlist(lapply(rho_parameters,function(x) x[[1]][i]))) #map over virtual groups; take 1st element (=median) and take individual i
    scaling = (unlist(lapply(rho_parameters,function(x) x[[2]][i]))) #map over virtual groups; take 2nd element (=mad) and take individual i


    location = median(location) #median over groups
    scaling = median(scaling) #median over groups

  }

  for(t in 1:TT) {

    #calculate the objective function, while making sure that NA's do not have any effect in the sum

    #calculate the estimationerror:
    E_prep = ERRORS_VIRTUAL[[k]][i,t]
    if(use_robust) { #define the rho-function (bisquare)

      E_prep_loc_scale = (E_prep - location) / scaling #this is a scalar
      E = Mpsi(E_prep_loc_scale, cc = 4.685, psi = "bisquare", deriv = -1) #rho-functie on scaled errors

    } else { #non-robuste version:
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
#' @param number_of_group_factors number of group factors
solveFG <- function(TT, number_of_groups, number_of_group_factors){
  solve_FG_FG_times_FG = list()
  for(group in 1:number_of_groups) {
    if(number_of_group_factors[group] > 0) {
      FG = factor_group[[group]]
      solve_FG_FG_times_FG[[group]] = solve(FG%*%t(FG))%*%FG #this is actually just FG/T (the code comes from AndoBai)
      rm(FG)
    } else {
      #make 0- matrix with 1 row
      solve_FG_FG_times_FG[[group]] = matrix(0, nrow = 1, ncol = TT)
    }
  }
  return(solve_FG_FG_times_FG)
}

#' Function that estimates group membership.
#'
#' @return list: 1st element contains group membership and second element contains the values which are used to determine group membership
#' @inheritParams estimate_theta
#' @examples
#' #This function needs several initial parameters to be initialized in order to work on itself
#' library(RobClustTimeSeries)
#' library(tidyverse)
#' library(robustbase)
#'
#' X = RobClustTimeSeries::X_dgp3
#' Y = RobClustTimeSeries::Y_dgp3
#' #Set estimations for group factors and its loadings, and group membership to the real value
#' lambda_group = RobClustTimeSeries::lambda_group_real_dgp3
#' factor_group = RobClustTimeSeries::factor_group_real_dgp3
#' g = RobClustTimeSeries::g_real_dgp3
#' #There are no common factors to be estimated  -> but needs placeholder
#' lambda = matrix(0,nrow=1,ncol=300)
#' comfactor = matrix(0,nrow=1,ncol=30)
#
#' use_robust = TRUE
#' grid = expand.grid(1:300,1:30)
#' #Choose how coefficients of the observable are estimated
#' homogeneous_coefficients = FALSE
#' heterogeneous_coefficients_groups = FALSE
#' heterogeneous_coefficients_individuals = TRUE
#' ABDGP1 = FALSE
#' ABintercept = TRUE
#' theta = estimate_theta(NN = 300, TT = 30,
#'   number_of_groups = 3, number_of_group_factors = c(3,3,3), number_of_common_factors = 0,
#'   number_of_variables = 3,number_vars_estimated=3,num_factors_may_vary=FALSE)[[1]]
#' grid = grid_add_variables(grid,theta, lambda, comfactor, NN = 300, TT = 30,
#'   number_of_variables = 3, number_vars_estimated = 3, number_of_groups = 3)
#' eclipz = FALSE
#' g_new = update_g(NN = 300, TT = 30, number_of_groups = 3, number_of_variables = 3,
#'   number_vars_estimated = 3,
#'   number_of_group_factors = c(3,3,3),
#'   number_of_common_factors = 0)[[1]]
#' @export
update_g <- function(NN = aantal_N, TT = aantal_T,
                     number_of_groups = aantalgroepen,
                     number_of_variables = aantalvars,
                     number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                     number_of_group_factors = aantalfactoren_groups,
                     number_of_common_factors = aantalfactoren_common) {



  if(do_we_estimate_group_factors(number_of_group_factors)) { #if there are groupfactors estimated
    solve_FG_FG_times_FG = solveFG(TT, number_of_groups, number_of_group_factors)

    #we calculate FgLg (groupfactors times grouploadings) for all the possible groups to which individual i could end up:

    virtual_grouped_factor_structure = lapply(1:number_of_groups, function(y) calculate_virtual_factor_and_lambda_group(y, solve_FG_FG_times_FG, NN, TT,
                                                                                                                        number_of_variables = number_of_variables,
                                                                                                                        number_vars_estimated = number_vars_estimated,
                                                                                                                        number_of_group_factors = number_of_group_factors,
                                                                                                                        number_of_common_factors = number_of_common_factors))


  } else {
    virtual_grouped_factor_structure = NA
  }


  LF = (t(lambda) %*% comfactor)

  #init matrix with objectivefunctionvalues for all groups
  matrix_obj_values = matrix(NA, nrow = NN, ncol = number_of_groups)


  #calculate errors for each (both virtual & real) group
  ERRORS_VIRTUAL = lapply(1:number_of_groups, function(x) calculate_errors_virtual_groups(x,LF,virtual_grouped_factor_structure, NN, TT,
                                                                                          number_of_variables,
                                                                                          number_of_common_factors,
                                                                                          number_of_group_factors))


  if(use_robust) {
    rho_parameters = lapply(1:number_of_groups,function(x) define_rho_parameters(ERRORS_VIRTUAL[[x]])) #(parameter object = NA -> returns median and madn of the calculated error term)
  } else {
    rho_parameters = NA
  }

  for(i in 1:NN) {
    obj_values = map_dbl(1:number_of_groups, function(x) calculate_obj_for_g(i, x, ERRORS_VIRTUAL, rho_parameters, TT = TT))
    g[i] = which.min(obj_values)
    matrix_obj_values[i,] = obj_values
  }

  #vectorizing is not faster: g = sapply(1:NN, function(z) which.min(map_dbl(1:number_of_groups, function(x) calculate_obj_for_g(z, x, ERRORS_VIRTUAL, rho_parameters))))
  return(list(g, matrix_obj_values))

}




#' Returns the estimated groupfactorstructure.
#'
#' This is the same function as calculate_FL_group_estimated(), but with adjustable parameters
#' @param f groupfactors
#' @param l grouploadings
#' @param j Number of groupfactors that are being used in the calculation.
#' @inheritParams calculate_FL_group_real
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @export
calculate_FL_group_estimated2 <- function(f,l,j,
                                          number_of_groups = aantalgroepen,
                                          number_of_common_factors = aantalfactoren_common,
                                          num_factors_may_vary = aantalfactors_verschillend_per_group,
                                          NN = aantal_N,
                                          TT = aantal_T) {

  temp = calculate_lgfg(l,f,number_of_groups, rep(j, number_of_groups), number_of_common_factors, num_factors_may_vary)
  FL_GROUP_geschat = lapply(1:NN, function(x) temp[[g[x]]][x,]) %>% unlist %>% matrix(nrow = TT) %>% t

  return(FL_GROUP_geschat)
}










#' Helpfunction in OF_vectorized3()
#'
#' @param i index of individual
#' @param t index of time
#' @param XTHETA matrixproduct of X and theta
#' @param LF matrixproduct of common factors and its loadings
#' @param group_memberships vector with group memberships
#' @param lgfg_list product of groupfactors and their loadings; list with length the number of groups
#' @param number_of_group_factors vector containing the number of group factors for all groups
OF_vectorized_helpfunction3 <- function(i,t,XTHETA,LF,
                                        group_memberships,
                                        lgfg_list,
                                        number_of_group_factors) {


  if(do_we_estimate_group_factors(number_of_group_factors) != 0) {
    result = as.numeric(Y[i,t] - XTHETA -
                          LF -
                          lgfg_list[[group_memberships[i]]][i,t]
    )^2

  } else {
    result = as.numeric(Y[i,t] - XTHETA -
                          LF
    )^2
  }

  return(result)
}

#' Calculates objective function: used in local_search + to determine "best_result".
#'
#' @param group_memberships Vector containing groupmembership for all individuals.
#' @param THETA theta
#' @param LAMBDA loadings of common factors
#' @param FACTOR common factors
#' @param LAMBDA_GROUP grouploadings
#' @param FACTOR_GROUP groupfactors
#' @inheritParams estimate_theta
#' @export
OF_vectorized3 <- function(group_memberships, THETA = theta,
                           LAMBDA = lambda, FACTOR = comfactor,
                           LAMBDA_GROUP = lambda_group, FACTOR_GROUP = factor_group,
                           NN = aantal_N,
                           number_of_groups = aantalgroepen,
                           number_of_common_factors = aantalfactoren_common,
                           number_of_group_factors = aantalfactoren_groups,
                           num_factors_may_vary = aantalfactors_verschillend_per_group) {
  #this is a list (length number of groups) of the product FgLg (which is the groupfactorstructure)
  lgfg_list = calculate_lgfg(LAMBDA_GROUP,FACTOR_GROUP, number_of_groups, number_of_group_factors, number_of_common_factors, num_factors_may_vary)

  if(homogeneous_coefficients | heterogeneous_coefficients_individuals) {
    return(sum(apply(grid,1,function(x) OF_vectorized_helpfunction3(x[1],x[2],x[3],x[4],group_memberships, lgfg_list, number_of_group_factors))))
  }
  if(heterogeneous_coefficients_groups) {
    #construct a vector with XTHETA-values depending on the group of the individuals:
    temp = grid %>% dplyr::select(starts_with("XTHETA"))
    XTHETA_parameter = sapply(1:NN, function(x) temp[x,g[x]])

    prep = grid %>%
      dplyr::select(-starts_with("XTHETA"))
    #used to put LF into the 3rd column -> x[3] in next line
    return(sum(apply(prep,1,function(x) OF_vectorized_helpfunction3(x[1],x[2],XTHETA_parameter,x[3],group_memberships, lgfg_list, number_of_group_factors))))
  }


}






#' Returns list (with as length the number of groups) with lgfg (product of grouploadings a&nd groupfactors).
#' Each element of the list with the assumption that all individuals are in the same group k.
#'
#' This function is used to speed up code.
#' @param lambda_group grouploadings
#' @param factor_group groupfactors
#' @inheritParams estimate_theta
#' @examples
#' library(tidyverse)
#' lambda_group = RobClustTimeSeries::lambda_group_real_dgp3
#' factor_group = RobClustTimeSeries::factor_group_real_dgp3
#' calculate_lgfg(lambda_group,factor_group,3,c(3,3,3),0,FALSE)
#' @export
calculate_lgfg <- function(lambda_group, factor_group, number_of_groups, number_of_group_factors, number_of_common_factors, num_factors_may_vary,
                           NN = aantal_N, TT = aantal_T) {
  lgfg_list = list()
  #define LgFg for each group (and later on select the correct element (correct group of individual i) of lgfg_list)
  for(k in 1:number_of_groups) {
    if(number_of_group_factors[k] > 0) {
      LG_clean = (as.matrix(lambda_group %>% arrange(id) %>% dplyr::select(-groep,-id)))[,1:number_of_group_factors[k]]

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
      lgfg_list[[q]] = matrix(0,nrow=NN,ncol=TT)
    }
  }


  return(lgfg_list)
}


#' Helpfunction in estimate_theta() for estimating theta.
#'
#' @param string can have values: "homogeen" or "heterogeen"
#' @param X_special preprocessed X (observable variables)
#' @param Y_special preprocessed Y
#' @param correct Defaults to TRUE. If TRUE, then My is not being used. My is an scaling thing that came out of Ando/Bai.
#' @param initialisatie boolean
#' @param indices individuals for which theta is being estimated
#' @param optimize_kappa indicates if kappa has to be optimized or not; defaults to FALSE
#' @inheritParams generate_Y
determine_theta <- function(string, X_special,Y_special, correct, initialisatie = FALSE, indices = NA, optimize_kappa = FALSE,
                           TT = aantal_T,
                           number_of_variables = aantalvars) {
  if(!(heterogeneous_coefficients_groups & initialisatie == TRUE)) {
    #calculate My, which is some sort of scaling thing, and subtract it from Y:
    if(string == "homogeen") {
      My <- t(rep(1,len=TT))%*%t(Y)/TT; My <- rep(1,len=TT)%*%My
    }
    Y_special = matrix(Y_special, nrow = length(indices), ncol = TT)
    if(string == "heterogeen") {
      My <- t(rep(1,len=TT))%*%t(Y_special)/TT; My <- rep(1,len=TT)%*%My
    }
    if(correct) My = My - My #set to zero

    if(initialisatie) {
      #initialisation of theta, so no factorstructure in Y_special
      Y_special = as.vector(t(Y)-My) #order: N1T1, N2T1,...N1T2,...N_endT_end
    } else {
      Y_special = as.vector(t(Y_special)-My) #order: N1T1, N2T1,...N1T2,...N_endT_end
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
    model <- LMROB(Y_special, X_special) #-> lmrob(Y_special ~ X_special, setting="KS2014")
  } else {
    #when "lmrob_and_classpca" exists, we use lmrob and classical pca
    if(exists("lmrob_and_classpca")) {
      model <- LMROB(Y_special, X_special)
    } else {

      #ncvreg, without weights
      model <- ncvreg(X_special, Y_special, family="gaussian", penalty="SCAD",lambda=(kappa_candidates))
      if(optimize_kappa) {
        minBIC <- 10^10
        best_lami = 0
        #print("indices")
        #print(indices)
        for(LAMi in 1:length(kappa_candidates)){
          theta_temp <- model$beta[,LAMi] #c(0,model$beta)[-2,LAMi]
          BIC <- sum( (Y[indices,]-cbind(1,X[indices,,]) %*% theta_temp )^2 )/(TT) + C * sigma2_max_model * log(TT)*sum(theta_temp!=0)/(TT)

          if(BIC<=minBIC){
            best_lami = LAMi
            minBIC <- BIC
            theta_temp <- model$beta[,LAMi] #c(0,model$beta)[-2,LAMi]
          }
        }
        #print(kappa_candidates[best_lami]) #-> so this is kappa[i]
        #Sys.sleep(0.1)
        return(theta_temp)
      }
    }
  }


  if(string == "homogeen") {
    if(class(model) == "lm" | class(model) == "lmrob") return(matrix(rep(model$coefficients, aantalgroepen), nrow = (number_of_variables + 1)))
    else return(matrix(rep(model$beta[,1], aantalgroepen), nrow = (number_of_variables + 1)))
  }
  if(string == "heterogeen") {
    if(class(model) == "lm" | class(model) == "lmrob") {
      if(ABDGP1 & ABintercept) { #=DGP1
        return(c(0,as.numeric(model$coefficients[-2])))
      } else {
        return(as.numeric(model$coefficients))
      }
    } else {
      if(ABDGP1 & ABintercept) { #=DGP1
        return(c(0,model$beta[-2,1]))
      } else {
        return(model$beta[,1])
      }
    }
  }

}



#' Estimates beta.
#'
#' Update step of algorithm to obtain new estimation for beta. Note that we call it theta because beta() exists in base R.
#' @inheritParams determine_theta
#' @param eclipz parameter to indicate using real world Eclipzdataset. Defaults to FALSE.
#' @param NN number of individuals
#' @param TT length of time series
#' @param number_of_groups number of groups estimated
#' @param number_of_group_factors number of groupfactors estimated
#' @param number_of_common_factors number of common factors estimated
#' @param number_of_variables number of observable variables
#' @param number_vars_estimated number of variables that are included in the algorithm and have their coefficient estimated. This is usually equal to number_of_variables.
#' @param num_factors_may_vary whether or not the number of groupfactors is constant over all groups or not
#' @return list: 1st element contains matrix (N columns: 1 for each element of the panel data) with estimated theta's.
#' @examples
#' #This function needs several initial parameters to be initialized in order to work on itself.
#' library(RobClustTimeSeries)
#' library(tidyverse)
#' library(robustbase)
#'
#' X = RobClustTimeSeries::X_dgp3
#' Y = RobClustTimeSeries::Y_dgp3
#' #Set estimations for group factors and its loadings, and group membership to the real value
#' lambda_group = RobClustTimeSeries::lambda_group_real_dgp3
#' factor_group = RobClustTimeSeries::factor_group_real_dgp3
#' g = RobClustTimeSeries::g_real_dgp3
#' #There are no common factors to be estimated  -> but needs placeholder
#' lambda = matrix(0,nrow=1,ncol=300)
#' comfactor = matrix(0,nrow=1,ncol=30)
#
#' use_robust = TRUE
#' #Choose how coefficients of the observable are estimated
#' homogeneous_coefficients = FALSE
#' heterogeneous_coefficients_groups = FALSE
#' heterogeneous_coefficients_individuals = TRUE
#' ABDGP1 = FALSE
#' ABintercept = TRUE
#' theta = estimate_theta(NN = 300, TT = 30,
#'   number_of_groups = 3, number_of_group_factors = c(3,3,3), number_of_common_factors = 0,
#'  number_of_variables = 3,number_vars_estimated=3,num_factors_may_vary=FALSE)[[1]]
#' @export
estimate_theta <- function(optimize_kappa = FALSE, eclipz = FALSE,
                                  NN = aantal_N,
                                  TT = aantal_T,
                                  number_of_groups = aantalgroepen,
                                  number_of_group_factors = aantalfactoren_groups,
                                  number_of_common_factors = aantalfactoren_common,
                                  number_of_variables = aantalvars,
                                  number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                                  num_factors_may_vary = aantalfactors_verschillend_per_group) {
  if(number_of_variables > 0) {


    if(homogeneous_coefficients) {
      #X needs to be in the form of (NN*TT x p matrix)
      X_special = X_restructured
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
        Y_special[i,] = Y_special[i,] -
          t(lambda[,i]) %*% comfactor[,]
        if(max(number_of_group_factors, na.rm=T) > 0) {
          #subtract group factorstructure
          Y_special[i,] = Y_special[i,] -
            temp[i,]
        }

      }

      theta = determine_theta("homogeen",X_special, Y_special, TRUE, indices = 1:NN, TT = TT, number_of_variables = number_of_variables)

    }
    if(heterogeneous_coefficients_groups) {
      theta = matrix(NA, nrow = (number_of_variables + 1), ncol = number_of_groups)
      for(group in 1:number_of_groups) {
        #select parts of X and Y of this group
        indices_group = which(g == group)
        if(length(indices_group) == 0) {
          message("Empty groep! This should never be seen.")
          Sys.sleep(15)
        }

        #X needs to be in the form of (NN*TT x p matrix)
        X_special = restructure_X_to_order_slowN_fastT(array(X[indices_group,,], dim = c(length(indices_group), TT, number_of_variables)),
                                                       eclipz,
                                                       number_of_variables = number_of_variables,
                                                       number_vars_estimated = number_vars_estimated)
        #define Y* as Y - FcLc - FgLg:
        Y_special = matrix(Y[indices_group,], nrow = length(indices_group)) %>% unlist


        for(i in 1:length(indices_group)) {
          index = indices_group[i]
          LAMBDAGROUP = as.matrix(lambda_group %>% filter(id %in% index) %>% dplyr::select(-groep,-id))
          Y_special[i,] = Y_special[i,] -
            t(lambda[,index]) %*% comfactor[,] -
            LAMBDAGROUP %*% factor_group[[g[index]]]
        }


        theta[,group] = determine_theta("heterogeen", X_special, Y_special, TRUE, indices = indices_group, TT = TT, number_of_variables = number_of_variables)

      }
    }
    if(heterogeneous_coefficients_individuals) {


      #Vectorized version: Microbenchmark: This takes only 75% of time for (1000,30)-system.
      X_local = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated, eclipz) #makes X smaller when number_vars_estimated < number_of_variables
      if(eclipz) {
        number_of_vars = number_of_variables
      } else {
        number_of_vars = number_vars_estimated
      }
      X_special_list = lapply(1:NN, function(x) restructure_X_to_order_slowN_fastT(matrix(X_local[x,,], ncol = number_of_vars),
                                                                                   eclipz,
                                                                                   number_of_variables = number_of_variables,
                                                                                   number_vars_estimated = number_vars_estimated))

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
        LAMBDAGROUP = LAMBDAGROUP[1:number_of_group_factors[g[i]]]


        Y_special = Y_special -
          t(lambda[,i]) %*% comfactor[,] -
          LAMBDAGROUP %*% factor_group[[g[i]]]


        return(Y_special)
      }

      Y_special_list = lapply(1:NN, function(x) make_Y_special(x) )

      if(optimize_kappa) {
        theta = pmap(list(X_special_list, Y_special_list, 1:NN),  function(x,y,z) determine_theta("heterogeen",x, y, TRUE, indices = z,  TT = TT, number_of_variables = number_of_variables) )
      } else {
        #note that mapply would be about 10% faster
        theta = map2(X_special_list, Y_special_list,  function(x,y) determine_theta("heterogeen",x, y, TRUE, indices = NA,  TT = TT, number_of_variables = number_of_variables) )
        #theta_new = mapply( function(x,y) { determine_theta("heterogeen",x, y, TRUE, indices = NA,  TT = TT, number_of_variables = number_of_variables) }, x = X_special_list, y = y_special_list )

      }
      ########################################################
      #Possible use of future::map2 instead of map2:
      #(is 2 to 3 times faster (microbenchmark) for (300,30)-system)
      #BUT:  ERRORS:
      #   Error: 'map2' is not an exported object from 'namespace:future'
      #theta = future::map2(X_special_list, Y_special_list,  function(x,y) determine_theta("heterogeen",x, y, TRUE, indices = NA,  TT = TT, number_of_variables = number_of_variables) )
      ########################################################

      theta = matrix(unlist(theta),ncol = NN)

    }


    return(list(theta))
  } else {
    return(list(NA))
  }
}



#' Calculates W = Y - X*THETA
#'
#' @param theta theta
#' @param g vector with group membership
#' @inheritParams estimate_theta
#' @return NxT matrix
#' @export
calculate_W <- function(theta, g ,
                        NN = aantal_N,
                        TT = aantal_T,
                        number_of_variables = aantalvars,
                        number_vars_estimated = SCHATTEN_MET_AANTALVARS, eclipz = FALSE) {
  W = matrix(0, nrow = NN, ncol = TT) #T x N-matrix in original paper , but I rather define as NxT

  #if number_vars_estimated < number_of_variables the obsoleterows in theta were already erased -> do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated, eclipz)

  if(number_of_variables > 0) {
    if(homogeneous_coefficients) {
      #non-dependence on g -> take first column
      for(i in 1:NN) W[i,] = Y[i,] - t(cbind(1,X[i,,]) %*% as.matrix(theta[,1]))
    }
    if(heterogeneous_coefficients_groups) {
      for(i in 1:NN) W[i,] = Y[i,] - t(cbind(1,X[i,,]) %*% as.matrix(theta[,g[i]]))
    }
    if(heterogeneous_coefficients_individuals) {
      for(i in 1:NN) {
        if(aantal_N_fulldata == 3112 & eclipz) { #this combination has issues with format of Y; not clear why
          W[i,] = matrix(unlist(Y[i,]), nrow = 1) - t(cbind(1,X[i,,]) %*% as.matrix(theta[,i]))
        } else {
          W[i,] = Y[i,] - t(cbind(1,X[i,,]) %*% as.matrix(theta[,i]))
        }
      }
    }
  } else {
    for(i in 1:NN) W = as.matrix(Y)
  }



  if(anyNA(W)) { #NA's veroorzaakt because an individual has zero weight


    #handle NA's: which currently means: remove the rows with NA's
    W = handleNA(W)[[1]]
  }
  return(W)
}

#' Calculates Z = Y - X*THETA - LgFg, to use in estimate of common factorstructure
#' @inheritParams calculate_W
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param initialise boolean
#' @inheritParams estimate_theta
#' @export
calculate_Z_common <- function(theta, g, lgfg_list, initialise = FALSE,
                               NN = aantal_N,
                               TT = aantal_T,
                               number_of_variables = aantalvars,
                               number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                               number_of_group_factors = aantalfactoren_groups, eclipz = FALSE) {
  Z = matrix(0, nrow = NN, ncol = TT) #T x N-matrix in paper , maar ik definieer liever als NxT

  #if number_vars_estimated < number_of_variables the obsoleterows in theta were already erased -> do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated, eclipz)

  for(i in 1:NN) {
    y = Y[i,] %>% as.numeric
    if(do_we_estimate_group_factors(number_of_group_factors)) { #when there are group factors to be estimated
      if(exists("TESTPERTMM") & class(lgfg_list) != "list") {
        #then lgfg_list is not a list, but a df of dimension NxT
        LF_GROUP = lgfg_list[i,]
      } else {
        LF_GROUP = lgfg_list[[g[i]]][i,]
      }
    } else {
      LF_GROUP = 0
    }
    if(number_of_variables > 0) {
      if(homogeneous_coefficients | heterogeneous_coefficients_groups) THETA = as.matrix(theta[,g[i]])
      if(heterogeneous_coefficients_individuals) THETA = as.matrix(theta[,i])
      Z[i,] = y - t(cbind(1,X[i,,]) %*% THETA) -
        LF_GROUP
    } else {
      Z[i,] = y -
        LF_GROUP
    }
  }


  if(anyNA(Z)) {
    #print("NA's in Z (calculate_Z_common) -> wegwerken, want eigenvectoren (schatten van F) kunnen anders niet berekend worden")

    #handle NA's: which currently means: remove the rows with NA's
    Z = handleNA(Z)[[1]]
  }
  return(Z)
}

#' Calculates Z = Y - X*THETA - LF), to use in estimate of groupfactorstructure.
#' @inheritParams calculate_W
#' @inheritParams estimate_theta
#' @param lambda common factor loadings
#' @param comfactor common factors
#' @param group the number of the group
#' @param initialise boolean
#' @export
calculate_Z_group <- function(theta, g, lambda, comfactor, group, initialise,
                              TT = aantal_T,
                              number_of_variables = aantalvars,
                              number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                              eclipz = FALSE,
                              number_of_common_factors = aantalfactoren_common) {
  #print("calc Z_group")
  indices_group = which(g == group)
  if(length(indices_group) == 0) {
    print("empty group (calculate_Z_group())")
    Sys.sleep(1)
  }

  Z = matrix(0, nrow = length(indices_group), ncol = TT) #Nj x T matrix
  #if number_vars_estimated < number_of_variables the obsoleterows in theta were already erased -> do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated, eclipz)

  for(i in 1:length(indices_group)) { #loop over aantal elementen in groep
    index = indices_group[i]

    #define XT
    if(number_of_variables > 0) {
      if(homogeneous_coefficients | heterogeneous_coefficients_groups)  THETA = as.matrix(theta[,g[index]])
      if(heterogeneous_coefficients_individuals) THETA = as.matrix(theta[,index])

      XT = t(cbind(1,X[index,,]) %*% THETA)
    } else {
      XT = rep(0,TT)
    }

    #calculate Z
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
    #print("NA's in Z (calculate_Z_group) -> wegwerken, want eigenvectoren (schatten van F) kunnen anders niet berekend worden")

    #handle NA's: which currently means: remove the rows with NA's
    Z = handleNA(Z)[[1]]
  }

  return(Z)
}

#' Solves a very specific issue with MacroPCA.
#'
#' MacroPCA crashes Rstudio with certain dimensions of the input. Solve this by doubling every row. No information is added by this, so no influence on end result,
#' but crashes are evaded.
#' @param object input
evade_crashes_macropca <- function(object) {
  #--------------------MacroPCA seems to make Rstudio crash when the dimension of object = (27,193) -----------------
  size_with_crashes1 = c( 27,  27,  25,  43,  66,  24)#,  44,  41)
  size_with_crashes2 = c(193, 310, 600, 560, 387, 374)#, 201, 201)
  changed_objectsize = FALSE
  message("Test size of object (for MacroPCA)")
  for(i in 1:length(size_with_crashes1)) {
    if(changed_objectsize == FALSE & dim(object)[1] == size_with_crashes1[i] & dim(object)[2] == size_with_crashes2[i]) {
      message("Rstudio would crash -> double amount of rows")
      print(i)
      object = rbind(object, object)
      message(dim(object))
      changed_objectsize = TRUE
    }
  }
  message("Test size: done")
  return(object)
}

#' Helpfunction in robustpca().
#'
#' It handles possible thrown errors in MacroPCA.
#' @param object input
#' @param temp this is the result of the trycatch block of using macropca on object
#' @param KMAX parameter kmax in MacroPCA
#' @param number_eigenvectors number of principal components that are needed
handle_macropca_errors <- function(object,temp,KMAX,number_eigenvectors) {

  if("error" %in% class(temp)) {
    message("*******************************")
    message(paste("MacroPCA with .... fails, -> use different amount of eigenvectors. Start with a couple more and decrease 1 by 1."))
    Sys.sleep(2)
    teller = 0
    while("error" %in% class(temp)) {
      teller = teller + 1
      print(paste("teller:",teller))


      temp =  tryCatch(
          cellWise::MacroPCA(object, k = max(12, number_eigenvectors) - teller, MacroPCApars = list(kmax=KMAX))$loadings[,1:number_eigenvectors],
          error = function(e) { message(e); return(e) }
      )

      if(teller >= 10) {
        message("-------------infinite loop (MacroPCA does not work with any k) -> use (classical) eigen(): has to be squared matrix -> take covariance matrix of object-------------")
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
#' KMAX: Different values for kmax give different factors, but the  product lambda*factor stays constant. Note that
#'  this number needs to be big enough, otherwise eigen() will be used.
#' Variation in k does give different results for lambda*factor
#'
#' NOTE: this function may (not certain) crash (in case of MACROPCA) when ncol(object) >> nrow(object)
#'   Actually it crashes with specific values of dim(object). For example when dim(object) = c(193,27).
#'   This is solved with evade_crashes_macropca()
#' @param object input
#' @param number_eigenvectors number of eigenvectors to extract
#' @param KMAX The maximal number of principal components to compute. This is a paramater in cellWise::MacroPCA()
#' @export
robustpca <- function(object, number_eigenvectors, KMAX = 20) {

  print(paste("*************************************************robust PCA with:",number_eigenvectors))
  print(dim(object))

  ######################
  # MacroPCA
  ######################

  object = evade_crashes_macropca(object)

    print(number_eigenvectors)
    if(number_eigenvectors > KMAX) {
      #Note that when k > kmax, k gets the value of kmax.
      message("MacroPCA is (through kmax) limited to 20 factors.")
      Sys.sleep(15)
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
    temp =  tryCatch(
      cellWise::MacroPCA(object, k = max(macropca_kmax, number_eigenvectors), MacroPCApars = list(kmax=KMAX))$loadings[,1:number_eigenvectors],
      error = function(e) { message(e); return(e) }
    )

    temp = handle_macropca_errors(object, temp, KMAX,number_eigenvectors)
    return(list(temp,NA))


}

#' Function to determine a robust covariance matrix, using MRCD.
#'
#' @param object input
#' @param ALPHA alpha of MRCD
#' @export
get_robust_covmatrix <- function(object, ALPHA = ALPHA_COVMRCD) {
  MRCD = CovMrcd(object, alpha = ALPHA)
  return(MRCD$cov)
}

#' Helpfunction in estimate_factor() for the robust approach.
#'
#' It has some options used for testing.


#'Estimates common factor F
#'
#'The estimator for F, see Anderson (1984), is equal to the first r eigenvectors (multiplied by sqrt(T) due to the restriction F'F/T = I)
#'associated with first r largest eigenvalues of the matrix WW' (which is TxT)
#'W is called Z in Ando/Bai (2017)
#' @inheritParams calculate_W
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param initialise boolean
#' @inheritParams estimate_theta
#' @return r x T matrix
#' @export
estimate_factor <- function(theta, g, lgfg_list, initialise = FALSE,
                            NN = aantal_N,
                            TT = aantal_T,
                            number_of_common_factors = aantalfactoren_common) {


  #initialisation: grouped structure NA
  if(initialise) {
    W = calculate_W(theta, g)
  } else {
    W = calculate_Z_common(theta, g, lgfg_list) #this is a "rows_without_NA x T - matrix"
  }

  #Define the object on which (robust or classical) PCA will be performed
  if(use_robust) {
    temp = prepare_for_robpca(W)
  } else {
    #Classical case
    #calculate 1/NT * W'W
    temp = t(W)%*%W / (NN * TT) #TxT matrix (this division by NT does not actually matter for the calculation of the eigenvectors)

  }



  #first r eigenvectors * sqrt(T)

  #take number_of_common_factors eigenvectors
  if(use_robust) {
    message(str_c("Estimate ",number_of_common_factors," common factors"))
    temp2 = robustpca(temp, number_of_common_factors) #robust pca
    schatterF = t(sqrt(TT) * temp2[[1]])
    scores = temp2[[2]]
    rm(temp2)
  } else {
    schatterF = t(sqrt(TT) * eigen(temp)$vectors[,1:number_of_common_factors])
  }
  if(number_of_common_factors == 0) {
    #then schatterF is TT, which is still an apropriate size to use in this case -> just set to zero
    schatterF = schatterF - schatterF
  }


  return(schatterF)
}

#' Helpfunction: prepares object to perform robust PCA on.
#'
#' It contains options to use the classical or robust covmatrix or no covariance matrix at all.
#' Is currently dependent on global variables.
#' @param object this is the object of which we may take the covariance matrix and then to perform robust PCA on
#' @param NN N
#' @param TT T
#' @param option 1 (robust covmatrix), 2 (classical covmatrix), 3 (no covmatrix)
prepare_for_robpca <- function(object, NN = aantal_N, TT = aantal_T, option = 3) {
  if(option == 1) temp = get_robust_covmatrix(object)

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
#' @inheritParams estimate_theta
#' @export
estimate_factor_group <- function(theta, g, lambda, comfactor, initialise = FALSE,
                                  NN = aantal_N,
                                  TT = aantal_T,
                                  number_of_groups = aantalgroepen,
                                  number_of_group_factors = aantalfactoren_groups) {
  schatterF = list()
  scores = list()

  if(eclipz & expert_based_initial_factors & step == 0) {
    schatterF = define_expert_based_initial_factors(number_of_groups)
  } else {
    for(group in 1:number_of_groups) {
      print(group)
      if(number_of_group_factors[group] > 0) {
        if(TT < number_of_group_factors[group]) {
          message("-> too many factor compared to TT")
        }

        Wj = calculate_Z_group(theta, g, lambda, comfactor, group, initialise)



        #Take a limit of 10 individuals per group (since get_robust_covmatrix() does not work at small sizes)
        if(use_robust & nrow(Wj) > 10) {


          temp = prepare_for_robpca(Wj)

          temp2 = robustpca(temp, number_of_group_factors[group])
          schatterF[[group]] = t(sqrt(TT) * temp2[[1]]) #robust pca
          scores[[group]] = temp2[[2]]

          rm(temp2)
        } else {
          if(use_robust) print("Dit is een zeer kleine groep -> CovMrcd werkt niet; robpca werkt ook niet -> gebruik niet-robuuste versie")
          temp = t(Wj)%*%Wj #delen door NT maakt geen verschil voor de eigenvectors
          scores[[group]] = NA
          schatterF[[group]] = t(sqrt(TT) * eigen(temp)$vectors[,1:number_of_group_factors[group]])
        }



      } else {
        schatterF[[group]] = matrix(0, 1, TT)
      }
    }
  }


  return(schatterF)
}

#' calculates factor loadings of common factors
#' @inheritParams calculate_W
#' @param comfactor common factors
#' @param lgfg_list This is a list (length number of groups) containing FgLg for every group.
#' @param initialise boolean
#' @inheritParams estimate_theta
#' @export
calculate_lambda <- function(theta, comfactor, g, lgfg_list, initialise = FALSE,
                             NN = aantal_N, TT = aantal_T,
                             number_of_common_factors = aantalfactoren_common) {


  if(initialise) {
    W = calculate_W(theta, g)
  }  else {
    W = calculate_Z_common(theta, g, lgfg_list)
  }


  if(use_robust) {
    lambda = return_robust_lambdaobject(W, NA, type = 2, FACTOR = comfactor, number_of_common_factors = nrow(comfactor))
  } else {
    lambda = t(W %*% t(comfactor) / TT)
  }

  if(number_of_common_factors == 0 & !initialise) {
    #then schatterF is of size 1 x TT, which is still an apropriate size to use in this case -> just set to zero
    lambda = lambda - lambda
  }

  #for income-series with NA's, the lambda's cannot be determined -> set to zero ...
  temp = data.frame(matrix(NA, nrow = max(number_of_common_factors,1), ncol = NN)) #new dataframe
  temp[,rows_without_NA] = lambda #put calculated lambda's in df
  if(sum(rows_with_NA,na.rm=T) > 0)  temp[,rows_with_NA] = 0 #set the rest to zero
  lambda = as.matrix(temp)


  return(lambda)
}

#'calculates factor loadings of groupfactors
#'
#'returns object which includes group and id of the individuals
#' @param theta theta
#' @param factor_group group factors
#' @param lambda loadings
#' @param comfactor common factors
#' @param initialise boolean
#' @param UPDATE1 option to indicate the number of groupfactors is updated during the algorithm; defults to FALSE
#' @param UPDATE2 option to indicate the number of common factors is updated during the algorithm; defults to FALSE
#' @inheritParams estimate_theta
#' @inheritParams calculate_W
#' @export
calculate_lambda_group <- function(theta, factor_group, g, lambda, comfactor, initialise = FALSE, UPDATE1 = update1, UPDATE2 = update2,
                                   TT = aantal_T,
                                   number_of_groups = aantalgroepen,
                                   number_of_group_factors = aantalfactoren_groups) {


  lambda_local = list()

  for(group in 1:number_of_groups) {
    if(number_of_group_factors[group] > 0) {
      Wj = calculate_Z_group(theta, g, lambda, comfactor, group, initialise)
      #robust things:
      if(use_robust) {


        #Since in the classical approach each lambda is the mean of a set of products of F and Y (Z),
        #   (lambda_N1 = (F_T1 * Y_N1T1 + F_T2 * Y_N1T2 + ...) / TT)
        #   we can replace this mean by an M-estimator in the robust approach
        if(UPDATE1) {
          lambda_local[[group]] = return_robust_lambdaobject(Wj, group, type = 3, FACTOR_GROUP = factor_group, number_of_group_factors = unlist(lapply(factor_group, nrow)))
        } else {
          lambda_local[[group]] = return_robust_lambdaobject(Wj, group, type = 3)
        }


      } else {

        FGG = factor_group[[group]]

        if(dim(t(FGG))[1] == 1) {
          FGG = matrix(FGG, nrow = 1)

        }
        lambda_local[[group]] = Wj %*% t(FGG) / TT
      }

    } else {
      lambda_local[[group]] = matrix(0,
                                     length(rows_without_NA[(rows_without_NA %in% which(g == group))]),
                                     1)
    }
  }

  #for income-series with NA's, the lambda's cannot be determined -> set to zero with this function
  add_lambdas_from_NArows <- function(DF, groep) {
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
  lambda_local2 = lambda_local2 %>% mutate(groep = 1, id = rows_without_NA[(rows_without_NA %in% which(g == 1))])


  #fill up with "rows_with_NA" of this groep:
  #lambda_group can not be defined for these rows, so put to 0
  lambda_local2 = lambda_local2 %>% rbind(add_lambdas_from_NArows(lambda_local2,1)) %>% arrange(id)

  if(number_of_groups > 1) {
    for(i in 2:number_of_groups) {
      df_to_add = data.frame(lambda_local[[i]])
      names(df_to_add) = str_c("X",1:ncol(df_to_add))

      df_to_add = df_to_add %>% mutate(groep = i, id = rows_without_NA[(rows_without_NA %in% which(g == i))])

      df_to_add = df_to_add %>% rbind(add_lambdas_from_NArows(df_to_add,i)) %>% arrange(id)

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

  return(lambda_local2)
}


#' Function which is used to have a df ("grid") with data (individualindex, timeindex, XT and LF) available.
#'
#' It is used in update_g().
#' @param grid dataframe containing values for X*theta and LF (product of common factor and its loadings)
#' @param theta theta
#' @param lambda loadings of the common factors
#' @param comfactor common factors
#' @inheritParams estimate_theta
#' @export
grid_add_variables <- function(grid, theta, lambda, comfactor,
                               NN = aantal_N,
                               TT = aantal_T,
                               number_of_variables = aantalvars,
                               number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                               number_of_groups = aantalgroepen) {
  stopifnot(number_of_variables > 0 & number_of_variables < 9) #code exists up to 8 groups
  if(number_of_variables > 0) {
    #for homogeneous theta (1 -> 4 at this moment), we only need 1 column as all columns are the same
    if(homogeneous_coefficients) {
      theta_used = as.matrix(theta[,1])
      #calculate matrix multiplications outside the OF-functions to speed up:
      grid$XTHETA =  (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta_used))
    } else {
      if(heterogeneous_coefficients_groups) {
        #for each group: define XT
        if(number_of_groups > 0) grid$XTHETA1 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,1]))
        if(number_of_groups > 1) grid$XTHETA2 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,2]))
        if(number_of_groups > 2) grid$XTHETA3 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,3]))
        if(number_of_groups > 3) grid$XTHETA4 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,4]))
        if(number_of_groups > 4) grid$XTHETA5 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,5]))
        if(number_of_groups > 5) grid$XTHETA6 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,6]))
        if(number_of_groups > 6) grid$XTHETA7 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,7]))
        if(number_of_groups > 7) grid$XTHETA8 = (apply(grid, 1, function(x) c(1, X[x[1],x[2],]) %*% theta[,8]))
        # if(number_of_groups > 8) {
        #   message("add code (grid_add_variables")
        #   Sys.sleep(900)
        # }
      }
      if(heterogeneous_coefficients_individuals) {
        grid$XTHETA = c(calculate_XT_estimated(NN = NN, TT = TT, number_of_variables = number_of_variables, number_vars_estimated = number_vars_estimated))
      }
    }
  } else {
    if(homogeneous_coefficients | heterogeneous_coefficients_individuals) {
      grid$XTHETA = 0
    }
    if(heterogeneous_coefficients_groups) {
      if(number_of_groups > 0) grid$XTHETA1 = 0
      if(number_of_groups > 1) grid$XTHETA2 = 0
      if(number_of_groups > 2) grid$XTHETA3 = 0
      if(number_of_groups > 3) grid$XTHETA4 = 0
      if(number_of_groups > 4) grid$XTHETA5 = 0
      if(number_of_groups > 5) grid$XTHETA6 = 0
      if(number_of_groups > 6) grid$XTHETA7 = 0
      if(number_of_groups > 7) grid$XTHETA8 = 0
      # if(number_of_groups > 8) {
      #   print("add code-")
      #   Sys.sleep(900)
      # }

    }

  }
  LF = t(lambda) %*% comfactor
  grid$LF = c(LF)

  return(grid)
}













#' Function with as input a dataframe. (this will be "Y" or "to_divide") It filters out rows with NA.
#'
#' @param df input
#' @export
handleNA <- function(df) {
  rows_with_NA = which(apply(df,1,function(x) sum(is.na(x)) != 0)) #rownumbers of rows containing NA
  #print(paste("There are ",length(rows_with_NA),"of",nrow(df), "rows with NA's."))

  return(list(na.omit(df), rows_with_NA))
}

#' Removes NA's in LG (in function calculate_virtual_factor_and_lambda_group() )
#' @param df input
#' @export
handleNA_LG <- function(df) {
  result = (as.matrix(data.frame(df) %>% mutate_if(is.numeric , replace_na, replace = 0)))
  return(result)
}



#' Helpunction to shorten code: use a and b
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
#' Helpunction to shorten code: use a and b
#' @param number_of_group_factors number of group factors
do_we_estimate_group_factors <- function(number_of_group_factors) {
  if(mean(number_of_group_factors, na.rm = T) == 0) {
    b = 0
  } else {
    b = 1
  }
  return(b)
}


#' Calculates the error term Y - X*THETA - LF - LgFg.
#' @return NxT matrix
#' @param no_common_factorstructure if there is a common factorstructure being estimated
#' @param no_group_factorstructure if there is a group factorstructure being estimated
#' @inheritParams grid_add_variables
#' @inheritParams estimate_theta
#' @export
calculate_error_term <- function(no_common_factorstructure = FALSE, no_group_factorstructure = FALSE,
                                 NN = aantal_N,
                                 TT = aantal_T,
                                 number_of_variables = aantalvars,
                                 number_of_groups = aantalgroepen,
                                 number_of_group_factors = aantalfactoren_groups,
                                 number_of_common_factors = aantalfactoren_common) {
  stopifnot(table(lambda_group$groep) == table(g)) #in this case lambda_groep would need to be updated

  u = matrix(NA,nrow = NN, ncol = TT)
  e = matrix(NA,nrow = NN, ncol = TT)
  lf = t(lambda) %*% comfactor
  if(number_of_variables > 0) {
    if(homogeneous_coefficients | heterogeneous_coefficients_groups) {
      xt = sapply(1:NN,
                  function(y) sapply(1:TT, function(x) c(1, X[y,x,]) %*% theta[,g[y]]))
    }
    if(heterogeneous_coefficients_individuals) {
      xt = t(calculate_XT_estimated(NN = NN, TT = TT, number_of_variables = number_of_variables)) #TxN matrix
    }
  } else {
    xt = matrix(0, nrow = TT, ncol = NN) #TxN
  }
  lf_group = list()
  group_membership = list()
  for(k in 1:number_of_groups) {
    LGclean = as.matrix(lambda_group %>% arrange(id) %>%
                          filter(groep == k) %>%
                          dplyr::select(starts_with("X")))
    lf_group[[k]] = LGclean[,1:number_of_group_factors[k]] %*% factor_group[[k]]
    group_membership[[k]] = data.frame(as.matrix(lambda_group %>% arrange(id) %>%
                                                   filter(groep == k) %>%
                                                   dplyr::select(groep, id)))
  }

  for(i in 1:NN) {
    lf_group_i = lf_group[[g[i]]]
    index = which(group_membership[[g[i]]]$id == i)
    for(t in 1:TT) {
      u[i,t] = Y[i,t] - xt[t,i]

      a = do_we_estimate_common_factors(number_of_common_factors)
      b = do_we_estimate_group_factors(number_of_group_factors)

      part2 = a * lf[i,t]
      part3 = b * lf_group_i[index,t]
      if(no_common_factorstructure) part2 = 0
      if(no_group_factorstructure) part3 = 0
      e[i,t] = u[i,t] - part2 - part3
    }
  }
  return(e)
}


#' Calculates the error term Y - X*THETA - LF - LgFg.
#'
#' Gives same value as calculate_error_term().
#' @inheritParams estimate_theta
#' @return NxT matrix
#' @export
calculate_error_term_individuals <- function(NN = aantal_N,
                                             TT = aantal_T,
                                             number_of_variables = aantalvars,
                                             number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                                             number_of_groups = aantalgroepen,
                                             number_of_group_factors = aantalfactoren_groups,
                                             number_of_common_factors = aantalfactoren_common,
                                             eclipz = FALSE) {
  u = matrix(NA,nrow = NN, ncol = TT)
  e = matrix(NA,nrow = NN, ncol = TT)
  lf = t(lambda) %*% comfactor

  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated, eclipz)
  if(number_of_variables > 0) {
    if(homogeneous_coefficients | heterogeneous_coefficients_groups) {
      xt = sapply(1:NN,
                  function(y) sapply(1:TT, function(x) c(1, X[y,x,]) %*% theta[,g[y]]))
    }
    if(heterogeneous_coefficients_individuals) {
      xt = sapply(1:NN,
                  function(y) sapply(1:TT, function(x) c(1, X[y,x,]) %*% theta[,y]))
    }

  } else {
    xt = matrix(0, nrow = TT, ncol = NN) #TxN
  }

  lf_group = list()
  group_membership = list()
  for(k in 1:number_of_groups) {
    LGclean = as.matrix(lambda_group %>% arrange(id) %>%
                          filter(groep == k) %>%
                          dplyr::select(starts_with("X")))
    lf_group[[k]] = LGclean[,1:number_of_group_factors[k]] %*% factor_group[[k]]
    group_membership[[k]] = data.frame(as.matrix(lambda_group %>% arrange(id) %>%
                                                   filter(groep == k) %>%
                                                   dplyr::select(groep, id)))

  }



  a = do_we_estimate_common_factors(number_of_common_factors)
  b = do_we_estimate_group_factors(number_of_group_factors)
  for(i in 1:NN) {
    lf_group_i = lf_group[[g[i]]]
    index = which(group_membership[[g[i]]]$id == i)
    u[i,] = (Y[i,] - xt[,i]) %>% as.numeric
    e[i,] = u[i,] - a * lf[i,] - b * lf_group_i[index,]


  }
  return(e)
}


#' Calculates sum of squared errors, divided by NT
#'
#' @param e matrix with error terms
#' @param NN N
#' @param TT T
#' @export
calculate_sigma2 <- function(e, NN = aantal_N, TT = aantal_T) {
  if(anyNA(e)) {
    message("There are NA's in e (calculate_sigma2(). This should not happen. ")
    Sys.sleep(3600)
  }
  sigma2 = sum(e*e)/(NN * TT)


  #check if equal to calculate_PIC_term1 (in case of non robust versions, because calculate_sigma2 gives only non-robust results)
  # if(!use_robust) {
  #   A = (abs(sigma2 - calculate_PIC_term1(e)))
  # }
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
#' @export
calculate_PIC_term1 <- function(e2) {
  if(use_robust) {
    #This replaces sum(z^2) by sum(rho(z)) with rho the bisquare function

    #RHO_E = Mpsi(e2, cc = 4.685, psi = "bisquare", deriv = -1)
    #NOTE: dividing by mad(e2) would intuitively make sense, but it scales the errors of configurations with bad estimations to the results of good estimations -> leads to wrong estimation of S & k
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
calculate_PIC_term4 <- function(temp, term4, Nj, NN = aantal_N) {


  if(exists("ALTERNATIVE_PIC_2")) {
    term4 = term4 + (temp * Nj/NN)
  } else {

    if(exists("ALTERNATIVE_PIC")) { #PIC of code AB2017
      if(ALTERNATIVE_PIC == TRUE) {
        message("We use an Alternative PIC")
        term4 = term4 + (temp * Nj/NN) #=weigh with "P/PP"
      } else {
        term4 = term4 + temp
      }
    } else {
      #=Standard PIC (from Ando/Bai2017-paper)
      term4 = term4 + temp
    }
  }
  return(term4)
}

#' function to test alternative PIC's
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
    term2 = term2 / aantal_T
    term3 = C * term3
    term4 = C * term4
  }
  return(list(term2,term3,term4))
}

#' Function to determine PIC (panel information criterium)
#'
#' This depends on kappa1 -> kappaN, through p (=number of nonzero elements of theta).
#' The parameter 'sigma2' is the non-robust sigma2. As it is only used in term 2 to 4, it does not actually matter what its value is (needs to be > 0).
#' Could be set to 1 as well.
#' @param C determines relative contribution of the penamlty terms
#' @inheritParams estimate_theta
#' @param e2 matrix with error terms
#' @param sigma2 scalar: sum of squared error terms, scaled by NT
#' @export
calculate_PIC <- function(C, number_of_common_factors, number_of_group_factors, e2, sigma2,
                          NN = aantal_N,
                          TT = aantal_T,
                          number_of_variables = aantalvars,
                          number_of_groups = aantalgroepen) {


  term1 = calculate_PIC_term1(e2)

  if(number_of_variables > 0) {
    if(homogeneous_coefficients | heterogeneous_coefficients_groups) {
      #p is number of nonzero elements of theta; we need the sum of all p's (inclusief intercept)
      p_sum = sum(sapply(1:NN, function(x) sum(theta[,g[x]] != 0)))
    } else{
      p_sum = sum(sapply(1:NN, function(x) sum(theta[,x] != 0)))
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

















#' Calculates the product of X*theta (the real value of this).
#' @inheritParams estimate_theta
#' @export
calculate_XT_real <- function(NN = aantal_N, TT = aantal_T, number_of_variables = aantalvars) {
  if(homogeneous_coefficients) { #relevant for DGP5 (BRamatiCroux)
    XT_real = theta_real[1,] + X[,,1] * theta_real[1,]
    stopifnot(number_of_variables == 1)
  }
  if(heterogeneous_coefficients_groups) {
    if(number_of_variables > 0 & !eclipz & !ABDGP1) {
      message("Note that this is not the default option. Extend code (calculate_XT_real()) if needed.")
    } else {
      XT_real = NA
    }
  }
  if(heterogeneous_coefficients_individuals) {
    if(number_of_variables > 0 & !eclipz & !ABDGP1) {
      XT_real = (t(sapply(1:NN,
                          function(y) sapply(1:TT, function(x) c(1, X[y,x,]) %*% theta_real[,g_real][,y]))))
    } else {
      XT_real = NA
    }
  }

  return(XT_real)
}

#' When we run the algorithm with a different number of observable variables then the number we actually have, we need to reformat X.
#'
#' @inheritParams create_theta_real
#' @inheritParams generate_Y
#' @inheritParams initialise_theta
adapt_X_estimating_less_variables <- function(number_of_variables,
                                              number_vars_estimated,
                                              eclipz = FALSE) {
  #if number_vars_estimated < number_of_variables, then the obsolete rows in theta are already erased -> do the same in X
  if(!eclipz & number_vars_estimated < number_of_variables) {
    X = X[,,1:number_vars_estimated]
  }
  return(X)
}

#' Calculates (the estimated value of) the matrix X*theta.
#'
#' @inheritParams estimate_theta
#' @export
calculate_XT_estimated <- function(NN = aantal_N, TT = aantal_T,
                                   number_of_variables = aantalvars,
                                   number_vars_estimated = SCHATTEN_MET_AANTALVARS,
                                   eclipz = FALSE) {
  #if number_vars_estimated < number_of_variables, then the obsolete rows in theta are already erased -> now do the same in X
  X = adapt_X_estimating_less_variables(number_of_variables, number_vars_estimated, eclipz)

  if(number_of_variables > 0) {
    if(homogeneous_coefficients) { #only designed for DGP 5 at the moment
      XT_geschat = theta[1] + X[,,1] * theta[1]
      stopifnot(number_of_variables == 1)
    }
    if(heterogeneous_coefficients_groups) {
      if(number_vars_estimated > 0) {
        XT_geschat = t(sapply(1:NN,
                 function(x) matrix(cbind(1, X[x,,]) %*% theta[,g[x]], nrow = 1)))
      } else {
        XT_geschat = NA
      }
    }
    if(heterogeneous_coefficients_individuals) {
      if(number_vars_estimated > 0) {
        XT_geschat = t(sapply(1:NN,
                 function(x) matrix(cbind(1, X[x,,]) %*% theta[,x], nrow = 1)))

      } else {
        XT_geschat = NA
      }
    }
  } else {
    XT_geschat = matrix(0, nrow = NN, ncol = TT)
  }
  return(XT_geschat)
}


#' Calculate real groupfactorstructure.
#' @return list with NjxT matrices
#' @param LAMBDA_GROUP_REAL real group factor loadings
#' @param FACTOR_GROUP_REAL real group factors
#' @param number_of_groups_real real number of groups
#' @param number_of_group_factors_real real number of group factors for each group
#' @param number_of_common_factors_real real number of common factors
#' @inheritParams estimate_theta
#' @inheritParams generate_Y
#' @param ABDGP1_local gives information about which DGP we use; TRUE of FALSE
#' @export
calculate_FL_group_real <- function(NN = aantal_N, TT = aantal_T,
                                    LAMBDA_GROUP_REAL = lambda_group_real, FACTOR_GROUP_REAL = factor_group_real,
                                    number_of_groups_real = aantalgroepen_real,
                                    number_of_group_factors_real = aantalfactoren_groups_real,
                                    number_of_common_factors_real = aantalfactoren_common_real,
                                    num_factors_may_vary = aantalfactors_verschillend_per_group,
                                    eclipz = FALSE,
                                    ABDGP1_local = ABDGP1) {

  if(exists("DGP_Bramati_Croux")) { #when using DGP5 theer are no factors -> return NA
    return(NA)
  }

  if(!eclipz & !ABDGP1_local) {
    temp = calculate_lgfg(LAMBDA_GROUP_REAL, FACTOR_GROUP_REAL,
                          number_of_groups_real,
                          number_of_group_factors_real,
                          number_of_common_factors_real,
                          num_factors_constant)
    FL_GROUP_real = lapply(1:NN, function(x) temp[[g_real[x]]][x,]) %>% unlist %>% matrix(nrow = TT) %>% t
  } else {
    FL_GROUP_real = NA
  }
  return(FL_GROUP_real)
}


#' Returns the estimated groupfactorstructure.
#' @param LAMBDA_GROUP loadings of group factors
#' @param FACTOR_GROUP group factors
#' @inheritParams calculate_FL_group_real
#' @inheritParams generate_Y
#' @inheritParams update_g
#' @return list with NjxT matrices
#' @export
calculate_FL_group_estimated <- function(LAMBDA_GROUP = lambda_group, FACTOR_GROUP = factor_group,
                                         NN = aantal_N,
                                         TT = aantal_T,
                                         number_of_groups = aantalgroepen,
                                         number_of_group_factors = aantalfactoren_groups,
                                         number_of_common_factors = aantalfactoren_common,
                                         num_factors_may_vary = aantalfactors_verschillend_per_group) {
  temp = calculate_lgfg(LAMBDA_GROUP, FACTOR_GROUP, number_of_groups, number_of_group_factors, number_of_common_factors, num_factors_may_vary)
  if(!is.na(temp)) {
    FL_group_est = lapply(1:NN, function(x) temp[[g[x]]][x,]) %>% unlist %>% matrix(nrow = TT) %>% t
  } else {
    FL_group_est = NA
  }
  return(FL_group_est)
}

#' Function to calculate MSE of theta.
#'
#' (calculated with formula of AndoBai2017: supplementary: p19).
#' Notes:
#' Heterogeneous_coefficients_individuals is the default case.
#' For DGP 1 & 2: When the number of variables in X is not equal to the standard of 3 it returns NA.
#' @param theta estimated values of theta
#' @param theta_real real values of theta
#' @param without_intercept boolean to remove the intercept in the calculation
#' @param ABintercept boolean to remove the de facto intercept in the calculation (this is only for DGP 1 & 2, due to how they are defined)
#' @inheritParams estimate_theta
#' @inheritParams calculate_FL_group_real
#' @export
calculate_mse_theta <- function(theta, theta_real, without_intercept = FALSE, ABintercept = FALSE,
                                NN = aantal_N,
                                number_of_variables = aantalvars,
                                ABDGP1_local = ABDGP1) {

  if(homogeneous_coefficients) { #relevant for DGP05 (Bramati-Croux)
    mse = mean( (theta - theta_real)^2 )
    if(without_intercept) mse = mean( (theta[-1,] - theta_real[-1,])^2 )
    return(mse)
  }

  #Heterogeneous_coefficients_groups is not the default case.
  #As mse_heterogeneous_groups() does not encompass every single option anymore we return NA.
  if(heterogeneous_coefficients_groups) {
    #pre_mse = mse_heterogeneous_groups()
    message("Heterogeneous_coefficients_groups is not the default case. Return NA.")
    return(NA)
  }

  if(heterogeneous_coefficients_individuals) { #Default case; In case of DGP 1 & 2 it returns NA when number of variables is not equal to 3.
    if(without_intercept) {
      if(ABintercept) { #->calculate MSE without real and without de facto intercept
        theta = matrix(theta[c(-1,-2),],ncol=NN)
        theta_real = theta_real[c(-1,-2),]
      } else { #calculate MSE without real intercept
        theta = theta[-1,]
        theta_real = theta_real[-1,]
      }
    }


    pre_mse = rep(NA,NN)
    for(i in 1:NN) {
      #Case of DGP01 and DGP02: These dgp's (can) contain added noise: v1, v2, ...
      if(ABDGP1_local) {
        if(number_of_variables == 3) { #if not TRUE, bvb if > 3, then v4,v5,... should be added
          if(without_intercept) {
            if(ABintercept) {
              #stopifnot(length(theta[,i]) == 2)
              afw = theta[,i] - (theta_real[,g_real[i]] + c(v2[i], v3[i])) #v1 not necessary (it belongs to the de facto intercept)
            } else {
              afw = theta[,i] - (theta_real[,g_real[i]] + c(v1[i], v2[i], v3[i])) #theta_real[,...] has number_of_variables elements
            }
          } else {
            afw = theta[,i] - (theta_real[,g_real[i]] + c(0, v1[i], v2[i], v3[i])) #theta_real[,...] has (number_of_variables + 1) elements
          }
        } else { #when number_of_variables != 3

          print("return NA")
          return(NA)

        }
      } else { #Case of DGP03 and DGP04:
        afw = theta[,i] - theta_real[,g_real[i]]
      }
      pre_mse[i] = mean(afw^2) #this is (theta_hat - theta)^2 (mean() goes over the number of variables)
    }
  }


  return(mean(pre_mse)) #this is E[(theta_hat - theta_real)^2]
}




#' Calculates VC², to determine the stability of the found number of groups and factors over the subsamples.
#'
#' VC² depends on C (this is the scale parameter in PIC). When VC² is equal to zero, the found number of groups and factors
#' are the same over the subsamples.
#' @param rcj dataframe containg the numer of groupfactors for all candidate C's and all subsamples
#' @param rc dataframe containg the numer of common factors for all candidate C's and all subsamples
#' @param C_candidates candidates for C (parameter in PIC)
#' @inheritParams calculate_lambda_group
#' @export
calculate_VCsquared <- function(rcj,rc,C_candidates, UPDATE1 = FALSE, UPDATE2 = FALSE,
                                number_of_group_factors = aantalfactoren_groups) {
  VC_squared = rep(NA, length(C_candidates))

  #vector with max number of groups
  Smax_local = rep(Smax,length(C_candidates))

  NUMBER_SUBSETS = length(a1) #number of subsamples used

  for(C_local in C_candidates) {
    C_index = which(C_candidates == C_local)
    part1 = 0
    for(i in 1:NUMBER_SUBSETS) {
      part1 = part1 + (rc[C_index,i] - sum(rc[C_index,])/NUMBER_SUBSETS)^2
    }
    part1 = part1 / NUMBER_SUBSETS

    part2 = 0
    for(j in 1:Smax_local[C_index]) {
      #print(paste("---",j))
      if(j <= length(number_of_group_factors) | UPDATE1 | UPDATE2) {
        part2_part = 0
        for(aA in 1:NUMBER_SUBSETS) {
          temp = 0
          for(bA in 1:NUMBER_SUBSETS) {
            if(j <= length(unlist(rcj[C_index, bA] %>% str_split("_")))) {
              temp = sum(temp, as.numeric(map(str_split(rcj[C_index, bA],"_"),j)), na.rm=T)
            }
          }
          if(j <= length(unlist(rcj[C_index, aA] %>% str_split("_")))) {
            rcjNT = as.numeric(map(str_split(rcj[C_index, aA], "_"), j))
          } else {
            rcjNT = 0
          }
          if(is.na(rcjNT)) rcjNT = 0
          part2_part = part2_part + (rcjNT - temp/NUMBER_SUBSETS)^2
        }
        part2_part = part2_part / NUMBER_SUBSETS
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
#' KS2014 is the recommended setting. It however does lead to increased computational time.
#' @param parameter_y dependent variable in regression
#' @param parameter_x independent variables in regression
#' @param nointercept if TRUE it performs regression without an intercept
#' @param nosetting option to remove the recommended setting in lmrob(). It is faster. Defaults to FALSE.
#' @export
LMROB <- function(parameter_y, parameter_x, nointercept = FALSE, nosetting = exists("run_lmrob_without_setting")) {

  if(is.na(parameter_x)) { #when there are no independent variables
    result2 = tryCatch(
      {
        if(nosetting) {
          result = lmrob(parameter_y ~ 1)
        } else {
          result = lmrob(parameter_y ~ 1, setting="KS2014")
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

          result2 = lmrob(parameter_y ~ parameter_x + 0, setting="KS2014")

      } else {
        #sometimes Error in if (init$scale == 0) happens. In that case: run without setting = KS2014.

          result2 = tryCatch(
            {
              lmrob(parameter_y ~ parameter_x, setting="KS2014")
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
#' Sigma2 is the sum of the squared errors, divided by NT. We need the sigma2 of the maxmodel to use (in term 2,3,4 of the PIC) instead of the configuration-dependent sigma2. (based on text in AndoBai 2016).
#' sigma2_max_model could actually be set to 1 as well, as it can be absorbed in parameter C of the PIC.
#' @inheritParams determine_theta
#' @inheritParams calculate_lambda_group
#' @param number_of_groups_candidates vector with candidate values for the number of groups
#' @param number_of_common_factors number of common factors
#' @export
calculate_sigma2maxmodel <- function(optimize_kappa = FALSE, UPDATE1 = update1, UPDATE2 = update2,
                                     number_of_groups = aantalgroepen,
                                     number_of_groups_candidates = aantalgroepen_candidates,
                                     number_of_group_factors = aantalfactoren_groups,
                                     number_of_common_factors = aantalfactoren_common
                                     ) {
  if(!optimize_kappa) {
    if(!UPDATE1 & !UPDATE2 &
       mean(number_of_group_factors,na.rm=T) == k_max &
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
#' The PIC was until now calculated with a sigma2 specific to the configuration (in terms of number of groups and factors).
#' Because the method requires sigma2 to be equal over all configurations (see proofs of different papers of Ando/Bai) we replace sigma2 by sigma2 of the configuration
#' with maximum number of groups and factors (this is the last one that is run).
#' @param all_pic contains PIC for all candidate C's and all subsamples
#' @param sigma2_max_model sigma2 of model with maximum number of groups and factors
#' @inheritParams calculate_lambda_group
#' @export
adapt_allpic_with_sigma2maxmodel <- function(all_pic, sigma2_max_model, UPDATE1 = update1, UPDATE2 = update2) {
  if(is.null(all_pic))  message("Warning in adapt_allpic_with_sigma2maxmodel(): all_pic is empty")
  if(!UPDATE1 & !UPDATE2) {
    message("Use sigma2_max_model in PIC (in object all_pic).")
    for(i in 1:nrow(all_pic)) {
      SIGMA2 = as.numeric(df_results$sigma2[i])
      all_pic[i,] = (all_pic[i,] - SIGMA2) / SIGMA2 * sigma2_max_model + SIGMA2
    }
    message("...done")
  }
  return(all_pic)
}
