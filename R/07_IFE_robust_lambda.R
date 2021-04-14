




#' Help-function for return_robust_lambdaobject().
#'
#' Uses the almost classical lambda to create a robust lambda by using M estimation.
#'
#' @param almost_classical_lambda matrix where the mean of each row is equal to the classical lambda
#' @param fastoption Uses nlm() instead of optim(). This is faster.
#' @param fastoption2 experimental parameter: can speed nlm() up (10%), but loses accuracy. May benefit from finetuning.
#' @param application_covid_DRL indicator of using the coviddataset
#' @importFrom stats optim
#' @importFrom stats nlm
#' @return M-estimator of location of the parameter, by minimizing sum of rho()
determine_robust_lambda <- function(almost_classical_lambda, fastoption = TRUE, fastoption2 = FALSE, application_covid_DRL = FALSE) {

  ############
  #speedtests:
  # test = t(apply(Y[1:500,], 1, function(x) x * factor_for_grouping[,1]))
  # micro = microbenchmark::microbenchmark(f1=apply(test, 1, function(x) determine_robust_lambda(x, fastoption = FALSE)),
  #                                        f2=apply(test, 1, function(x) determine_robust_lambda(x, fastoption = TRUE, fastoption2 = FALSE)),
  #                                        f3=apply(test, 1, function(x) determine_robust_lambda(x, fastoption = TRUE, fastoption2 = TRUE)),times = 50)
  # time gain of using nlm():
  # when N = 500, time gain is on average 117/173 (68%) when using fastoption.
  # when N = 2000, time gain is on average 475/713 (67%) when using fastoption.

  #the option fastoption2 wins about 10% extra, but would not be accurate enough.
  ############
  if(application_covid_DRL) {
    #nlm() does not find minimum, so use optim()
    fastoption = FALSE
  }

  almost_classical_lambda = unlist(almost_classical_lambda) #because sometimes is not a vector (bvb with eclipz with N=3112)
  MADCL =  mad(almost_classical_lambda) #this is a value for one person and one factor (depends on i & r)



  #special case where number of factors are updated during algorithm.
  #Then sometimes groups are empty leading to no estimation of groupfactors
  #and to 0/0 errors in this optimisation.
  #Solve this by setting the denominator > 0
  if(min(almost_classical_lambda) == 0 & max(almost_classical_lambda) == 0) {
    MADCL = 0.000001
  }
  #do the same when almost_classical_lambda has too many zero's
  if(MADCL == 0) {
    MADCL = 0.000001
  }
  #m-estimate: minimize sum of rho's in scaled observations
  sum_of_rho <- function(x) {
    sum( Mpsi((almost_classical_lambda - x) / MADCL , cc = 4.685, psi = "bisquare", deriv = -1) )
  }

  if(fastoption) {
    if(fastoption2) {
      robust_lambda = nlm(f = sum_of_rho, p = 0, steptol=1e-3, gradtol=1e-3)$estimate
    } else {
      robust_lambda = nlm(f = sum_of_rho, p = 0)$estimate
    }

  } else {
    robust_lambda = optim(par = 0, fn = sum_of_rho, method = "L-BFGS-B")$par
  }

  return(robust_lambda)
}







#' Calculates robust loadings
#'
#' Uses the almost classical lambda (this is an object of which the mean equals to the classical lambda) to create a robust lambda by using M estimation
#' @param Y_like_object this is Y_ster or W or W_j
#' @param group index of group
#' @param type scalar which shows in which setting this function is used
#' @param factor_group_RRN estimation of groupfactors
#' @param number_of_group_factors_RRN number of group factors
#' @param comfactor_RRN estimation of common factors
#' @param number_of_common_factors_RRN number of common factors
#' @param NN N
#' @param eclipz_RRN indicator of using the eclipzdataset
#' @param application_covid indicator of using the coviddataset
#' @param verbose when TRUE, it prints messages
#' @return Nxk dataframe
#' @export
return_robust_lambdaobject <- function(Y_like_object, group, type,
                                       NN = aantal_N,
                                       factor_group_RRN = factor_group,
                                       number_of_group_factors_RRN = aantalfactoren_groups,
                                       comfactor_RRN = comfactor,
                                       number_of_common_factors_RRN = aantalfactoren_common,
                                       eclipz_RRN = FALSE,
                                       application_covid = FALSE,
                                       verbose = FALSE) {

  if(verbose) message(paste("type =",type))
  if(type == 1) {  #used in calculate_virtual_factor_and_lambda_group()
    if(number_of_group_factors_RRN[group] > 0) { #this can be zero when updating the number of factors
      LG_local = data.frame(matrix(NA, nrow = NN, ncol = number_of_group_factors_RRN[group]))


      for(ii in 1:NN) {
        #if(verbose) message(ii)
        for(rr in 1:number_of_group_factors_RRN[group]) {
          #if(verbose) message(rr)

          almost_classical_lambda = (Y_like_object[ii,] * t(factor_group_RRN[[group]])[,rr])  #the mean of this is equal to classical lambda

          if(eclipz_RRN) almost_classical_lambda = unlist(almost_classical_lambda) #to evade errors
          if(verbose) message(paste("start determine_robust_lambda()"))
          robust_lambda = determine_robust_lambda(almost_classical_lambda, application_covid_DRL = application_covid)
          if(verbose) message(paste("end determine_robust_lambda()"))
          LG_local[ii,rr] = robust_lambda
        }
      }
    } else { #otherwise set all to zero
      LG_local = data.frame(matrix(0, nrow = NN, ncol = 1))
    }


    return(LG_local)
  }


  if(type == 2) { #used in calculate_lambda()
    if(number_of_common_factors_RRN > 0) {
      lambda = matrix(NA, nrow = number_of_common_factors_RRN, ncol = NN)

      for(ii in 1:NN) {
        for(rr in 1:number_of_common_factors_RRN) {

          almost_classical_lambda = (Y_like_object[ii,] * t(comfactor_RRN)[,rr]) #the mean of this is equal to the classical lambda

          robust_lambda = determine_robust_lambda(almost_classical_lambda, application_covid_DRL = application_covid)
          lambda[rr,ii] = robust_lambda
        }
      }
    } else {
      lambda = matrix(rep(0,NN), nrow = 1) #all zero's
    }
    return(lambda)
  }



  if(type == 3) { #used in calculate_lambda_group()
    lambda_local = data.frame(matrix(NA,nrow = length(which(g == group)), ncol = number_of_group_factors_RRN[group]))

    for(ii in 1:length(which(g == group))) {
      for(rr in 1:number_of_group_factors_RRN[group]) {
        almost_classical_lambda = (Y_like_object[ii,] * t(factor_group_RRN[[group]])[,rr])  #the mean of this is equal to classical lambda

        robust_lambda = determine_robust_lambda(almost_classical_lambda, application_covid_DRL = application_covid)

        lambda_local[ii,rr] = robust_lambda
      }
    }
    return(lambda_local) #this is the result for 1 group

  }


  if(type == 4) { #used in initialisation of grouping (which works by estimating a lambdalike object)
    #-> uses 'factor_for_grouping'
    lambda_local = data.frame(matrix(NA,
                                     nrow = nrow(Y_like_object),
                                     ncol = ncol(factor_for_grouping)))




    for(ii in 1:nrow(Y_like_object)) {
      for(rr in 1:ncol(factor_for_grouping)) {
        stopifnot(length(Y_like_object[ii,]) == length(factor_for_grouping[,rr]))
        #reason can be that input in macropca() contained columns with zero or tiny median absolute deviation.
        #These are ignored in analysis.

        almost_classical_lambda = (Y_like_object[ii,] * factor_for_grouping[,rr]) #Do not use the transpose here, because we use factor_for_grouping here.

        lambda_local[ii,rr] = determine_robust_lambda(almost_classical_lambda, application_covid_DRL = application_covid)
      }
    }

    #############################################
    #this is equally fast/slow for (N=500,T=100), so keep the double loop from above
    # for(rr in 1:ncol(factor_for_grouping)) {
    #   test = t(apply(Y_like_object, 1, function(x) x * factor_for_grouping[,1]))
    #   lambda_local[,rr] = apply(test, 1, function(x) determine_robust_lambda(x))
    # }
    #############################################

    # micro = microbenchmark::microbenchmark(f1(),f2(),times = 20)
    # summary(micro)



    return(lambda_local)

  }
}
