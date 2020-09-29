


#' Calculates vector of which the mean is equal to classical lambda -> function should be obsolete
#'
#' @param Y_like_object NxT matrix
#' @param factor_like_object
#' @param NUMBER_FACTORS number of factors (common or group, depending on context)
#' @param group index of group
#' @param NN N
#' @return matrix where the mean of each row is equal to the classical lambda.
all_almost_classical_lambda <- function(Y_like_object, factor_like_object, NUMBER_FACTORS, group, NN) {
  ding = matrix(NA, nrow = NN, ncol = aantal_T)
  for(ii in 1:NN) {
    for(rr in 1:NUMBER_FACTORS) {

      almost_classical_lambda = (Y_like_object[ii,] * t(factor_like_object)[,rr])  #the mean of this is equal to classical lambda
      ding[ii,] = almost_classical_lambda
    }
  }
  return(ding)
}



#' Help-function for return_robust_lambdaobject().
#'
#' Uses the almost classical lambda to create a robust lambda by using M estimation.
#' @param almost_classical_lambda matrix where the mean of each row is equal to the classical lambda
#' @return M-estimator of location of the parameter, by minimizing sum of rho()
determine_robust_lambda <- function(almost_classical_lambda) {

  almost_classical_lambda = unlist(almost_classical_lambda) #because sometimes is not a vector (bvb with eclipz with N=3112)
  MADCL =  mad(almost_classical_lambda) #this is a value for one person and one factor (depends on i & r)

  #special case where number of factors are updated during algorithm.
  #Then sometimes groups are empty leading to no estimation of groupfactors
  #and to 0/0 errors in this optimisation.
  #Solve this by setting the denominator > 0
  if(min(almost_classical_lambda) == 0 & max(almost_classical_lambda) == 0) {
    MADCL = 0.001
  }
  #m-estimate: minimize sum of rho's in scaled observations
  sum_of_rho <- function(x) {
    sum( Mpsi((almost_classical_lambda - x) / MADCL , cc = 4.685, psi = "bisquare", deriv = -1) )
  }

  robust_lambda = optim(par = 0, fn = sum_of_rho,method = "L-BFGS-B")$par
  return(robust_lambda)
}







#' Calculates robust loadings
#'
#' Uses the almost classical lambda (this is an object of which the mean equals to the classical lambda) to create a robust lambda by using M estimation
#' @param Y_like_object this is Y_ster or W or W_j
#' @param group index of group
#' @param type scalar which shows in which setting this function is used
#' @param FACTOR_GROUP estimation of groupfactors
#' @param number_of_group_factors number of group factors
#' @param FACTOR estimation of common factors
#' @param number_of_common_factors number of common factors
#' @param NN N
#' @return Nxk dataframe
#' @export
return_robust_lambdaobject <- function(Y_like_object, group, type,
                                       FACTOR_GROUP = factor_group,
                                       number_of_group_factors = aantalfactoren_groups,
                                       FACTOR = comfactor,
                                       number_of_common_factors = aantalfactoren_common,
                                       NN = aantal_N) {


  if(type == 1) {  #used in calculate_virtual_factor_and_lambda_group()
    LG_local = data.frame(matrix(NA,nrow = NN, ncol = number_of_group_factors[group]))

    almost_classical_lambda = lapply(1:NN, function(x) lapply(1:number_of_group_factors[group], function(y) (Y_like_object[x,] * t(FACTOR_GROUP[[group]])[,y])))

    for(ii in 1:NN) {
      for(rr in 1:number_of_group_factors[group]) {

        almost_classical_lambda = (Y_like_object[ii,] * t(FACTOR_GROUP[[group]])[,rr])  #the mean of this is equal to classical lambda

        if(eclipz) almost_classical_lambda = unlist(almost_classical_lambda) #to evade errors
        #print(almost_classical_lambda)
        robust_lambda = determine_robust_lambda(almost_classical_lambda)
        LG_local[ii,rr] = robust_lambda
      }
    }



    return(LG_local)
  }


  if(type == 2) { #used in calculate_lambda()
    if(number_of_common_factors > 0) {
      lambda = matrix(NA, nrow = number_of_common_factors, ncol = NN)

      for(ii in 1:NN) {
        for(rr in 1:number_of_common_factors) {

          almost_classical_lambda = (Y_like_object[ii,] * t(FACTOR)[,rr]) #the mean of this is equal to the classical lambda

          robust_lambda = determine_robust_lambda(almost_classical_lambda)
          lambda[rr,ii] = robust_lambda
        }
      }
    } else {
      lambda = matrix(rep(0,NN), nrow = 1) #all zero's
    }
    return(lambda)
  }



  if(type == 3) { #used in calculate_lambda_group()
    lambda_local = data.frame(matrix(NA,nrow = length(which(g == group)), ncol = number_of_group_factors[group]))

    for(ii in 1:length(which(g == group))) {
      for(rr in 1:number_of_group_factors[group]) {
        almost_classical_lambda = (Y_like_object[ii,] * t(FACTOR_GROUP[[group]])[,rr])  #the mean of this is equal to classical lambda
        #print(mean(almost_classical_lambda))

        robust_lambda = determine_robust_lambda(almost_classical_lambda)

        lambda_local[ii,rr] = robust_lambda
      }
    }
    return(lambda_local) #this is the result for 1 group

  }


  if(type == 4) { #used in initialisation of grouping (which works by estiamting a lambdalike object)
    #-> uses 'factor_for_grouping'
    lambda_local = data.frame(matrix(NA,
                                     nrow = nrow(Y_like_object),
                                     ncol = nrow(factor_for_grouping)))
    for(ii in 1:nrow(Y_like_object)) {
      for(rr in 1:nrow(factor_for_grouping)) {
        almost_classical_lambda = (Y_like_object[ii,] * t(factor_for_grouping)[,rr])
        lambda_local[ii,rr] = determine_robust_lambda(almost_classical_lambda)
      }
    }
    return(lambda_local)

  }
}
