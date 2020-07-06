
#Script contains functions to estimate robust lambda's


#' calculates vector of which the mean is equal to classical lambda
#'
#' @return matrix where the mean of each row is equal to the classical lambda.
#' @export
all_almost_classical_lambda <- function(Y_like_object, factor_like_object, NUMBER_FACTORS, group, NUMBER_PERSONS) {
  ding = matrix(NA, nrow = NUMBER_PERSONS, ncol = aantal_T)
  for(ii in 1:NUMBER_PERSONS) {
    for(rr in 1:NUMBER_FACTORS) {

      almost_classical_lambda = (Y_like_object[ii,] * t(factor_like_object)[,rr])  #the mean of this is equal to classical lambda
      ding[ii,] = almost_classical_lambda
    }
  }
  return(ding)
}



#' Help-function for return_robust_lambdaobject().
#'
#' Uses the almost classical lambda to create a robust lambda by using M estimation
#' @return M-estimator of location of the parameter, by minimizing sum of rho()
determine_robust_lambda <- function(almost_classical_lambda) {
  #print("determine_robust_lambda")
  almost_classical_lambda = unlist(almost_classical_lambda) #because sometimes is not a vector (bvb bij eclipz met N= 3112)
  MADCL =  mad(almost_classical_lambda) #this is a value for one person and one factor (depends on i & r)
  #print(MADCL)

  if(exists("robust_lambda_MM")) {
    # print("testing stuff...")
    # print(almost_classical_lambda)
    #This is the scale belonging to an S-estimator
    MADCL = LMROB(almost_classical_lambda, NA)$scale #dit is: lmrob(almost_classical_lambda ~ 1, setting="KS2014")$scale
  }

  #m-estimate: minimize sum of rho's in scaled observations
  sum_of_rho <- function(x) {
    sum( Mpsi((almost_classical_lambda - x) / MADCL , cc = 4.685, psi = "bisquare", deriv = -1) )
  }
  # if(eclipz) {
  #   print(almost_classical_lambda)
  #   print(MADCL)
  # }

  #check if rho-function = x^2 -> then this should be the same as the classical lambda (= mean(almost_classical_lambda))
  # testje <- function(x) {
  #   return((almost_classical_lambda-x)^2 %>% sum)
  # }
  # optim(par = 0, fn = testje,method = "L-BFGS-B")$par -  mean(almost_classical_lambda) < 1e13 #moet TRUE zijn

  robust_lambda = optim(par = 0, fn = sum_of_rho,method = "L-BFGS-B")$par
  return(robust_lambda)
}

#One might think that using a predefined function would be faster then using optim().
#Microbenchmark says it is not true: time of determine_robust_lambda2() is slightly bigger even.
#Although I had not checked how close the results are to each other. (Also this function does not work (yet) for > 1 factor.)
# determine_robust_lambda2 <- function(almost_classical_lambda) {
#   robust_lambda = rlm(almost_classical_lambda~1,psi="psi.bisquare")
#   return(robust_lambda)
# }









#' Calculates robust loadings
#'
#' Uses the almost classical lambda to create a robust lambda by using M estimation
#' @param Y_like_object this is Y_ster or W or W_j
#' @param group group
#' @param type scalar which shows in which setting this function is used
#' @return Nxk dataframe
#' @export
return_robust_lambdaobject <- function(Y_like_object, group, type,
                                       FACTOR_GROUP = factor_group, AANTALFACTOREN_GROUPS = aantalfactoren_groups,
                                       FACTOR = comfactor, AANTALFACTOREN_COMMON = aantalfactoren_common) {
  #print("return_robust_lambdaobject")
  #print(type)
  #used in calculate_virtual_factor_and_lambda_group
  if(type == 1) {
    LG_local = data.frame(matrix(NA,nrow = aantal_N, ncol = AANTALFACTOREN_GROUPS[group]))
    #Deze MADCL diende om de sigma in the M-estimator onafhankelijk van i te maken
    #MADCL = mad(all_almost_classical_lambda(Y_like_object, FACTOR_GROUP[[group]], AANTALFACTOREN_GROUPS[group], group, aantal_N))

    almost_classical_lambda = lapply(1:aantal_N, function(x) lapply(1:AANTALFACTOREN_GROUPS[group], function(y) (Y_like_object[x,] * t(FACTOR_GROUP[[group]])[,y])))

    for(ii in 1:aantal_N) {
      for(rr in 1:AANTALFACTOREN_GROUPS[group]) {

        almost_classical_lambda = (Y_like_object[ii,] * t(FACTOR_GROUP[[group]])[,rr])  #the mean of this is equal to classical lambda

        if(eclipz) almost_classical_lambda = unlist(almost_classical_lambda) #anders krijgen we foutmeldingen

        #print(mean(almost_classical_lambda))
        # print(ii)
        # print(rr)
        robust_lambda = determine_robust_lambda(almost_classical_lambda)
        LG_local[ii,rr] = robust_lambda
      }
    }

#Testing vectorizing with this:
# Is not faster.

    # hulpf <- function(x,y) {
    #   almost_classical_lambda = (Y[x,] * t(factor_group[[group]])[,y])
    #   if(eclipz) almost_classical_lambda = unlist(almost_classical_lambda) #anders krijgen we foutmeldingen
    #   return(determine_robust_lambda(almost_classical_lambda))
    # }
    # f2<-function() {
    #   robust_lambda = lapply(1:aantal_N, function(x) lapply(1:aantalfactoren_groups[group], function(y) hulpf(x,y)))
    #   LG_local2 = data.frame(t(matrix(map(robust_lambda, unlist) %>% unlist, nrow = aantalfactoren_groups[group])))
    #   return(LG_local2)
    # }



    return(LG_local)
  }

  #used in calculate_lambda()
  if(type == 2) {
    if(AANTALFACTOREN_COMMON > 0) {
      lambda = matrix(NA, nrow = AANTALFACTOREN_COMMON, ncol = aantal_N)
      #MADCL = mad(all_almost_classical_lambda(Y_like_object, FACTOR, AANTALFACTOREN_COMMON, NA, aantal_N))
      for(ii in 1:aantal_N) {
        for(rr in 1:AANTALFACTOREN_COMMON) {

          almost_classical_lambda = (Y_like_object[ii,] * t(FACTOR)[,rr]) #the mean of this is equal to the classical lambda

          robust_lambda = determine_robust_lambda(almost_classical_lambda)
          lambda[rr,ii] = robust_lambda
        }
      }
    } else {
      lambda = matrix(rep(0,aantal_N), nrow = 1) #all zero's
    }
    return(lambda)
  }


  #used in calculate_lambda_group()
  if(type == 3) {
    lambda_local = data.frame(matrix(NA,nrow = length(which(g == group)), ncol = AANTALFACTOREN_GROUPS[group]))

    #length(which(g == group)) == nrow(Y_like_object)
    #MADCL = mad(all_almost_classical_lambda(Y_like_object, FACTOR_GROUP[[group]], AANTALFACTOREN_GROUPS[group], group, length(which(g == group))  ))
    for(ii in 1:length(which(g == group))) {
      for(rr in 1:AANTALFACTOREN_GROUPS[group]) {
        almost_classical_lambda = (Y_like_object[ii,] * t(FACTOR_GROUP[[group]])[,rr])  #the mean of this is equal to classical lambda
        #print(mean(almost_classical_lambda))

        robust_lambda = determine_robust_lambda(almost_classical_lambda)


        lambda_local[ii,rr] = robust_lambda
      }
    }
    return(lambda_local) #this is the result for 1 group

  }

  #used in initialisation of grouping (which works by estiamting a lambdalike object)
  #-> uses 'factor_for_grouping'
  if(type == 4) {
    lambda_local = data.frame(matrix(NA,
                                     nrow = nrow(Y_like_object),
                                     ncol = nrow(factor_for_grouping)))
    for(ii in 1:nrow(Y_like_object)) {
      for(rr in 1:nrow(factor_for_grouping)) {
        #print(ii)
        #print(rr)
        almost_classical_lambda = (Y_like_object[ii,] * t(factor_for_grouping)[,rr])
        #print(almost_classical_lambda[1:3])
        #print(mean(almost_classical_lambda))
        #print(almost_classical_lambda)
        lambda_local[ii,rr] = determine_robust_lambda(almost_classical_lambda)
        #print("***")
      }
    }
    return(lambda_local)

  }
}
