
# Prints a dataframe
#
# make_df_simulationresults <- function() {
#   #Results:
#   df_sim = df_sim %>% filter(!is.na(S))
#   df_sim_special1 = df_sim_special1 %>% filter(!is.na(S))
#   df_sim_special2 = df_sim_special2 %>% filter(!is.na(S))
#   df_sim_special3 = df_sim_special3 %>% filter(!is.na(S))
#   df_sim_special4 = df_sim_special4 %>% filter(!is.na(S))
#   df_sim_special5 = df_sim_special5 %>% filter(!is.na(S))
#
#   simresults = matrix(NA,6,3)
#   simresults[1,1] = (df_sim %>% filter(S == 3) %>% nrow) / nrow(df_sim) * 100
#   simresults[1,2] = (df_sim %>% filter(k1 == 3) %>% nrow) / nrow(df_sim) * 100
#   simresults[1,3] = (df_sim %>% filter(S == 3,k1 == 3) %>% nrow) / nrow(df_sim %>% filter(S == 3)) * 100
#
#   simresults[2,1] = (df_sim_special1 %>% filter(S == 3) %>% nrow) / nrow(df_sim_special1) * 100
#   simresults[2,2] = (df_sim_special1 %>% filter(k1 == 3) %>% nrow) / nrow(df_sim_special1) * 100
#   simresults[2,3] = (df_sim_special1 %>% filter(S == 3,k1 == 3) %>% nrow) / nrow(df_sim_special1 %>% filter(S == 3)) * 100
#
#   simresults[3,1] = (df_sim_special2 %>% filter(S == 3) %>% nrow) / nrow(df_sim_special2) * 100
#   simresults[3,2] = (df_sim_special2 %>% filter(k1 == 3) %>% nrow) / nrow(df_sim_special2) * 100
#   simresults[3,3] = (df_sim_special2 %>% filter(S == 3,k1 == 3) %>% nrow) / nrow(df_sim_special2 %>% filter(S == 3)) * 100
#
#
#   simresults[4,1] = (df_sim_special3 %>% filter(S == 3) %>% nrow) / nrow(df_sim_special3) * 100
#   simresults[4,2] = (df_sim_special3 %>% filter(k1 == 3) %>% nrow) / nrow(df_sim_special3) * 100
#   simresults[4,3] = (df_sim_special3 %>% filter(S == 3,k1 == 3) %>% nrow) / nrow(df_sim_special3 %>% filter(S == 3)) * 100
#
#   simresults[5,1] = (df_sim_special4 %>% filter(S == 3) %>% nrow) / nrow(df_sim_special4) * 100
#   simresults[5,2] = (df_sim_special4 %>% filter(k1 == 3) %>% nrow) / nrow(df_sim_special4) * 100
#   simresults[5,3] = (df_sim_special4 %>% filter(S == 3,k1 == 3) %>% nrow) / nrow(df_sim_special4 %>% filter(S == 3)) * 100
#
#   simresults[6,1] = (df_sim_special5 %>% filter(S == 3) %>% nrow) / nrow(df_sim_special5) * 100
#   simresults[6,2] = (df_sim_special5 %>% filter(k1 == 3) %>% nrow) / nrow(df_sim_special5) * 100
#   simresults[6,3] = (df_sim_special5 %>% filter(S == 3,k1 == 3) %>% nrow) / nrow(df_sim_special5 %>% filter(S == 3)) * 100
#
#   print(simresults)
# }


#' Makes df with summarystats about % S and % k correct
#'
#' @param conditional_k_correct calculates percentage that number of groupfactors is correct under assumption that number of groups is correct
#' @export
prepare_dfQ <- function(df, cutsize, conditional_k_correct = FALSE) {
  df2 = df %>%
    mutate(grens = cut(sd_epsilon,breaks = c(-Inf,seq(1e-14,2,cutsize),Inf)),
           grens2 = as.character(grens),
           grens2 = as.numeric(str_sub(grens2, str_locate(grens2,",")[,1] + 1, nchar(grens2)-1))) %>%
    group_by(grens2) %>% mutate(S2 = S == 2,S3 = S == 3,S4 = S == 4,
                                #groups: before filtering
                                meanS = mean(S),
                                S2 = mean(S2)*100, S3 = mean(S3)*100, S4 = mean(S4)*100) %>%
    ungroup()

  if(conditional_k_correct) {
    df2 = df2 %>%
      filter(S == 3)
  }
  df2 = df2 %>% group_by(grens2) %>% mutate(gf2 = k1 == 2,gf3 = k1 == 3,gf4 = k1 == 4) %>%
    summarize(n = n(),
              #factors: after filtering
              meank1 = mean(k1),
              gf2 = mean(gf2)*100, gf3 = mean(gf3)*100, gf4 = mean(gf4)*100,

              meanS = mean(meanS),S2 = mean(S2),S3 = mean(S3),S4 = mean(S4), #-> this is to make sure meanS, S2,S3 and S4 are in the final data frame


              r2_factorgroupstructure = mean(r2_factorgroupstructure, na.rm = T),
              r2_factorgroup = mean(r2_factorgroup, na.rm = T),
              r2_lambdagroup = mean(r2_lambdagroup, na.rm = T))
  return(df2)
}

# Filters inputdataframe
#'
df_filter <- function(df, aantal_N_filter = 300, FGR_filter = 0.001, LGR_filter = 0.001, error_term_special_filter = 0) {
  return(df %>% filter(aantal_N == aantal_N_filter & aantal_T == 30,
                       FGR_FACTOR == FGR_filter & LGR_FACTOR == LGR_filter & S_real == 3 & k1_real == 3 &
                         T_FACTOR == 1 & error_term_special == error_term_special_filter))
}

#' Makes a plot of ...
#'
#' @export
simplot1 <- function(aantal_N_filter, FGR_filter, LGR_filter, error_term_special_filter = 0, cutsize = 0.00005,
                     n_limit = 50, conditional_k_correct = FALSE, xlimiet = 0.02) {
  Q_sim = prepare_dfQ(df_sim %>% df_filter(aantal_N_filter = aantal_N_filter, FGR_filter = FGR_filter, LGR_filter = LGR_filter, error_term_special_filter = error_term_special_filter),cutsize = cutsize, conditional_k_correct = conditional_k_correct) %>%
    filter(n > n_limit)
  Q_sim_special1 = prepare_dfQ(df_sim_special1 %>% df_filter(aantal_N_filter = aantal_N_filter, FGR_filter = FGR_filter, LGR_filter = LGR_filter, error_term_special_filter = error_term_special_filter),cutsize = cutsize, conditional_k_correct = conditional_k_correct) %>%
    filter(n > n_limit)
  Q_sim_special2 = prepare_dfQ(df_sim_special2 %>% df_filter(aantal_N_filter = aantal_N_filter, FGR_filter = FGR_filter, LGR_filter = LGR_filter, error_term_special_filter = error_term_special_filter),cutsize = cutsize, conditional_k_correct = conditional_k_correct) %>%
    filter(n > n_limit)
  Q_sim_special3 = prepare_dfQ(df_sim_special3 %>% df_filter(aantal_N_filter = aantal_N_filter, FGR_filter = FGR_filter, LGR_filter = LGR_filter, error_term_special_filter = error_term_special_filter),cutsize = cutsize, conditional_k_correct = conditional_k_correct) %>%
    filter(n > n_limit)
  Q_sim_special4 = prepare_dfQ(df_sim_special4 %>% df_filter(aantal_N_filter = aantal_N_filter, FGR_filter = FGR_filter, LGR_filter = LGR_filter, error_term_special_filter = error_term_special_filter),cutsize = cutsize, conditional_k_correct = conditional_k_correct) %>%
    filter(n > n_limit)
  Q_sim_special5 = prepare_dfQ(df_sim_special5 %>% df_filter(aantal_N_filter = aantal_N_filter, FGR_filter = FGR_filter, LGR_filter = LGR_filter, error_term_special_filter = error_term_special_filter),cutsize = cutsize, conditional_k_correct = conditional_k_correct) %>%
    filter(n > n_limit)
  Q_sim_special6 = prepare_dfQ(df_sim_special6 %>% df_filter(aantal_N_filter = aantal_N_filter, FGR_filter = FGR_filter, LGR_filter = LGR_filter, error_term_special_filter = error_term_special_filter),cutsize = cutsize, conditional_k_correct = conditional_k_correct) %>%
    filter(n > n_limit)

  print(Q_sim_special6)
  print(Q_sim_special5[1:3,])
  print(Q_sim_special4[1:3,])
  print(Q_sim_special3[1:3,])
  print(Q_sim_special2[1:3,])
  print(Q_sim_special1[1:3,])
  print(Q_sim[1:3,])
  print(nrow(Q_sim))
  print(nrow(Q_sim_special6))

  if(nrow(Q_sim) > 0) {
    p3 = ggplot(Q_sim_special1) + ggtitle("Aangepaste PIC 1: Aantal") +
      geom_point(aes(grens2,meanS,group=1)) + geom_line(aes(grens2,meanS,group=1)) +
      geom_point(aes(grens2,meank1,group=1), color = "blue") + geom_line(aes(grens2,meank1,group=1), color = "blue") +
      geom_line(data = Q_sim,aes(grens2,meanS,group=1, alpha = 0.2)) +
      geom_line(data = Q_sim,aes(grens2,meank1,group=1, alpha = 0.2), color = "blue") +
      geom_line(data = Q_sim_special2,aes(grens2,meank1,group=1, alpha = 0.2), color = "red") +
      geom_line(data = Q_sim_special3,aes(grens2,meank1,group=1, alpha = 0.2), color = "green") +
      guides(alpha = FALSE)
    if(nrow(Q_sim_special5) > 0) {
      p3 = p3 +
        geom_line(data = Q_sim_special4,aes(grens2,meank1,group=1, alpha = 0.2), color = "orange") +
        geom_line(data = Q_sim_special5,aes(grens2,meank1,group=1, alpha = 0.2), color = "yellow") +
        geom_line(data = Q_sim_special6,aes(grens2,meank1,group=1, alpha = 0.2), color = "purple")
    }
    p4 = ggplot(Q_sim_special1) + ggtitle("Aangepaste PIC 1: Percent correct") +
      geom_point(aes(grens2,S3,group=1)) + geom_line(aes(grens2,S3,group=1)) +
      geom_point(aes(grens2,gf3,group=1), color = "blue") + geom_line(aes(grens2,gf3,group=1), color = "blue") +
      geom_line(data = Q_sim,aes(grens2,S3,group=1, alpha = 0.2)) +
      geom_line(data = Q_sim,aes(grens2,gf3,group=1, alpha = 0.2), color = "blue") +
      geom_line(data = Q_sim_special2,aes(grens2,gf3,group=1, alpha = 0.2), color = "red") +
      geom_line(data = Q_sim_special3,aes(grens2,gf3,group=1, alpha = 0.2), color = "green")+
      ylab("Accuracy") +
      guides(alpha = FALSE)
    if(nrow(Q_sim_special5) > 0) {
      p4 = p4 +
        geom_line(data = Q_sim_special4,aes(grens2,gf3,group=1, alpha = 0.2), color = "orange") +
        geom_line(data = Q_sim_special5,aes(grens2,gf3,group=1, alpha = 0.2), color = "yellow") +
        geom_line(data = Q_sim_special6,aes(grens2,gf3,group=1, alpha = 0.2), color = "purple")
    }

    hoofdtitel = str_c("FGR = ", FGR_filter, " & LGR = ",LGR_filter)
    if(conditional_k_correct) {
      hoofdtitel = str_c(hoofdtitel," (+ conditional on correct S)")
    }
    grid.arrange(
      ggplot(Q_sim) + ggtitle("Standaard PIC: Aantal")  +
        geom_point(aes(grens2,meanS,group=1)) + geom_line(aes(grens2,meanS,group=1)) +
        geom_point(aes(grens2,meank1,group=1), color = "blue") +
        geom_line(aes(grens2,meank1,group=1), color = "blue") +
        coord_cartesian(xlim = c(0,xlimiet)),
      ggplot(Q_sim) + ggtitle("Standaard PIC: Percent correct")  +
        geom_point(aes(grens2,S3,group=1)) + geom_line(aes(grens2,S3,group=1)) +
        geom_point(aes(grens2,gf3,group=1), color = "blue") +
        geom_line(aes(grens2,gf3,group=1), color = "blue") +
        ylab("Accuracy") +
        coord_cartesian(xlim = c(0,xlimiet)),
      p3 + coord_cartesian(xlim = c(0,xlimiet)),
      p4 + coord_cartesian(xlim = c(0,xlimiet)),
      nrow = 2, ncol = 2,
      top = hoofdtitel
    )
  } else {
    print("too little data")
  }
}
