#' Makes a plot of the (weighted) mean and median income for each group
#'
#' Weights are based on the relative number of people that retire in a given year ("minpensioenjaar")
#' @param band adds sd-bands
#' @param groei if FALSE then it uses incomelevels instead of incomegrowth
#' @param add_grenzen adds grenzen van de pensioenbijdragen (enkel voor groei == FALSE)
#' @param Y_norm: normaliseren van de inkomens(groei)
#' @export
plot_typegevallen <- function(band = FALSE, groei = TRUE, add_grenzen = TRUE, Y_norm = TRUE) {


  print("Y_normalized(i,t) = Y(i,t) - mean(Y(t)) + mean(Y)")
  print("weighted mean and median values with pensioenjaar (currently this is implemented minpensioenjaar; maxpensioenjaar can be different)")
  print("sd en mad zijn niet gewogen (NaN-issues bij sd)")
  df = matrices_to_df(g, factor, lambda, use_Y_beforeMGandMZ = TRUE)

  grid2 = grid %>%
    bind_cols(df %>% dplyr::select(Y, starts_with("X"), group))

  normalise_Y <- function(grid2) {
    print("normalize Y with gt and g_bar (uses Y_beforeMGandMZ)")
    if(Y_norm) {
      grid2 = grid2 %>%
        group_by(Var2) %>% #group over time
        mutate(gt = mean(Y)) %>% #define gt = mean growth over all individuals
        ungroup() %>%
        mutate(g_bar = mean(Y)) %>% #define g_bar = mean growth over all individuals and all time
        mutate(Y = Y - gt + g_bar) #normalize Y (values are of Y_beforeMGandMZ (see parameter of matrices_to_df())) with gt and g_bar
    }
    return(grid2)
  }

  if(eclipz) {
    #collection of minpensioenjaar_aantal:
    temp = X_other %>%
      mutate(group = g) %>%
      group_by(group,minpensioenjaar) %>%
      mutate(minpensioenjaar_aantal = n()) %>%
      ungroup %>%
      dplyr::select(minpensioenjaar_aantal)

    #mean and median income b in 1984 by group; weighted with minpensioenjaar_aantal
    starting_Value = X_other %>%
      dplyr::select(starts_with("inkomen")) %>%
      mutate(mean = rowMeans(.),
             median = apply(.,1,median),
             group = g,
             minpensioenjaar_aantal = temp$minpensioenjaar_aantal) %>%
      group_by(group) %>%
      summarise(start_mean = weighted.mean(mean, 1/minpensioenjaar_aantal),
                start_median = spatstat::weighted.median(median, 1/minpensioenjaar_aantal))


    grid2_summary = normalise_Y(grid2) %>%
      mutate(minpensioenjaar = rep(X_other$minpensioenjaar, aantal_T)) %>% #add minimumpensioenjaar
      group_by(group,minpensioenjaar) %>%
      mutate(minpensioenjaar_aantal = n()) %>%
      ungroup() %>%
      group_by(group,Var2) %>% #group by group and T
      summarise(mean = weighted.mean(Y, 1/minpensioenjaar_aantal), #Y has the values of Y_beforeMGandMZ
                median = spatstat::weighted.median(Y, 1/minpensioenjaar_aantal), #Y has the values of Y_beforeMGandMZ
                sd = sd(Y), #geeft NaN: sqrt(Hmisc::wtd.var(Y, 1/minpensioenjaar_aantal))  #Y has the values of Y_beforeMGandMZ
                mad = mad(Y) #Y has the values of Y_beforeMGandMZ
      ) %>%
      left_join(starting_Value) %>%
      mutate(inkomen_mean = ifelse(Var2 == 1, start_mean*(1 + mean),NA),
             inkomen_median = ifelse(Var2 == 1, start_median*(1 + median),NA)) %>%
      ungroup()
    for(t in 2:aantal_T) {
      grid2_summary = grid2_summary %>%
        mutate(inkomen_mean = ifelse(Var2 == t, lag(inkomen_mean)*(1 + mean),inkomen_mean)) %>%
        mutate(inkomen_median = ifelse(Var2 == t, lag(inkomen_median)*(1 + median),inkomen_median))
    }


    data = import_inkomensgrenzen()
    data = data %>% filter(jaar >= 1985 & jaar <= eclipz_eindjaar)
    grid2_summary$mingrens = rep(data$BD_MIN , aantalgroepen)
    grid2_summary$maxgrens = rep(data$BD_MAX..maximum.bedrijfsinkomsten.in.pensioenberekening. , aantalgroepen)
    grid2_summary$tussengrens = rep(data$tussenplafond2 , aantalgroepen)

  } else { #if not eclipz
    starting_Value = data.frame(group = c(1,2,3), start_mean = 1000, start_median = 1000)

    grid2_summary = normalise_Y(grid2) %>%
      group_by(group,Var2) %>% #group by group and T
      summarise(mean = mean(Y_beforeMGandMZ),
                median = median(Y_beforeMGandMZ),
                sd = sd(Y_beforeMGandMZ) ,
                mad = mad(Y_beforeMGandMZ)
      ) %>%
      left_join(starting_Value) %>%
      mutate(inkomen_mean = ifelse(Var2 == 1, start_mean*(1 + mean),NA),
             inkomen_median = ifelse(Var2 == 1, start_median*(1 + median),NA)) %>%
      ungroup()
    for(t in 2:aantal_T) {
      grid2_summary = grid2_summary %>%
        mutate(inkomen_mean = ifelse(Var2 == t, lag(inkomen_mean)*(1 + mean),inkomen_mean)) %>%
        mutate(inkomen_median = ifelse(Var2 == t, lag(inkomen_median)*(1 + median),inkomen_median))
    }
  }




  if(groei) {
    titel_mean = "Mean income growth"
    if(Y_norm) titel_mean = "Mean (normalised) income growth"
    plot1 = ggplot(grid2_summary, aes(Var2, mean, color = as.factor(group))) +
      geom_point() +
      geom_line() +
      ggtitle(titel_mean) +
      xlab("Years") +
      ylab("Income growth") +
      scale_color_discrete(name="Group")
    plot2 = ggplot(grid2_summary, aes(Var2, mean, color = as.factor(group))) +
      geom_line(aes(Var2,median)) +
      ggtitle(str_replace(titel_mean,"Mean","Median")) +
      xlab("Years") +
      ylab("Income growth") +
      scale_color_discrete(name="Group")
    if(band) {
      plot1 = plot1 + geom_ribbon(aes(ymin = pmax(0, mean - sd/2), ymax = mean + sd/2, fill = "band", group = as.factor(group)), alpha = 0.2)
      plot2 = plot2 + geom_ribbon(aes(ymin = pmax(0, median - mad/2), ymax = median + mad/2, fill = "band", group = as.factor(group)), alpha = 0.2)
    }
  } else {
    titel_mean = "Mean income, by group"
    if(Y_norm) titel_mean = "Mean (normalised) income, by group"
    plot1 = ggplot(grid2_summary, aes(Var2, inkomen_mean, color = as.factor(group))) +
      geom_point() +
      geom_line() +
      #geom_line(aes(Var2,inkomen_median)) +
      ggtitle(titel_mean) +
      xlab("Years") +
      ylab("Income in EUR (x1000)") +
      scale_color_discrete(name="Group") +
      scale_y_continuous(labels = function(x) format(x/1000, scientific = FALSE))
    plot2 = ggplot(grid2_summary, aes(Var2, inkomen_mean, color = as.factor(group))) +
      geom_line(aes(Var2,inkomen_median)) +
      ggtitle(str_replace(titel_mean,"Mean","Median")) +
      xlab("Years") +
      ylab("Income in EUR") +
      scale_color_discrete(name="Group")


    if(add_grenzen & eclipz) {
      plot1 = plot1 +
        geom_line(aes(Var2, mingrens), color = "black",alpha = 0.3) +
        geom_line(aes(Var2, tussengrens), color = "black",alpha = 0.3) +
        geom_line(aes(Var2, maxgrens), color = "black",alpha = 0.3)
      plot2 = plot2 +
        geom_line(aes(Var2, mingrens), color = "black",alpha = 0.3) +
        geom_line(aes(Var2, tussengrens), color = "black",alpha = 0.3) +
        geom_line(aes(Var2, maxgrens), color = "black",alpha = 0.3)
    }

  }


  grid.arrange(plot1, plot2, nrow = 2, ncol = 1)

}
