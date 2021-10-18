
library(cowplot)
library(dplyr)
library(ggplot2)
library(nimue)
library(shiny)
library(squire)

countries <- sort(unique(squire::population$country))
iso3cs <- squire::population$iso3c[match(countries, squire::population$country)]

#' Run a nimue model using online fits
#'
#' @title nimue run
#' @param iso3c ISO3C country code
#' @param inf_eff Vaccine efficacy against infection. Default 0.9
#' @param dis_eff Vaccine efficacy against infection. Default 0.96
#' @param forecast How long into future to simulate. Default = 0
#' @param json_path Path of json for fits. Default = NULL, will grab online fit
#' @param dose_factor Scaling of OWID in data total vaccinations per day. I.e
#'   Default = 1, means the efficacy specified for the vaccine is achieved after
#'   1 dose. If 2, we halve the vaccinatioons given out per day to approximate
#'   say needing two doses.
#' @param max_vaccine Vaccine doses per week. Default = NULL which will use
#'   2.5% of pop size
#' @param vaccine_uptake Max vaccine coverage for all age groups. Default = 0.8
#' @param vaccine_available Max vaccine availability. Default = 0.2
#' @param vaccine_durability Vaccine protection duration. Default = 1095
#' @param risk_proportion Proportion high risk and thus vaccine prioritised.
#'   Default = 0.1
#' @param future_Rt Changes to Rt in the future.
#'   Default = numeric(0), which is no change. 1.5 would cause Rt to be 1.5
#' @param future_Rt_change Proportional changes to Rt in the future.
#'   Default = numeric(0), which is no change. 1.5 would cause the final Rt
#'   value to be multiplied by 1.5. If future_Rt is given, future_Rt takes
#'   priority
#' @param tt_Rt_changes Timing of changes to Rt in the future.
#'   Default = numeric(0), which is no change
#' @param future_vaccines Daily vaccines in the future.
#'   Default = numeric(0), which is no vaccine changes. 5000 would be 5000 vaccines
#'   delivered per day.
#' @param tt_future_vaccines Timing of vaccines in the future.
#'   Default = numeric(0), which is no vaccine changes from current roll out.
#' @param final_coverage Final vaccine coverage by end of forecast period. Default
#'   is numeric(0) which does nothing. If set, this will overright future vaccines
#'   such that coverage is reached linearly
#' @param strategy Rollout strategy for covidsim. Must be one of:
#'   "HCW and Elderly" (default), "HCW, Elderly and "High-Risk", "Elderly", "All"
#'
#' @export
create_vacc_fit <- function(iso3c,
                            json_path = NULL,
                            country = NULL,
                            inf_eff = 0.9,
                            dis_eff = 0.96,
                            forecast = 120,
                            dose_factor = 1,
                            max_vaccine = NULL,
                            vaccine_uptake = 0.8,
                            vaccine_available = 0.95,
                            vaccine_durability = 5000,
                            risk_proportion = 0.1,
                            future_Rt = numeric(0),
                            future_Rt_changes = numeric(0),
                            tt_Rt_changes = numeric(0),
                            future_vaccines = numeric(0),
                            final_coverage = numeric(0),
                            tt_future_vaccines = numeric(0),
                            strategy = "HCW, Elderly and High-Risk") {

  # get what country this is
  country <- squire::population$country[squire::population$iso3c==iso3c][1]

  # grab the json from the data exports
  if(is.null(json_path)) {
    file_path <- "https://raw.githubusercontent.com/mrc-ide/global-lmic-reports/master/"
    json_path <- file.path(file_path,iso3c,"input_params.json")
  }
  json <- jsonlite::read_json(json_path)

  ## get country specific params
  contact_matrix = squire::get_mixing_matrix(country)
  population = squire::get_population(country)

  # get inputs from json
  betas <- unlist(lapply(json, "[[", "beta_set"))
  betas_min <- unlist(lapply(json, "[[", "beta_set_min"))
  betas_max <- unlist(lapply(json, "[[", "beta_set_max"))
  betas <- unlist(lapply(json, "[[", "beta_set"))
  tt_R0 <- unlist(lapply(json, "[[", "tt_beta"))
  dates <- unlist(lapply(json, "[[", "date"))
  deaths <- unlist(lapply(json, "[[", "deaths"))

  # trim to input data
  date_deaths <- unlist(lapply(json, function(x){
    if("deaths" %in% names(x)) {
      x["date"]
    } else {
      NULL
    }
  }))
  Rts <- unlist(lapply(json, "[[", "Rt"))
  Rts <- Rts[which(dates <= max(date_deaths))]
  betas <- betas[which(dates <= max(date_deaths))]
  betas_min <- betas_min[which(dates <= max(date_deaths))]
  betas_max <- betas_max[which(dates <= max(date_deaths))]
  tt_R0 <- tt_R0[which(dates <= max(date_deaths))]
  dates <- dates[which(dates <= max(date_deaths))]

  # get vaccine data if in json (if statement only here as previously this did not exist)
  if("max_vaccine" %in% names(json[[1]])) {
    max_vaccine <- unlist(lapply(json, "[[", "max_vaccine"))
    max_vaccine <- max_vaccine[which(dates <= max(date_deaths))]
  } else {

    # if null use pop size
    if (is.null(max_vaccine)) {
      # default in covidsim in 2.5% of the population to recieve per week so /7 here
      max_vaccine <- as.integer(sum(population$n)*0.025/7)
    }
    max_vaccine <- c(rep(0, length(tt_R0) - length(max_vaccine)), max_vaccine)

  }

  if("vaccine_efficacy_infection" %in% names(json[[1]])) {
    vaccine_efficacy_infection <- lapply(json, "[[", "vaccine_efficacy_infection")
    vaccine_efficacy_infection <- vaccine_efficacy_infection[which(dates <= max(date_deaths))]
  } else {

    # if not use defaults
    vaccine_efficacy_infection <- inf_eff

  }

  if("vaccine_efficacy_disease" %in% names(json[[1]])) {
    vaccine_efficacy_disease <- lapply(json, "[[", "vaccine_efficacy_disease")
    vaccine_efficacy_disease <- vaccine_efficacy_disease[which(dates <= max(date_deaths))]
  } else {

    # if not use defaults
    vaccine_efficacy_disease <- dis_eff

  }

  # future betas based on user inputs
  durs <- nimue:::default_durations()
  probs <- nimue:::default_probs()
  if(length(future_Rt_changes) == 0) {
    future_beta <- nimue::beta_est_infectiousness(dur_IMild = durs$dur_IMild,
                                                  dur_ICase = durs$dur_ICase,
                                                  prob_hosp = probs$prob_hosp,
                                                  rel_infectiousness = rep(1, 17),
                                                  mixing_matrix = squire:::process_contact_matrix_scaled_age(
                                                    squire:::get_mixing_matrix(iso3c = iso3c), population$n),
                                                  R0 = future_Rt)
    future_beta_changes <- future_beta / tail(betas,1)
  } else {
    future_beta_changes <- future_Rt_changes
  }
  if(length(future_beta_changes) != length(tt_Rt_changes)) {
    stop("future_Rt or future_Rt_changes must be same length as tt_Rt_changes")
  }

  new_betas <- c(betas, tail(betas,1) * future_beta_changes)
  tt_s <- c(tt_R0+1, tt_Rt_changes + tail(tt_R0+1,1))

  # future vaccines
  if(length(final_coverage) > 0) {
    current_coverage <- sum(max_vaccine) / sum(squire::population$n[squire::population$iso3c==iso3c])
    final_coverage <- max(final_coverage, current_coverage)
    to_give <- (final_coverage - current_coverage) * sum(squire::population$n[squire::population$iso3c==iso3c])
    to_give <- as.integer(to_give / forecast)
    future_vaccines <- to_give
    tt_future_vaccines <- 1
  }
  new_vaccines <- c(max_vaccine, future_vaccines)
  tt_vacc <- c(tt_R0+1, tt_future_vaccines + tail(tt_R0+1,1))

  # Vaccine strategy
  vacc_json <- paste0("https://github.com/mrc-ide/nimue_js/releases/download/v1.0.10/", iso3c, ".json")
  vacc_strat_json <- jsonlite::read_json(vacc_json)

  # get what was used in the json
  if("vaccine_coverage" %in% names(json[[1]])) {
    vaccine_uptake <- json[[1]]$vaccine_coverage
  }
  if("vaccines_available" %in% names(json[[1]])) {
    vaccine_available <- json[[1]]$vaccines_available
  }
  if("vaccine_strategy" %in% names(json[[1]])) {
    strategy <- json[[1]]$vaccine_strategy
  }

  # get cov_mat for strategy
  if(strategy == "HCW and Elderly") {
    cov_mat <- matrix(unlist(vacc_strat_json$whoPriority), ncol = 17) * vaccine_uptake
  } else if (strategy == "HCW, Elderly and High-Risk") {
    cov_mat <- matrix(unlist(vacc_strat_json$etagePriority), ncol = 17)  * vaccine_uptake
  } else if (strategy == "Elderly") {
    cov_mat <- nimue::strategy_matrix("Elderly", max_coverage = vaccine_uptake, 0)
  } else if (strategy == "All") {
    cov_mat <- nimue::strategy_matrix("All", max_coverage = vaccine_uptake, 0)
  } else {
    stop('Incorrect strategy. Must be one of "HCW and Elderly", "HCW, Elderly and High-Risk", "Elderly", "All"')
  }

  # scale vaccine coverage for availability function
  scale_cov_mat <- function(cov_mat, vaccine_available, pop) {

    # total vaccs available
    tot_vaccines <- sum(pop*vaccine_available)

    # step 1, find when max allocation exceeds capacity
    step <- 1
    step_found <- FALSE
    tot_vaccs_steps <- 0
    cov_mat_dup_ex <- rbind(0, cov_mat)

    while(!step_found && step <= nrow(cov_mat)) {

      if(nrow(cov_mat) == 1) {
        step_found <- TRUE
      }

      vaccs_in_step <- sum((cov_mat_dup_ex[step+1, ] - cov_mat_dup_ex[step, ]) * pop)
      tot_vaccs_steps <- tot_vaccs_steps + vaccs_in_step
      if(tot_vaccs_steps > tot_vaccines) {
        step_found <- TRUE
      } else {
        step <- step+1
      }
    }

    # if we have enough vaccine return now
    if(step > nrow(cov_mat)) {
      return(cov_mat)
    }

    # set steps after max available reached to 0
    if(step < nrow(cov_mat)) {
      cov_mat[(step+1):nrow(cov_mat),] <- 0
    }

    # now set this step to be correct for available
    tots_given <- sum(cov_mat[step-1,] %*% pop)
    tots_tried <- sum(cov_mat[step,] %*% pop)
    remaining <- tot_vaccines - tots_given

    # next_group
    next_group <- cov_mat[step,]-cov_mat[step-1,]
    poss_to_vacc <- (next_group[which(next_group > 0)] * pop[which(next_group > 0)])
    new_cov <- (remaining/sum(poss_to_vacc)) * cov_mat[step, which(next_group > 0)]
    cov_mat[step, which(next_group > 0)] <- new_cov
    return(cov_mat)
  }
  cov_mat <- scale_cov_mat(cov_mat, vaccine_available, population$n)

  # format vaccine efficacies correctly
  for(i in seq_along(vaccine_efficacy_infection)) {
    if(vaccine_efficacy_disease[[i]] < vaccine_efficacy_infection[[i]]) {
      vaccine_efficacy_disease[[i]] <- vaccine_efficacy_infection[[i]]
    }
    vaccine_efficacy_disease[[i]] <- (vaccine_efficacy_disease[[i]] - vaccine_efficacy_infection[[i]]) / (1 - vaccine_efficacy_infection[[i]])
  }


  # NIMUE RUN
  det_out_vac <- nimue::run(
    country = country,
    dur_R = 365,
    use_dde = TRUE,
    # all changes to Rt are set here using beta and not R0 and R0 being set below is ignored internally
    beta_set = new_betas,
    seeding_cases = 5,
    seeding_age_order = 6:10,
    tt_R0 = tt_s,
    R0 = new_betas,
    max_vaccine = new_vaccines,
    tt_vaccine = tt_vacc,
    time_period = length(betas)+forecast,
    dur_V = vaccine_durability,
    vaccine_efficacy_infection = vaccine_efficacy_infection,
    tt_vaccine_efficacy_infection = seq_along(vaccine_efficacy_infection),
    vaccine_efficacy_disease = vaccine_efficacy_disease,
    tt_vaccine_efficacy_disease = seq_along(vaccine_efficacy_disease),
    rel_infectiousness_vaccinated = 0.5,
    vaccine_coverage_mat = cov_mat)

  # get results
  index <- squire:::odin_index(det_out_vac$model)
  D_index <- index$D
  inf_cumu_index <- index$infections_cumu
  hosp_demand_index <- index$hospital_demand
  icu_demand_index <- index$ICU_demand
  vacc_cumu_index <- index$vaccines_cumu

  # build data frame for main plot outputs
  df <- data.frame(date = as.Date(date_deaths), real = deaths)
  df2 <- data.frame(deaths = diff(rowSums(det_out_vac$output[,D_index,1])),
                    infections = diff(rowSums(det_out_vac$output[,inf_cumu_index,1])),
                    hospitilisations = rowSums(det_out_vac$output[-1,hosp_demand_index,1]),
                    critical = rowSums(det_out_vac$output[-1,icu_demand_index,1]))
  df2$date <- seq.Date(as.Date(dates[2]), as.Date(dates[1]) + length(df2$deaths), 1)
  df <- dplyr::left_join(df2, df, by = "date")

  # make simple plot for checking deaths
  plot <- ggplot(df, aes(date, real)) +
    geom_bar(aes(x = as.Date(date), y = real, fill = "Reported"),
             stat = "identity",
             fill = "#c59e96") +
    geom_line(aes(date, zoo::rollmean(real, 7, na.pad = TRUE), color = "7-day Weekly Mean"), lwd = 1) +
    geom_line(aes(date, deaths, color = "Deaths"), lwd = 1) +
    ylab("Deaths") +
    scale_fill_manual(values = "#c59e96") +
    scale_color_manual(values = c("black","#3f8ea7")) +
    xlab("") +
    scale_y_continuous(expand = c(0,0)) +
    ggpubr::theme_pubclean() +
    theme(axis.line = element_line(), legend.title = element_blank(), legend.key = element_blank()) +
    ggtitle(iso3c)

  # And now get the reff out
  # to do this we just need to first take our input beta and interpolate it
  # for all time points explored
  get_reff <- function(out, beta) {

    # mixing_matrix is already the mixing matrix that we pass to you in the country json files
    mixing_matrix <- squire:::process_contact_matrix_scaled_age(
      out$parameters$contact_matrix_set[[1]],
      out$parameters$population
    )

    t_now <- length(beta)

    # these parameters are found in pars_0.json that is imported in index.js
    dur_ICase <- out$parameters$dur_ICase
    dur_IMild <- out$parameters$dur_IMild
    rel_infectiousness <- out$odin_parameters$rel_infectiousness_vaccinated

    # vaccine efficacy is now time changing
    # so we make a list of all the arrays at each time point
    vei_list <- lapply(
      seq_len(nrow(out$odin_parameters$vaccine_efficacy_infection)),
      function(x) {
        out$odin_parameters$vaccine_efficacy_infection[x,,]
      })
    t_vei <- diff(c(out$odin_parameters$tt_vaccine_efficacy_infection, t_now))
    t_vei[1] <- t_vei[1]+1
    vei_list_long <- purrr::flatten(
      lapply(seq_along(t_vei), function(x) {
        rep(list(vei_list[[x]]), t_vei[x])
      }))

    prob_hosp_list <- lapply(
      seq_len(nrow(out$odin_parameters$prob_hosp)),
      function(x) {
        out$odin_parameters$prob_hosp[x,,]
      })
    t_prob_hosp <- diff(c(out$odin_parameters$tt_vaccine_efficacy_disease, t_now))
    t_prob_hosp[1] <- t_prob_hosp[1]+1
    prob_hosp_list_long <- purrr::flatten(
      lapply(seq_along(t_prob_hosp), function(x) {
        rep(list(prob_hosp_list[[x]]), t_prob_hosp[x])
      }))

    index <- nimue:::odin_index(out$model)

    # pop is a 17 length with population sizes in each age category
    pop <- out$parameters$population

    # in here we work out each time point the number of individuals in each age category in
    # the S compartment at each time point.
    susceptible <- array(
      out$output[seq(t_now),index$S,],
      dim=c(t_now, dim(index$S))
    )

    # We divide by the total population
    prop_susc <- sweep(susceptible, 2, pop, FUN='/')

    # We multiply by the effect of vaccines on onward infectiousness at each time point
    prop_susc <- vapply(
      seq_len(nrow(prop_susc)),
      FUN = function(i){ prop_susc[i,,]*vei_list_long[[i]]},
      FUN.VALUE = prop_susc[1,,]
    )

    # back into shape for next step
    prop_susc <- aperm(prop_susc, c(3,1,2))

    # Length 17 with relative R0 in each age category
    relative_R0_by_age <- lapply(prob_hosp_list_long, function(x) {
      x*dur_ICase + (1-x)*dur_IMild
    })

    # here we are looping over each time point to calculate the adjusted eigen
    # incorporating the proportion of the susceptible population in each age group
    adjusted_eigens <- vapply(
      seq(t_now),
      function(t) {
        Re(eigen(mixing_matrix * rowSums(prop_susc[t,,] * relative_R0_by_age[[t]]))$values[1])
      },
      numeric(1)
    )

    # multiply beta by the adjusted eigen at each time point to get Reff
    beta * adjusted_eigens
  }

  get_rts <- function(out, beta) {
  dur_ICase <- out$parameters$dur_ICase
  dur_IMild <- out$parameters$dur_IMild
  prob_hosp <- out$parameters$prob_hosp
  mixing_matrix <- squire:::process_contact_matrix_scaled_age(
    out$parameters$contact_matrix_set[[1]],
    out$parameters$population
  )
  beta * squire::adjusted_eigen(dur_IMild, dur_ICase, prob_hosp, mixing_matrix)
  }

  reff_beta <- approx(x = tt_s, y = new_betas, xout = seq_len(length(betas)+forecast), rule = 2, method = "constant")
  reff <- get_reff(det_out_vac, beta = reff_beta$y)
  rts <- get_rts(det_out_vac, beta = reff_beta$y)

  # again we have everything indexing from t = 1 so remove the initial value and combine with dated data.frame
  df$reff <- reff[-1]
  df$rt <- rts[-1]
  reff_plot <- ggplot(df, aes(date, reff)) +
    geom_line(color = "green") +
    geom_line(aes(y = rt), color = "black") +
    ylab("Reff (green), Rt (black)") +
    scale_color_manual(values = c("green", "black")) +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    ggpubr::theme_pubclean() +
    theme(axis.line = element_line(), legend.title = element_blank(), legend.key = element_blank()) +
    ggtitle(iso3c)

  hosp_plot <- ggplot(df, aes(date, hospitilisations)) +
    geom_line() +
    ylab("Hospitalisations") +
    scale_color_manual(values = c("black")) +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    ggpubr::theme_pubclean() +
    theme(axis.line = element_line(), legend.title = element_blank(), legend.key = element_blank()) +
    ggtitle(iso3c)

  icu_plot <- ggplot(df, aes(date, critical)) +
    geom_line() +
    ylab("Critical Care") +
    scale_color_manual(values = c("black")) +
    xlab("") +
    scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
    ggpubr::theme_pubclean() +
    theme(axis.line = element_line(), legend.title = element_blank(), legend.key = element_blank()) +
    ggtitle(iso3c)


  return(list("plot" = plot, "reff_plot" = reff_plot,
              "hosp_plot" = hosp_plot, "icu_plot" = icu_plot,
              "df" = df, "det_out_vac"=det_out_vac))

}

# ----------------------------
# default plots
# ----------------------------

get_out <- function(iso3c,
                vaccine_uptake = 0.8,
                vaccine_available = 0.95,
                vaccine_durability = 5000,
                future_Rt = numeric(0),
                tt_Rt_changes = numeric(0),
                strategy = "HCW, Elderly and High-Risk") {

  out <- create_vacc_fit(iso3c = iso3c,
                         forecast = as.integer(as.Date("2022-01-01") - Sys.Date())+180,
                         vaccine_uptake = vaccine_uptake,
                         vaccine_available = vaccine_available,
                         vaccine_durability = vaccine_durability,
                         future_Rt = future_Rt,
                         tt_Rt_changes = tt_Rt_changes,
                         strategy = strategy)

  out
}

epi_plot <- function(out) {

  cowplot::plot_grid(out$plot, out$reff_plot, ncol = 1, align = "v")

}

health_plot <- function(out) {

  cowplot::plot_grid(out$icu_plot, out$hosp_plot, ncol = 1, align = "v")

}

# ----------------------------
# UI
# ----------------------------

ui <- shinyUI(
  fluidPage(
    includeCSS("style.css"),
    br(),
    sidebarPanel(width = 3,
                 h3('COVID-19 Model Fits and Forecasts'),
                 br(),
                 uiOutput("country_selection"),
                 br(),
                 uiOutput("rt_selection"),
                 uiOutput("rt_date_selection"),
                 hr(),
                 uiOutput("vaccine_uptake"),
                 uiOutput("vaccine_available"),
                 uiOutput("vaccine_durability"),
                 hr(),
                 uiOutput("strategy"),
                 hr(),
                 actionButton("run", "Run Model")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Epidemic Trajectory",
                 br(),
                 plotOutput('plot1', height = "600px")
        ),
        tabPanel("Healthcare Demand",
                 br(),
                 plotOutput('plot2', height = "600px")
        )
      )
    )
  )
)

# ----------------------------
# SERVER
# ----------------------------

server <- function(input, output) {

  output$country_selection <- renderUI({
    selectInput('country', 'Country', countries,
                selected="Angola")
  })

  output$rt_selection <- renderUI({
    textInput("rt", "Future Rt. Enter future values of Rt separated by ','",
              value = "No Change")
  })

  output$rt_date_selection <- renderUI({
    textInput("rt_dates", "Timing of Future Rt. Enter dates as YYYY-MM-DD for dates of future changes to Rt separated by ','",
              value = "No Change")
  })

  output$vaccine_uptake <- renderUI({
    numericInput('vaccine_uptake', 'Vaccine Uptake (% of group targeted)', 80,
                 min = 10, max = 99)
  })

  output$vaccine_available <- renderUI({
    numericInput('vaccine_available', 'Vaccine Courses Available (% of population):', 95,
                 min = 10, max = 99)
  })

  output$vaccine_durability <- renderUI({
    numericInput('vaccine_durability', 'Vaccine Durability (Days)', 5000,
                 min = 100, max = 10000)
  })

  output$strategy <- renderUI({
    selectInput("strategy", "Prioritisation & Allocation",
                c("HCW and Elderly", "HCW, Elderly and High-Risk","Elderly","All"),
                selected="HCW, Elderly and High-Risk")
  })


  # PLOTS
  res <- eventReactive(input$run, {

    # simple inputs
    iso3c <- iso3cs[match(input$country, countries)]
    vaccine_uptake <- as.numeric(input$vaccine_uptake)/100
    vaccine_available <- as.numeric(input$vaccine_available)/100
    vaccine_durability <- as.numeric(input$vaccine_durability)
    strategy <- as.character(input$strategy)

    # rt inputs
    if(input$rt_dates == "No Change" || input$rt == "No Change") {
      future_Rt = numeric(0)
      tt_Rt_changes = numeric(0)
    } else {
      dates <- trimws(strsplit(input$rt_dates, ",", fixed = TRUE)[[1]])
      dates <- as.Date(dates)
      tt_Rt_changes <- as.integer(dates - Sys.Date())
    }

    if(input$rt_dates == "No Change" || input$rt == "No Change") {
      future_Rt = numeric(0)
      tt_Rt_changes = numeric(0)
    } else {
      rts <- trimws(strsplit(input$rt, ",", fixed = TRUE)[[1]])
      future_Rt <- as.numeric(rts)
    }

    get_out(iso3c = iso3c,
        vaccine_uptake = vaccine_uptake,
        vaccine_available = vaccine_available,
        vaccine_durability = vaccine_durability,
        future_Rt = future_Rt,
        tt_Rt_changes = tt_Rt_changes,
        strategy = strategy)


  })

  output$plot1 <- renderPlot({epi_plot(res())})
  output$plot2 <- renderPlot({health_plot(res())})

}

# ----------------------------
# DEPLOY
# ----------------------------

shinyApp(ui, server)
