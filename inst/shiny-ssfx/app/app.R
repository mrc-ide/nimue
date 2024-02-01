
library(cowplot)
library(dplyr)
library(ggplot2)
library(nimue)
library(shiny)
library(squire)
library(tidyr)

# Set up global defaults
# -------------------------------------------------------------------
date_origin <- as.Date("2024-01-21")

date_0 <- as.Date("2024-02-06")
date_1 <- as.Date("2024-03-08")
date_2 <- as.Date("2024-07-10")
date_3 <- as.Date("2024-11-20")
date_4 <- as.Date("2026-01-01")

time_period <- as.integer(date_4 - date_0)
R0 <- 3.5

quick_format <- function(x, var_select, date_0) {

  d <- nimue:::odin_index(x$model)

  do.call(rbind,lapply(var_select, function(i){
    do.call(rbind, lapply(seq_len(dim(x$output)[3]), function(y) {
      df <- data.frame(y = rowSums(x$output[,d[[i]],y]), compartment = i, t = x$output[,"t",y])
      #df$t <- seq_len(nrow(df)) - nrow(df)
      df$replicate <- y
      df$date <- df$t + date_0
      return(df)
    }))
  }))

}

seed_options <- c(12, 50, 80)
tt_R0_options <- c(3, 10, 17)
R0_open_options <- c(1.25, 1.3, 1.4)
tt_vacc_options <- c(150, 164, 171)
max_vaccine_options <- c(150000, 150000, 150000)

# all ui combs
grid <- expand.grid(ui1 = 1:3, ui2 = 1:3, ui3 = 1:3)

# remove those not available
grid <- grid[!(grid$ui1 %in% c(2,3) &  grid$ui2 == 1),]

# set up grid
grid$D1 <- 0
grid$D2 <- 0
grid$D3 <- 0

# current economy
econ_df <- data.frame(date = c(as.Date("2024-01-01"), date_origin),
                      compartment = "ROS Economy",
                      y = c(0, -1))

global <- list()

# default plot
default_plot <- function(r1, econ_df, message, date_0, filter) {
  quick_format(r1, c("D"), date_0 = date_0) %>%
  select(date, y, compartment) %>%
  filter(date <= filter) %>%
  complete(date = seq.Date(as.Date("2024-01-01"), max(date), 1)) %>%
  mutate(y = replace_na(y, 0)) %>%
  mutate(compartment = replace_na(compartment, "D")) %>%
  rbind(
    quick_format(r1, c("infections_cumu"), date_0 = date_0) %>%
      select(date, y, compartment) %>%
      filter(date <= filter) %>%
      complete(date = seq.Date(as.Date("2024-01-01"), max(date), 1)) %>%
      mutate(y = replace_na(y, 0)) %>%
      mutate(compartment = replace_na(compartment, "infections_cumu"))) %>%
  rbind(econ_df) %>%
  mutate(y = round(y)) %>%
  mutate(y = replace(y, compartment == "infections_cumu", round(y[compartment == "infections_cumu"]/80))) %>%
  mutate(compartment = replace(compartment, compartment == "infections_cumu", "SSFX Cases")) %>%
  mutate(compartment = replace(compartment, compartment == "D", "SSFX Deaths")) %>%
  ggplot(aes(x = date, y = y)) +
  geom_line(linewidth = 1.5) +
    ggtitle(filter,
            subtitle = message)+
  ylab("") +
  xlab("Date") +
  facet_wrap(~compartment, scales = "free_y", ncol = 1, switch = "y", strip.position = "left") +
  ggpubr::theme_pubclean(base_size = 18) +
  theme(axis.line = element_line(), legend.title = element_blank(), legend.key = element_blank(), panel.grid.major.x = element_blank(),
        strip.placement = "outside", strip.background = element_rect(fill = NA),
        panel.spacing = unit(1, "cm", data = NULL),
        plot.title = element_text(size = 18, face = "bold"),
        plot.subtitle = element_text(size = 14))
}

#' Run a phase 1 model
#'
#' @param ui_1 UI_1 option
phase_1 <- function(ui_1, econ_df) {

  seed <- seed_options[as.integer(ui_1)]

  r1 <- nimue::run("United Kingdom", R0 = R0, time_period = time_period,
                   seeding_cases = seed, ICU_bed_capacity = 1e5,
                   max_vaccine = 0)

  # Key outputs
  r1$output[,"time",1] <- r1$output[,"time",1] - as.integer(date_1 - date_0)
  r1$output[,"t",1] <- r1$output[,"time",1]
  D1 <- quick_format(r1, "D", date_0 = date_1) %>% filter(date == date_1) %>% pull(y)
  I1 <- quick_format(r1, "infections_cumu", date_0 = date_1) %>% filter(date == date_1) %>% pull(y)

  # Do the economics
  if(ui_1 == 1) {
    econ_df <- rbind(econ_df, data.frame(date = date_1, compartment = "ROS Economy", y = -1))
  } else {
    econ_df <- rbind(econ_df, data.frame(date = date_1, compartment = "ROS Economy", y = 4))
  }

  message <- paste0("\nSSFX Deaths: ", round(D1), "\n\nSSFX Cases: ", round(I1/80),
                    "\n\nROS Economy: ", paste0(tail(econ_df$y,1), "%\n"))
    plot <- default_plot(r1, econ_df, message, date_1, date_1)

  global$r1 <<- r1
  global$econ_df <<- econ_df

  return(list("plot" = plot, "r1" = r1, "econ_df" = econ_df))

}

#' Run a phase 2 model
phase_2 <- function(ui_2) {


  r1 <- global$r1
  econ_df <- global$econ_df

  tt_lock <- tt_R0_options[as.integer(ui_2)]
  r_lock <- 0.8
  r1 <- squire::projections(r1, R0 = c(R0, 1.1, r_lock), tt_R0 = c(0, tt_lock, tt_lock+7))

  # Key outputs
  r1$output[,"time",1] <- r1$output[,"time",1] - as.integer(date_2 - date_1)
  r1$output[,"t",1] <- r1$output[,"time",1]
  D2 <- quick_format(r1, "D", date_0 = date_2) %>% filter(date == date_2) %>% pull(y)
  I2 <- quick_format(r1, "infections_cumu", date_0 = date_2) %>% filter(date == date_2) %>% pull(y)


  # Do the economics
  if(as.integer(ui_2) == 1) {
    econ_df <- rbind(econ_df, data.frame(date = date_2, compartment = "ROS Economy", y = -5))
  } else if (as.integer(ui_2) == 2) {
    econ_df <- rbind(econ_df, data.frame(date = date_2, compartment = "ROS Economy", y = -10))
  } else {
    econ_df <- rbind(econ_df, data.frame(date = date_2, compartment = "ROS Economy", y = -20))
  }

  message <- paste0("\nSSFX Deaths: ", round(D2), "\n\nSSFX Cases: ", round(I2/80),
                    "\n\nEconomy: ", paste0(tail(econ_df$y,1), "%\n"))
  plot <- default_plot(r1, econ_df, message, date_2, date_2)

  global$r2 <<- r1
  global$econ_df_2 <<- econ_df
  global$D2 <<- D2

  return(list("plot" = plot, "r1" = r1, "econ_df" = econ_df, "D2" = D2))

}

#' Run a phase 3 model
phase_3 <- function(ui_1, ui_2, ui_3) {


  r1 <- global$r2
  econ_df <- global$econ_df_2
  D2 <- global$D2

  # Phase 3
  R0_1 <- R0_open_options[as.integer(ui_3)]
  tt_vaccine <- tt_vacc_options[as.integer(ui_3)]
  max_vaccine <- max_vaccine_options[as.integer(ui_3)]

  # did a variant arise
  if(ui_1 >1 && ui_2 == 3) {
    var <- 1.5
  } else {
    var <- 1
  }

  variant_R0 <- R0_1 * var
  r1 <- squire::projections(r1, R0 = c(0.8, R0_1, variant_R0), tt_R0 = c(0, 30, 130),
                            model_user_args = list(list("max_vaccine" = c(0, max_vaccine),
                                                        "tt_vaccine" = c(0, tt_vaccine))))

  # Key outputs
  D3 <- quick_format(r1, "D", date_0 = date_2) %>% filter(date == date_4) %>% pull(y)
  I3 <- quick_format(r1, "infections_cumu", date_0 = date_2) %>% filter(date == date_4) %>% pull(y)



  # Do the economics
  if(as.integer(ui_3) == 1) {
    econ_gain <- 5
  } else if (as.integer(ui_3) == 2) {
    econ_gain <- 15
  } else {
    econ_gain <- 25
  }
  if(D3 > 100000){
    econ_gain <- -1
  }

  econ_df <- rbind(econ_df, data.frame(date = date_4, compartment = "ROS Economy", y = econ_gain))

  message <- paste0("\nSSFX Deaths: ", round(D3), "\n\nSSFX Cases: ", round(I3/80),
                    "\n\nEconomy: ", paste0(tail(econ_df$y,1), "%\n"))
  plot <- default_plot(r1, econ_df, message, date_2, date_4)

  global$D3 <<- D3

  return(list("plot" = plot, "r1" = r1, "econ_df" = econ_df))


}

results <- function(ui_1, ui_2, ui_3) {

  science <- global$D3 < 50000
  economist <- (global$D3 < 25000) | (global$D3 < 100000 & as.integer(ui_1) > 1 & as.integer(ui_3) > 1)
  government <- (global$D3 < 100000 & as.integer(ui_3) > 1)

  textp <- c()
  if(!science) {textp <- c("Scientists")}
  if(!economist) {textp <- c(textp, "Economists")}
  if(!government) {textp <- c(textp, "Government")}
  textp <- paste0("Unfortunately, the following departments did not achieve their aim: \n\n", paste(textp, collapse = ", "))

  if(science && economist && government) {
    list("plot" = ggplot() + ggtitle("Congratulations. The scientists, economists and government all achieved their aim") + theme_void(base_size = 18))
  } else {
    list("plot" = ggplot() + ggtitle(textp) + theme_void(base_size = 18))
  }
}



# ----------------------------
# UI
# ----------------------------

ui <- shinyUI(
  fluidPage(
    includeCSS("style.css"),
    br(),
    sidebarPanel(width = 3,
                 img(src="ssfx_logo.png",height=48,width=120),
                 br(),
                 br(),
                 h3('Government Policy Options'),
                 hr(),
                 uiOutput("dummy"),
                 uiOutput("ui_1"),
                 actionButton("run1", "Implement Phase 1 Policy"),
                 hr(),
                 uiOutput("ui_2"),
                 actionButton("run2", "Implement Phase 2 Policy"),
                 hr(),
                 uiOutput("ui_3"),
                 actionButton("run3", "Implement Phase 3 Policy"),
                 hr(),
                 actionButton("run4", "Result")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Phase 1",
                 br(),
                 plotOutput('plot1', height = "800px")
        ),
        tabPanel("Phase 2",
                 br(),
                 plotOutput('plot2', height = "800px")
        ),
        tabPanel("Phase 3",
                 br(),
                 plotOutput('plot3', height = "800px")
        ),
        tabPanel("Result",
                 br(),
                 plotOutput('plot4', height = "800px")
        )
      )
    )
  )
)

# ----------------------------
# SERVER
# ----------------------------

server <- function(input, output) {

  output$ui_1 <- renderUI({
    radioButtons('ui_1', 'Phase 1', selected = NA, choiceNames = paste("Option", 1:3), choiceValues = 1:3)
  })

  output$ui_2 <- renderUI({
    radioButtons('ui_2', 'Phase 2', selected = NA, choiceNames = paste("Option", 1:3), choiceValues = 1:3)
  })

  output$ui_3 <- renderUI({
    radioButtons('ui_3', 'Phase 3', selected = NA, choiceNames = paste("Option", 1:3), choiceValues = 1:3)
  })

  output$dummy <- renderUI({
    radioButtons('dummy', 'Department:', selected = NA, choices = c("Economy", "Science"))
  })

  # PLOTS
  res1 <- eventReactive(input$run1, {

    phase_1(input$ui_1, econ_df)

  })


  res2 <- eventReactive(input$run2, {

    phase_2(input$ui_2)

  })

  res3 <- eventReactive(input$run3, {

    phase_3(input$ui_1, input$ui_2, input$ui_3)

  })

  res4 <- eventReactive(input$run4, {

    results(input$ui_1, input$ui_2, input$ui_3)

  })

  output$plot1 <- renderPlot({res1()$plot})
  output$plot2 <- renderPlot({res2()$plot})
  output$plot3 <- renderPlot({res3()$plot})
  output$plot4 <- renderPlot({res4()$plot})

}

# ----------------------------
# DEPLOY
# ----------------------------

shinyApp(ui, server)
