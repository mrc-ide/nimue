
#' Prioritisation matrix
#' 
#' Choose from one of four preset prioritisation matrices.
#'
#' @param strategy Prioritisation strategy
#' \itemize{
#'  \item{"All"}{"Target all age groups"}
#'  \item{"Elderly"}{"Prioritise by most elderly first, working to younger age groups at each step"}
#'  \item{"Working Elderly Children"}{"Prioritise working ages followed by elderly followed by children"}
#'  \item{"Risk Elderly Working Children"}{"Prioritise a sub-population of thoe working age, followed by elderly, followed by remiaining working age followed by children"}  
#' }
#' @param max_coverage Maximum coverage in any one age group
#' @param risk_proportion Proportion of working age individuals to prioritise if "Risk Elderly Working Children" is choosen
#'
#' @return Prioritisation coverage matrix
#' @export
strategy_matrix <- function(strategy, max_coverage = 0.8, risk_proportion = 0){
  if(!strategy %in% c("All", "Elderly", "Working Elderly Children", "Risk Elderly Working Children")){
    stop("Strategy must be one of: All, Elderly, Working Elderly Children, Risk Elderly Working Children")
  }
  if(max_coverage < 0 | max_coverage > 1){
    stop("max_coverage must be between 0 and 1")
  }
  if(risk_proportion < 0 | risk_proportion > 1){
    stop("risk_proportion must be between 0 and 1")
  }

  if(strategy == "All"){
    m <- p_all(max_coverage)
  }
  if(strategy == "Elderly"){
    m <- p_elderly(max_coverage)
  }
  if(strategy == "Working Elderly Children"){
    m <- p_working_elderly_children(max_coverage)
  }
  if(strategy == "Risk Elderly Working Children"){
    m <- p_risk_elderly_working_children(max_coverage, risk_proportion)
  }

  return(m)
}

#' Prioritise all
#'
#' @inheritParams strategy_matrix
#'
#' @return Prioritisation coverage matrix
p_all <- function(max_coverage){
  matrix(max_coverage, ncol = 17, nrow = 1)
}

#' Prioritise Elderly
#'
#' @inheritParams strategy_matrix
#'
#' @return Prioritisation coverage matrix
p_elderly <- function(max_coverage){
  m <- matrix(0, ncol = 17, nrow = 17)
  m[upper.tri(m, diag = TRUE)] <- max_coverage
  m <- t(apply(m, 2, rev))
  return(m)
}

#' Prioritise Working > Elderly > Children
#'
#' @inheritParams strategy_matrix
#'
#' @return Prioritisation coverage matrix
p_working_elderly_children <- function(max_coverage){
  m <- matrix(0, ncol = 17, nrow = 3)
  m[1, 4:12] <- max_coverage
  m[2, 4:17] <- max_coverage
  m[3, ] <- max_coverage
  return(m)
}

#' Prioritise Working age risk group > Elderly > working age > Children
#'
#' @inheritParams strategy_matrix
#'
#' @return Prioritisation coverage matrix
p_risk_elderly_working_children <- function(max_coverage, risk_proportion){
  m <- matrix(0, ncol = 17, nrow = 4)
  m[, 4:12] <- max_coverage * risk_proportion
  m[2, 13:17] <- max_coverage
  m[3, 4:17] <- max_coverage
  m[4,] <- max_coverage
  return(m)
}
