#' Filter Data for Multiple Cancer Sites
#'
#' Filters a dataset to include only records for specified cancer sites.
#'
#' @param data A data frame containing cancer data with a 'cancer_site' column.
#' @param cancer_sites A character vector of cancer site names to filter for.
#'
#' @return A filtered data frame containing only records where cancer_site
#'   matches one of the specified sites.
#'
#' @examples
#' \dontrun{
#' # Filter for lung and breast cancers
#' filtered_data <- filter_multiple_cancers(
#'   data = cancer_data,
#'   cancer_sites = c("Lung", "Breast")
#' )
#' }
#'
#' @export
#' @importFrom dplyr filter
filter_multiple_cancers <- function(data, cancer_sites) {
  filtered_data <-  data %>%
    filter(cancer_site %in% cancer_sites)
  return(filtered_data)
}

#' Calculate Other-Cause Mortality Rates
#'
#' This function calculates mortality rates for causes other than the selected (MCED) cancers.
#'
#' @param all_cause_data A data frame containing all-cause mortality data with
#'   columns: Notes, Age, Gender, Year, Deaths, Population, all_crude_Rate.
#' @param MCED_cause_data A data frame containing MCED cancer-specific mortality
#'   data with columns: Year, Age, Gender, cancer_site, Deaths, Population.
#' @param selected_cancers A character vector specifying which cancer sites to
#'   include in the MCED calculation.
#'
#' @return A data frame with other-cause death rates containing columns:
#'   age, sex, year, other_cause_rate.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters MCED cancer data for selected cancer sites
#'   \item Converts all cause and MCED cancer-specific central death rates to probabilities using the formula
#'         of Rosenberg (2006)
#'   \item Calculates other-cause death probability as the difference
#'   \item Convert other cause death probabilities to rates using inverse formula of Rosenberg (2006).
#' }
#'
#' @examples
#' \dontrun{
#' other_rates <- get_other_cause_mortality(
#'   all_cause_data = mortality_data,
#'   MCED_cause_data = cancer_data,
#'   selected_cancers = c("Lung", "Colorectal", "Breast")
#' )
#' }
#'
#' @export
#' @importFrom dplyr filter select mutate rename group_by summarize left_join if_else first
#' @importFrom magrittr %>%
get_other_cause_mortality <- function(all_cause_data, MCED_cause_data, selected_cancers) {

  # Filter MCED cancers
  filtered_cancer_data <- filter_multiple_cancers(data = MCED_cause_data, cancer_sites = selected_cancers)

  # Adjust all-cause crude rate for ages > 84 using merged population
  all_cause_death_rate <- all_cause_data %>%
    select(Notes, Age, Gender, Year, Deaths, Population, all_crude_Rate) %>%
    mutate(all_crude_Rate = ifelse(Age > 84, (Deaths / Population) * 100000, all_crude_Rate)) %>%
    rename(sex = Gender, age = Age, year = Year) %>%
    mutate(Type = "All-cause death rates") %>%
    filter(!is.na(all_crude_Rate), !is.na(sex), !is.na(year))

  # Convert all-cause death rates to probabilities
  all_cause_death_probability <- all_cause_death_rate %>%
    mutate(all_death_probability = (all_crude_Rate / 100000) / (1 + 0.5 * all_crude_Rate / 100000)) %>%
    mutate(Type = "All-cause death probability")

  # Process MCED cancer death rates
  MCED_cancer_death_rate <- filtered_cancer_data %>%
    select(Year, Age, Gender, cancer_site, Deaths, Population) %>%
    mutate(crude_Rate = (Deaths / Population) * 100000) %>%
    rename(sex = Gender, age = Age, year = Year) %>%
    group_by(year, age, sex) %>%
    summarize(
      MCED_Deaths = sum(Deaths, na.rm = TRUE),
      Population = first(Population),
      MCED_crude_Rate = sum(crude_Rate, na.rm = TRUE),
      .groups = 'drop'
    )

  # Convert MCED cancer death rates to probabilities
  MCED_cause_death_probability <- all_cause_death_rate %>%
    left_join(MCED_cancer_death_rate %>% select(year, age, sex, MCED_crude_Rate), by = c("year", "age", "sex")) %>%
    mutate(MCED_death_probability = (MCED_crude_Rate / 100000) / (1 + 0.5 * all_crude_Rate / 100000)) %>%
    select(age, sex, year, all_crude_Rate, MCED_crude_Rate, MCED_death_probability)

  # Compute probability of death from other causes
  other_cause_death_probability <- all_cause_death_probability %>%
    left_join(MCED_cause_death_probability %>% select(year, age, sex, MCED_death_probability), by = c("year", "age", "sex")) %>%
    mutate(other_cause_probability = if_else(is.na(MCED_death_probability),
                                             all_death_probability,
                                             all_death_probability - MCED_death_probability)) %>%
    select(age, sex, year, all_death_probability, MCED_death_probability, other_cause_probability) %>%
    mutate(Type = "Other-cause death probability")

  # Convert other-cause probabilities to rates
  other_cause_death_rate <- other_cause_death_probability %>%
    mutate(other_cause_rate = (other_cause_probability / (1 - 0.5 * other_cause_probability)) * 100000) %>%
    select(age, sex, year, other_cause_rate) %>%
    mutate(Type = "Other-cause death rate")

  return(other_cause_death_rate)
}

#' Create Other-Cause Death Survival Table
#'
#' Generates a survival table for other-cause mortality (non-MCED cancers) for
#' a specific sex, year, and starting age, including cumulative survival probabilities.
#'
#' @param all_cause_data A data frame containing all-cause mortality data.
#' @param MCED_cause_data A data frame containing MCED cancer-specific mortality data.
#' @param selected_cancers A character vector of cancer sites to include in MCED calculations.
#' @param the_sex A character string specifying the sex ("Male" or "Female").
#' @param the_year A numeric value specifying the year.
#' @param the_starting_age A numeric value specifying the starting age for survival calculations.
#'
#' @return A data frame containing survival probabilities with columns:
#'   age, sex, year, other_cause_rate, Type, surv.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Calls \code{get_other_cause_mortality} to get other-cause death rates
#'   \item Filters for the specified sex and year
#'   \item Calculates cumulative survival probabilities using \code{cumprod}
#'   \item Normalizes survival to 1.0 at the starting age
#' }
#'
#' @examples
#' \dontrun{
#' surv_table <- make_othercause_death_table(
#'   all_cause_data = mortality_data,
#'   MCED_cause_data = cancer_data,
#'   selected_cancers = c("Lung", "Breast"),
#'   the_sex = "Female",
#'   the_year = 2020,
#'   the_starting_age = 50
#' )
#' }
#'
#' @export
#' @importFrom dplyr filter mutate
#' @importFrom magrittr %>%
make_othercause_death_table <- function(all_cause_data, MCED_cause_data, selected_cancers,
                                        the_sex, the_year, the_starting_age) {

  othercause_surv=get_other_cause_mortality(all_cause_data = all_cause_data, MCED_cause_data = MCED_cause_data,
                                            selected_cancers = selected_cancers) %>%
    filter(sex == the_sex, year == the_year)%>%
    mutate(surv=cumprod(1-other_cause_rate/100000))%>%
    mutate(surv=ifelse(age<=the_starting_age, 1, surv/surv[age==the_starting_age]))

  return(othercause_surv)

}

#' Simulate Other-Cause Death Time
#'
#' Simulates the time to death from other causes (non-MCED cancers) based on
#' a survival table, with optional seed setting for reproducibility.
#'
#' @param othercause_death_table A data frame containing survival data with
#'   columns 'age' and 'surv', typically output from \code{make_othercause_death_table}.
#' @param ID A numeric value used as a seed for random number generation.
#'
#' @return A numeric value representing the simulated time to death from other causes.
#'
#' @examples
#' \dontrun{
#' # Simulate without setting seed
#' death_time <- sim_othercause_death(survival_table)
#'
#' # Simulate with reproducible seed
#' death_time <- sim_othercause_death(survival_table, ID = 123)
#' }
#'
#' @export
sim_othercause_death <- function(othercause_death_table, ID = NA) {
  if(!is.na(ID)){
    set.seed(ID)
  }
  the_time=gettime(time = othercause_death_table$age, surv = othercause_death_table$surv)

  return(the_time)

}
