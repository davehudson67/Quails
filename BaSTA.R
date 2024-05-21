library(tidyverse)
library(BaSTA)

## load data
quails <- read_csv("Survival.FINAL.csv") 

quails <- quails %>%
  mutate(sex = ifelse(sex == "Male", "m", "f")) %>%
  mutate(line = as.factor(line)) %>%
  mutate(time = time + 1)

quails$ID <- 1:nrow(quails)

max_day <- max(quails$time) + 5
days <- 1:max_day

create_capture_history <- function(df) {
  max_day <- max(df$time) + 5
  days <- 1:max_day
  
  # Initialize the capture history matrix
  capture_history <- data.frame(matrix(0, nrow = nrow(df), ncol = length(days) + 3))
  colnames(capture_history) <- c("ID", "Birth_Day", "Death_Day", as.numeric(days))
  
  for (i in 1:nrow(df)) {
    id <- df$ID[i]
    death_day <- ifelse(df$status[i] == 1, df$time[i], 0)
    capture_history[i, "ID"] <- id
    capture_history[i, "Birth_Day"] <- 1
    capture_history[i, "Death_Day"] <- death_day
    
    if (death_day != 0) {
      capture_history[i, as.character(2:(death_day - 1))] <- 1
    } else {
      capture_history[i, as.character(2:max_day)] <- 1
    }
  }
  
  return(capture_history)
}

CH <- create_capture_history(quails)

Ch <- CH[1:10,]

newData <- DataCheck(Ch, studyStart = 1 , studyEnd = 1093, autofix = rep(0, 7), silent = FALSE)

out <- basta(CH, studyStart = 1, studyEnd = 1093, nsim = 4,
             parallel = TRUE, ncpus = 4, model = "GO", shape = "simple")
