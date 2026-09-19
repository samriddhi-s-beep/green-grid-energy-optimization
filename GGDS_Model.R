# ------------------------------------------------------------
# D3M Assignment - GreenGrid Energy Optimization Model
# This script builds and solves MILP models for different datasets
# using HiGHS and GLPK solvers.
# ------------------------------------------------------------
# Loading necessary library
library(magrittr)
library(ompr)
library(ompr.roi)
library(ROI.plugin.glpk)
library(ROI.plugin.highs)

# START LOGGING: using sink() function to get solver message in text file
sink("40529135_solverlogs.txt")

############################################################
# FUNCTION: SOLVE MODEL
GGE_model <- function(file_path) {

# Each row corresponds to a different parameter (as per assignment structure)
  
  ############################################################
  # READ DATA
  ############################################################
  
  lines <- readLines(file_path)
  
  # Helper function to split comma-separated values and convert to numeric
  split_vec <- function(x) as.numeric(strsplit(x, ",")[[1]])
  
  D <- split_vec(lines[1])
  T <- length(D)
  
  alpha <- split_vec(lines[2])
  beta  <- split_vec(lines[3])
  
  m <- split_vec(lines[4])
  S <- length(m)
  
  total_lines <- length(lines)
  A <- total_lines - 4 - S
  
  ############################################################
  # ALLIANCES
  ############################################################
  
  alliances <- list()
  for (k in 1:A) {
    alliances[[k]] <- split_vec(lines[4 + k])
  }
  
  ############################################################
  # COST MATRIX
  ############################################################
  
  C <- matrix(0, nrow = S, ncol = T)
  for (s in 1:S) {
    C[s, ] <- split_vec(lines[4 + A + s])
  }
  
  ############################################################
  # BUILD BASE MODEL (PART B)
  ############################################################
  
  M <- max(D)
  
  model <- MIPModel() %>%
    
    add_variable(x[s, t], s = 1:S, t = 1:T, lb = 0) %>%
    add_variable(y[s], s = 1:S, type = "binary") %>%
    
    set_objective(
      sum_expr(C[s, t] * x[s, t], s = 1:S, t = 1:T),
      sense = "min"
    ) %>%
    
    add_constraint(
      sum_expr(x[s, t], s = 1:S) == D[t],
      t = 1:T
    ) %>%
    
    add_constraint(
      x[s, t] <= alpha[t] * D[t],
      s = 1:S, t = 1:T
    ) %>%
    
    add_constraint(
      sum_expr(x[s, t], t = 1:T) >= m[s] * y[s],
      s = 1:S
    ) %>%
    
    # LINKING CONSTRAINT
    add_constraint(
      x[s, t] <= D[t] * y[s],
      s = 1:S, t = 1:T
    )
  
  ############################################################
  # ALLIANCE CONSTRAINTS
  ############################################################
  
  for (k in 1:A) {
    model <- model %>%
      add_constraint(
        sum_expr(x[s, t], s = alliances[[k]]) <= beta[t] * D[t],
        t = 1:T
      )
  }
  
  ############################################################
  # SOLVE MILP (PART B)
  ############################################################
  
  time_highs <- system.time(
    res_highs <- solve_model(
      model,
      with_ROI(solver = "highs", control = list(time_limit = 3600, verbose = TRUE))
    )
  )
  
  time_glpk <- system.time(
    res_glpk <- solve_model(
      model,
      with_ROI(solver = "glpk", control = list(tm_limit = 3600000, verbose = TRUE))
    )
  )
  
  obj_highs <- objective_value(res_highs)
  obj_glpk  <- objective_value(res_glpk)
  
  runtime_highs <- time_highs["elapsed"]
  runtime_glpk  <- time_glpk["elapsed"]
  
  ############################################################
  # LP RELAXATION (PART B)
  ############################################################
  
  model_lp <- MIPModel() %>%
    
    add_variable(x[s, t], s = 1:S, t = 1:T, lb = 0) %>%
    add_variable(y[s], s = 1:S, lb = 0, ub = 1) %>%
    
    set_objective(
      sum_expr(C[s, t] * x[s, t], s = 1:S, t = 1:T),
      sense = "min"
    ) %>%
    
    add_constraint(sum_expr(x[s, t], s = 1:S) == D[t], t = 1:T) %>%
    
    add_constraint(x[s, t] <= alpha[t] * D[t], s = 1:S, t = 1:T) %>%
    
    add_constraint(sum_expr(x[s, t], t = 1:T) >= m[s] * y[s], s = 1:S) %>%
    
    add_constraint(x[s, t] <= D[t] * y[s], s = 1:S, t = 1:T)
  
  for (k in 1:A) {
    model_lp <- model_lp %>%
      add_constraint(
        sum_expr(x[s, t], s = alliances[[k]]) <= beta[t] * D[t],
        t = 1:T
      )
  }
  
  res_lp <- solve_model(model_lp, with_ROI(solver = "highs", verbose = TRUE))
  obj_lp <- objective_value(res_lp)
  
  ############################################################
  # PART D: FAIRNESS MODEL
  ############################################################
  
  model_d <- model %>%
    
    add_variable(z[t], t = 1:T, lb = 0) %>%
    
    add_constraint(
      z[t] >= C[s, t] - M * (1 - y[s]),
      s = 1:S, t = 1:T
    ) %>%
    
    add_constraint(
      sum_expr(C[s, t] * x[s, t], s = 1:S) <= 1.2 * z[t] * D[t],
      t = 1:T
    )
  
  ############################################################
  # SOLVE PART D
  ############################################################
  
  time_highs_d <- system.time(
    res_highs_d <- solve_model(
      model_d,
      with_ROI(solver = "highs", control = list(time_limit = 3600, verbose = TRUE))
    )
  )
  
  time_glpk_d <- system.time(
    res_glpk_d <- solve_model(
      model_d,
      with_ROI(solver = "glpk", control = list(tm_limit = 3600000, verbose = TRUE))
    )
  )
  
  obj_highs_d <- objective_value(res_highs_d)
  obj_glpk_d  <- objective_value(res_glpk_d)
  
  runtime_highs_d <- time_highs_d["elapsed"]
  runtime_glpk_d  <- time_glpk_d["elapsed"]
  
  # LP RELAXATION (PART D)
  
  model_lp_d <- MIPModel() %>%
    
    add_variable(x[s, t], s = 1:S, t = 1:T, lb = 0) %>%
    add_variable(y[s], s = 1:S, lb = 0, ub = 1) %>%
    add_variable(z[t], t = 1:T, lb = 0) %>%
    
    set_objective(
      sum_expr(C[s, t] * x[s, t], s = 1:S, t = 1:T),
      sense = "min"
    ) %>%
    
    add_constraint(sum_expr(x[s, t], s = 1:S) == D[t], t = 1:T) %>%
    
    add_constraint(x[s, t] <= alpha[t] * D[t], s = 1:S, t = 1:T) %>%
    
    add_constraint(sum_expr(x[s, t], t = 1:T) >= m[s] * y[s], s = 1:S) %>%
    
    add_constraint(x[s, t] <= D[t] * y[s], s = 1:S, t = 1:T) %>%
    
    add_constraint(
      z[t] >= C[s, t] - M * (1 - y[s]),
      s = 1:S, t = 1:T
    ) %>%
    
    add_constraint(
      sum_expr(C[s, t] * x[s, t], s = 1:S) <= 1.2 * z[t] * D[t],
      t = 1:T
    )
  
  for (k in 1:A) {
    model_lp_d <- model_lp_d %>%
      add_constraint(
        sum_expr(x[s, t], s = alliances[[k]]) <= beta[t] * D[t],
        t = 1:T
      )
  }
  
  res_lp_d <- solve_model(model_lp_d, with_ROI(solver = "highs", verbose = TRUE))
  obj_lp_d <- objective_value(res_lp_d)
  # RETURN RESULTS
  
  return(list(
    PART_B = data.frame(
      MILP_HiGHS = obj_highs,
      Runtime_HiGHS = runtime_highs,
      MILP_GLPK = obj_glpk,
      Runtime_GLPK = runtime_glpk,
      LP_Relaxation = obj_lp
    ),
    PART_D = data.frame(
      MILP_HiGHS = obj_highs_d,
      Runtime_HiGHS = runtime_highs_d,
      MILP_GLPK = obj_glpk_d,
      Runtime_GLPK = runtime_glpk_d,
      LP_Relaxation = obj_lp_d
    )
  ))
}

# RUN DATASETS

small_dataset  <- GGE_model("40529135_small.csv")
medium_dataset <- GGE_model("40529135_medium.csv")
large_dataset  <- GGE_model("40529135_large.csv")

sink()

# VIEW RESULTS

small_dataset$PART_B
small_dataset$PART_D
medium_dataset$PART_B
medium_dataset$PART_D
large_dataset$PART_B
large_dataset$PART_D

sessionInfo()







