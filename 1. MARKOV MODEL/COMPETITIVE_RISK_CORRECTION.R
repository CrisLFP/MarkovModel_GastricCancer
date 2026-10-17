## ------- Competitive risk correction function ------- ##

## competitive risk for negative transitions probabilites

competitive_risk_correction <- function (
    
  transition_array,              ## original transition array
  df_neg_probs,                  ## dfs with negative probs present in transition array
  dead_state_name = "Dead"){     ## name of the dead state
  
  
  ## re-sacled array at the end of the competitive_risk_correction function.
  correct_array = transition_array
  
  ## 1. action when no negative prob
  if (nrow(df_neg_probs) == 0) {
    
    message("There's no negative probabilities in transition array. Good job! Keep moving")
    return(correct_array)
    
  } ## end if()
  
  
  ## 2. action when negative probs exists 
  for(i in 1:nrow(df_neg_probs)){
    cycle     = df_neg_probs$cycle[i]
    org_state = df_neg_probs$inc_state[i]
    
    
    ## 2.1 extract dead probabilities
    p_dead <- correct_array[org_state,        ## strat state at i cycle
                            dead_state_name,  ## dead name
                            cycle]            ## i cycle
    
    ## id fixed probs (exclude stay in one state and death)
    actual_row  <- correct_array[org_state, , cycle]
    state_name <- names(actual_row)
    
    progression_state <- state_name[state_name != org_state &
                                      state_name != dead_state_name &
                                      actual_row > 0]
    
    
    ## 2.3 sum progression probabilities
    p_fixed <- sum(actual_row[progression_state])
    
    
    ## 2.4 apply re-sclaing for competitive risks
    p_survive    <- 1 - p_dead
    p_available  <- pmax(0, p_survive)
    
    scale_factor <- ifelse(p_fixed == 0,
                           0,
                           pmin(1, p_available / p_fixed))
    
    
    ## 2.5 Adjust and reassign scaled probabilities in transsition array
    for (place in progression_state) {
      correct_array[org_state, place, cycle] <- actual_row[place] * scale_factor
      
    } ## end adjust and reassign for() loop
    
    
    ## 2.6 put new scale probability in original place
    p_scale_new <- pmax(0, 1 - sum(correct_array[org_state, progression_state, cycle])- p_dead)
    
    correct_array[org_state, org_state, cycle] <- p_scale_new
    
  } ## end for() loop
  
  
  ## 2.6 function outputs
  message("Complete re-scaling: ", nrow(df_neg_probs), " cycles adjusted.")
  
  return(correct_array) ## return re-scales array 
  
}   ## end competitive_risk_correction function
