#' Implement Drop-out Mechanism 
#' 
#' 
#'
#' Pass the dataframe of data and implement different types of drop-out to yield missing data
#' that corresponds to different missing data categories
#' TODO: Add checks, like if 'cmcar' is selected, and there's no covariate, kick error
#' cannot implement dropout at baseline currently - should that be an option? 
#' @param data pass a dataframe generated using function `sim_dat()`
#' @param type_dropout pass character or character vector specifying the type of drop-out you want, 
#' options are c('mcar', 'cmcar', 'mar', 'mnar')
#' @param prop.miss  proportion missing, either a numeric scalar or vector, if vector, must be equal to number of timepoints
#' if it's a scalar it will be the proportion at the last timepoint, amount at other timepoints will be equally spread out
#' @param stochastic.component  control how deterministic drop-out is by adding noise to the governing parameter; 
#' this is just the variance of a standard normal
#' @return returns a dataframe containing the simulated data, with the original scores and the new scores with missingness 
#' @export
#' 




dropout <- function(dat, 
                    type_dropout = NULL, #c('mcar', 'cmcar', 'mar', 'mnar'),
                    prop.miss = 0.3,  # numeric scalar or vector, if vector, must be equal to number of timepoints
                    stochastic.component = 0.2 # control how deterministic drop-out is by adding noise to the governing parameter
){
  
  
  if (is.null(type_dropout)) stop('Specify the type of dropout: `mcar`, `cmcar`, `mar`, `mnar` ')
  
  N <- length(unique(dat$USUBJID))
  number.timepoints <- length(unique(dat$Time))
  
  if (length(prop.miss) == 1) { 
    
    prop.miss <- seq(0, prop.miss, length.out = number.timepoints)
    
  } else {
    
    if (length(prop.miss) != number.timepoints) stop('Vector of proportion missing does not match number of timepoints')
    
  }
  
#----------------------------------------------------------------------------------------------------------- 

  for (score in type_dropout) {
    
    score <- paste0('Y_', score)
    # Add the variable:
    dat$Y_new <- dat$Y_comp
    colnames(dat)[colnames(dat) == 'Y_new'] <- score
    
    
    ###############
    # Create Drop-out:
    #######################
    # Theta is the parameter governing drop-out     
    
    for(tt in 2:number.timepoints){ # start at timepoint 2, no drop-out at timepoint 1!
      
      
      # previous drop out?
      prev.score <- dat[ dat$Time == paste0('Time_', tt - 1) , c('USUBJID', score) ] 
      ### Monotone drop-out: drop out at earlier timepoints -> still dropped out
      ids <- prev.score[ is.na(prev.score[ , score]) , 'USUBJID'] #subject IDs
      ids_t <- dat[ , 'USUBJID'] %in% ids & dat[ , 'Time'] == paste0('Time_', tt) #subject IDs at this timepoint
      dat[ ids_t, score] <- NA
      
      
      # MNAR
      if(score == 'Y_mnar'){
        
        theta <- dat[ dat$Time == paste0('Time_', tt) , c('USUBJID', score)  ] # current score
        colnames(theta)[colnames(theta) == score] <- 'theta'
        
      }
      
      
      # MAR
      if(score == 'Y_mar'){
        
        theta <- dat[ dat$Time == paste0('Time_', tt - 1) , c('USUBJID', score) ]  # previous score
        colnames(theta)[colnames(theta) == score] <- 'theta'

      }
      
      
      
      # Conditional MCAR
      if(score == 'Y_cmcar'){
        
        theta <- dat[ dat$Time == paste0('Time_', tt) , c('USUBJID', 'Covariate', score) , drop = F ] # Covariate
        #drop.cov <- dat[ dat$Time == paste0('Time_', tt) , c('USUBJID', score) , drop = F ]
        drop.cov <- is.na(theta[ , score, drop = T]) 
        theta[ drop.cov, 'Covariate'] <- NA
        theta <- theta[, c('USUBJID', 'Covariate')]
        colnames(theta)[colnames(theta) == 'Covariate'] <- 'theta'
        
      }  
      
      
      # MCAR - fixed 1.28.21
      if(score == 'Y_mcar'){
        
        theta <- dat[ dat$Time == paste0('Time_', tt) , c('USUBJID', score)  ] # current score
        rand.var <-  rnorm(n = length(theta[ , score]), mean = 0, sd = 1)   # replace with a random score
        rand.var[ is.na(theta[ , score]) ] <- NA
        theta[ , score] <- rand.var
        colnames(theta)[colnames(theta) == score] <- 'theta'
        
      }
      
      
      ### Number dropping out at timepoint t
      prop.t <- diff(prop.miss)[tt-1]
      nt <- N*prop.t # number of subjects that will drop out at timepoint t
      
      # Add stochastic component:
      theta[ , 'theta'] <- theta[ , 'theta'] + rnorm(n = nrow(theta), mean = 0, sd = stochastic.component)
      
      # Order theta:
      tmp <- order(theta[ , 'theta'], decreasing = T, na.last = T)
      out <- theta[tmp, 'USUBJID']  # tmp and out overlap by chance! 
      out <- out[1:nt] # pull subject IDs, n = nt
      
      # Make scores NA for subjects in "out" at timepoint "tt" in the dataframe:
      drop.it <- dat[ , 'USUBJID'] %in% out & dat[ , 'Time'] == paste0('Time_', tt)
      dat[ drop.it , score] <- NA
      
      
    }# end loop over timepoints
    
    
  } # end loop over drop-out types
  
  
  
  
  return(dat)
  
} # end Function to implement drop-out


