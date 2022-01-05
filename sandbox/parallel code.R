

rm(list = ls())
gc()

# libraries that i use:
library(foreach)
library(doParallel)

# Initialize parallelization
cl<-makeCluster(8, outfile="useless.txt") # How many cores do you have?
print(cl)
registerDoParallel(cl)
getDoParWorkers()
getDoParName()



  # Parallelized code:
  out <- foreach(repl= 1:number.repl, .combine='rbind',
                 .errorhandling='pass',
                 .packages=c('foreach', 'doParallel'), .verbose = TRUE)  %dopar% {


                   # Here you put the main function that you want to parallelize


                 } #end parallelization

  cat('Simulation Replication: ', repl, '\n')








