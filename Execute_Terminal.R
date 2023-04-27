############################################################################
## This script implement the designated R file with parallel computation
###############################################################################
## For terminal using Args
## Check lines 39-42 & 57 for progress bar set-up
## Mute above lines to remove the progress bars
## Check line 50 for core usage - make changes if needed
###################### Get Args ##############################
args = commandArgs(trailingOnly = TRUE)
start_condition = args[1]
end_condition = args[2]

###################### Install and Load packages ############################## 
.libPaths()
dir.create("dependencies")# Dependencies
.libPaths(c("dependencies", .libPaths()))
install.packages("pacman", repos='http://ftp.ussg.iu.edu/CRAN/')
pacman::p_load(dplyr, tidyselect, tidyr, mirt, lavaan, stringi, foreach, Hmisc, doParallel, progress)

######################### Setup directories ########################## 
path = setwd(".")
project.dir = paste(path, '/PolyLD', sep = '') # Project directory
records.dir = paste(project.dir, '/Records/', sep="")# Record directory for replication records
output.dir = paste(project.dir, '/Output/', sep = "") # Output directory for condition-specific outputs

foreach (dir = c(project.dir, records.dir, output.dir)) %do% dir.create(dir)
foreach (dir = c(project.dir, records.dir, output.dir)) %do% list.files(dir)

###############################################################################
###############################################################################
## Read the script file
## Setup the correct folder path before running the script
script = paste(project.dir, "/_Final_Condition&Function.R", sep = "")
source(script)

iteration = 500
condct = end_condition - start_condition + 1

pbi = progress_bar$new(total = condct, 
                       clear = F, width = 100,
                       format = "(:spin) [:bar] :percent 
                         [Elapsed: :elapsedfull || Estimated remaining: :eta]")

##########################################
############## Parallel computing ########
##########################################
#### Load the package
library(doParallel)
#detectCores() # 8 cores
core = 10 # how many cores you want to use to parallel
#### obtain analysis results using parallel computing ####
############
CLUSTER = makeCluster(core, type = "FORK") # establish clusters
registerDoParallel(CLUSTER)
strtime = Sys.time()

foreach (cd = start_condition : end_condition) %dopar% {
  pbi$tick()# Adding a progress bar
  Results(cd, iteration)
}

endtime = Sys.time()
difftime(endtime,strtime)
############
stopCluster(CLUSTER) # free-up clusters

### There should be #.iteration records in condition file and #.condition files in record folder ###

######## Identify cd and it for results #########




