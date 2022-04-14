# info necessary:
# R-Version, Windowsversion
# Ram, CPU, cores used
# Packages
# packages: benchmarkme
# 
##
# takes the number of cores used as input and returns
# detailed session info.
##
session_all <- function(no_of_cores_used = "unknown") {
  
  list(sessioninfo = sessionInfo(),
       CPU_used = benchmarkme::get_cpu(),
       RAM_available = benchmarkme::get_ram(),
       no_of_cores_used = no_of_cores_used,
       no_of_threads_available = parallel::detectCores(logical = TRUE),
       no_of_physical_cores_available = parallel::detectCores(logical = FALSE))
}