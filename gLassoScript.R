library(glassoFast)
options <- commandArgs(trailingOnly = TRUE)

cov.file <- options[1]
l.scle  <- as.integer(options[2])
out.file <- options[3]


S <- readRDS(cov.file)

lambda <- 1/l.scale
print(lambda)

RHO <- matrix(lambda, nrow = nrow(S), ncol = ncol(S))
diag(RHO ) <- 0

cat('starting glasso\n')
system.time(gl <- glassoFast(S, rho = RHO))

cat('\nprocessing done\n')
#saveRDS(gl, out.file)
