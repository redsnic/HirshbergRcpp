library(Rcpp)
library(purrr)
library(pryr)
library(dplyr)
library(profmem)

sourceCpp("Utils.cpp", cacheDir = getwd())
sourceCpp("NWSW.cpp", cacheDir = getwd())
sourceCpp("Hirshberg.cpp", cacheDir = getwd())

source("RUtils.R")

# ESEMPIO

seq1 <- "ACA"
seq2 <- "CCC"


a_NWSW <- align_NWSW(seq1, seq2, FALSE,
           nucleotides.order,
           DNA.base.matrix,
           -1)

a_HIRSH <- align_Hirshberg(seq1, seq2, FALSE,
           nucleotides.order,
           DNA.base.matrix,
           -1)

print_alignment(seq1, seq2, a_NWSW)

print(
score_alignment(seq1, seq2, FALSE,
           nucleotides.order,
           DNA.base.matrix,
           -1, a_NWSW)
)

print_alignment(seq1, seq2, a_HIRSH)

print(
score_alignment(seq1, seq2, FALSE,
           nucleotides.order,
           DNA.base.matrix,
           -1, a_HIRSH)
)

# BENCHMARKING

# Ns <- c(1000,3000,5000,7000, 10000, 20000, 30000, 50000, 70000, 100000)
# benchs_global <- map(Ns, ~ profile_time(TRUE,. , nrep = ifelse(. <= 20000, 100, 10))) 
# benchs_local <- map(Ns, ~ profile_time(FALSE,. , nrep = ifelse(. <= 20000, 100, 10))) 
# 
# save(benchs_global, benchs_local, file = "benchmark.Rds")

# # TESTING LOCAL 
# print("TESTING LOCAL")
# 
# good <- 0
# bad <- 0
# for(i in 1:100){
#   seq1 <- paste(collapse = "", sample(c("A","C","G","T"), sample(100:1000,1,replace = TRUE), replace=TRUE))
#   seq2 <- paste(collapse = "", sample(c("A","C","G","T"), sample(100:1000,1,replace = TRUE), replace=TRUE))
#   if(test.score(seq1,seq2,FALSE)==TRUE)
#     good<-good+1
#   else
#     bad<-bad+1
# }
# 
# print(good/(good+bad))
# 
# # TESTING GLOBAL
# print("TESTING GLOBAL")
# 
# good <- 0
# bad <- 0
# for(i in 1:100){
#   seq1 <- paste(collapse = "", sample(c("A","C","G","T"), sample(100:1000,1,replace = TRUE), replace=TRUE))
#   seq2 <- paste(collapse = "", sample(c("A","C","G","T"), sample(100:1000,1,replace = TRUE), replace=TRUE))
#   if(test.score(seq1,seq2,TRUE)==TRUE)
#     good<-good+1
#   else
#     bad<-bad+1
# }
# 
# print(good/(good+bad))
