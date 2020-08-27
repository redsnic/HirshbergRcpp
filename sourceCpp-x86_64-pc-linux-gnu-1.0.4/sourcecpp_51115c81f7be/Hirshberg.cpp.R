`.sourceCpp_5_DLLInfo` <- dyn.load('/home/redsnic/Documenti/RStuff/POLICRITI/sourceCpp-x86_64-pc-linux-gnu-1.0.4/sourcecpp_51115c81f7be/sourceCpp_6.so')

align_Hirshberg <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_align_Hirshberg')
print_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, alignment) {}, TRUE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_print_alignment')
score_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost, alignment) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_score_alignment')
read_fasta <- Rcpp:::sourceCppFunction(function(path) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_read_fasta')
align_NWSW <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost) {}, FALSE, `.sourceCpp_5_DLLInfo`, 'sourceCpp_5_align_NWSW')

rm(`.sourceCpp_5_DLLInfo`)
