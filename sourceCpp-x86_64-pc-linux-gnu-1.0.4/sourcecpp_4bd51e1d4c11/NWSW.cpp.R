`.sourceCpp_9_DLLInfo` <- dyn.load('/home/redsnic/git/HirshbergRcpp/sourceCpp-x86_64-pc-linux-gnu-1.0.4/sourcecpp_4bd51e1d4c11/sourceCpp_10.so')

align_NWSW <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost) {}, FALSE, `.sourceCpp_9_DLLInfo`, 'sourceCpp_9_align_NWSW')
print_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, alignment) {}, TRUE, `.sourceCpp_9_DLLInfo`, 'sourceCpp_9_print_alignment')
score_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost, alignment) {}, FALSE, `.sourceCpp_9_DLLInfo`, 'sourceCpp_9_score_alignment')
read_fasta <- Rcpp:::sourceCppFunction(function(path) {}, FALSE, `.sourceCpp_9_DLLInfo`, 'sourceCpp_9_read_fasta')

rm(`.sourceCpp_9_DLLInfo`)
