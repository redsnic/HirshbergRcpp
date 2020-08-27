`.sourceCpp_11_DLLInfo` <- dyn.load('/home/redsnic/git/HirshbergRcpp/sourceCpp-x86_64-pc-linux-gnu-1.0.4/sourcecpp_4bd52a17ac9a/sourceCpp_12.so')

align_Hirshberg <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_align_Hirshberg')
print_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, alignment) {}, TRUE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_print_alignment')
score_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost, alignment) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_score_alignment')
read_fasta <- Rcpp:::sourceCppFunction(function(path) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_read_fasta')
align_NWSW <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost) {}, FALSE, `.sourceCpp_11_DLLInfo`, 'sourceCpp_11_align_NWSW')

rm(`.sourceCpp_11_DLLInfo`)
