`.sourceCpp_7_DLLInfo` <- dyn.load('/home/redsnic/git/HirshbergRcpp/sourceCpp-x86_64-pc-linux-gnu-1.0.4/sourcecpp_4bd57fce9c51/sourceCpp_8.so')

print_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, alignment) {}, TRUE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_print_alignment')
score_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost, alignment) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_score_alignment')
read_fasta <- Rcpp:::sourceCppFunction(function(path) {}, FALSE, `.sourceCpp_7_DLLInfo`, 'sourceCpp_7_read_fasta')

rm(`.sourceCpp_7_DLLInfo`)
