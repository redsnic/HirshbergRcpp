`.sourceCpp_3_DLLInfo` <- dyn.load('/home/redsnic/Documenti/RStuff/POLICRITI/sourceCpp-x86_64-pc-linux-gnu-1.0.4/sourcecpp_511177a5c94b/sourceCpp_4.so')

align_NWSW <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_align_NWSW')
print_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, alignment) {}, TRUE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_print_alignment')
score_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost, alignment) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_score_alignment')
read_fasta <- Rcpp:::sourceCppFunction(function(path) {}, FALSE, `.sourceCpp_3_DLLInfo`, 'sourceCpp_3_read_fasta')

rm(`.sourceCpp_3_DLLInfo`)
