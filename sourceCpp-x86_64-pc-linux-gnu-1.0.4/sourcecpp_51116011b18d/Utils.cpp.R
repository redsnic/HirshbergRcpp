`.sourceCpp_1_DLLInfo` <- dyn.load('/home/redsnic/Documenti/RStuff/POLICRITI/sourceCpp-x86_64-pc-linux-gnu-1.0.4/sourcecpp_51116011b18d/sourceCpp_2.so')

print_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, alignment) {}, TRUE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_print_alignment')
score_alignment <- Rcpp:::sourceCppFunction(function(s1, s2, is_global, alphabet, score_matrix, gap_cost, alignment) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_score_alignment')
read_fasta <- Rcpp:::sourceCppFunction(function(path) {}, FALSE, `.sourceCpp_1_DLLInfo`, 'sourceCpp_1_read_fasta')

rm(`.sourceCpp_1_DLLInfo`)
