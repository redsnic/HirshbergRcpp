#ifndef NWSW
#define NWSW

#include <Rcpp.h>
#include "Utils.h"
#include "vector_reader.h"

using namespace Rcpp; 

void init_global(std::vector<std::vector<double>*>* alignment_matrix, double gap_cost);
void init_bktr_borders(std::vector<std::vector<edit_string>*>* backtracking_matrix);

CharacterVector align_NWSW(
        CharacterVector s1,
        CharacterVector s2,
        bool is_global, 
        CharacterVector alphabet,
        NumericMatrix score_matrix,
        double gap_cost);

std::vector<edit_string>* alignment_naive( 
        vector_reader<int>* nseq1,
        vector_reader<int>* nseq2,
        bool is_global,
        std::string* seq_alphabet,
        NumericMatrix score_matrix,
        double gap_cost);

#endif 