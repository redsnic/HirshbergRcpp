#ifndef HIRSHBERG
#define HIRSHBERG

#include <Rcpp.h>
#include "Utils.h"
#include "NWSW.h"

using namespace Rcpp;

CharacterVector align_Hirshberg(
        CharacterVector s1, CharacterVector s2,
        bool is_global,
        CharacterVector alphabet,
        NumericMatrix score_matrix,
        double gap_cost);


std::vector<edit_string>* alignment_Hirshberg( 
        vector_reader<int>* nseq1, vector_reader<int>* nseq2, bool is_global,
        std::string* seq_alphabet,
        NumericMatrix score_matrix,
        double gap_cost);

typedef struct {
    std::vector<double>* line;
    std::vector<edit_string>* bktr;
    double score;   
} ALGB_OUT;

#endif