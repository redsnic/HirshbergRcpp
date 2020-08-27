#include "NWSW.h"
using namespace Rcpp;
 
// prepara la matrice di programmazione dinamica con le informazioni per il backtracking sui 
// bordi (in quanto queste sono ovvie)
void init_bktr_borders(std::vector<std::vector<edit_string>*>* backtracking_matrix){
    for(int i=1; i<nrow(backtracking_matrix); i++)
        backtracking_matrix->at(i)->at(0) = edit_string::insertion;
    for(int j=1; j<ncol(backtracking_matrix); j++)
        backtracking_matrix->at(0)->at(j) = edit_string::deletion;
}

// wrapper per R
// [[Rcpp::export]]
CharacterVector align_NWSW(CharacterVector s1, CharacterVector s2,
                           bool is_global,
                           CharacterVector alphabet,
                           NumericMatrix score_matrix,
                           double gap_cost){

    std::string* seq1 = new std::string(as<std::string>(s1));
    std::string* seq2 = new std::string(as<std::string>(s2));
    std::string* seq_alphabet = new std::string(as<std::string>(alphabet));
    
    vector_reader<int>* nseq1 = new vector_reader<int>(encode_seq(seq1, seq_alphabet),0,seq1->length()-1);
    vector_reader<int>* nseq2 = new vector_reader<int>(encode_seq(seq2, seq_alphabet),0,seq2->length()-1);
    
    CharacterVector out = readable_alignment(
        alignment_naive(nseq1, nseq2, is_global,
                        seq_alphabet, score_matrix, gap_cost)
    );
    
    // cleanup
    delete(seq1);
    delete(seq2);
    delete(nseq1);
    delete(seq_alphabet);
    delete(nseq2);
    
    return(out);

}



/*
 * Schema della disposizione delle sequenze:
 * 
 -- s e q 1  --
 |
 s
 e
 q
 2
 |
*/
// algoritmo di Needelman-Wunsch e Smith-Waterson
std::vector<edit_string>* alignment_naive( vector_reader<int>* nseq1, vector_reader<int>* nseq2, bool is_global,
                                           std::string* seq_alphabet,
                                           NumericMatrix score_matrix,
                                           double gap_cost){
    
    std::vector<std::vector<double>*>* alignment_matrix =
        make_matrix(nseq2->size()+1,nseq1->size()+1, 0.);

    std::vector<std::vector<edit_string>*>* backtracking_matrix =
        make_matrix(nseq2->size()+1,nseq1->size()+1, edit_string::deletion);
    
    
    /* init */
    if(is_global == 1 || !is_on_borders(0,nseq1)){
        for(int i=1; i<nrow(alignment_matrix); i++)
            alignment_matrix->at(i)->at(0) = a(alignment_matrix, i-1, 0) + gap_cost;
    }
    if(is_global == 1 || !is_on_borders(0,nseq2)){
        for(int j=1; j<ncol(alignment_matrix); j++)
            alignment_matrix->at(0)->at(j) = a(alignment_matrix, 0, j-1) + gap_cost;
    }
    
    // print_matrix(alignment_matrix);

    init_bktr_borders(backtracking_matrix);
     
    // PD
    for(int i=1; i<nrow(alignment_matrix); i++){
        for(int j=1; j<ncol(alignment_matrix); j++){
            int el1 = nseq1->at(j-1);
            int el2 = nseq2->at(i-1);
            double up   = a(alignment_matrix, i-1, j) + 
                ((is_global==0 && (is_on_borders_end(j-1,nseq1)))?0:gap_cost);
            // l'inserimento locale non deve penalizzare i gap in coda
            double left = a(alignment_matrix, i, j-1) + 
                ((is_global==0 && (is_on_borders_end(i-1,nseq2)))?0:gap_cost);
            // l'inserimento locale non deve penalizzare i gap in coda
            double diag = a(alignment_matrix, i-1, j-1) + score_matrix(el1,el2);
            if( (up > left) && (up > diag)){ //up
                alignment_matrix->at(i)->at(j) = up;
                backtracking_matrix->at(i)->at(j) = edit_string::insertion;
            }else if(left > diag){ //left
                alignment_matrix->at(i)->at(j) = left;
                backtracking_matrix->at(i)->at(j) = edit_string::deletion;
            }else{ //diag
                alignment_matrix->at(i)->at(j) = diag;
                backtracking_matrix->at(i)->at(j) = edit_string::match;
            }
        }
    }

    
    // BACKTRACKING
    std::vector<edit_string>* out = new std::vector<edit_string>();
    int i = nrow(alignment_matrix)-1;
    int j = ncol(alignment_matrix)-1;

    
    while(i>0 || j>0){
        if(a(backtracking_matrix,i,j) == edit_string::insertion){ //up
            out->push_back(edit_string::insertion);
            i--;
        }else if(a(backtracking_matrix,i,j) == edit_string::deletion){ //left
            out->push_back(edit_string::deletion);
            j--; 
        }else{ //diag
            out->push_back(edit_string::match);
            i--;
            j--;
        }
    }
    
    /* l'output è stato inserito dalla coda alla testa */
    std::reverse(out->begin(),out->end());
    
    // cleanup 
    del(alignment_matrix);
    del(backtracking_matrix);
    
    return out;

}






#include <Rcpp.h>
// align_NWSW
CharacterVector align_NWSW(CharacterVector s1, CharacterVector s2, bool is_global, CharacterVector alphabet, NumericMatrix score_matrix, double gap_cost);
RcppExport SEXP sourceCpp_9_align_NWSW(SEXP s1SEXP, SEXP s2SEXP, SEXP is_globalSEXP, SEXP alphabetSEXP, SEXP score_matrixSEXP, SEXP gap_costSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< bool >::type is_global(is_globalSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type score_matrix(score_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type gap_cost(gap_costSEXP);
    rcpp_result_gen = Rcpp::wrap(align_NWSW(s1, s2, is_global, alphabet, score_matrix, gap_cost));
    return rcpp_result_gen;
END_RCPP
}
// print_alignment
void print_alignment(CharacterVector s1, CharacterVector s2, CharacterVector alignment);
RcppExport SEXP sourceCpp_9_print_alignment(SEXP s1SEXP, SEXP s2SEXP, SEXP alignmentSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type alignment(alignmentSEXP);
    print_alignment(s1, s2, alignment);
    return R_NilValue;
END_RCPP
}
// score_alignment
double score_alignment(CharacterVector s1, CharacterVector s2, bool is_global, CharacterVector alphabet, NumericMatrix score_matrix, double gap_cost, CharacterVector alignment);
RcppExport SEXP sourceCpp_9_score_alignment(SEXP s1SEXP, SEXP s2SEXP, SEXP is_globalSEXP, SEXP alphabetSEXP, SEXP score_matrixSEXP, SEXP gap_costSEXP, SEXP alignmentSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< bool >::type is_global(is_globalSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type score_matrix(score_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type gap_cost(gap_costSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type alignment(alignmentSEXP);
    rcpp_result_gen = Rcpp::wrap(score_alignment(s1, s2, is_global, alphabet, score_matrix, gap_cost, alignment));
    return rcpp_result_gen;
END_RCPP
}
// read_fasta
std::string read_fasta(String path);
RcppExport SEXP sourceCpp_9_read_fasta(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(read_fasta(path));
    return rcpp_result_gen;
END_RCPP
}
