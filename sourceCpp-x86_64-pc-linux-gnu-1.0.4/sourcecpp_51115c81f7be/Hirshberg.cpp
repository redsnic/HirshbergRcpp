#include "Hirshberg.h"

using namespace Rcpp;

// Wrapper per la chiamata con R
//[[Rcpp::export]]
CharacterVector align_Hirshberg(CharacterVector s1, CharacterVector s2,
                                bool is_global,
                                CharacterVector alphabet,
                                NumericMatrix score_matrix,
                                double gap_cost){
    
    std::string* seq1 = new std::string(as<std::string>(s1));
    std::string* seq2 = new std::string(as<std::string>(s2));
    std::string* seq_alphabet = new std::string(as<std::string>(alphabet));
     
    vector_reader<int>* nseq1 = new vector_reader<int>(encode_seq(seq1, seq_alphabet),0,seq1->size()-1);
    vector_reader<int>* nseq2 = new vector_reader<int>(encode_seq(seq2, seq_alphabet),0,seq2->size()-1);
    
    CharacterVector out = readable_alignment(
        alignment_Hirshberg(nseq1, nseq2, is_global,
                        seq_alphabet, score_matrix, gap_cost)
    );
    
    // cleanup
    delete(seq1);
    delete(seq2);
    delete(seq_alphabet);
    
    return out;
    
} 

// calcolo dello score secondo gli algoritmi di NW/SW in spazio lineare
// il nome è preso dal paper di Hishberg
ALGB_OUT* algorithm_B(
        vector_reader<int>* nseq1, vector_reader<int>* nseq2, bool is_global,
        std::string* seq_alphabet,
        NumericMatrix score_matrix,
        double gap_cost){
    
    std::vector<double>* line1 = make_vector(nseq1->size()+1, 0.);
    std::vector<double>* line2 = make_vector(nseq1->size()+1, 0.);
    
    std::vector<double>* upper_border = make_vector(nseq1->size()+1, 0.);
    std::vector<double>* left_border = make_vector(nseq2->size()+1, 0.);
    
    std::vector<edit_string>* bktr1 = make_vector(nseq1->size()+1, edit_string::insertion);
    std::vector<edit_string>* bktr2 = make_vector(nseq1->size()+1, edit_string::insertion);
    
    // precalcolo dei bordi della matrice (per semplicità implementativa)
    for(int i = 0; i<upper_border->size(); i++){
        if(is_global == 0 && is_on_borders(0,nseq2)){
            upper_border->at(i) = 0;
        } else {
            upper_border->at(i) = gap_cost*i;    
        }
        line1->at(i) = upper_border->at(i);
        bktr1->at(i) = edit_string::deletion;
    }
    
    for(int i = 0; i<left_border->size(); i++){
        if(is_global == 0 && is_on_borders(0,nseq1)){
            left_border->at(i) = 0;
        } else { 
            left_border->at(i) = gap_cost*i;    
        }
    }
    
    // algoritmo di programmazione dinamica
    for(int depth = 1; depth<=nseq2->size(); depth++){
        
        line2->at(0) = left_border->at(depth);
        bktr2->at(0) = edit_string::insertion;
        
        for(int i=1; i<line1->size(); i++){
        
            int el1 = nseq1->at(i-1);
            int el2 = nseq2->at(depth-1);
            
            // gestione delle diverse direzioni
            double up = 0.;
            double left = 0.;
            double diag = 0.;
            
            if(is_global == 0 && is_on_borders_end(i-1,nseq1)){
                up = line1->at(i);
            }else{
                up = line1->at(i) + gap_cost;
            }
            
            if(is_global == 0 && is_on_borders_end(depth-1,nseq2)){
                left = line2->at(i-1);
            }else{
                left = line2->at(i-1) + gap_cost;
            }
            
            diag = line1->at(i-1) + score_matrix(el1,el2);
            
            if(up > left && up > diag){
                line2->at(i) = up;
                bktr2->at(i) = edit_string::insertion;
            }else if(left>diag){
                line2->at(i) = left;
                bktr2->at(i) = edit_string::deletion;
            }else{
                line2->at(i) = diag;
                bktr2->at(i) = edit_string::match;
            }
        }
        
        // passaggio alla linea successiva della matrice
        // tramite swap
        std::vector<double>* swap = line1;
        line1 = line2;
        line2 = swap;
        
        std::vector<edit_string>* swapbktr = bktr1;
        bktr1 = bktr2;
        bktr2 = swapbktr;
        
        
    }
    
    ALGB_OUT* out = new ALGB_OUT;
    out->line = line1;
    out->bktr = bktr1;
    out->score = line1->at(line1->size()-1);
    // cleanup bordi
    delete(line2);
    delete(bktr2);
    delete(upper_border);
    delete(left_border);
    return out;
     
} 

// esecuzione ricorsiva dell'algoritmo di Hirshberg
std::vector<edit_string>* alignment_Hirshberg( 
        vector_reader<int>* nseq1, vector_reader<int>* nseq2, bool is_global,
        std::string* seq_alphabet,
        NumericMatrix score_matrix,
        double gap_cost){
    
    int l1 = nseq1->size();
    int l2 = nseq2->size();
    
    // casi base
    if(l1 == 0 && l2 == 0){
        return new std::vector<edit_string>();
    }
    if(l1 <= 1 || l2 <= 1){
        return alignment_naive(nseq1, nseq2, is_global,
                               seq_alphabet, score_matrix, gap_cost); 
    }
    
    // calcolo degli score sulle due metà
    int split = (int) ((l2-1)/2);
    
    vector_reader<int>* nseq2_low = new vector_reader<int>(nseq2,0,split);
    vector_reader<int>* nseq2_hi = new vector_reader<int>(nseq2,split+1,nseq2->size()-1);
    
    ALGB_OUT* res_low = algorithm_B(nseq1,nseq2_low,is_global,
                                   seq_alphabet, score_matrix, gap_cost);
    
    nseq1->reverse();
    nseq2_hi->reverse();
    
    ALGB_OUT* res_hi = algorithm_B(nseq1,nseq2_hi,is_global,
                                  seq_alphabet, score_matrix, gap_cost);
    
    nseq1->reverse(); // reset 
    
    // cleanup indici ausiliari
    delete(nseq2_hi);
    delete(nseq2_low);
    
    // join degli allineamenti
    int max_pos = 0;
    int max = INT_MIN;
    
    vector_reader<double>* ans_hi = new vector_reader<double>(res_hi->line,0,res_hi->line->size()-1);
    vector_reader<double>* ans_low = new vector_reader<double>(res_low->line,0,res_low->line->size()-1);
    ans_hi->reverse();
    
    for(int i=0; i<ans_hi->size(); i++){
        int val = ans_low->at(i) + ans_hi->at(i);
        if(val > max){
            max_pos = i-1;
            max = val;
        }    
    }
    
    // cleanup elementi per il join
    delete(ans_hi);
    delete(ans_low);
    delete(res_hi->bktr);
    delete(res_hi->line);
    delete(res_low->bktr);
    delete(res_low->line);
    delete(res_hi);
    delete(res_low);
    
    // ricorsione
    std::vector<edit_string>* rec1 = alignment_Hirshberg(
        new vector_reader<int>(nseq1,0,max_pos),
        new vector_reader<int>(nseq2,0,split),
        is_global,
        seq_alphabet,
        score_matrix,
        gap_cost
    );
    
    std::vector<edit_string>* rec2 = alignment_Hirshberg(
        new vector_reader<int>(nseq1,max_pos+1,nseq1->size()-1),
        new vector_reader<int>(nseq2,split+1,nseq2->size()-1),
        is_global,
        seq_alphabet,
        score_matrix,
        gap_cost
    );

    // concatena i due risultati
    for(int i=0; i<rec2->size();i++)
        rec1->push_back(rec2->at(i));
    
    // cleanup finale
    delete(rec2);
    // vanno bene qui per come operano le chiamate ricorsive
    delete(nseq1);
    delete(nseq2);
    
    return rec1;
    
}



 








#include <Rcpp.h>
// align_Hirshberg
CharacterVector align_Hirshberg(CharacterVector s1, CharacterVector s2, bool is_global, CharacterVector alphabet, NumericMatrix score_matrix, double gap_cost);
RcppExport SEXP sourceCpp_5_align_Hirshberg(SEXP s1SEXP, SEXP s2SEXP, SEXP is_globalSEXP, SEXP alphabetSEXP, SEXP score_matrixSEXP, SEXP gap_costSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type s1(s1SEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type s2(s2SEXP);
    Rcpp::traits::input_parameter< bool >::type is_global(is_globalSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type alphabet(alphabetSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type score_matrix(score_matrixSEXP);
    Rcpp::traits::input_parameter< double >::type gap_cost(gap_costSEXP);
    rcpp_result_gen = Rcpp::wrap(align_Hirshberg(s1, s2, is_global, alphabet, score_matrix, gap_cost));
    return rcpp_result_gen;
END_RCPP
}
// print_alignment
void print_alignment(CharacterVector s1, CharacterVector s2, CharacterVector alignment);
RcppExport SEXP sourceCpp_5_print_alignment(SEXP s1SEXP, SEXP s2SEXP, SEXP alignmentSEXP) {
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
RcppExport SEXP sourceCpp_5_score_alignment(SEXP s1SEXP, SEXP s2SEXP, SEXP is_globalSEXP, SEXP alphabetSEXP, SEXP score_matrixSEXP, SEXP gap_costSEXP, SEXP alignmentSEXP) {
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
RcppExport SEXP sourceCpp_5_read_fasta(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(read_fasta(path));
    return rcpp_result_gen;
END_RCPP
}
// align_NWSW
CharacterVector align_NWSW(CharacterVector s1, CharacterVector s2, bool is_global, CharacterVector alphabet, NumericMatrix score_matrix, double gap_cost);
RcppExport SEXP sourceCpp_5_align_NWSW(SEXP s1SEXP, SEXP s2SEXP, SEXP is_globalSEXP, SEXP alphabetSEXP, SEXP score_matrixSEXP, SEXP gap_costSEXP) {
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
