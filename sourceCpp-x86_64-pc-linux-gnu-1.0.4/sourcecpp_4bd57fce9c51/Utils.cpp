#include "Utils.h"
using namespace Rcpp;

// stampa l'allineamento in formato leggibile
// [[Rcpp::export]]
void print_alignment(CharacterVector s1, CharacterVector s2,
                       CharacterVector alignment){
    
    // conversione edit string testuale
    std::string* alignment_str = new std::string(as<std::string>(alignment));
    std::vector<edit_string>* algn = new std::vector<edit_string>();
    for(int i = 0; i<alignment_str->length(); i++){
        auto c = alignment_str->at(i);
        if(c == 'i'){
            algn->push_back(edit_string::insertion);    
        }else if(c == 'd'){
            algn->push_back(edit_string::deletion);
        }else{
            algn->push_back(edit_string::match);
        }
    }
    
    std::string* in_seq1 = new std::string(as<std::string>(s1));
    std::string* in_seq2 = new std::string(as<std::string>(s2));
    
    std::vector<char>* seq1 = new std::vector<char>();
    std::vector<char>* printAlgn = new std::vector<char>();
    std::vector<char>* seq2 = new std::vector<char>();
    
    // scorro la edit string e le sequenze con 3 puntatori
    int p1 = 0;
    int p2 = 0;
    int pA = 0;
    
    while(pA < algn->size()){
        edit_string op = algn->at(pA);
        if(op == edit_string::insertion){
            seq1->push_back('_');
            printAlgn->push_back(' ');
            seq2->push_back(in_seq2->at(p2));
            p2++;
        }else if(op == edit_string::deletion){
            seq1->push_back(in_seq1->at(p1));
            printAlgn->push_back(' ');
            seq2->push_back('_');
            p1++;
        }else{ 
            seq1->push_back(in_seq1->at(p1));
            if(in_seq1->at(p1) == in_seq2->at(p2)){
                printAlgn->push_back('|');  
            }else{
                printAlgn->push_back('.');
            }
            seq2->push_back(in_seq2->at(p2));
            p1++;
            p2++;
        }
        pA++;
    }
    
    print_vector(seq1);
    print_vector(printAlgn);
    print_vector(seq2);
    
    // cleanup
    delete(alignment_str);
    delete(algn);
    delete(in_seq1);
    delete(in_seq2);
    delete(seq1);
    delete(printAlgn);
    delete(seq2);
    
}

// calcola in modo indipendente lo score dell'allineamento (testing)
// [[Rcpp::export]]
double score_alignment(CharacterVector s1, CharacterVector s2,
                       bool is_global,
                       CharacterVector alphabet,
                       NumericMatrix score_matrix,
                       double gap_cost, 
                       CharacterVector alignment){
    // conversione edit string
    std::string* alignment_str = new std::string(as<std::string>(alignment));
    std::vector<edit_string>* algn = new std::vector<edit_string>();
    for(int i = 0; i<alignment_str->length(); i++){
        auto c = alignment_str->at(i); 
        if(c == 'i'){
            algn->push_back(edit_string::insertion);    
        }else if(c == 'd'){
            algn->push_back(edit_string::deletion);
        }else{
            algn->push_back(edit_string::match);
        }
    }
    
    // caso base
    if(algn->size() == 0){
        return 0;
    }
    
    std::string* seq1 = new std::string(as<std::string>(s1));
    std::string* seq2 = new std::string(as<std::string>(s2));
    std::string* seq_alphabet = new std::string(as<std::string>(alphabet));
    
    // scorrimento con puntatori multipli
    int p1 = 0;
    int p2 = 0;
    int pA = 0;
    double score = 0;
    
    // gestione dei gap iniziali e finali a costo 0
    // se allineamento locale 
    if(is_global == 0){
        // FWD
        int pA_rev = algn->size()-1;
        edit_string choosen_side = edit_string::match;
        while(pA_rev >= 0 && algn->at(pA_rev) != edit_string::match){
            if(choosen_side != edit_string::match && choosen_side != algn->at(pA_rev)){
                break;
            }
            choosen_side = algn->at(pA_rev);
            algn->pop_back();
            pA_rev--;
        }
        // REV
        choosen_side = edit_string::match;
        while(pA < algn->size() && algn->at(pA) != edit_string::match){
            if(choosen_side != edit_string::match && choosen_side != algn->at(pA)){
                break;
            }
            choosen_side = algn->at(pA);
            if(algn->at(pA) != edit_string::insertion){
                p1++;    
            }else{
                p2++;
            }
            pA++;
        }
        
    } 
    
    // valutazione del punteggio
    while(pA < algn->size()){
        edit_string op = algn->at(pA);
        if(op == edit_string::insertion){
            p2++;
            score+=gap_cost;
        }else if(op == edit_string::deletion){
            p1++;
            score+=gap_cost;
        }else{ 
            score+=score_matrix(encode(seq1->at(p1), seq_alphabet), encode(seq2->at(p2), seq_alphabet));
            p1++;
            p2++;
        }
        pA++;
    }
    
    // cleanup
    delete(alignment_str);
    delete(algn);
    delete(seq1);
    delete(seq2);
    delete(seq_alphabet);
    
    return score;
}


// legge un file fasta e ne
// estrae la sequenza contenuta
// in una stringa
// [[Rcpp::export]]
std::string read_fasta(String path){
    std::ifstream input_file;
    input_file.open(path.get_cstring());
    std::string sequence("");
    if(input_file.is_open()){
        std::string line;
        while(getline(input_file, line)){
            if(line.length() > 0 && line.at(0) != '>'){
                sequence.append(line);
            }
        }
        input_file.close();
    } else {
        std::cout << "Error opening file: " << path.get_cstring() << "\n";
    }
    return sequence;
}

// codifica in modo numerico un carattere dell'alfabeto
// usato per la specifica sequenza da allineare
int encode(char nt, std::string* seq_alphabet){
    for(int i=0; i<seq_alphabet->length(); i++){
        if(nt == seq_alphabet->at(i)){
            return i;
        }
    }
    return -1;
}

// decodifica a carattere l'identificativo numerico associato all'alfabeto
// usato per la specifica sequenza da allineare
char decode(int i, std::string* seq_alphabet){
    return seq_alphabet->at(i);
}

// encode su stringa
std::vector<int>* encode_seq(std::string* seq, std::string* seq_alphabet){
    std::vector<int>* out = new std::vector<int>();
    for(int i = 0; i<seq->length(); i++){
        out->push_back(encode(seq->at(i), seq_alphabet));
    }
    return out;
}

// decode su stringa
std::string* decode_seq(std::vector<int>* seq, std::string* seq_alphabet){
    std::vector<char>* out = new std::vector<char>();
    for(int i = 0; i<seq->size(); i++){
        out->push_back(decode(seq->at(i), seq_alphabet));
    }
    std::string* out_str = new std::string(out->begin(), out->end());
    delete(out);
    return out_str;
}

// trasforma un vettore di allineamento in formato testuale leggibile
// IMPORTANTE: cancella il vettore in input!
CharacterVector readable_alignment(std::vector<edit_string>* algn){
    std::vector<char>* out_str = new std::vector<char>();
    for(int i = 0; i<algn->size(); i++){
        edit_string p = algn->at(i);
        if(p == insertion){
            out_str->push_back('i');
        }else if(p == deletion){ 
            out_str->push_back('d');
        }else{
            out_str->push_back('.');
        }
    }
    // cleanup
    delete(algn);
    return wrap(std::string(out_str->begin(), out_str->end()));
}

// funzioni ausiliarie per l'allineamento locale, indicano quando
// è necessario valutare 0 il costo dei gap

// indica se la posizione sul vector_reader indicata è quella di 
// inizio per la sequenza considerata
bool is_on_borders(int pos, vector_reader<int>* seq){
    return (( seq->is_true_first(pos) && !seq->is_rev())  ||
            ( seq->is_true_last(pos) && seq->is_rev()));
}

// indica se la posizione sul vector_reader indicata è quella di 
// fine per la sequenza considerata
bool is_on_borders_end(int pos, vector_reader<int>* seq){
    return (( seq->is_true_last(pos) && !seq->is_rev())  ||
            ( seq->is_true_first(pos) && seq->is_rev()));
}




#include <Rcpp.h>
// print_alignment
void print_alignment(CharacterVector s1, CharacterVector s2, CharacterVector alignment);
RcppExport SEXP sourceCpp_7_print_alignment(SEXP s1SEXP, SEXP s2SEXP, SEXP alignmentSEXP) {
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
RcppExport SEXP sourceCpp_7_score_alignment(SEXP s1SEXP, SEXP s2SEXP, SEXP is_globalSEXP, SEXP alphabetSEXP, SEXP score_matrixSEXP, SEXP gap_costSEXP, SEXP alignmentSEXP) {
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
RcppExport SEXP sourceCpp_7_read_fasta(SEXP pathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type path(pathSEXP);
    rcpp_result_gen = Rcpp::wrap(read_fasta(path));
    return rcpp_result_gen;
END_RCPP
}
