
#ifndef UTILS
#define UTILS

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include "vector_reader.h"

using namespace Rcpp;

// componenti edit string
enum edit_string
{
    insertion,  
    deletion,
    match
};
 
std::string read_fasta(String path);

// gestione alfabeto
int encode(char nt, std::string* seq_alphabet);
char decode(int i, std::string* seq_alphabet);
std::vector<int>* encode_seq(std::string* seq, std::string* seq_alphabet);
std::string* decode_seq(std::vector<int>* seq, std::string* seq_alphabet);
CharacterVector readable_alignment(std::vector<edit_string>* algn);

// TEMPLATES (devono essere scritti direttemente nell'header)

// Crea una matrice inizializzata a 0
template <typename T>
std::vector<std::vector<T>*>* make_matrix(int nrow, int ncol, T init){
    std::vector<std::vector<T>*>* v = new  std::vector<std::vector<T>*>();
    for(int i=0; i<nrow; i++){
        v->push_back(new std::vector<T>());
        for(int j=0; j<ncol; j++){
            v->at(i)->push_back(init);
        }
    }
    return v;
}

// Crea una matrice inizializzata a 0 (vettore di vettori)
template <typename T>
std::vector<T>* make_vector(int size, T init){
    std::vector<T>* v = new  std::vector<T>();
    for(int i=0; i<size; i++){
        v->push_back(init);
    }
    return v;
}

// numero di righe di una matrice
template <typename T>
int nrow(std::vector<std::vector<T>*>* M){
    return M->size();
}

// numero di colonne di una matrice
template <typename T>
int ncol(std::vector<std::vector<T>*>* M){
    if (nrow(M) == 0) return 0;
    return M->at(0)->size();
}

// accede alla posizione i,j della matrice M
template <typename T>
T a(std::vector<std::vector<T>*>* M, int i, int j){
    return M->at(i)->at(j);
}

// cancella una matrice
template <typename T>
void del(std::vector<std::vector<T>*>* M){
    for(int i=0; i<M->size(); i++){
        delete(M->at(i));
    }
    delete(M);
}

// stampa a schermo un vettore
template <typename T>
void print_vector(T* v){
    for(int i=0; i<v->size(); i++){
        Rcout << v->at(i) << " "; 
    }
    Rcout << std::endl;
}

// stampa a schermo una matrice
template <typename T>
void print_matrix(T* M){
    for(int i=0; i<nrow(M); i++){
        print_vector(M->at(i));
    }
}

// altre funzioni

double score_alignment(CharacterVector s1, CharacterVector s2,
                       bool is_global,
                       CharacterVector alphabet,
                       NumericMatrix score_matrix,
                       double gap_cost, 
                       CharacterVector alignment);

void print_alignment(CharacterVector s1, CharacterVector s2,
                     CharacterVector alignment);


bool is_on_borders(int pos, vector_reader<int>* seq);
bool is_on_borders_end(int pos, vector_reader<int>* seq);

#endif
