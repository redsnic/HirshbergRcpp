#ifndef VECTOR_READER
#define VECTOR_READER

#include <Rcpp.h>
using namespace Rcpp;

template <class T> class vector_reader
{
public:
    vector_reader(std::vector<T>*, int, int);
    vector_reader(vector_reader<T>*, int, int);
    int true_position(int);
    bool is_true_last(int pos);
    bool is_true_first(int pos);
    T& at(int);
    int size();
    int offset(int pos);
    void reverse();
    std::vector<T>* get_buffer();
    bool is_rev();
private:
    int begin,end;
    std::vector<T>* buffer;
    bool is_reverse;
};

// indica se la direzione di lettura è rovescia
template <class T>
bool vector_reader<T>::is_rev(){
    return this->is_reverse;
}

// creazione a partire da vettore
template <class T>
vector_reader<T>::vector_reader(std::vector<T>* buffer, int begin, int end){
    this->begin = begin;
    this->end = end;
    this->buffer = buffer;
    this->is_reverse = false; 
} 

// creazione a partire da altro vector_reader
template <class T>
vector_reader<T>::vector_reader(vector_reader<T>* other, int begin, int end){
    if(other->is_reverse == 0){
        this->begin = begin+other->begin;
        this->end = end+other->begin;
    }else{
        //   ..B..e....b...E..
        //        B    E
        this->begin = other->end-end;
        this->end = other->end-begin;
    }
    this->buffer = other->buffer;
    this->is_reverse = other->is_reverse; 
} 

// dimensione considerata da questo vector_reader
template <class T>
int vector_reader<T>::size(){
    return std::max(0, this->end - this->begin + 1);
}

// getter per l'elemento alla posizione indicata
template <class T>
T& vector_reader<T>::at(int pos){
    // con verifica dei limiti
    if(this->true_position(pos) < 0 || this->true_position(pos) >= this->buffer->size()){
        Rcerr << "Out of bound vector_reader BUFFER access at position " << pos << std::endl;
        Rcerr << begin << ":" << end << " with buffer->size() = " <<
            buffer->size() << std::endl;
        throw 1;
    }
    if(is_reverse == 1){
        if(this->end - pos < this->begin || this->size() == 0){ // check bounds
            Rcerr << "Out of bound vector_reader access at position " << pos << std::endl;
            Rcerr << begin << ":" << end << " with buffer->size() = " <<
                buffer->size() << std::endl;
            throw 1;  
        }else{
            return this->buffer->at((this->end)-pos);
        }
        
    }else{
        if(this->begin + pos > this->end || this->size() == 0){ // check bounds
            Rcerr << "Out of bound vector_reader access at position " << pos << std::endl;
            Rcerr << begin << ":" << end << " with buffer->size() = " <<
                buffer->size() << std::endl;
            throw 1;  
        }else{
            return this->buffer->at((this->begin)+pos); 
        }
    }
}

// posizione tradotta sul vettore originale (senza controllo)
template <class T>
int vector_reader<T>::true_position(int pos){
    if(this->is_reverse)
        return this->end - pos;
    return this->begin + pos;
}

// inverte il senso di lettura
template <class T>
void vector_reader<T>::reverse(){
    this->is_reverse = !this->is_reverse;
}

// restituisce il vettore originale letto da questo vector_reader
template <class T>
std::vector<T>* vector_reader<T>::get_buffer(){
    return this->buffer;
}

// indica se la posizione indicata punta in coda al vettore originale
// utilizzato per gestire il matching locale
template <class T>
bool vector_reader<T>::is_true_last(int pos){
    return ( (this->true_position(pos) == this->get_buffer()->size()-1) && this->size()>0 && pos<this->size() );
}

// indica se la posizione indicata punta in testa al vettore originale
template <class T>
bool vector_reader<T>::is_true_first(int pos){
    return (this->true_position(pos) == 0 && this->size()>0 && pos<this->size());
}

// offset sul vettore originale 
template <class T>
int vector_reader<T>::offset(int pos){
    if(this->is_reverse){
        return (this->buffer->size() - this->true_position(pos) - 1); 
    }
    return (this->true_position(pos));
}

#endif