//
// Created by Josh Pughe-Sanford on 1/25/16.
//

#include <cstring>
#include "QuantumWalk.h"

double QuantumWalk::A_r[];
double QuantumWalk::A_i[];
double QuantumWalk::B_r[];
double QuantumWalk::B_i[];
double QuantumWalk::mat[];

QuantumWalk::QuantumWalk(){
    double* IC = this->IC;
    A_r[2] = IC[0];
    A_r[3] = IC[1];
    A_i[2] = IC[2];
    A_i[3] = IC[3];
}
QuantumWalk::QuantumWalk(double ur, double ui, double lr, double li) {
    double* IC = this->IC;
    IC[0] = A_r[2] = ur;
    IC[1] = A_r[3] = lr;
    IC[2] = A_i[2] = ui;
    IC[3] = A_i[3] = li;
}
double powerOfTwo = 0;
bool QuantumWalk::log_cond(double t){
    if(t >= pow(2, powerOfTwo)){
        powerOfTwo++;
        return true;
    } else return false;
}
void QuantumWalk::generate_matrix_array(){
    long iter = this->iterations;
    mat[2] = 0;
    mat[3] = 0;
    for(long i = 0; i < iter+3; i++){
        double theta;
        if (i<1) theta = M_PI_2;
        else theta = ((double) std::rand() / (RAND_MAX))*M_PI_2;
//        else theta = 0; //M_PI_4;
        //std::cout << i << " " <<theta <<" " <<cos(theta) <<" " <<sin(theta) << " ________ " << 4*i+6 << "  " << 4*i+7 << " " << 4*i << " " << 4*i+1 << std::endl;
        mat[4*i+6] = cos(theta);   //index 11 of matrix i
        mat[4*i+7] = sin(theta);   //index 12 of matrix i
        mat[4*i]   =-sin(theta);   //index 21 of matrix i
        mat[4*i+1] = cos(theta);   //index 22 of matrix i
    }
}

void QuantumWalk::flush() {
    memset(A_r, 0, sizeof(A_r));
    memset(A_i, 0, sizeof(A_i));
    memset(B_r, 0, sizeof(B_r));
    memset(B_i, 0, sizeof(B_i));
    double* IC = this->IC;
    A_r[2] = IC[0];
    A_r[3] = IC[1];
    A_i[2] = IC[2];
    A_i[3] = IC[3];
    this->_t = 0;
    powerOfTwo = 0;
    this->log_x.clear();
    this->log_t.clear();
}

void QuantumWalk::run_simulation() {
    long iter = this->iterations;
    long t = this->_t;
    double* mat = this->mat;
    double * src_r;
    double * src_i;
    double * dest_r;
    double * dest_i;
//    if(this->log_cond(t)) {
//        this->generate_statistics(A_r, A_i);
//    }
    for(long j = 0; j<iter; j++) {
//        this->print_norm_array();
//        this->print_amplitude_array();
        if(j%2==0){
            src_r  = A_r; src_i  = A_i;
            dest_r = B_r; dest_i = B_i;
            for(long i=0; i < t+3; i+=2){
                dest_r[2*i]   = (src_r[2*i-2]*mat[4*i+2]+src_r[2*i-1]*mat[4*i+3]);   //real part of upper component
                dest_r[2*i+1] = (src_r[2*i+2]*mat[4*i+4]+src_r[2*i+3]*mat[4*i+5]);   //real part of lower component
                dest_i[2*i]   = (src_i[2*i-2]*mat[4*i+2]+src_i[2*i-1]*mat[4*i+3]);   //imag part of upper component
                dest_i[2*i+1] = (src_i[2*i+2]*mat[4*i+4]+src_i[2*i+3]*mat[4*i+5]);   //imag part of lower component
            }
        } else {
            src_r  = B_r; src_i  = B_i;
            dest_r = A_r; dest_i = A_i;
            for(long i=1; i < t+3; i+=2){
                dest_r[2*i]   = (src_r[2*i-2]*mat[4*i+2]+src_r[2*i-1]*mat[4*i+3]);   //real part of upper component
                dest_r[2*i+1] = (src_r[2*i+2]*mat[4*i+4]+src_r[2*i+3]*mat[4*i+5]);   //real part of lower component
                dest_i[2*i]   = (src_i[2*i-2]*mat[4*i+2]+src_i[2*i-1]*mat[4*i+3]);   //imag part of upper component
                dest_i[2*i+1] = (src_i[2*i+2]*mat[4*i+4]+src_i[2*i+3]*mat[4*i+5]);   //imag part of lower component
            }
        }
        //update variables
        this->_t = t = t + this->_time_step;
        if(this->log_cond(t)) {
//            this->print_norm_array(dest_r, dest_i);
            this->generate_statistics(dest_r, dest_i);
        }
    }
}

void QuantumWalk::generate_statistics(double src_r[], double src_i[]){
    double x = 0.0;
    double y,t;
    double c = 0.0;
    long T = this->_t;
    double pos_step = this->_pos_step;
    double pos_offset =-1;
    double norm;
    double pos;
    long i = 0;
    if(T%2==0) i = 1;
    for(; i<T+3; i+=2) {
        pos = pos_offset+(i*pos_step);
        norm = (src_r[2*i]*src_r[2*i]+src_i[2*i]*src_i[2*i])
               +(src_r[2*i+1]*src_r[2*i+1]+src_i[2*i+1]*src_i[2*i+1]);
        // Kahan Summation Algorithm: reduces rounding errors in low order bits
        y = (norm * pos * pos) - c;
        t = x + y;
        c = (t - x) - y;
        x = t;
    }
    this->log_x.push_back(x);
    this->log_t.push_back(T);
    std::cout << "\tgenerating statistics\t(x2,t)->\t(" << x << ", " << T <<")"<< std::endl;
}

void QuantumWalk::write_statistics(std::ofstream &file){
    std::vector<double> x = this->log_x;
    std::vector<double> t = this->log_t;
    for (std::vector<int>::size_type i = 0; i != t.size(); i++) file << t[i] << ",";
    file << std::endl;
    for(std::vector<int>::size_type i = 0; i != x.size(); i++) file << x[i] << ",";
    file << std::endl;
}

void QuantumWalk::write_statistics(std::ofstream &file, bool log_time){
    std::vector<double> x = this->log_x;
    if(log_time) {
        std::vector<double> t = this->log_t;
        for (std::vector<int>::size_type i = 0; i != t.size(); i++) file << t[i] << ",";
        file << std::endl;
    }
    for(std::vector<int>::size_type i = 0; i != x.size(); i++) file << x[i] << ",";
    file << std::endl;
}

void QuantumWalk::print_statistics(){
    std::vector<double> log_x = this->log_x;
    std::vector<double> log_t = this->log_t;
    for (std::vector<double>::const_iterator i = log_x.begin(); i != log_x.end(); ++i)
        std::cout << *i << ' ';
    std::cout<<std::endl;
    for (std::vector<double>::const_iterator i = log_t.begin(); i != log_t.end(); ++i)
        std::cout << *i << ' ';
}

void QuantumWalk::print_norm_array(double src_r[], double src_i[]){
    long T = this->_t;
    double norm;
    std::cout << "[";
    for(long i = 0; i<T+3; i++) {
        if((T+i)%2 == 0) norm = 0;
        else norm = (src_r[2*i]*src_r[2*i]+src_i[2*i]*src_i[2*i])
                    +(src_r[2*i+1]*src_r[2*i+1]+src_i[2*i+1]*src_i[2*i+1]);
        std::cout << norm << " ";
    }
    std::cout << "]" << std::endl;
}

void QuantumWalk::print_norm_array(){
    long T = this->_t;
    double * src_r; double * src_i; double norm; long i;
    if(T%2==0){ src_r  = A_r; src_i  = A_i; i = 1; } else { src_r  = B_r; src_i  = B_i; i = 0;}
    std::cout << "[";
    for(; i<T+3; i++) {
        if((T+i)%2 == 0) norm = 0;
        else norm = (src_r[2*i]*src_r[2*i]+src_i[2*i]*src_i[2*i])
                    +(src_r[2*i+1]*src_r[2*i+1]+src_i[2*i+1]*src_i[2*i+1]);
        std::cout << norm << " ";
    }
    std::cout << "]" << std::endl;
}

void QuantumWalk::print_amplitude_array(){
    long T = this->_t;
    double * src_r; double * src_i; double norm; long i;
    if(T%2==0){ src_r  = A_r; src_i  = A_i; i = 1; } else { src_r  = B_r; src_i  = B_i; i = 0;}
    std::cout << "[";
    for(; i<T+3; i++) {
        if((T+i)%2 == 0) norm = 0;
        else norm = (src_r[2*i]*src_r[2*i]+src_i[2*i]*src_i[2*i])
                    +(src_r[2*i+1]*src_r[2*i+1]+src_i[2*i+1]*src_i[2*i+1]);
        std::cout << "{ " << src_r[2*i] << "+ i" << src_i[2*i] << ", " << src_r[2*i+1] << "+ i" << src_i[2*i+1] << " } ";
    }
    std::cout << "]" << std::endl;
}

void QuantumWalk::write_norm_array(std::ofstream &file){
    long T = this->_t;
    double * src_r; double * src_i; double norm; long i;
    if(T%2==0){ src_r  = A_r; src_i  = A_i; i = 1; } else { src_r  = B_r; src_i  = B_i; i = 0;}
    for(; i<T+3; i++) {
        if((T+i)%2 == 0) norm = 0;
        else norm = (src_r[2*i]*src_r[2*i]+src_i[2*i]*src_i[2*i])
                    +(src_r[2*i+1]*src_r[2*i+1]+src_i[2*i+1]*src_i[2*i+1]);
        file << norm << "\t";
    }
    file << std::endl;
}
