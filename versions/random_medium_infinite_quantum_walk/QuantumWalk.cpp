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
    long iter = this->iterations;
    A_r[2*iter] = IC[0];
    A_r[2*iter+1] = IC[1];
    A_i[2*iter] = IC[2];
    A_i[2*iter+1] = IC[3];
}
QuantumWalk::QuantumWalk(double ur, double ui, double lr, double li) {
    double* IC = this->IC;
    long iter = this->iterations;
    IC[0] = A_r[2*iter] = ur;
    IC[1] = A_r[2*iter+1] = lr;
    IC[2] = A_i[2*iter] = ui;
    IC[3] = A_i[2*iter+1] = li;
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
    long i;
    double theta;
    for(long X = -iter; X<=iter; X++) {
        i = 2 * (X + iter);
        if (X == 0) theta = M_PI_4;
        else theta = ((double) std::rand() / (RAND_MAX)) * M_PI_2;
        mat[i  ] = sin(theta);   //index 11 and -22 of matrix i
        mat[i+1] = cos(theta);   //index 12 and 21 of matrix i
    }
}

//void QuantumWalk::generate_ultrametric_matrix_array(double e){
//    long iter = this->iterations;
//    long i;
//    double theta;
//    mat[2*iter  ] = sin(M_PI_4);
//    mat[2*iter+1] = cos(M_PI_4);
//    int m = int(log(iter)/log(2));
//    int j, k;
//    for(int i = 1; i<=m; i++){
//        int interval = int(pow(2, i-1));
//        double theta_i = M_PI_4*pow(e, i);
//        j=0;
//        while(interval*(2*j+1)<=iter){
//            k = interval*(2*j+1);
//            mat[2*(iter+k)]   = sin(theta_i);
//            mat[2*(iter+k)+1] = cos(theta_i);
//            mat[2*(iter-k)]   = sin(theta_i);
//            mat[2*(iter-k)+1] = cos(theta_i);
////            std::cout << "X= " << k << " : " << theta_i << std::endl;
//            j++;
//        }
//    }
//}

void QuantumWalk::flush() {
    memset(A_r, 0, sizeof(A_r));
    memset(A_i, 0, sizeof(A_i));
    memset(B_r, 0, sizeof(B_r));
    memset(B_i, 0, sizeof(B_i));
    long iter = this->iterations;
    double* IC = this->IC;
    A_r[2*iter] = IC[0];
    A_r[2*iter+1] = IC[1];
    A_i[2*iter] = IC[2];
    A_i[2*iter+1] = IC[3];
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
    long i;
    for(long j = 0; j<iter; j++) {
//        this->print_norm_array();
//        this->print_amplitude_array();
        if(j%2==0){
            src_r  = A_r; src_i  = A_i;
            dest_r = B_r; dest_i = B_i;
            for(long x=-(t+1); x <= (t+1); x+=2){ // want to go from x=-(t+1) to (t+1) --> -N+i/2 = ±(t+1) --> i = 2N ± 2(t+1)
                i = 2*(iter+x);
                //      psiU_{x} = M(x-1)_11*psiU_{x-1}+M(x-1)_12*psiL_{x-1}+M(x+1)_21*psiU_{x+1}+M(x+1)_22*psiL_{x+1}
                // -->  psi_{i} = psiU(x) = M(i/2+t)*psi_{i-2}+M(i/2+t+1)*psi_{i-1}+M(i/2+t+1)*psi_{i+2}-M(i/2+t)*psi_{i+3}
                dest_r[i]   = src_r[i-2]*mat[i-2]+src_r[i-1]*mat[i-1];   //real part of upper component
                dest_r[i+1] = src_r[i+2]*mat[i+3]-src_r[i+3]*mat[i+2];   //real part of lower component
                dest_i[i]   = src_i[i-2]*mat[i-2]+src_i[i-1]*mat[i-1];   //real part of upper component
                dest_i[i+1] = src_i[i+2]*mat[i+3]-src_i[i+3]*mat[i+2];   //imag part of lower component
            }
        } else {
            src_r  = B_r; src_i  = B_i;
            dest_r = A_r; dest_i = A_i;
            for(long x=-(t+1); x <= (t+1); x+=2){ // want to go from x=-(t+1) to (t+1) --> -N+i/2 = ±(t+1) --> i = 2N ± 2(t+1)
                i = 2*(iter+x);
                dest_r[i]   = src_r[i-2]*mat[i-2]+src_r[i-1]*mat[i-1];   //real part of upper component
                dest_r[i+1] = src_r[i+2]*mat[i+3]-src_r[i+3]*mat[i+2];   //real part of lower component
                dest_i[i]   = src_i[i-2]*mat[i-2]+src_i[i-1]*mat[i-1];   //real part of upper component
                dest_i[i+1] = src_i[i+2]*mat[i+3]-src_i[i+3]*mat[i+2];   //imag part of lower component
            }
        }
        //update variables
        this->_t = t = t + this->_time_step;
        if(this->log_cond(t)) {
            //this->print_norm_array(); //calling this method is very slow and freezes clion for large t due to the lag of printing to console
            std::cout << "\tnsum:\t"<< this->norm_sum() << std::endl;
            this->generate_statistics(dest_r, dest_i);
        }
    }
}

double QuantumWalk::norm_sum(){
    double sum = 0.0;
    double y,t;
    double c = 0.0;
    long T = this->_t;
    long iter = this->iterations;
    double norm;
    double * src_r; double * src_i;
    long X=-T;
    long i;
    if(T%2==0){src_r  = A_r; src_i  = A_i;} else { src_r  = B_r; src_i  = B_i;}
    for(; X <= T; X+=2){
        i = 2*(iter+X);
        norm = (src_r[i]*src_r[i]+src_i[i]*src_i[i])
               +(src_r[i+1]*src_r[i+1]+src_i[i+1]*src_i[i+1]);
        // Kahan Summation Algorithm: reduces rounding errors in low order bits
        y = norm - c;
        t = sum + y;
        c = (t - sum) - y;
        sum = t;
    }
    return sum;
}

void QuantumWalk::generate_statistics(double src_r[], double src_i[]){
    double x = 0.0;
    double y,t;
    double c = 0.0;
    long T = this->_t;
    long iter = this->iterations;
    double norm;
    double pos;
    long X=-(T+1)+1;
    long i;
    for(; X <= (T+1); X+=2){
        i = 2*(iter+X);
        pos = X;
        norm = (src_r[i]*src_r[i]+src_i[i]*src_i[i])
               +(src_r[i+1]*src_r[i+1]+src_i[i+1]*src_i[i+1]);
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

void QuantumWalk::write_time(std::ofstream &file, bool endl){
    std::vector<double> t = this->log_t;
    for (std::vector<int>::size_type i = 0; i != t.size(); i++) file << t[i] << ",";
    if(endl) file << std::endl;
}

void QuantumWalk::write_mds(std::ofstream &file, bool endl){
    std::vector<double> x = this->log_x;
    for(std::vector<int>::size_type i = 0; i != x.size(); i++) file << x[i] << ",";
    if(endl) file << std::endl;
}

void QuantumWalk::write_norm_array(std::ofstream &file, bool endl){
    long T = this->_t;
    long iter = this->iterations;
    double norm;
    double * src_r; double * src_i;
    long i;
    if(T%2==0){src_r  = A_r; src_i  = A_i;} else { src_r  = B_r; src_i  = B_i;}
    for(long X=-T; X < (T+1); X+=2){
        i = 2*(iter+X);
        norm = (src_r[i]*src_r[i]+src_i[i]*src_i[i])+(src_r[i+1]*src_r[i+1]+src_i[i+1]*src_i[i+1]);
        file << norm << ",";
    }
    if(endl) file << std::endl;
}

void QuantumWalk::print_norm_array(){
    long T = this->_t;
    long iter = this->iterations;
    double norm;
    double * src_r; double * src_i;
    long i;
    if(T%2==0){src_r  = A_r; src_i  = A_i;} else { src_r  = B_r; src_i  = B_i;}
    std::cout << "[";
    for(long X=-T; X < (T+1); X+=2){
        i = 2*(iter+X);
        norm = (src_r[i]*src_r[i]+src_i[i]*src_i[i])+(src_r[i+1]*src_r[i+1]+src_i[i+1]*src_i[i+1]);
        std::cout << norm << " ";
    }
    std::cout << "] for t = " << T << std::endl;
}
