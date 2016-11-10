//
// Created by Josh Pughe-Sanford on 1/25/16.
//

#ifndef QUANTUMDRIFT_QUANTUMWALK_H
#define QUANTUMDRIFT_QUANTUMWALK_H

#include <vector>
#include <cstdlib>
#include <math.h>
#include <iostream>
#include <fstream>

const long iter = 10;

class QuantumWalk {
public:
    QuantumWalk();
    QuantumWalk(double ur, double ui, double lr, double li);
    unsigned long iterations = iter;
    void run_simulation();
    void generate_statistics(double src_r[], double src_i[]);
    void generate_matrix_array();
    void flush();
    std::vector<double> log_x;
    std::vector<double> log_t;
    void print_norm_array(double src_r[], double src_i[]);
    void print_norm_array();
    void print_amplitude_array();
    void write_norm_array(std::ofstream & file);
    void write_statistics(std::ofstream & file);
    void write_statistics(std::ofstream & file, bool log_time);
    void print_statistics();
private:
    long _t = 0;
    double _time_step = 1;
    double _pos_step = 1;
    bool log_cond(double t);
    double IC[4] = {1/sqrt(2),0,0,1/sqrt(2)};
    static double A_r[2*(iter+3)];
    static double A_i [2*(iter+3)];
    static double B_r [2*(iter+3)];
    static double B_i [2*(iter+3)];
    static double mat [4*(iter+3)];
};


#endif //QUANTUMDRIFT_QUANTUMWALK_H
