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

const long iter = 32768; //1048576; //32768;

class PersistantRandomWalk {
public:
    PersistantRandomWalk();
    PersistantRandomWalk(double u, double l);
    unsigned long iterations = iter;
    void run_simulation(std::ofstream & file);
    void generate_statistics(std::ofstream & file);
    void generate_matrix_array();
//    void generate_ultrametric_matrix_array(double e);
    void flush();
    std::vector<double> log_x;
    std::vector<double> log_t;
    void print_norm_array();
    void write_time(std::ofstream & file, bool endl);
    void write_norm_array(std::ofstream & file, bool endl);
    void write_mds(std::ofstream & file, bool endl);
    double norm_sum();
private:
    long _t = 0;
    double _time_step = 1;
    double _pos_step = 1;
    bool log_cond(double t);
    double IC[2] = {1/sqrt(2),1/sqrt(2)};
    static double A_r [4*(iter+1)];
    static double B_r [4*(iter+1)];
    static double mat [4*(iter+1)];
};


#endif //QUANTUMDRIFT_QUANTUMWALK_H
