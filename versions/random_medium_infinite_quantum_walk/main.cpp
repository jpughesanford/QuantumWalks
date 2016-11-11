#include "QuantumWalk.h"
#include <sys/time.h>
#include <iostream>
#include <fstream>

/* ASSUMPTIONS:
 * - Coin matrix elements are real numbers -> however, this algorithm can be easily expanded to complex matrices using mat_i and mat_r
 * - Initial Conditions are only non-zero at x=0 -->  not doing so will require a change in algorithm and 2x loss in efficiecy
 */
//todo make arrays 2x shorter, using constance of zeros

std::ofstream mdsfile;
std::ofstream nfile;
static unsigned long long get_timestamp () {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (unsigned long long)now.tv_sec * 1000000;
}
 
int main() {
    //create output log
    time_t rawtime = time(0);
    struct tm * now = localtime(&rawtime);
    char filename[80];
    strftime(filename, 80, "%y_%m_%d_at_%I_%M_%S_%p", now);
    std::string description = "";
    std::cout << "Please enter description:\n> ";
    getline(std::cin, description);
    std::string filepath = "./outlog/randinf_";
    std::cout << filepath + filename + ".csv" << std::endl;
    mdsfile.open((filepath + filename + ".csv").c_str());
    mdsfile << "description\t" << description << std::endl;
    nfile.open((filepath + filename + "_norm.csv").c_str());
    nfile << "description\t" << description << std::endl;
    nfile << "header, t, norm array" << std::endl;

//    //check runtime of simulation
//    unsigned long long t0 = get_timestamp();
//    QuantumWalk quantumWalk;
//    quantumWalk.generate_matrix_array();
//    quantumWalk.run_simulation();
//    unsigned long long t1 = get_timestamp();
//
//    long double secs = (t1 - t0) / 1000000.0L;
//    std::cout << "Sec: " << secs<< std::endl;

//    //write statistics and close file
//    quantumWalk.write_statistics(mdsfile);
//    mdsfile.close();
//    return 0;

//    //check runtime of simulation
//    unsigned long long t0 = get_timestamp();
//    QuantumWalk quantumWalk;
//    quantumWalk.generate_matrix_array();
//    quantumWalk.run_simulation();
//    unsigned long long t1 = get_timestamp();
//
//    long double secs = (t1 - t0) / 1000000.0L;
//    std::cout << "Sec: " << secs<< std::endl;
//
//    //write statistics and close file
//    quantumWalk.write_norm_array(mdsfile, true);
//    mdsfile.close();
//    return 0;

    //routine for logging
    QuantumWalk quantumWalk;
    for(int i=0; i<1; i++){
        std::cout << "run: " << i << std::endl;
        quantumWalk.generate_matrix_array();
        quantumWalk.run_simulation(nfile);
//        if(i==0) quantumWalk.write_statistics(mdsfile, true);
//        else quantumWalk.write_statistics(mdsfile, false);
        if(i==0) quantumWalk.write_time(mdsfile, true);
        quantumWalk.write_mds(mdsfile, true);
        quantumWalk.flush();
    }
    mdsfile.close();
    return 0;
}
