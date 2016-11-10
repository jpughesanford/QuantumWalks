#include "UltraQuantumWalk.h"
#include <sys/time.h>
#include <iostream>
#include <fstream>

/* ASSUMPTIONS:
 * - Coin matrix elements are real numbers -> however, this algorithm can be easily expanded to complex matrices using mat_i and mat_r
 * - Initial Conditions are only non-zero at x=0 -->  not doing so will require a change in algorithm and 2x loss in efficiecy
 */

std::ofstream outputfile;
static unsigned long long get_timestamp () {
    struct timeval now;
    gettimeofday (&now, NULL);
    return  now.tv_usec + (unsigned long long)now.tv_sec * 1000000;
}

//void generate_matrix(){
//    long iter = 32;
//    double mat[iter+1];
//    int m = log(iter)/log(2);
//    std::cout << "m: " << m << std::endl;
//    std::cout << "0: " << 0 << std::endl;
//    for(int i = 1; i<=m; i++){
//        int interval = pow(2, i-1);
//        int j=0;
//        while(interval*(2*j+1)<=iter){
//            mat[interval*(2*j+1)] = i;
//            j++;
//        }
//    }
//    for(int j=1; j<iter; j++){
//        std::cout << j << ": " << mat[j] << std::endl;
//    }
//}
 
int main() {
    //create output log
    time_t rawtime = time(0);
    struct tm * now = localtime(&rawtime);
    char filename[80];
    strftime(filename, 80, "%y_%m_%d_at_%I_%M_%S_%p", now);
    std::string description = "";
    std::cout << "Please enter description:\n> ";
    getline(std::cin, description);
    std::string filepath = "./outlog/ultra_";
    std::cout << filepath + filename + ".csv" << std::endl;
    outputfile.open ((filepath + filename + ".csv").c_str());
    outputfile << "description\t" << description << std::endl;

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
//    quantumWalk.write_statistics(outputfile);
//    outputfile.close();
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
//    quantumWalk.write_norm_array(outputfile);
//    outputfile.close();
//    return 0;
    //routine for logging
    //todo store norm somehow as a check of accuracy
    //todo account for non infinite axis
    //todo change output methods to n, e, t
    UltraQuantumWalk quantumWalk;

    //for norm array mode at various t values
    double e = .75;
    quantumWalk.generate_matrix_array(e);
    outputfile << "t, norm array" << std::endl;
    for (int i = 1; i <= quantumWalk.iterations; i = i*2) {
        std::cout << "starting t->" << i << "...\t";
        quantumWalk.step_simulation(i);
        outputfile << i << ",";
        quantumWalk.write_norm_array(outputfile, true);
        std::cout << "finished with norm = " << quantumWalk.norm_sum() << std::endl;
    }
    quantumWalk.flush();
    // end norm array mode

//    double e0 = 0;
//    double de = .1;
//    for(double e=e0; e<=1.; e+=de){
//        std::cout << "run: episilon = " << e << std::endl;
//        quantumWalk.generate_matrix_array(e);
//        quantumWalk.run_simulation();
//        if(e==e0) quantumWalk.write_time(outputfile, true);
//        quantumWalk.write_mds(outputfile, false);
//        outputfile << e << "," << quantumWalk.norm_sum() << std::endl;
//        quantumWalk.flush();
//    }
    outputfile.close();
    return 0;
}

