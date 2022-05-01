#ifndef CONSTANTS_h
#define CONSTANTS_h
#include <chrono>
#include <vector>
int K;
int R;
double eps;
double Gamma = 0.1;
int BUCKET_SIZE;
int POOL_SIZE;
int nbuckets;
std::string file_name;
bool sampling_maxsat;
int F;
double heparam = 0;

unsigned TIMEOUT = 9990;
unsigned SMALL_TIMEOUT = 1000;
bool verbose = false;
vector<int> b;
int beta = 1000; // some random values
int iter = 100; // some random values
int variant = 1; // 1 for hoa and 2 for prob


std::chrono::high_resolution_clock::time_point start_time;
#endif