#ifndef CONSTANTS_h
#define CONSTANTS_h
#include <chrono>
int K;
int R;
double eps;
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

std::chrono::high_resolution_clock::time_point start_time;
#endif