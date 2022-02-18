#ifndef CONSTANTS_h
#define CONSTANTS_h
#include <chrono>
static int K;
static int R;
static double eps;
static const double alpha = 0.1;
static int BUCKET_SIZE;
static int POOL_SIZE;
static int nbuckets;
static std::string file_name;
static bool sampling_maxsat;
static int F;
static double heparam = 0;
static int clause_policy = 2;

static unsigned TIMEOUT = 9990;
static unsigned SMALL_TIMEOUT = 1000;
static bool verbose = false;

static std::chrono::high_resolution_clock::time_point start_time;
#endif