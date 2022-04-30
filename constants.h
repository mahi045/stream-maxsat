#ifndef CONSTANTS_h
#define CONSTANTS_h
#include <chrono>
static int K;
static int R;
static double eps;
static double alpha = 0.1;
static int BUCKET_SIZE;
static int POOL_SIZE;
static int nbuckets;
static std::string file_name;
static bool sampling_maxsat;
static int F;
static double heparam = 0;
static int clause_policy = 2;
static bool median_heu = true;
static bool decision_heu = true;
static double npercentile = 0.75;

static unsigned TIMEOUT = 9990;
static unsigned SMALL_TIMEOUT = 1000;
static bool verbose = false;
static bool use_pool = true;
static bool use_hard = true;

static std::chrono::high_resolution_clock::time_point start_time;
static int number_of_no_assignment = 0;
static bool use_fixed_memory = true;
static int total_memory = 4000;
static bool use_filtering_condition = true;
#endif