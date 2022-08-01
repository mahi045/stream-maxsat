#ifndef CONSTANTS_h
#define CONSTANTS_h
#include <chrono>
#include <cstdint>
static int K;
static int R;
static double eps;
static double alpha = 0.1;
static int BUCKET_SIZE;
static uint64_t POOL_SIZE;
static int nbuckets;
static std::string file_name;
static bool sampling_maxsat;
static int F;
static double heparam = 0;
static int clause_policy = 2;
static bool median_heu = true;
static bool decision_heu = true;
static double npercentile = 0.5;
static bool L_1 = true;
// hoa means uniform sampling
static bool hoa = false;
static bool postprocessing = false;
static bool store_all = false;

static unsigned TIMEOUT = 9990;
static unsigned SMALL_TIMEOUT = 1000;
static bool verbose = false;
static bool use_pool = true;

static std::chrono::high_resolution_clock::time_point start_time;
static int number_of_no_assignment = 0;
static bool use_fixed_memory = true;
static int total_memory = 4000;
static int iter = 10;

static bool log_of_beta = true;
static bool random_sat_of_beta = true;
static bool expectation_of_clause = true;
static double fraction_of_memory = 0.1;
static double fraction_of_memory_bucket = 0.02;
static double Gamma = 0.1;
static int beta_factor = 3;
static bool default_variable = false;
#endif