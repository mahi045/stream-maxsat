/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
 * Open-WBO, Copyright (c) 2013-2017, Ruben Martins, Vasco Manquinho, Ines Lynce
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include "utils/Options.h"
#include "utils/ParseUtils.h"
#include "utils/System.h"
#include <errno.h>
#include <signal.h>
#include <zlib.h>

#include <fstream>
#include <iostream>
#include <map>
#include <stdlib.h>
#include <string>
#include <vector>

#ifdef SIMP
#include "simp/SimpSolver.h"
#else
#include "core/Solver.h"
#endif

#include "MaxSAT.h"
#include "MaxTypes.h"
#include "ParserMaxSAT.h"
#include "constants.h"
#include "sampling.h"
#include "streaming.h"

#define VER1_(x) #x
#define VER_(x) VER1_(x)
#define SATVER VER_(SOLVERNAME)
#define VER VER_(VERSION)

using NSPACE::cpuTime;
using NSPACE::OutOfMemoryException;
using NSPACE::IntOption;
using NSPACE::BoolOption;
using NSPACE::StringOption;
using NSPACE::IntRange;
using NSPACE::parseOptions;
using namespace openwbo;

//=================================================================================================

static MaxSAT *mxsolver;
MaxSATFormula *maxsat_formula;
static void SIGINT_exit(int signum) {
  mxsolver->printAnswer(_UNKNOWN_);
  exit(_UNKNOWN_);
}

//=================================================================================================
#if !defined(_MSC_VER) && !defined(__MINGW32__)
void limitMemory(uint64_t max_mem_mb)
{
// FIXME: OpenBSD does not support RLIMIT_AS. Not sure how well RLIMIT_DATA works instead.
#if defined(__OpenBSD__)
#define RLIMIT_AS RLIMIT_DATA
#endif

    // Set limit on virtual memory:
    if (max_mem_mb != 0){
        rlim_t new_mem_lim = (rlim_t)max_mem_mb * 1024*1024;
        rlimit rl;
        getrlimit(RLIMIT_AS, &rl);
        if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
            rl.rlim_cur = new_mem_lim;
            if (setrlimit(RLIMIT_AS, &rl) == -1)
                printf("c WARNING! Could not set resource limit: Virtual memory.\n");
        }
    }

#if defined(__OpenBSD__)
#undef RLIMIT_AS
#endif
}
#else
void limitMemory(uint64_t /*max_mem_mb*/)
{
    printf("c WARNING! Memory limit not supported on this architecture.\n");
}
#endif


#if !defined(_MSC_VER) && !defined(__MINGW32__)
void limitTime(uint32_t max_cpu_time)
{
    if (max_cpu_time != 0){
        rlimit rl;
        getrlimit(RLIMIT_CPU, &rl);
        if (rl.rlim_max == RLIM_INFINITY || (rlim_t)max_cpu_time < rl.rlim_max){
            rl.rlim_cur = max_cpu_time;
            if (setrlimit(RLIMIT_CPU, &rl) == -1)
                printf("c WARNING! Could not set resource limit: CPU-time.\n");
        }
    }
}
#else
void limitTime(uint32_t /*max_cpu_time*/)
{
    printf("c WARNING! CPU-time limit not supported on this architecture.\n");
}
#endif

//=================================================================================================
// Main:

int main(int argc, char **argv) {
  printf(
      "c\nc Open-WBO:\t a Modular MaxSAT Solver -- based on %s (%s version)\n",
      SATVER, VER);
  printf("c Version:\t September 2018 -- Release: 2.1\n");
  printf("c Authors:\t Ruben Martins, Vasco Manquinho, Ines Lynce\n");
  printf("c Contributors:\t Miguel Neves, Saurabh Joshi, Norbert Manthey, Mikolas Janota\n");
  printf("c Contact:\t open-wbo@sat.inesc-id.pt -- "
         "http://sat.inesc-id.pt/open-wbo/\nc\n");
  try {
    NSPACE::setUsageHelp("c USAGE: %s [options] <input-file>\n\n");

#if defined(__linux__)
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw);
    newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(newcw);
    printf(
        "c WARNING: for repeatability, setting FPU to use double precision\n");
#endif

    // BoolOption printmodel("Open-WBO", "print-model", "Print model.\n", true);

    // StringOption printsoft("Open-WBO", "print-unsat-soft", "Print unsatisfied soft claues in the optimal assignment.\n", NULL);

    // IntOption verbosity("Open-WBO", "verbosity",
    //                     "Verbosity level (0=minimal, 1=more).\n", 0,
    //                     IntRange(0, 1));

    // IntOption cpu_lim("Open-WBO", "cpu-lim",
    //                   "Limit on CPU time allowed in seconds.\n", 0,
    //                   IntRange(0, INT32_MAX));

    // IntOption mem_lim("Open-WBO", "mem-lim",
    //                   "Limit on memory usage in megabytes.\n", 0,
    //                   IntRange(0, INT32_MAX));

    // IntOption algorithm("Open-WBO", "algorithm",
    //                     "Search algorithm "
    //                     "(0=wbo,1=linear-su,2=msu3,3=part-msu3,4=oll,5=best)."
    //                     "\n",
    //                     5, IntRange(0, 5));

    // IntOption partition_strategy("PartMSU3", "partition-strategy",
    //                              "Partition strategy (0=sequential, "
    //                              "1=sequential-sorted, 2=binary)"
    //                              "(only for unsat-based partition algorithms).",
    //                              2, IntRange(0, 2));

    // IntOption graph_type("PartMSU3", "graph-type",
    //                      "Graph type (0=vig, 1=cvig, 2=res) (only for unsat-"
    //                      "based partition algorithms).",
    //                      2, IntRange(0, 2));

    // BoolOption bmo("Open-WBO", "bmo", "BMO search.\n", true);

    // IntOption cardinality("Encodings", "cardinality",
    //                       "Cardinality encoding (0=cardinality networks, "
    //                       "1=totalizer, 2=modulo totalizer).\n",
    //                       1, IntRange(0, 2));

    // IntOption amo("Encodings", "amo", "AMO encoding (0=Ladder).\n", 0,
    //               IntRange(0, 0));

    // IntOption pb("Encodings", "pb", "PB encoding (0=SWC,1=GTE,2=Adder).\n", 1,
    //              IntRange(0, 2));

    // IntOption formula("Open-WBO", "formula",
    //                   "Type of formula (0=WCNF, 1=OPB).\n", 0, IntRange(0, 1));

    // IntOption weight(
    //     "WBO", "weight-strategy",
    //     "Weight strategy (0=none, 1=weight-based, 2=diversity-based).\n", 2,
    //     IntRange(0, 2));

    // BoolOption symmetry("WBO", "symmetry", "Symmetry breaking.\n", true);

    // IntOption symmetry_lim(
    //     "WBO", "symmetry-limit",
    //     "Limit on the number of symmetry breaking clauses.\n", 500000,
    //     IntRange(0, INT32_MAX));
    IntOption R_value("Open-WBO", "Rvalue",
                        "the ratio of bucket size and pool size "
                        "\n",
                        5, IntRange(1, 100));

    IntOption K_value("Open-WBO", "Kvalue",
                        "the value of K "
                        "\n",
                        1, IntRange(1, 100));
    IntOption F_value("Open-WBO", "Fvalue",
                        "the value of F "
                        "\n",
                        2, IntRange(0, 10000));
       IntOption Timeout_value("Open-WBO", "timeout",
                            "the value of timeout "
                            "\n",
                            10, IntRange(1, 20000));
       IntOption Immediate_timeout_value("Open-WBO", "small-timeout",
                            "the value of small timeout "
                            "\n",
                            10, IntRange(1, 10000));
    Glucose::DoubleOption epsilon("Open-WBO", "epsilon",
                        "the value of ep "
                        "\n",
                        0.25);
    Glucose::DoubleOption hp("Open-WBO", "heparam",
                        "the value of heparam "
                        "\n",
                        0.05);
    Glucose::DoubleOption al("Open-WBO", "alpha",
                        "the value of alpha "
                        "\n",
                        0.1);
    Glucose::DoubleOption percentile("Open-WBO", "percentile",
                        "the value of npercentile "
                        "\n",
                        0.75);                    
    IntOption c_p("Open-WBO", "clause_policy",
                        "the selection of clause policy "
                        "\n",
                        2, IntRange(0, 2));                    
    BoolOption sampling("WBO", "sampling", "Symmetry breaking.\n", false);
    BoolOption print_verbose("WBO", "print-verbose", "Printing the verbose.\n", false);
    BoolOption pool_c("WBO", "pool", "use pool.\n", true);
    BoolOption hard_c("WBO", "conflict", "adding the non-conflict variables.\n", true);
    BoolOption decision_c("WBO", "decision", "enable decision heuristic.\n", false);
    BoolOption use_median("WBO", "median", "use median as F.\n", true);
    parseOptions(argc, argv, true);
    R = (int)R_value;
    K = (int)K_value;
    F = (int)F_value;
    eps = (double)epsilon;
    alpha = (double)al;
    clause_policy = (int)c_p; 
    use_hard = (bool) hard_c;
    use_pool = (bool) pool_c;
    decision_heu = (bool) decision_c;

    heparam = (double) hp;
    npercentile = (double) percentile;
    median_heu = (bool) use_median;

    TIMEOUT = (int) Timeout_value;
    SMALL_TIMEOUT = (int) Immediate_timeout_value;

    printf("R_value: %d\n", R);
    printf("K_value: %d\n", K);
    printf("Epsilon: %f\n", (double)epsilon);
//     printf("Streaming: %d\n", (int)sampling);
    printf("Timeout for complete maxsat: %d\n", (int)TIMEOUT);
    printf("Timeout for partition maxsat: %d\n", (int)SMALL_TIMEOUT);

    
    // // Try to set resource limits:
    // if (cpu_lim != 0) limitTime(cpu_lim);
    // if (mem_lim != 0) limitMemory(mem_lim);

    double initial_time = cpuTime();
    MaxSAT *S = NULL;

    // switch ((int)algorithm) {
    // case _ALGORITHM_WBO_:
    //   S = new WBO(verbosity, weight, symmetry, symmetry_lim);
    //   break;

    // case _ALGORITHM_LINEAR_SU_:
    //   S = new LinearSU(verbosity, bmo, cardinality, pb);
    //   break;

    // case _ALGORITHM_PART_MSU3_:
    //   S = new PartMSU3(verbosity, partition_strategy, graph_type, cardinality);
    //   break;

    // case _ALGORITHM_MSU3_:
    //   S = new MSU3(verbosity);
    //   break;

    // case _ALGORITHM_OLL_:
    //   S = new OLL(verbosity, cardinality);
    //   break;

    // case _ALGORITHM_BEST_:
    //   break;

    // default:
    //   printf("c Error: Invalid MaxSAT algorithm.\n");
    //   printf("s UNKNOWN\n");
    //   exit(_ERROR_);
    // }

    // signal(SIGXCPU, SIGINT_exit);
    // signal(SIGTERM, SIGINT_exit);

    if (argc == 1) {
      printf("c Warning: no filename.\n");
    }

    gzFile in = (argc == 1) ? gzdopen(0, "rb") : gzopen(argv[1], "rb");
    if (in == NULL)
      printf("c ERROR! Could not open file: %s\n",
             argc == 1 ? "<stdin>" : argv[1]),
          printf("s UNKNOWN\n"), exit(_ERROR_);

    maxsat_formula = new MaxSATFormula();
    file_name = std::string(argv[1]);   
    sampling_maxsat = sampling;
    verbose = print_verbose;
    // if ((int)formula == _FORMAT_MAXSAT_) {
      parseMaxSATFormula(in, maxsat_formula);
      maxsat_formula->setFormat(_FORMAT_MAXSAT_);
    // } else {
    //   ParserPB *parser_pb = new ParserPB();
    //   parser_pb->parsePBFormula(argv[1], maxsat_formula);
    //   maxsat_formula->setFormat(_FORMAT_PB_);
    // }
    gzclose(in);
    if (sampling) {
       printf("Running sampling version of maxsat!!!\n");
       sample_clauses(maxsat_formula);
    }
    else {
       printf("Running streaming version of maxsat!!!\n");
      //  streaming_maxsat(maxsat_formula);
    }

    printf("c |                                                                "
           "                                       |\n");
    printf("c ========================================[ Problem Statistics "
           "]===========================================\n");
    printf("c |                                                                "
           "                                       |\n");

    if (maxsat_formula->getFormat() == _FORMAT_MAXSAT_)
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n",
          "MaxSAT");
    else
      printf(
          "c |  Problem Format:  %17s                                         "
          "                          |\n",
          "PB");

    if (maxsat_formula->getProblemType() == _UNWEIGHTED_)
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Unweighted");
    else
      printf("c |  Problem Type:  %19s                                         "
             "                          |\n",
             "Weighted");

    printf("c |  Number of variables:  %12d                                    "
           "                               |\n",
           maxsat_formula->nVars());
    printf("c |  Number of hard clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nHard());
    printf("c |  Number of soft clauses:    %7d                                "
           "                                   |\n",
           maxsat_formula->nSoft());
    printf("c |  Number of cardinality:     %7d                                "
           "                                   |\n",
           maxsat_formula->nCard());
    printf("c |  Number of PB :             %7d                                "
           "                                   |\n",
           maxsat_formula->nPB());
    double parsed_time = cpuTime();

    printf("c |  Parse time:           %12.2f s                                "
           "                                 |\n",
           parsed_time - initial_time);
    printf("c |                                                                "
           "                                       |\n");
     
    delete S;

  } catch (OutOfMemoryException &) {
    sleep(1);
    printf("c Error: Out of memory.\n");
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  } catch(MaxSATException &e) {
    sleep(1);
    printf("c Error: MaxSAT Exception: %s\n", e.getMsg());
    printf("s UNKNOWN\n");
    exit(_ERROR_);
  }
}
// ./open-wbo -epsilon=0.25  maxcut-hamming.wcnf
