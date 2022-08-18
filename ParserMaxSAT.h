/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
 * MiniSat,  Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
 *           Copyright (c) 2007-2010, Niklas Sorensson
 * Open-WBO, Copyright (c) 2013-2015, Ruben Martins, Vasco Manquinho, Ines Lynce
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

#ifndef ParserMaxSAT_h
#define ParserMaxSAT_h

#include <cstdint>
#include <cstdio>
#include <stdio.h>
#include <cmath>
#include "MaxSATFormula.h"
#include "core/SolverTypes.h"
#include "utils/ParseUtils.h"
#include "constants.h"
#include "streaming.h"
#include "sampling.h"

#ifdef HAS_EXTRA_STREAMBUFFER
#include "utils/StreamBuffer.h"
#endif

using NSPACE::mkLit;
using NSPACE::StreamBuffer;

namespace openwbo {

//=================================================================================================
// DIMACS Parser:

template <class B> static uint64_t parseWeight(B &in) {
  uint64_t val = 0;
  while ((*in >= 9 && *in <= 13) || *in == 32)
    ++in;
  if (*in < '0' || *in > '9')
    fprintf(stderr, "PARSE ERROR! Unexpected char: %c\n", *in), exit(3);
  while (*in >= '0' && *in <= '9')
    val = val * 10 + (*in - '0'), ++in;
  return val;
}

template <class B, class MaxSATFormula>
static uint64_t readClause(B &in, MaxSATFormula *maxsat_formula,
                           vec<Lit> &lits) {
  int parsed_lit, var;
  int64_t weight = 1;
  lits.clear();
  if (maxsat_formula->getProblemType() == _WEIGHTED_)
    weight = parseWeight(in);
  assert(weight > 0);

  for (;;) {
    parsed_lit = parseInt(in);
    if (parsed_lit == 0)
      break;
    var = abs(parsed_lit) - 1;
    while (var >= maxsat_formula->nVars())
      maxsat_formula->newVar();
    lits.push((parsed_lit > 0) ? mkLit(var) : ~mkLit(var));
  }
  return weight;
}

vector<uint32_t> sample_k_vars(int n, double k)
{
  int t = ceil(n * k);
  // cout << "k " << k << endl;
  vector<uint32_t> b(t);
  for (std::size_t i = 0; i != n; ++i)
  {
    std::uniform_int_distribution<> dis(0, i);
    std::size_t j = dis(gen);
    if (j < b.size())
    {
      if (i < b.size())
      {
        b[i] = b[j];
      }
      b[j] = i;
    }
  }
  return b;
}
int get_assignment(MaxSATFormula *maxsat_formula) {
  string line, variable;
  string delim = " ";
  int lit;
  ifstream resultfile2(result_file_name);
  int assigned_variable = 0;
  if (resultfile2.good()) {
    bool no_assign = false;
    while (getline(resultfile2, line)) {
      if (line.rfind("v ") == 0) {
        no_assign = true;
        auto start = 0U;
        auto end = line.find(delim);
        while (end != std::string::npos) {
          variable = line.substr(start, end - start);
          start = end + delim.length();
          end = line.find(delim, start);
          if (variable != "v") {
            lit = stoi(variable);
            if (lit < 0) {
              maxsat_formula->assignment[abs(lit)] = l_False;
              assigned_variable++;
            } else {
              maxsat_formula->assignment[lit] = l_True;
              assigned_variable++;
            }
          }
        }
      }
    }
    std::cout << "Total number of assigned variables: " << assigned_variable << std::endl;
  }
  else {
    std::cout << "No assignment is given !!!" << std::endl;
  }
  if (assigned_variable > 0) {
    random_device rd;
    mt19937 rng{12345};
    uniform_real_distribution<> dis(0,1.0);
    bool default_value = dis(maxsat_formula->getRNG()) > 0.5;
    for (auto index = 1; index < maxsat_formula->assignment_default.size(); index++) {
      if (maxsat_formula->assignment[index] == l_Undef) {
        if (dis(maxsat_formula->getRNG()) > 0.5)
          maxsat_formula->assignment_default[index] = l_True;
        else 
          maxsat_formula->assignment_default[index] = l_False;
      }
    }
  }
  return assigned_variable;
}
template <class B, class MaxSATFormula>
static void parseMaxSAT(B &in, MaxSATFormula *maxsat_formula) {
  vec<Lit> lits;
  ofstream assignfile;
  uint64_t hard_weight = UINT64_MAX;
  uint64_t num_var, num_cla;
  vector<bool> occ;
  uint32_t variables_not_found = 0;
  mpz_init_set_ui(maxsat_formula->clause_weight_sum, 0);
  mpz_init_set_ui(maxsat_formula->bucket_clause_weight, 0);
  mpz_init_set_ui(maxsat_formula->unsat_weight, 0);
  mpz_init_set_ui(maxsat_formula->unsat_weight_false, 0);
  maxsat_formula->weight_sampler.clear();
  double w;
  maxsat_formula->weight_map.clear();
  maxsat_formula->unsat_clauses.clear();
  // random_device rd;
  // mt19937 rng{12345};
  // uniform_real_distribution<> dis(0,1.0);
  int assigned_variable = 0;
  // bool default_value = dis(maxsat_formula->getRNG()) > 0.5;
  for (;;) {
    skipWhitespace(in);
    if (*in == EOF)
      break;
    else if (*in == 'p') {
      if (eagerMatch(in, "p cnf")) {
        num_var = parseInt(in); // Variables
        num_cla = parseInt(in); // Clauses
        init_stream(maxsat_formula, num_var, num_cla);
      } else if (eagerMatch(in, "wcnf")) {
        maxsat_formula->setProblemType(_WEIGHTED_);
        num_var = parseInt(in); // Variables
        num_cla = parseInt(in); // Clauses
        occ.assign(num_var+1, false);
        variables_not_found = num_var;
        if (*in != '\r' && *in != '\n') {
          hard_weight = parseWeight(in);
          maxsat_formula->setHardWeight(hard_weight);
        }
        init_stream(maxsat_formula, num_var, num_cla);
      } else
        printf("c PARSE ERROR! Unexpected char: %c\n", *in),
            printf("s UNKNOWN\n"), exit(_ERROR_);
      assigned_variable = get_assignment(maxsat_formula);
    } else if (*in == 'c' || *in == 'p')
      skipLine(in);
    else {
      uint64_t weight = readClause(in, maxsat_formula, lits);
      if (!ignore_hard || weight < hard_weight - 1)
        mpz_add_ui(maxsat_formula->clause_weight_sum, maxsat_formula->clause_weight_sum, weight);
      if (maxsat_formula->m.size() < lits.size() + 1)
      {
        maxsat_formula->m.resize(lits.size() + 1, 0);
        maxsat_formula->m[lits.size()] = weight;
      }
      else
      {
        maxsat_formula->m[lits.size()] += weight;
      }
      int len = lits.size();
      if (maxsat_formula->weight_map.find(std::make_pair(weight,len)) == maxsat_formula->weight_map.end()) {
        maxsat_formula->weight_map[std::make_pair(weight,len)] = 1;
      } 
      else {
        maxsat_formula->weight_map[std::make_pair(weight,len)] = maxsat_formula->weight_map[std::make_pair(weight,len)] + 1;
      }
      bool unsat = true;
      bool unsat_d_false = true;
      for (int j = 0; j < lits.size(); j++)
      {
        w = (double) weight / pow(2, lits.size() - 1);
        if (sign(lits[j])) w *= -1;
        maxsat_formula->bias += (abs(maxsat_formula->var_bias[var(lits[j])] + w) - abs(maxsat_formula->var_bias[var(lits[j])]));
        maxsat_formula->var_bias[var(lits[j])] += w;
        if (occ[var(lits[j])] == false) {
          variables_not_found--;
          occ[var(lits[j])] = true;
        }
        
        if (sign(lits[j])) {
          if (maxsat_formula->assignment[var(lits[j]) + 1] == l_False) {
            unsat = false;
            unsat_d_false = true;
          }
          else if (unsat && maxsat_formula->assignment[var(lits[j]) + 1] == l_Undef && maxsat_formula->assignment_default[var(lits[j]) + 1] == l_False) {
            unsat_d_false = false;
          }
        }
        else {
          if (maxsat_formula->assignment[var(lits[j]) + 1] == l_True) {
            unsat = false;
            unsat_d_false = true;
          }
          else if (unsat && maxsat_formula->assignment[var(lits[j]) + 1] == l_Undef && maxsat_formula->assignment_default[var(lits[j]) + 1] == l_True) {
            unsat_d_false = false;
          }
        }
        if (!unsat) 
          break;
      }
      if (unsat) {
        len = lits.size();
        if (maxsat_formula->unsat_clauses.find(std::make_pair(weight, len)) ==
            maxsat_formula->unsat_clauses.end()) {
          maxsat_formula->unsat_clauses[std::make_pair(weight, len)] = 1;
        } else {
          maxsat_formula->unsat_clauses[std::make_pair(weight, len)] =
              maxsat_formula->unsat_clauses[std::make_pair(weight, len)] + 1;
        }
        if (!ignore_hard || weight < hard_weight - 1)
          mpz_add_ui(maxsat_formula->unsat_weight, maxsat_formula->unsat_weight,
                    weight);
      }
      if (unsat && unsat_d_false && assigned_variable > 0) {
        mpz_add_ui(maxsat_formula->unsat_weight_false, maxsat_formula->unsat_weight_false,
                    weight);
      }
      // if (unsat && unsat_d_true && assigned_variable > 0) {
      //   mpz_add_ui(maxsat_formula->unsat_weight_true, maxsat_formula->unsat_weight_true,
      //               weight);
      // }
      
      // if (weight < hard_weight ||
      //     maxsat_formula->getProblemType() == _UNWEIGHTED_) {
      //   assert(weight > 0);
      //   // Updates the maximum weight of soft clauses.
      //   maxsat_formula->setMaximumWeight(weight);
      //   // Updates the sum of the weights of soft clauses.
      //   maxsat_formula->updateSumWeights(weight);
      //   maxsat_formula->addSoftClause(weight, lits);
      // } else
      //   maxsat_formula->addHardClause(lits);
        
      // if ((maxsat_formula->nSoft() > 0) && (maxsat_formula->nSoft() % BUCKET_SIZE == 0)) {
      //   printf("%d-th bucket !! \n", maxsat_formula->nSoft() / BUCKET_SIZE);
      //   streaming_maxsat(maxsat_formula);
      //   maxsat_formula->clearBucket();
      // }
    }
  }
  double bias_thre = bias_threshold(maxsat_formula);
  double coff;
  uint64_t sum = 0;
  if (maxsat_formula->bias > bias_thre) {
    for (int k = 1; k <= maxsat_formula->m.size(); k++) {
      if (maxsat_formula->var_bias[k] >= 0) {
        continue;
      }
      sum += floor(((double) (k) / pow(2, k)) * maxsat_formula->m[k]);
    }
    sum += ceil(maxsat_formula->bias / 2);
  }
  else {
    for (int k = 1; k <= maxsat_formula->m.size(); k++) {
      if (maxsat_formula->m[k] == 0) {
        continue;
      }
      sum += floor((1 - (double) (1) / pow(2, k)) * maxsat_formula->m[k]);
    }
    sum += ceil(((double) (maxsat_formula->bias * maxsat_formula->bias) / (4 * bias_thre)));
  }
  // printf("Lower bound of MaxSAT: %ju\n", sum);
  assignfile.open("result_k_maxsat_" + file_name);
  assignfile << "v ";
  printf("Sum of weight: %s\n", mpz_get_str (NULL, 10, maxsat_formula->clause_weight_sum));
  printf("Sum of unsat weight: %s\n", mpz_get_str (NULL, 10, maxsat_formula->unsat_weight));
  if (assigned_variable > 0) {
    printf("Sum of unsat (d: default) weight: %s\n", mpz_get_str (NULL, 10, maxsat_formula->unsat_weight_false));
    // printf("Sum of unsat (d: true) weight: %s\n", mpz_get_str (NULL, 10, maxsat_formula->unsat_weight_true));
  }
  // printf("Number of variables not found in the formula: %u\n", variables_not_found);
  // printf("v");
  if (maxsat_formula->bias > bias_thre) {
    for (uint32_t k = 1; k <= maxsat_formula->nVars(); k++) {
      if (maxsat_formula->var_bias[k] >= 0) {
        assignfile << k << " ";
      }
      else {
        assignfile <<  "-" << k  << " ";
      }
    }
  }
  else {
    vector<uint32_t> var = sample_k_vars(maxsat_formula->nVars(), 0.5 - (double) (maxsat_formula->bias) 
        / (2 * bias_thre));
    std::sort(var.begin(), var.end());
    int in = 0;
    for (uint32_t k = 1; k <= maxsat_formula->nVars(); k++) {
      if (maxsat_formula->var_bias[k] >= 0) {
        if(in < var.size() && var[in] == k) {
          assignfile << "-" << k << " ";
          in++;
        }
        else 
          assignfile << k << " ";
      }
      else {
        if(in < var.size() && var[in] == k) {
          assignfile << k << " ";
          in++;
        }
        else 
          assignfile << "-" << k << " ";
      }
    }
  }
  // assignfile << "0" << endl;
  // if (maxsat_formula->nSoft() % BUCKET_SIZE > 0) {
  //   printf("%d-th bucket !! \n", (maxsat_formula->nSoft() / BUCKET_SIZE) + 1);
  //   streaming_maxsat(maxsat_formula);
  // }
  
  std::cout << "Here is the clause statistic: " << std::endl;
  for (auto &x : maxsat_formula->weight_map) {
    std::cout << "(weight:" << x.first.first << ", len:" << x.first.second
              << "):" << x.second << ", ";
  }
  std::cout << "\nHere is the unsat clause statistic: " << std::endl;
  for (auto &x : maxsat_formula->unsat_clauses) {
    std::cout << "(weight:" << x.first.first << ", len:" << x.first.second
              << "):" << x.second << ", ";
  }
  // assert(maxsat_formula->nSoft() == maxsat_formula->weight_sampler.size());
}

// Inserts problem into solver.
//
template <class MaxSATFormula>
static void parseMaxSATFormula(gzFile input_stream,
                               MaxSATFormula *maxsat_formula) {
  StreamBuffer in(input_stream);
  parseMaxSAT(in, maxsat_formula);
  if (maxsat_formula->getMaximumWeight() == 1)
    maxsat_formula->setProblemType(_UNWEIGHTED_);
  else
    maxsat_formula->setProblemType(_WEIGHTED_);

  // maxsat_formula->setInitialVars(maxsat_formula->nVars());
}

//=================================================================================================
} // namespace openwbo

#endif
