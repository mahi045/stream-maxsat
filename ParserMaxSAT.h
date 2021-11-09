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

template <class B, class MaxSATFormula>
static void parseMaxSAT(B &in, MaxSATFormula *maxsat_formula) {
  vec<Lit> lits;
  uint64_t hard_weight = UINT64_MAX;
  uint64_t num_var, num_cla;
  mpz_init_set_ui(maxsat_formula->clause_weight_sum, 0);
  mpz_init_set_ui(maxsat_formula->bucket_clause_weight, 0);
  maxsat_formula->weight_sampler.clear();
  double w;
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
        if (*in != '\r' && *in != '\n') {
          hard_weight = parseWeight(in);
          maxsat_formula->setHardWeight(hard_weight);
        }
        init_stream(maxsat_formula, num_var, num_cla);
      } else
        printf("c PARSE ERROR! Unexpected char: %c\n", *in),
            printf("s UNKNOWN\n"), exit(_ERROR_);
    } else if (*in == 'c' || *in == 'p')
      skipLine(in);
    else {
      uint64_t weight = readClause(in, maxsat_formula, lits);
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
      
      for (int j = 0; j < lits.size(); j++)
      {
        w = (double) weight / pow(2, lits.size() - 1);
        if (sign(lits[j])) w *= -1;
        maxsat_formula->bias += (abs(maxsat_formula->var_bias[var(lits[j])] + w) - abs(maxsat_formula->var_bias[var(lits[j])]));
        maxsat_formula->var_bias[var(lits[j])] += w;
      }
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
  printf("Lower bound of MaxSAT: %ju\n", sum);
  printf("v");
  if (maxsat_formula->bias > bias_thre) {
    for (uint32_t k = 1; k <= min(maxsat_formula->nVars(), 10); k++) {
      if (maxsat_formula->var_bias[k] >= 0) {
        printf(" %ju", k);
      }
      else {
        printf(" -%ju", k);
      }
    }
  }
  else {
    for (uint32_t k = 1; k <= min(maxsat_formula->nVars(), 10); k++) {
      if (maxsat_formula->var_bias[k] >= 0) {
        printf(" %ju", k);
      }
      else {
        printf(" -%ju", k);
      }
    }
  }
  printf(" 0\n");
  // if (maxsat_formula->nSoft() % BUCKET_SIZE > 0) {
  //   printf("%d-th bucket !! \n", (maxsat_formula->nSoft() / BUCKET_SIZE) + 1);
  //   streaming_maxsat(maxsat_formula);
  // }
  printf("Sum of weight: %s\n", mpz_get_str (NULL, 10, maxsat_formula->clause_weight_sum));
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
