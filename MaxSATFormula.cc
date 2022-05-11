/*!
 * \author Ruben Martins - ruben@sat.inesc-id.pt
 *
 * @section LICENSE
 *
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

#include <cstdint>
#include <iostream>

#include "MaxSATFormula.h"
#include <cmath>
#define pow2(n) ( 1 << (n) ) /* https://stackoverflow.com/questions/101439/the-most-efficient-way-to-implement-an-integer-based-power-function-powint-int */

using namespace std;
using namespace openwbo;

MaxSATFormula *MaxSATFormula::copyMaxSATFormula() {
  assert(format == _FORMAT_MAXSAT_);

  MaxSATFormula *copymx = new MaxSATFormula();
  copymx->setInitialVars(nVars());

  for (int i = 0; i < nVars(); i++)
    copymx->newVar();

  for (int i = 0; i < nSoft(); i++)
    copymx->addSoftClause(getSoftClause(i).weight, getSoftClause(i).clause);

  for (int i = 0; i < nHard(); i++)
    copymx->addHardClause(getHardClause(i).clause);

  copymx->setProblemType(getProblemType());
  copymx->updateSumWeights(getSumWeights());
  copymx->setMaximumWeight(getMaximumWeight());
  copymx->setHardWeight(getHardWeight());

  return copymx;
}

// Adds a new hard clause to the hard clause database.
void MaxSATFormula::addHardClause(vec<Lit> &lits) {
  hard_clauses.push();
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);
  new (&hard_clauses[hard_clauses.size() - 1]) Hard(copy_lits);
  n_hard++;
}

// Adds a new soft clause to the hard clause database.
void MaxSATFormula::addSoftClause(uint64_t weight, vec<Lit> &lits) {
  soft_clauses.push();
  vec<Lit> vars;
  Lit assump = lit_Undef;
  mpz_add_ui(clause_weight_sum, clause_weight_sum, weight); // update the weight sum
  uint64_t w = weight / pow(1 + alpha, lits.size() - 1);
  if (w < 1) {
    w = 1;
  }
  mpz_add_ui(bucket_clause_weight, bucket_clause_weight, w); // update the weight sum
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);
  if (clause_policy == 0) {
    if (use_pool) weight_sampler.push_back(1);
  }
  else if (clause_policy == 1) {
    if (use_pool) weight_sampler.push_back(weight);
  }
  else {
    if (use_pool) weight_sampler.push_back(w);
  }
  new (&soft_clauses[soft_clauses.size() - 1])
      Soft(copy_lits, weight, assump, vars);
  n_soft++;
}

void MaxSATFormula::clearBucket() {
  soft_clauses.clear();
}
void MaxSATFormula::addPoolClause(uint64_t weight, vec<Lit> &lits) {
  pool_clauses.push();
  vec<Lit> vars;
  Lit assump = lit_Undef;
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);
  uint64_t w = weight / pow(1 + alpha, lits.size() - 1);
  w = (w < 1) ? 1 : w;
  weight_pool.push_back(w);
  max_weight_pool = max(max_weight_pool, weight);
  new (&pool_clauses[pool_clauses.size() - 1])
      Soft(copy_lits, weight, assump, vars);
  n_pool++;
}

void MaxSATFormula::PoolCapacity() {
  cout << "pool_clauses.capacity(): " << pool_clauses.capacity() * sizeof(pool_clauses) << endl;
}

void MaxSATFormula::PrintPoolClause(int index) {
  cout << "Printing pool clause" << endl;
  cout << getPoolClause(index).weight << " ";
  for (int j = 0; j < getPoolClause(index).clause.size();
       j++) {
    if (sign(getPoolClause(index).clause[j])) {
      cout << "-";
    }
    cout << var(getPoolClause(index).clause[j]) + 1
             << " ";
  }
  cout << "0" << endl;
}

// void test_update_function(MaxSATFormula *maxsat_formula) {
//     cout << "Before update" << endl;
//     cout << maxsat_formula->getSoftClause(11).weight << " ";
//     for (int j = 0; j < maxsat_formula->getSoftClause(11).clause.size(); j++) {
//         cout << var(maxsat_formula->getSoftClause(11).clause[j]) + 1<< " ";
//     }
//     cout << endl;
//     cout << maxsat_formula->getSoftClause(13).weight << " ";
//     for (int j = 0; j < maxsat_formula->getSoftClause(13).clause.size(); j++) {
//         cout << var(maxsat_formula->getSoftClause(13).clause[j]) + 1<< " ";
//     }
//     cout << endl;
//     maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(11).weight, maxsat_formula->getSoftClause(11).clause, 13);
//     cout << "After update" << endl;
//     cout << maxsat_formula->getSoftClause(11).weight << " ";
//     for (int j = 0; j < maxsat_formula->getSoftClause(11).clause.size(); j++) {
//         cout << var(maxsat_formula->getSoftClause(11).clause[j]) + 1<< " ";
//     }
//     cout << endl;
//     cout << maxsat_formula->getSoftClause(13).weight << " ";
//     for (int j = 0; j < maxsat_formula->getSoftClause(13).clause.size(); j++) {
//         cout << var(maxsat_formula->getSoftClause(13).clause[j]) + 1 << " ";
//     }
//     cout << endl;
// }

void MaxSATFormula::updatePoolClause(uint64_t weight, vec<Lit> &lits, int pos) {
  // pool_clauses.push();
  vec<Lit> vars;
  Lit assump = lit_Undef;
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);
  uint64_t w = weight / pow(1 + alpha, lits.size() - 1);
  w = (w < 1) ? 1 : w;
  weight_pool[pos] = w;
  max_weight_pool = max(max_weight_pool, weight);
  new (&pool_clauses[pos])
      Soft(copy_lits, weight, assump, vars);
}

// Adds a new soft clause to the hard clause database with predefined relaxation
// variables.
void MaxSATFormula::addSoftClause(uint64_t weight, vec<Lit> &lits,
                                  vec<Lit> &vars) {
  soft_clauses.push();
  Lit assump = lit_Undef;
  mpz_add_ui(clause_weight_sum, clause_weight_sum, weight); // update the weight sum
  uint64_t w = weight / pow(1 + alpha, lits.size() - 1);
  w = (w < 1) ? 1 : w;
  mpz_add_ui(bucket_clause_weight, bucket_clause_weight, w); // update the weight sum
  vec<Lit> copy_lits;
  lits.copyTo(copy_lits);
  if (clause_policy == 0) {
    if (use_pool) weight_sampler.push_back(1);
  }
  else if (clause_policy == 1) {
    if (use_pool) weight_sampler.push_back(weight);
  }
  else {
    if (use_pool) weight_sampler.push_back(w);
  }
  new (&soft_clauses[soft_clauses.size() - 1])
      Soft(copy_lits, weight, assump, vars);
  n_soft++;
}

int MaxSATFormula::nInitialVars() {
  return n_initial_vars;
} // Returns the number of variables in the working MaxSAT formula.

void MaxSATFormula::setInitialVars(int vars) { n_initial_vars = vars; }

int MaxSATFormula::nVars() {
  return n_vars;
} // Returns the number of variables in the working MaxSAT formula.

int MaxSATFormula::nSoft() {
  return n_soft;
} // Returns the number of soft clauses in the working MaxSAT formula.

int MaxSATFormula::nPool() {
  return n_pool;
} // Returns the number of pool clauses in the working MaxSAT formula.

int MaxSATFormula::nHard() {
  return n_hard;
} // Returns the number of hard clauses in the working MaxSAT formula.

void MaxSATFormula::newVar(int v) {
  if(v == -1) n_vars++;
  else if(v > n_vars) n_vars = v;
} // Increases the number of variables in the working MaxSAT formula.

// Makes a new literal to be used in the working MaxSAT formula.
Lit MaxSATFormula::newLiteral(bool sign) {
  Lit p = mkLit(nVars(), sign);
  newVar();
  return p;
}

void MaxSATFormula::setProblemType(int type) {
  problem_type = type;
} // Sets the problem type.

int MaxSATFormula::getProblemType() {
  return problem_type; // Return the problem type.
}

// 'ubCost' is initialized to the sum of weights of the soft clauses.
void MaxSATFormula::updateSumWeights(uint64_t weight) {
  if (weight != hard_weight)
    sum_soft_weight += weight;
}

// The initial 'currentWeight' corresponds to the maximum weight of the soft
// clauses.
void MaxSATFormula::setMaximumWeight(uint64_t weight) {
  if (weight > max_soft_weight && weight != hard_weight)
    max_soft_weight = weight;
}

uint64_t MaxSATFormula::getMaximumWeight() { return max_soft_weight; }

void MaxSATFormula::setHardWeight(uint64_t weight) {
  hard_weight = weight;
} // Sets the weight of hard clauses.

Soft &MaxSATFormula::getSoftClause(int pos) {
  assert(pos < nSoft());
  return soft_clauses[pos];
}

Hard &MaxSATFormula::getHardClause(int pos) {
  assert(pos < nHard());
  return hard_clauses[pos];
}

Soft &MaxSATFormula::getPoolClause(int pos) {
  assert(pos < nPool());
  return pool_clauses[pos];
}

void MaxSATFormula::addPBConstraint(PB *p) {

  // Add constraint to formula data structure.
  if (p->isClause()) {
    addHardClause(p->_lits);
  } else if (p->isCardinality()) {
    if (!p->_sign) {
      p->changeSign();
    }
    cardinality_constraints.push(new Card(p->_lits, p->_rhs));

  } else {
    if (!p->_sign) {
      p->changeSign();
    }

    pb_constraints.push(new PB(p->_lits, p->_coeffs, p->_rhs, p->_sign));
  }
}

int MaxSATFormula::newVarName(char *varName) {
  int id = varID(varName);
  if (id == var_Undef) {
    id = nVars();
    newVar();
    std::string s(varName);
    std::pair<std::string, int> nv(s, id);
    std::pair<int, std::string> ni(id, s);
    _nameToIndex.insert(nv);
    _indexToName.insert(ni);
  }
  return id;
}

int MaxSATFormula::varID(char *varName) {
  std::string s(varName);

  nameMap::const_iterator iter = _nameToIndex.find(s);
  if (iter != _nameToIndex.end()) {
    return iter->second;
  }
  return var_Undef;
}

void MaxSATFormula::convertPBtoMaxSAT() {
  assert(objective_function != NULL);
  vec<Lit> unit_soft(1);

  // Convert objective function to soft clauses
  for (int i = 0; i < objective_function->_lits.size(); i++) {
    assert(objective_function->_coeffs[i] > 0);
    unit_soft[0] = ~objective_function->_lits[i];
    addSoftClause(objective_function->_coeffs[i], unit_soft);

    // Updates the maximum weight of soft clauses.
    setMaximumWeight(objective_function->_coeffs[i]);
    // Updates the sum of the weights of soft clauses.
    updateSumWeights(objective_function->_coeffs[i]);
  }

  if (getMaximumWeight() == 1)
    setProblemType(_UNWEIGHTED_);
  else
    setProblemType(_WEIGHTED_);
}
unordered_set<uint32_t> MaxSATFormula::pick_k_clauses(int k, bool reverse = false) {
  int rnd_max = weight_sampler.size();
    int ntake = k;

    if (reverse == false) {
      for (int index = 0; index < rnd_max; index++) {
        weight_sampler[index] = hard_clause_identifier - weight_sampler[index];
      }
    }

    /* determine smallest power of two that is larger than N */
    int tree_levels = ceil(log2((double) rnd_max));

    /* initialize vector with place-holders for perfectly-balanced tree */
    std::vector<double> tree_weights(pow2(tree_levels + 1));

    /* compute sums for the tree leaves at each node */
    int offset = pow2(tree_levels) - 1;
    for (int ix = 0; ix < rnd_max; ix++) {
        tree_weights[ix + offset] = weight_sampler[ix];
    }
    for (int ix = pow2(tree_levels+1) - 1; ix > 0; ix--) {
        tree_weights[(ix - 1) / 2] += tree_weights[ix];
    }

    /* sample according to uniform distribution */
    double rnd_subrange, w_left;
    double curr_subrange;
    int curr_ix;
    std::unordered_set<uint32_t> sampled(ntake);
    for (int el = 0; el < ntake; el++) {

        /* go down the tree by drawing a random number and
           checking if it falls in the left or right sub-ranges */
        curr_ix = 0;
        curr_subrange = tree_weights[0];
        for (int lev = 0; lev < tree_levels; lev++) {
            rnd_subrange = std::uniform_real_distribution<double>(0, curr_subrange)(rng);
            w_left = tree_weights[2 * curr_ix + 1];
            curr_ix = 2 * curr_ix + 1 + (rnd_subrange >= w_left);
            curr_subrange = tree_weights[curr_ix];
        }

        /* finally, add element from this iteration */
        sampled.insert(curr_ix - offset);

        /* now remove the weight of the chosen element */
        tree_weights[curr_ix] = 0;
        for (int lev = 0; lev < tree_levels; lev++) {
            curr_ix = (curr_ix - 1) / 2;
            tree_weights[curr_ix] =   tree_weights[2 * curr_ix + 1]
                                    + tree_weights[2 * curr_ix + 2];
        }
    }
  return sampled;
}
unordered_set<uint32_t> MaxSATFormula::pick_k_clauses_from_pool(int k) {
  int rnd_max = weight_pool.size();
    int ntake = k;

    /* determine smallest power of two that is larger than N */
    int tree_levels = ceil(log2((double) rnd_max));

    /* initialize vector with place-holders for perfectly-balanced tree */
    std::vector<double> tree_weights(pow2(tree_levels + 1));

    /* compute sums for the tree leaves at each node */
    int offset = pow2(tree_levels) - 1;
    for (int ix = 0; ix < rnd_max; ix++) {
        tree_weights[ix + offset] = max_weight_pool + 1 - weight_pool[ix];
        assert(tree_weights[ix + offset] >= 0);
    }
    for (int ix = pow2(tree_levels+1) - 1; ix > 0; ix--) {
        tree_weights[(ix - 1) / 2] += tree_weights[ix];
    }

    /* sample according to uniform distribution */
    double rnd_subrange, w_left;
    double curr_subrange;
    int curr_ix;
    std::unordered_set<uint32_t> sampled(ntake);
    for (int el = 0; el < ntake; el++) {

        /* go down the tree by drawing a random number and
           checking if it falls in the left or right sub-ranges */
        curr_ix = 0;
        curr_subrange = tree_weights[0];
        for (int lev = 0; lev < tree_levels; lev++) {
            rnd_subrange = std::uniform_real_distribution<double>(0, curr_subrange)(rng);
            w_left = tree_weights[2 * curr_ix + 1];
            curr_ix = 2 * curr_ix + 1 + (rnd_subrange >= w_left);
            curr_subrange = tree_weights[curr_ix];
        }

        /* finally, add element from this iteration */
        sampled.insert(curr_ix - offset);

        /* now remove the weight of the chosen element */
        tree_weights[curr_ix] = 0;
        for (int lev = 0; lev < tree_levels; lev++) {
            curr_ix = (curr_ix - 1) / 2;
            tree_weights[curr_ix] =   tree_weights[2 * curr_ix + 1]
                                    + tree_weights[2 * curr_ix + 2];
        }
    }
  return sampled;
}

void MaxSATFormula::status_pool() {
  // std::unordered_map<int, int> weight_map;
  std::map<std::pair<int,int>, int> weight_map;
  weight_map.clear();
  for (int cla_index = 0; cla_index < nPool(); cla_index++) {
    int weight = getPoolClause(cla_index).weight;
    int len = getPoolClause(cla_index).clause.size();
    if (weight_map.find(std::make_pair(weight,len)) == weight_map.end()) {
      weight_map[std::make_pair(weight,len)] = 1;
    }
    else {
      weight_map[std::make_pair(weight,len)] = weight_map[std::make_pair(weight,len)] + 1;
    }
  }
  std::cout << "Here is the clause pool: " << std::endl;
  for (auto& x: weight_map)
    std::cout << "(" <<x.first.first << ", " << x.first.second << "):" << x.second << ", ";
  std::cout << "Pool size: " << nPool() << std::endl;
}
