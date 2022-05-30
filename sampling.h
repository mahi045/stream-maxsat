#ifndef STREAMING_h
#define STREAMING_h
#include "constants.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <istream>
#include <random>
#include <stdio.h>
#include <sstream>
#include <string>
#include <unordered_set>
#include <algorithm>

using NSPACE::vec;
using namespace openwbo;
using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int currentUsedSizeinVM(){ //Note: this value is in KB!
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != NULL){
        if (strncmp(line, "VmSize:", 7) == 0){
            result = parseLine(line);
            break;
        }
    }
    fclose(file);
    return result;
}

vector<int> sample_k_items(int n, int k) {
    vector<int> b(k);
    for(std::size_t i = 0; i != n; ++i) {
        std::uniform_int_distribution<> dis(0, i);
        std::size_t j = dis(gen);
        if(j < b.size()) {
            if (i < b.size()) {
                b[i] = b[j];
            }
            b[j] = i;
        }
    }
    return b;
}

void init_sampling(MaxSATFormula *maxsat_formula, uint64_t var, uint64_t cla) {
    // POOL_SIZE = min((uint64_t) (K * var / (eps * eps)), cla);
    int random_k = 0;
    if (random_sat_of_beta) {
        int k = 1;
        for (; ; k++) {
            if (var * pow(2, k) > k * cla) {
                break;
            }
        }
        random_k = k;
        // cout << "The k of random k sat is: " << k << endl;
    }
    int log_of_clause = 0;
    if (log_of_beta) {
        log_of_clause = ceil(log2(cla));
        // cout << "The log clauses is: log2(" << cla << ") = " << maxsat_formula->beta << endl;
    }
    int expectation = 0;
    if (expectation_of_clause) {
        double exp = 0;
        for (auto &x : maxsat_formula->clause_map)
        {
            exp += ((double) x.first.second * x.second / cla);
        }
        expectation =  4 * ceil(exp);
        // cout << "The expected clause lenght is: E[clause_lenght] = " << maxsat_formula->beta << endl;
    }
    int minimum = min(expectation, min(log_of_clause, random_k));
    int maximum = max(expectation, max(log_of_clause, random_k));
    // expectation + log_of_clause + random_k - minimum -
    maxsat_formula->beta =  expectation + log_of_clause + random_k - minimum - maximum;
    cout << "median(" << expectation << "," << random_k << "," << log_of_clause << ")= " << maxsat_formula->beta << endl;
    POOL_SIZE = (3000 * fraction_of_memory * 1000 * 1000) / (4 * (maxsat_formula->beta) + sizeof(Soft));
    POOL_SIZE = min(POOL_SIZE, cla);
    BUCKET_SIZE = ceil((double) POOL_SIZE / R);
    if (POOL_SIZE == cla) {
        store_all = true;
    }
    cout << "The pool size is: " << POOL_SIZE << ", which is " << (double) POOL_SIZE / var << " factor of n" << endl;
    cout << "The number of clauses is " <<  (double) cla / var << " factor of n" << endl;
    // setting the capacity of pool
    maxsat_formula->createPool(POOL_SIZE);
    if (hoa && POOL_SIZE < cla) {
        postprocessing = true;
        iter = ceil(log(cla) / (log (10) * 0.25));
        cout << "the value of iterations: " << iter << endl;
    }
    // maxsat_formula->assignment.growTo(var + 1, l_Undef);
    // maxsat_formula->var_bias.growTo(var + 1, 0);
    // maxsat_formula->weight_pool.clear();
}

void modify_pool(MaxSATFormula *maxsat_formula) { 
    // POOL_SIZE = min((int) (K * maxsat_formula->nVars() / (eps * eps)), maxsat_formula->nSoft());
    // BUCKET_SIZE = POOL_SIZE / R;
    // test_update_function(maxsat_formula);
    int bucket_start;
    int index_bucket, index_pool;
    double positive_phase, negative_phase;
    // maxsat_formula->weight_sampler.clear();
    int bucket_index = maxsat_formula->nSoft() / BUCKET_SIZE - 1;
    int bound = (maxsat_formula->nSoft() % BUCKET_SIZE) ? maxsat_formula->nSoft() % BUCKET_SIZE : BUCKET_SIZE;
    // cout << "sizeof(maxsat_formula->getSoftClause(0).weight) => " << sizeof(Soft) << endl;
    // cout << "sizeof(maxsat_formula->getSoftClause(0).weight) => " << maxsat_formula->getSoftClause(0).clause.size() << endl;
    for (int i = 0; i < bound; i++) {
        if (!(i % BUCKET_SIZE)) {
            bucket_start = i;
        }
        if (!((i + 1) % BUCKET_SIZE) || i + 1 == bound) {

            if (!store_all && maxsat_formula->clause_seen_so_far + BUCKET_SIZE > POOL_SIZE) {
                // this clauses will be replaced
                // int clause_need_replace = ((double) (BUCKET_SIZE) / (BUCKET_SIZE + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
                int clause_need_replace = ceil(((double) (mpz_get_d(maxsat_formula->bucket_clause_weight)) / (mpz_get_d(maxsat_formula->clause_weight_sum))) 
                    * maxsat_formula->nPool());
                int remaining_clause = min(bound, BUCKET_SIZE);
                // cout << "clause_need_replace: " << clause_need_replace << endl;
                clause_need_replace = min(clause_need_replace, remaining_clause);
                // cout << "clause_need_replace: " << clause_need_replace << endl;
                // if (false && i + 1 == bound) {
                //     cout << "Now it is executed" << endl;
                //     remaining_clause = i + 1;
                //     clause_need_replace = ((double) (remaining_clause) / (remaining_clause + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
                // }
                if (!hoa && !L_1) {
                    // it is special case for random sampling 
                    clause_need_replace = ((double) (remaining_clause) / (remaining_clause + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
                }
                clause_need_replace = min(clause_need_replace, remaining_clause);
                unordered_set<uint32_t> replaced_clause_pool = maxsat_formula->pick_k_clauses_from_pool(clause_need_replace);
                unordered_set<uint32_t> replaced_clause_bucket = maxsat_formula->pick_k_clauses(clause_need_replace, true);
                cout << " From bucket index " << bucket_start + bucket_index * BUCKET_SIZE << " to " << i + bucket_index * BUCKET_SIZE << " => ";
                auto start_itr1 = replaced_clause_pool.begin();
                auto start_itr2 = replaced_clause_bucket.begin();
                // cout << "replaced_clause_pool.size(): " << replaced_clause_pool.size() << endl;
                // cout << "replaced_clause_bucket.size(): " << replaced_clause_bucket.size() << endl;
                // cout << "clause_need_replace: " << clause_need_replace << endl;
                assert(replaced_clause_pool.size() == replaced_clause_bucket.size());
                for (;start_itr1 != replaced_clause_pool.end(); start_itr1++, start_itr2++) {
                    // index_bucket = *start_itr2 + maxsat_formula->clause_seen_so_far - 1;
                    index_bucket = *start_itr2;
                    index_pool = *start_itr1;
                    // debugfile << maxsat_formula->getSoftClause(index_bucket).weight << " is replacing " << maxsat_formula->getPoolClause(index_pool).weight << endl;
                    // debugfile << "Before update: " << maxsat_formula->getPoolClause(index_pool).weight << " ";
                    maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(index_bucket).weight, 
                        maxsat_formula->getSoftClause(index_bucket).clause, index_pool);
                    // debugfile << "After update: " << maxsat_formula->getPoolClause(index_pool+1).weight << " ";
                }
                cout << "Total " << clause_need_replace << " clauses replaced from pool !!!" << endl;
            }
            else {
                for (auto cla_index = 0; cla_index != i + 1; cla_index++) {
                    maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(cla_index).weight, 
                        maxsat_formula->getSoftClause(cla_index).clause, maxsat_formula->last_index_in_pool);
                    maxsat_formula->last_index_in_pool++; // increment the next position in pool (ESSENTIAL)
                }
            }
            maxsat_formula->clause_seen_so_far += (i + 1);
            mpz_init(maxsat_formula->bucket_clause_weight);
            maxsat_formula->weight_sampler.clear();
            maxsat_formula->weight_sampler.shrink_to_fit();
            if (verbose)
                maxsat_formula->status_pool();
        }
    }
}

void sample_clauses(MaxSATFormula *maxsat_formula) {
    ofstream myfile, assignfile;
    string line, variable;
    int lit;
    string delim = " ";
    std::string sampled_maxsat_file = "sampled_" + file_name;
    if (hoa) {
        sampled_maxsat_file = "hoa_sampled_" + file_name;
    }
    myfile.open(sampled_maxsat_file);
    myfile << "p wcnf " + to_string(maxsat_formula->nVars()) + " " +
                  to_string(POOL_SIZE) + " " +
                  to_string(maxsat_formula->getHardWeight())
           << endl;

    for (int cla_index = 0; cla_index < maxsat_formula->nPool(); cla_index++) {
      myfile << maxsat_formula->getPoolClause(cla_index).weight << " ";
      for (int j = 0;
           j < maxsat_formula->getPoolClause(cla_index).clause.size(); j++) {
        if (sign(maxsat_formula->getPoolClause(cla_index).clause[j])) {
          myfile << "-";
        }
        myfile << var(maxsat_formula->getPoolClause(cla_index).clause[j]) + 1
                 << " ";
      }
      myfile << "0" << endl;
    }

    myfile.close();
    int available_memory = total_memory;
    if (use_fixed_memory) {
      int used_memory = currentUsedSizeinVM() / 1024;
      available_memory = (available_memory > used_memory)
                             ? (available_memory - used_memory)
                             : available_memory;
    }
    cout << "The available memory: " << available_memory
         << endl;
    std::ostringstream stringStream;
    int timeout = TIMEOUT - 100;
    timeout = (timeout < 10) ? 10 : timeout;
    stringStream << "./open-wbo_static -print-model -cpu-lim="<< to_string(timeout) << " -mem-lim=" << available_memory << " " << sampled_maxsat_file + " > " + "result_" + sampled_maxsat_file;
    
    system(stringStream.str().c_str());

    
    maxsat_formula->assignment.growTo(maxsat_formula->nVars() + 1, l_Undef);
    string result_file_name = "result_" + sampled_maxsat_file;
    ifstream resultfile(result_file_name);
    while (getline(resultfile, line))
    {
        if (line.rfind("v ") == 0)
        {
            auto start = 0U;
            auto end = line.find(delim);
            while (end != std::string::npos)
            {
                variable = line.substr(start, end - start);
                start = end + delim.length();
                end = line.find(delim, start);
                if (variable != "v")
                {
                    lit = stoi(variable);
                    if (lit < 0)
                    {
                        maxsat_formula->assignment[abs(lit)] = l_False;
                    }
                    else
                    {
                        maxsat_formula->assignment[lit] = l_True;
                    }
                }
            }
        }
    }
    resultfile.close();
    mpz_t cost, least_cost;
    vector<int> best_per;
    mpz_init_set_ui(least_cost, 0);
    vector<int> per_var;
    if (postprocessing) {
        cout << "The post-processing algorithm is running ...";
        bool unsat = false;
        for (int itr = 0; itr < iter; itr++)
        {
            mpz_init_set_ui(cost, 0);
            for (int cla_index = 0; cla_index < maxsat_formula->nPool(); cla_index++)
            {
                unsat = true;
                for (int j = 0; j < maxsat_formula->getPoolClause(cla_index).clause.size(); j++)
                {
                    Lit p = maxsat_formula->getPoolClause(cla_index).clause[j];
                    if (sign(p)) {
                        if (maxsat_formula->assignment[var(p) + 1] == l_False) {
                            unsat = false;
                        }
                    }
                    else {
                        if (maxsat_formula->assignment[var(p) + 1] == l_True) {
                            unsat = false;
                        }
                    }
                }
                if (unsat) {
                    mpz_add_ui(cost, cost,
                        maxsat_formula->getPoolClause(cla_index).weight);
                }
            }
            if (itr == 0) {
                mpz_set(least_cost, cost);
            }
            else {
                if (mpz_cmp(least_cost, cost) > 0) {
                    mpz_set(least_cost, cost);
                    best_per.clear();
                    best_per.assign(per_var.begin(), per_var.end());
                }
            }
            // undo the last perturb
            for (int per_atom = 0; per_atom < per_var.size(); per_atom++) {
                if (maxsat_formula->assignment[per_var[per_atom]] == l_True) {
                    maxsat_formula->assignment[per_var[per_atom]] == l_False;
                }
                else if (maxsat_formula->assignment[per_var[per_atom]] == l_False)
                {
                    maxsat_formula->assignment[per_var[per_atom]] == l_True;
                }
            }
            per_var.clear();
            // the new perturb
            per_var = sample_k_items(maxsat_formula->nVars(), (int) (maxsat_formula->nVars() * Gamma));
            // cout << "The size of per_var: " << per_var.size() << endl;
            for (int per_atom = 0; per_atom < per_var.size(); per_atom++) {
                if (maxsat_formula->assignment[per_var[per_atom]] == l_True) {
                    maxsat_formula->assignment[per_var[per_atom]] == l_False;
                }
                else if (maxsat_formula->assignment[per_var[per_atom]] == l_False)
                {
                    maxsat_formula->assignment[per_var[per_atom]] == l_True;
                }
            }
        }
        // undo the last perturb
        for (int per_atom = 0; per_atom < per_var.size(); per_atom++)
        {
            if (maxsat_formula->assignment[per_var[per_atom]] == l_True) {
                maxsat_formula->assignment[per_var[per_atom]] == l_False;
            }
            else if (maxsat_formula->assignment[per_var[per_atom]] == l_False) {
                maxsat_formula->assignment[per_var[per_atom]] == l_True;
            }
        }
        per_var.clear();
        // assigning the value to best perturbtion
        for (int per_atom = 0; per_atom < best_per.size(); per_atom++)
        {
            if (maxsat_formula->assignment[best_per[per_atom]] == l_True) {
                maxsat_formula->assignment[best_per[per_atom]] == l_False;
            }
            else if (maxsat_formula->assignment[best_per[per_atom]] == l_False) {
                maxsat_formula->assignment[best_per[per_atom]] == l_True;
            }
        }
    }
    // finally write the final assignment
    assignfile.open("result_" + sampled_maxsat_file);
    assignfile << "v ";
    for (int variable = 1; variable <= maxsat_formula->nVars(); variable++)
    {
        if (maxsat_formula->assignment[variable] == l_True) {
            assignfile << variable << " ";
        }
        else if (maxsat_formula->assignment[variable] == l_False) {
            assignfile << -variable << " ";
        }
    }
    assignfile.close();
    
    return;
}

#endif