#ifndef SAMPLING_h
#define SAMPLING_h

#include <cstdio>
#include <iostream>
#include <stdio.h>
#include <string>
#include <unordered_set>
#include <vector>
#include <random>
#include <fstream>
#include "constants.h"
#include <sstream>


using NSPACE::vec;
using namespace openwbo;
using namespace std;
std::random_device rd;
std::mt19937 gen(rd());

/* the function chooses k elements from n elements */
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

void sample_clauses(MaxSATFormula *maxsat_formula) {
    POOL_SIZE = K * maxsat_formula->nVars() / (eps * eps * R); //CAUTIOUS: this special case for big stream of clauses
    if (POOL_SIZE >= maxsat_formula->nSoft()) {
        POOL_SIZE = maxsat_formula->nSoft();
        cout << "the number of sampled clauses: " << POOL_SIZE << " (same as original)" << endl;
    }
    else {
        cout << "the number of sampled clauses: " << POOL_SIZE << " (less than original clauses " << maxsat_formula->nSoft() << ")" << endl;
    }
    b.resize(maxsat_formula->nSoft());
    beta = ceil(log(maxsat_formula->nSoft()) / Gamma);
    cout << "the value of beta: " << beta << endl;
    iter = ceil(log(maxsat_formula->nSoft()) / (log (10) * eps));
    cout << "the value of iterations: " << iter << endl;
    uint64_t nlong_clause = 0;
    bool postprocess = false;
    string line, variable;
    int lit;
    string delim = " ";
    unordered_set<uint32_t> in_pool;
    // cout << POOL_SIZE << " " <<  maxsat_formula->nSoft() << endl;
    if (POOL_SIZE == maxsat_formula->nSoft()) {
        std::iota(b.begin(), b.end(), 0);
    }
    // else if (maxsat_formula->nSoft() - POOL_SIZE < POOL_SIZE) {
    //     std::iota(b.begin(), b.end(), 0);
    //     vector<int> temp_b = vector<int>();
    //     in_pool = maxsat_formula->pick_k_clauses(maxsat_formula->nSoft() - POOL_SIZE, false);
    //     // cout << in_pool.size() << " " << b.size() << " ";
    //     for (auto itr = b.begin(); itr != b.end(); ++itr) {
    //         if (in_pool.find(*itr) == in_pool.end()) {
    //             temp_b.push_back(*itr);
    //         }
    //     }
    //     b.clear();
    //     b = temp_b;
    //     assert(b.size() == POOL_SIZE);
    // }
    else {
      b.clear();
      postprocess = true;
      for (int j = 0; j < maxsat_formula->nSoft(); j++) {
        if (maxsat_formula->getSoftClause(j).clause.size() >= beta) {
          maxsat_formula->weight_sampler[j] = 0;
          nlong_clause++;
        }
      }
      in_pool = maxsat_formula->pick_k_clauses(POOL_SIZE, true);
      for (auto itr = in_pool.begin(); itr != in_pool.end(); ++itr) {
        b.push_back(*itr);
      }
      assert(b.size() == POOL_SIZE);
    }
    sort(b.begin(), b.end());
    
    // if (nlong_clause == 0) 
    //     postprocess = false;
    // if (postprocess)
    //     cout << "We need postprocessing algorithm" << endl;
}

void run_sampled_maxsat(MaxSATFormula *maxsat_formula) {
    string line, variable;
    int lit;
    string delim = " ";
    ofstream myfile, assignfile;
    std::string sampled_maxsat_file = "hoa_sampled_" + file_name;
    if (variant == 2) {
        sampled_maxsat_file = "sampled_" + file_name;
    }
    myfile.open(sampled_maxsat_file);
    myfile << "p wcnf " + to_string(maxsat_formula->nVars()) + " " + to_string(POOL_SIZE) + " " + to_string(maxsat_formula->getHardWeight()) << endl;
    for (auto index = 0; index < b.size(); index++) {
        maxsat_formula->addPoolClause(maxsat_formula->getSoftClause(index).weight, 
                        maxsat_formula->getSoftClause(index).clause);
        myfile << maxsat_formula->getSoftClause(index).weight << " ";
        for (int j = 0; j < maxsat_formula->getSoftClause(index).clause.size(); j++) {
            if (sign(maxsat_formula->getSoftClause(index).clause[j])) {
                myfile << "-";
            }
            myfile << var(maxsat_formula->getSoftClause(index).clause[j]) + 1 << " ";
        }
        myfile << "0" << endl;
    }
    myfile.close();

    std::ostringstream stringStream;
    stringStream << "./open-wbo_static -print-model -cpu-lim="<< TIMEOUT << " " << sampled_maxsat_file + " > " + "result_" + sampled_maxsat_file;
    // calling the sampled maxsat query
    // cout << stringStream.str() << endl;
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