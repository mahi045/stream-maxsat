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
    vector<int> b(POOL_SIZE);
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
    POOL_SIZE = K * maxsat_formula->nVars() / (eps * eps);
    if (POOL_SIZE >= maxsat_formula->nSoft()) {
        POOL_SIZE = maxsat_formula->nSoft();
        cout << "the number of sampled clauses: " << POOL_SIZE << " (same as original)" << endl;
    }
    else {
        cout << "the number of sampled clauses: " << POOL_SIZE << " (less than original clauses " << maxsat_formula->nSoft() << ")" << endl;
    }
    vector<int> b = vector<int>(maxsat_formula->nSoft());
    unordered_set<uint32_t> in_pool;
    if (POOL_SIZE == maxsat_formula->nSoft()) {
        std::iota(b.begin(), b.end(), 0);
    }
    else if (maxsat_formula->nSoft() - POOL_SIZE < POOL_SIZE) {
        std::iota(b.begin(), b.end(), 0);
        vector<int> temp_b = vector<int>();
        in_pool = maxsat_formula->pick_k_clauses(maxsat_formula->nSoft() - POOL_SIZE, false);
        // cout << in_pool.size() << " " << b.size() << " ";
        for (auto itr = b.begin(); itr != b.end(); ++itr) {
            if (in_pool.find(*itr) == in_pool.end()) {
                temp_b.push_back(*itr);
            }
        }
        b.clear();
        b = temp_b;
        assert(b.size() == POOL_SIZE);
    }
    else {
        b.clear();
        in_pool = maxsat_formula->pick_k_clauses(POOL_SIZE, true);
        for (auto itr = in_pool.begin(); itr != in_pool.end(); ++itr) {
            b.push_back(*itr);
        }
        assert(b.size() == POOL_SIZE);
    }
    // b = sample_k_items(maxsat_formula->nSoft(), POOL_SIZE);
    ofstream myfile;
    std::string sampled_maxsat_file = "sampled_" + file_name;
    myfile.open(sampled_maxsat_file);
    myfile << "p wcnf " + to_string(maxsat_formula->nVars()) + " " + to_string(POOL_SIZE) + " " + to_string(maxsat_formula->getHardWeight()) << endl;
    for (auto index: b) {
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
    return;
}

#endif