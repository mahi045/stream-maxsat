#ifndef SAMPLING_h
#define SAMPLING_h

#include <cstdio>
#include <iostream>
#include <stdio.h>
#include <string>
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
    vector<int> b = sample_k_items(maxsat_formula->nSoft(), POOL_SIZE);
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
    stringStream << "./open-wbo_static -print-model -cpu-lim=10 " + sampled_maxsat_file + " > " + "result_" + sampled_maxsat_file;
    // calling the smapled maxsat query
    system(stringStream.str().c_str());
    return;
}

#endif