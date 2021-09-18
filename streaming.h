#ifndef STREAMING_h
#define STREAMING_h
#include "constants.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <istream>
#include <random>
#include <stdio.h>
#include <sstream>

using NSPACE::vec;
using namespace openwbo;
using namespace std;

void streaming_maxsat(MaxSATFormula *maxsat_formula) { 
    POOL_SIZE = K * maxsat_formula->nVars() / (eps * eps);
    BUCKET_SIZE = POOL_SIZE / R;
    // test_update_function(maxsat_formula);
    std::ostringstream stringStream;
    string line, variable;
    string delim = " ";
    int lit;
    maxsat_formula->occurance_list.growTo(2 * maxsat_formula->nVars() + 1, 0);
    maxsat_formula->assignment.growTo(maxsat_formula->nVars(), l_Undef);
    printf("Size of occurance list: %d\n", maxsat_formula->occurance_list.size());
    printf("Size of assignment list: %d\n", maxsat_formula->assignment.size());
    ofstream myfile;
    ifstream resultfile;
    string result_file_name;
    vec<int> incompatible;
    vec<int> agreed;
    int bucket_start;
    int index_bucket, index_pool;
    std::string stream_maxsat_file = "streaming_" + file_name;
    double w;
    for (int i = 0; i < maxsat_formula->nSoft(); i++) {
        for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
            w = maxsat_formula->getSoftClause(i).weight / maxsat_formula->getSoftClause(i).clause.size();
            maxsat_formula->occurance_list[var(maxsat_formula->getSoftClause(i).clause[j])] += w; 
        }
        if (!(i % BUCKET_SIZE)) {
            bucket_start = i;
            myfile.open(stream_maxsat_file);
            myfile << "p wcnf " + to_string(maxsat_formula->nVars()) + " " + to_string(BUCKET_SIZE) + " " + to_string(maxsat_formula->getHardWeight()) << endl;
        }
        myfile << maxsat_formula->getSoftClause(i).weight << " ";
        for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
            if (sign(maxsat_formula->getSoftClause(i).clause[j])) {
                myfile << "-";
            }
            myfile << var(maxsat_formula->getSoftClause(i).clause[j]) + 1 << " ";
        }
        myfile << "0" << endl;
        if (!((i + 1) % BUCKET_SIZE) || i + 1 == maxsat_formula->nSoft()) {
            for (int variable = 1; variable <= maxsat_formula->nVars(); variable++) {
                if (maxsat_formula->assignment[variable] == l_True) {
                    if (ceil(maxsat_formula->occurance_list[2 * variable]) >= ceil(maxsat_formula->occurance_list[2 * variable + 1])) {
                        myfile << ceil(maxsat_formula->occurance_list[2 * variable]) << " " << variable << " " << 0 << endl;
                        myfile << ceil(maxsat_formula->occurance_list[2 * variable + 1]) << " " << -variable << " " << 0 << endl;
                    }
                }
                else if (maxsat_formula->assignment[variable] == l_False) {
                    if (ceil(maxsat_formula->occurance_list[2 * variable]) <= ceil(maxsat_formula->occurance_list[2 * variable + 1])) {
                        myfile << ceil(maxsat_formula->occurance_list[2 * variable]) << " " << variable << " " << 0 << endl;
                        myfile << ceil(maxsat_formula->occurance_list[2 * variable + 1]) << " " << -variable << " " << 0 << endl;
                    }
                }
            }
            myfile.close();
            stringStream.str("");
            stringStream << "./open-wbo_static -print-model -cpu-lim=10 " + stream_maxsat_file + " > " + "result_" + stream_maxsat_file;
            // calling the smapled maxsat query
            system(stringStream.str().c_str());

            // reading the recent maxsat call
            incompatible.clear();
            agreed.clear();
            result_file_name = "result_" + stream_maxsat_file;
            ifstream resultfile2(result_file_name);
            while (getline(resultfile2, line)) {
                if (line.rfind("v ") == 0) {
                    auto start = 0U;
                    auto end = line.find(delim);
                    while (end != std::string::npos) {
                        variable = line.substr(start, end - start);
                        start = end + delim.length();
                        end = line.find(delim, start);
                        if (variable != "v") {
                            lit = stoi(variable);
                            if (lit < 0) {
                                if (maxsat_formula->assignment[abs(lit)] == l_True) {
                                    incompatible.push(lit);
                                }
                                else {
                                    maxsat_formula->assignment[abs(lit)] = l_False;
                                    agreed.push(lit);
                                }
                            }
                            else {
                                if (maxsat_formula->assignment[lit] == l_False) {
                                    incompatible.push(lit);
                                }
                                else {
                                    maxsat_formula->assignment[lit] = l_True;
                                    agreed.push(lit);
                                }
                            }
                        }
                    }
                }
            }
            cout << "Total " << incompatible.size() << " literals are incompatible" << endl;
            cout << "Total " << agreed.size() << " literals are compatibles" << endl;
            if (incompatible.size() > 0) {
                // now invoking maxsat query again
                myfile.open(stream_maxsat_file, std::ios_base::app);
                for (int cla_index = 0; cla_index < maxsat_formula->nPool();
                     cla_index++) {
                    myfile << maxsat_formula->getPoolClause(cla_index).weight << " ";
                    for (int j = 0; j < maxsat_formula->getPoolClause(cla_index).clause.size(); j++) {
                        if (sign(maxsat_formula->getPoolClause(cla_index).clause[j])) {
                            myfile << "-";
                        }
                        myfile << var(maxsat_formula->getPoolClause(cla_index).clause[j]) + 1 << " ";
                    }
                    myfile << "0" << endl;
                }
            }
            if (agreed.size() > 0) {
                // adding agreed literals as hard clauses
                for (auto lit_index = 0; lit_index < agreed.size(); lit_index++) {
                    myfile << maxsat_formula->getHardWeight() << " " << agreed[lit_index] << " 0" << endl;
                }
            }
            myfile.close();
            if (incompatible.size() > 0) {
                stringStream.str("");
                stringStream << "./open-wbo_static -print-model -cpu-lim=10 " +
                                  stream_maxsat_file + " > " + "result_" + stream_maxsat_file;
                // calling the new maxsat query
                system(stringStream.str().c_str());
                result_file_name = "result_" + stream_maxsat_file;
                ifstream resultfile1(result_file_name);
                while (getline(resultfile1, line)) {
                    if (line.rfind("v ") == 0) {
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
                                } else {
                                    maxsat_formula->assignment[lit] = l_True;
                                }
                            }
                        }
                    }
                }
            }


            if (maxsat_formula->clause_seen_so_far + BUCKET_SIZE > POOL_SIZE) {
                // this clauses will be replaced
                int clause_need_replace = ((double) (BUCKET_SIZE) / (BUCKET_SIZE + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
                int remaining_clause = BUCKET_SIZE;
                if (i + 1 == maxsat_formula->nSoft()) {
                    remaining_clause = i + 1 - bucket_start;
                    clause_need_replace = ((double) (remaining_clause) / (remaining_clause + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
                }
                vector<int> replaced_clause_pool = sample_k_items(maxsat_formula->nPool(), clause_need_replace);
                vector<int> replaced_clause_bucket = sample_k_items(remaining_clause, clause_need_replace);
                cout << " From bucket index " << bucket_start << " to " << i << endl;
                for (auto cla_index = 0; cla_index < clause_need_replace; cla_index++) {
                    cout << " Bucket " << index_bucket << "'th clause replace pool " << index_pool << "'th clause" << endl;
                    index_bucket = replaced_clause_bucket[cla_index] + maxsat_formula->clause_seen_so_far - 1;
                    index_pool = replaced_clause_pool[cla_index];
                    maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(index_bucket).weight, 
                        maxsat_formula->getSoftClause(index_bucket).clause, index_pool);

                    
                }
                cout << "Total " << clause_need_replace << " clauses replaced from pool !!!" << endl;
            }
            else {
                for (auto cla_index = maxsat_formula->clause_seen_so_far; cla_index != i + 1; cla_index++) {
                    maxsat_formula->addPoolClause(maxsat_formula->getSoftClause(cla_index).weight, 
                        maxsat_formula->getSoftClause(cla_index).clause);
                }
            }
            maxsat_formula->clause_seen_so_far = i + 1;
        }
    }
}
#endif