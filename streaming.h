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
#include <unordered_set>

using NSPACE::vec;
using namespace openwbo;
using namespace std;

void init_stream(MaxSATFormula *maxsat_formula, uint64_t var, uint64_t cla) {
    POOL_SIZE = min((uint64_t) (K * var / (eps * eps)), cla);
    BUCKET_SIZE = POOL_SIZE / R;
    maxsat_formula->occurance_list.growTo(2 * var + 1, 0.0);
    maxsat_formula->unit.growTo(2 * var + 1, 0.0);
    maxsat_formula->assignment.growTo(var + 1, l_Undef);
    printf("Size of occurance list: %d\n", maxsat_formula->occurance_list.size());
    printf("Size of assignment list: %d\n", maxsat_formula->assignment.size());
    maxsat_formula->weight_pool.clear();
}

void streaming_maxsat(MaxSATFormula *maxsat_formula) { 
    // POOL_SIZE = min((int) (K * maxsat_formula->nVars() / (eps * eps)), maxsat_formula->nSoft());
    // BUCKET_SIZE = POOL_SIZE / R;
    // test_update_function(maxsat_formula);
    std::ostringstream stringStream;
    string line, variable;
    string delim = " ";
    int lit;
    maxsat_formula->temp_occurance_list.growTo(2 * maxsat_formula->nVars() + 1, 0.0);
    maxsat_formula->unit_last.growTo(2 * maxsat_formula->nVars() + 1, 0.0);
    ofstream myfile, assignfile;
    ifstream resultfile;
    string result_file_name;
    vec<int> incompatible;
    vec<int> agreed;
    int bucket_start;
    int index_bucket, index_pool;
    double positive_phase, negative_phase, p1, p0;
    double positive_unit_phase, negative_unit_phase;
    maxsat_formula->weight_sampler.clear();
    std::string stream_maxsat_file = "streaming_" + file_name;
    double w;
    int var_ind = 0;
    int bucket_index = maxsat_formula->nSoft() / BUCKET_SIZE - 1;
    int bound = (maxsat_formula->nSoft() % BUCKET_SIZE) ? maxsat_formula->nSoft() % BUCKET_SIZE : BUCKET_SIZE;
    for (int i = 0; i < bound; i++) {
        for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
            w = (double) maxsat_formula->getSoftClause(i).weight / pow(2, maxsat_formula->getSoftClause(i).clause.size() - 1);
            var_ind = var(maxsat_formula->getSoftClause(i).clause[j]) * 2;
            if (sign(maxsat_formula->getSoftClause(i).clause[j])) {
                var_ind += 1;
            }
            if (maxsat_formula->getSoftClause(i).clause.size() == 1) {
                maxsat_formula->unit[var_ind] += w; 
                maxsat_formula->unit_last[var_ind] += w; 
            }
            maxsat_formula->occurance_list[var_ind] += w; 
            maxsat_formula->temp_occurance_list[var_ind] += w; 
            if (maxsat_formula->hard_clause_identifier <= static_cast<uint64_t>(ceil(maxsat_formula->occurance_list[var_ind]))) {
                maxsat_formula->hard_clause_identifier = static_cast<uint64_t>(ceil(maxsat_formula->occurance_list[var_ind]) + 2);
            }
            // if (maxsat_formula->hard_clause_identifier <= ceil(maxsat_formula->occurance_list[var(maxsat_formula->getSoftClause(i).clause[j])])) {
            //     assert(false);
            // }
        }
        if (!(i % BUCKET_SIZE)) {
            bucket_start = i;
            myfile.open(stream_maxsat_file);

        }
        if (!((i + 1) % BUCKET_SIZE) || i + 1 == bound) {
            myfile << "p wcnf " + to_string(maxsat_formula->nVars()) + " " + to_string(BUCKET_SIZE) + " " + to_string(maxsat_formula->hard_clause_identifier) << endl;
            for (auto start_index = bucket_start; start_index <= i;
                 start_index++) {
                myfile << maxsat_formula->getSoftClause(start_index).weight << " ";
                for (int j = 0;
                    j < maxsat_formula->getSoftClause(start_index).clause.size(); j++) {
                        if (sign(maxsat_formula->getSoftClause(start_index).clause[j])) {
                        myfile << "-";
                        }
                        myfile << var(maxsat_formula->getSoftClause(start_index).clause[j]) + 1
                            << " ";
                }
                myfile << "0" << endl;
            }
            for (int variable = 1; variable <= maxsat_formula->nVars(); variable++) {
                positive_phase = ceil(maxsat_formula->occurance_list[2 * (variable - 1)] - maxsat_formula->temp_occurance_list[2 * (variable - 1)]);
                negative_phase = ceil(maxsat_formula->occurance_list[2 * (variable - 1) + 1] - maxsat_formula->temp_occurance_list[2 * (variable - 1) + 1]);
                positive_unit_phase = ceil(maxsat_formula->unit[2 * (variable - 1)] - maxsat_formula->unit_last[2 * (variable - 1)]);
                negative_unit_phase = ceil(maxsat_formula->unit[2 * (variable - 1) + 1] - maxsat_formula->unit_last[2 * (variable - 1) + 1]);
                
                if (maxsat_formula->temp_occurance_list[2 * (variable - 1)] <= 1 && maxsat_formula->temp_occurance_list[2 * (variable - 1) + 1] <= 1) {
                    if (maxsat_formula->assignment[variable] == l_True) {
                        // if (positive_phase > negative_phase) {
                            myfile << static_cast<uint64_t>(positive_phase+positive_unit_phase) << " " << variable << " " << 0 << endl;
                            // myfile << static_cast<uint64_t>(negative_phase) << " " << -variable << " " << 0 << endl;
                        // }
                    }
                    else if (maxsat_formula->assignment[variable] == l_False) {
                        // if (positive_phase < negative_phase) {
                        //     myfile << static_cast<uint64_t>(positive_phase) << " " << variable << " " << 0 << endl;
                            myfile << static_cast<uint64_t>(negative_phase+negative_unit_phase) << " " << -variable << " " << 0 << endl;
                        // }
                    }
                }
                else if (positive_phase > 0 || negative_phase > 0) {
                    double slack = abs(positive_phase+positive_unit_phase - negative_phase - negative_unit_phase);
                    p1 = positive_phase + positive_unit_phase;
                    p0 = negative_phase + negative_unit_phase;
                    if (slack <= 1e-9 || slack >= positive_unit_phase + negative_unit_phase) {
                        if (maxsat_formula->assignment[variable] == l_True) {
                            if (p1 > p0) {
                                myfile << static_cast<uint64_t>(p1) << " " << variable << " " << 0 << endl;
                                myfile << static_cast<uint64_t>(p0) << " " << -variable << " " << 0 << endl;
                            }
                        }
                        else if (maxsat_formula->assignment[variable] == l_False) {
                            if (p1 < p0) {
                                myfile << static_cast<uint64_t>(p1) << " " << variable << " " << 0 << endl;
                                myfile << static_cast<uint64_t>(p0) << " " << -variable << " " << 0 << endl;
                            }
                        }
                    }
                    else {
                        double epsilon_1 = (double) (slack * (positive_unit_phase + negative_unit_phase) - slack * slack) /
                        (slack + slack + positive_phase + negative_phase - positive_unit_phase - negative_unit_phase);
                        if (maxsat_formula->assignment[variable] == l_True) {
                            if (p1 + epsilon_1 > p0 - epsilon_1) {
                                myfile << static_cast<uint64_t>(p1 + epsilon_1) << " " << variable << " " << 0 << endl;
                                myfile << static_cast<uint64_t>(p0 - epsilon_1) << " " << -variable << " " << 0 << endl;
                            }
                        }
                        else if (maxsat_formula->assignment[variable] == l_False) {
                            if (p1 - epsilon_1 < p0 + epsilon_1) {
                                myfile << static_cast<uint64_t>(p1 - epsilon_1) << " " << variable << " " << 0 << endl;
                                myfile << static_cast<uint64_t>(p0 + epsilon_1) << " " << -variable << " " << 0 << endl;
                            }
                        }
                    }
                }
                // reset the state of temp occurence list
                maxsat_formula->temp_occurance_list[2 * (variable - 1)] = 0;
                maxsat_formula->temp_occurance_list[2 * (variable - 1) + 1] = 0;
            }
            myfile.close();
            stringStream.str("");
            cout << "Calling maxsat query from clause = " << bucket_start + bucket_index * BUCKET_SIZE << " to clause = " << i + bucket_index * BUCKET_SIZE << endl;
            stringStream << "./open-wbo_static -print-model -cpu-lim=" << SMALL_TIMEOUT << " " + stream_maxsat_file + " > " + "result_" + stream_maxsat_file;
            // calling the smapled maxsat query
            system(stringStream.str().c_str());

            // reading the recent maxsat call
            incompatible.clear();
            agreed.clear();
            result_file_name = "result_" + stream_maxsat_file;
            ifstream resultfile2(result_file_name);
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
            resultfile2.close();
            if (!no_assign) {
                cout << " I found no assignment";
                exit(1);
            }
            cout << "Total " << incompatible.size() << " (" << agreed.size() << ") literals are incompatible (compatibles)" << endl;
            
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
                    myfile << maxsat_formula->hard_clause_identifier << " " << agreed[lit_index] << " 0" << endl;
                }
            }
            myfile.close();
            if (incompatible.size() > 0) {
                stringStream.str("");
                stringStream << "./open-wbo_static -print-model -cpu-lim=" << SMALL_TIMEOUT << " " <<
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
                resultfile1.close();
            }
            assignfile.open("result_" + stream_maxsat_file);
            assignfile << "v ";
            for (int variable = 1; variable <= maxsat_formula->nVars(); variable++) { 
                if (maxsat_formula->assignment[variable] == l_True) {
                    assignfile << variable << " ";
                }
                else if (maxsat_formula->assignment[variable] == l_False) {
                    assignfile << -variable << " ";
                }
            } 
            assignfile.close();

            if (maxsat_formula->clause_seen_so_far + BUCKET_SIZE > POOL_SIZE) {
                // this clauses will be replaced
                // int clause_need_replace = ((double) (BUCKET_SIZE) / (BUCKET_SIZE + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
                int clause_need_replace = ceil(((double) (mpz_get_d(maxsat_formula->bucket_clause_weight)) / (mpz_get_d(maxsat_formula->clause_weight_sum))) 
                    * maxsat_formula->nPool());
                int remaining_clause = BUCKET_SIZE;
                // cout << "clause_need_replace: " << clause_need_replace << endl;
                clause_need_replace = min(clause_need_replace, bound);
                // cout << "clause_need_replace: " << clause_need_replace << endl;
                if (false && i + 1 == bound) {
                    remaining_clause = i + 1;
                    clause_need_replace = ((double) (remaining_clause) / (remaining_clause + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
                }
                unordered_set<uint32_t> replaced_clause_pool = maxsat_formula->pick_k_clauses_from_pool(clause_need_replace);
                unordered_set<uint32_t> replaced_clause_bucket = maxsat_formula->pick_k_clauses(clause_need_replace, true);
                cout << " From bucket index " << bucket_start + bucket_index * BUCKET_SIZE << " to " << i + bucket_index * BUCKET_SIZE << " => ";
                auto start_itr1 = replaced_clause_pool.begin();
                auto start_itr2 = replaced_clause_bucket.begin();
                assert(replaced_clause_pool.size() == replaced_clause_bucket.size());
                for (;start_itr1 != replaced_clause_pool.end(); start_itr1++, start_itr2++) {
                    // index_bucket = *start_itr2 + maxsat_formula->clause_seen_so_far - 1;
                    index_bucket = *start_itr2;
                    index_pool = *start_itr1;
                    maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(index_bucket).weight, 
                        maxsat_formula->getSoftClause(index_bucket).clause, index_pool);
                }
                // for (auto cla_index = 0; cla_index < clause_need_replace; cla_index++) {
                //     // cout << " Bucket " << index_bucket << "'th clause replace pool " << index_pool << "'th clause" << endl;
                //     index_bucket = replaced_clause_bucket[cla_index] + maxsat_formula->clause_seen_so_far - 1;
                //     index_pool = replaced_clause_pool[cla_index];
                //     maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(index_bucket).weight, 
                //         maxsat_formula->getSoftClause(index_bucket).clause, index_pool);

                    
                // }
                cout << "Total " << clause_need_replace << " clauses replaced from pool !!!" << endl;
            }
            else {
                for (auto cla_index = 0; cla_index != i + 1; cla_index++) {
                    maxsat_formula->addPoolClause(maxsat_formula->getSoftClause(cla_index).weight, 
                        maxsat_formula->getSoftClause(cla_index).clause);
                }
            }
            maxsat_formula->clause_seen_so_far += (i + 1);
            mpz_init(maxsat_formula->bucket_clause_weight);
            maxsat_formula->weight_sampler.clear();
        }
    }
}
#endif