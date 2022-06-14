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
#include <algorithm>

using NSPACE::vec;
using namespace openwbo;
using namespace std;

int parseLine(char* line){
    // This assumes that a digit will be found and the line ends in " Kb".
    int i = strlen(line);
    const char* p = line;
    while (*p <'0' || *p > '9') p++;
    line[i-3] = '\0';
    i = atoi(p);
    return i;
}

int currentUsedSizeinVM(MaxSATFormula *maxsat_formula){ //Note: this value is in KB!
    // FILE* file = fopen("/proc/self/status", "r");
    // int result = -1;
    // char line[128];

    // while (fgets(line, 128, file) != NULL){
    //     if (strncmp(line, "VmSize:", 7) == 0){
    //         result = parseLine(line);
    //         break;
    //     }
    // }
    // fclose(file);
    // return result;
    uint64_t memory_consumed = maxsat_formula->effective_pool_size;
    memory_consumed += maxsat_formula->last_index_in_pool * sizeof(Soft);
    memory_consumed += maxsat_formula->memory_consumed_by_bucket;
    memory_consumed += sizeof(maxsat_formula);
    memory_consumed += sizeof(maxsat_formula->assignment[0]) * maxsat_formula->assignment.size();
    memory_consumed += sizeof(maxsat_formula->occurance_list[0]) * maxsat_formula->occurance_list.size();
    memory_consumed += sizeof(maxsat_formula->seen[0]) * maxsat_formula->seen.size();
    return memory_consumed; // the value is in bytes
}

bool init_stream(MaxSATFormula *maxsat_formula, uint64_t var, uint64_t cla) {
    // POOL_SIZE = min((uint64_t) (K * var / (eps * eps)), cla);
    
    // maxsat_formula->occurance_list.growTo(2 * var + 1, 0.0);
    // if (median_heu)
    //     maxsat_formula->occurance_F.resize(var + 1, 0.0);
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
    POOL_SIZE = (total_memory * fraction_of_memory * 1000 * 1000) / (4 * (maxsat_formula->beta) + sizeof(Soft));
    POOL_SIZE = min(POOL_SIZE, cla);
    BUCKET_SIZE = (total_memory * fraction_of_memory_bucket * 1000 * 1000); // it is size in terms of MB
    cout << "BUCKET_SIZE: (in bytes)" << BUCKET_SIZE << endl;
    cout << "The pool size is: " << POOL_SIZE << ", which is " << (double) POOL_SIZE / var << " factor of n" << endl;
    cout << "The number of clauses is " <<  (double) cla / var << " factor of n" << endl;
    // setting the capacity of pool
    maxsat_formula->assignment.growTo(var + 1, l_Undef);
    if (POOL_SIZE == cla) {
        cout << "No streaming algorithm !!!" << endl;
        return true; // no streaming algorithm
    }
    cout << "Run streaming algorithm !!!" << endl;
    maxsat_formula->occurance_list.growTo(2 * var + 1, 0.0);
    // if (median_heu)
    //     maxsat_formula->occurance_F.resize(var + 1, 0.0);
    maxsat_formula->createPool(POOL_SIZE);
    maxsat_formula->seen.assign(var + 1, false);

    // maxsat_formula->var_bias.growTo(var + 1, 0);
    printf("Size of occurance list: %d\n", maxsat_formula->occurance_list.size());
    printf("Size of assignment list: %d\n", maxsat_formula->assignment.size());
    // maxsat_formula->weight_pool.clear();
    return false;
}

// double bias_threshold(MaxSATFormula *maxsat_formula) {
//     double sum = 0;
//     double coff;
//     for (int k = 2; k <= maxsat_formula->m.size(); k++) {
//         if (maxsat_formula->m[k] == 0) {
//             continue;
//         }
//         coff = (double) (pow(2, k) - k - 1) / (pow(2, k-2));
//         sum += coff * maxsat_formula->m[k];
//     }
//     return sum;
// }

void run_maxsat_solver(MaxSATFormula *maxsat_formula) {
    string result_file_name;
    std::ostringstream stringStream;
    stringStream.str("");
    ofstream assignfile;
    string delim = " ", line, variable;
    int lit;
    std::chrono::duration<double> remaining_time;
    uint32_t remaining_time_second, timeout;
    std::chrono::high_resolution_clock::time_point current_time;
    current_time = std::chrono::high_resolution_clock::now();
    std::string stream_maxsat_file = "streaming_" + file_name;
    remaining_time =  current_time - start_time;
    remaining_time_second = floor((TIMEOUT - remaining_time.count() - 5));
    timeout = min(TIMEOUT, remaining_time_second);
    timeout = (timeout < 10) ? 10 : timeout;
    int available_memory = total_memory;
    if (use_fixed_memory)
    {
        int used_memory = 0; // it is non-stream case
        available_memory = (available_memory > used_memory) ? (available_memory - used_memory) : available_memory;
    }
    cout << "The available memory and time limit (only maxsat call): " << available_memory << " and " << timeout << endl;
    stringStream << "./open-wbo_static -print-model -cpu-lim=" << timeout << " -mem-lim=" << available_memory << " " << file_name + " > " + "result_" + stream_maxsat_file;
    // calling the smapled maxsat query
    system(stringStream.str().c_str());
    // renaming the output file
    stringStream.str("");
    std::string open_wbo_maxsat_file = "result_open_wbo_" + file_name;
    stringStream << "mv " << open_wbo_maxsat_file << " result_" + stream_maxsat_file;
    system(stringStream.str().c_str());

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
    ofstream myfile, poolfile, assignfile, debugfile;
    ifstream resultfile;
    string result_file_name;
    vec<int> incompatible;
    vec<int> agreed;
    int bucket_start;
    int index_bucket, index_pool;
    double positive_phase, negative_phase;
    // maxsat_formula->weight_sampler.clear();
    std::string stream_maxsat_file = "streaming_" + file_name;
    std::string pool_stream_maxsat_file = "pool_streaming_" + file_name;
    std::chrono::high_resolution_clock::time_point current_time;
    uint32_t remaining_buckets;
    std::chrono::duration<double> remaining_time;
    uint32_t remaining_time_second, timeout;
    debugfile.open("debug_" + file_name);
    double w, w_adj;
    int var_ind = 0;
    // int bucket_index = maxsat_formula->nSoft() / BUCKET_SIZE - 1;
    int bound = maxsat_formula->nSoft();
    maxsat_formula->in_bucket.assign(maxsat_formula->nVars() + 1, false);
    // cout << "sizeof(maxsat_formula->getSoftClause(0).weight) => " << sizeof(Soft) << endl;
    // cout << "sizeof(maxsat_formula->getSoftClause(0).weight) => " << maxsat_formula->getSoftClause(0).clause.size() << endl;
    for (int i = 0; i < bound; i++) {
        for (int j = 0; j < maxsat_formula->getSoftClause(i).clause.size(); j++) {
            w = (double) maxsat_formula->getSoftClause(i).weight / pow(2, maxsat_formula->getSoftClause(i).clause.size() - 1);
            // if (maxsat_formula->m.size() < maxsat_formula->getSoftClause(i).clause.size() + 1) {
            //     maxsat_formula->m.resize(maxsat_formula->getSoftClause(i).clause.size() + 1, 0);
            //     maxsat_formula->m[maxsat_formula->getSoftClause(i).clause.size()] = maxsat_formula->getSoftClause(i).weight;
            // }
            // else {
            //     maxsat_formula->m[maxsat_formula->getSoftClause(i).clause.size()] += maxsat_formula->getSoftClause(i).weight;
            // }
            var_ind = var(maxsat_formula->getSoftClause(i).clause[j]) * 2;
            if (sign(maxsat_formula->getSoftClause(i).clause[j])) {
                var_ind += 1;
            }
            maxsat_formula->occurance_list[var_ind] += w; 
            maxsat_formula->temp_occurance_list[var_ind] += w; 
            maxsat_formula->seen[var_ind / 2] = true; 
            maxsat_formula->in_bucket[var_ind / 2] = true; 
            // if (sign(maxsat_formula->getSoftClause(i).clause[j])) {
            //     w *= -1;
            // }
            // maxsat_formula->bias += (abs(maxsat_formula->var_bias[var_ind/2] + w) - abs(maxsat_formula->var_bias[var_ind/2]));
            // maxsat_formula->var_bias[var_ind/2] += w;
            // if (maxsat_formula->hard_clause_identifier <= static_cast<uint64_t>(ceil(maxsat_formula->occurance_list[var_ind]))) {
            //     maxsat_formula->hard_clause_identifier = static_cast<uint64_t>(ceil(maxsat_formula->occurance_list[var_ind]) + 2);
            // }
            // the median heuristic was here
            // here is the previous occurance_F
            // if (median_heu) {
            //     if (var_ind % 2 == 0) {
            //         double f = maxsat_formula->occurance_list[var_ind] >= maxsat_formula->occurance_list[var_ind + 1] ? 
            //             (maxsat_formula->occurance_list[var_ind] / maxsat_formula->occurance_list[var_ind + 1]) :
            //             maxsat_formula->occurance_list[var_ind + 1] / maxsat_formula->occurance_list[var_ind];
            //         if (isinf(f))  {
            //             maxsat_formula->occurance_F[var_ind/2] = 
            //             maxsat_formula->occurance_list[var_ind] >= maxsat_formula->occurance_list[var_ind + 1] ? 
            //             maxsat_formula->occurance_list[var_ind] : maxsat_formula->occurance_list[var_ind + 1];
            //         }
            //         else {
            //             maxsat_formula->occurance_F[var_ind/2] = f;
            //         }
            //     }
            //     else {
            //         double f = maxsat_formula->occurance_list[var_ind] >= maxsat_formula->occurance_list[var_ind - 1] ? 
            //             (maxsat_formula->occurance_list[var_ind] / maxsat_formula->occurance_list[var_ind - 1]) :
            //             maxsat_formula->occurance_list[var_ind - 1] / maxsat_formula->occurance_list[var_ind];
            //         if (isinf(f))  {
            //             maxsat_formula->occurance_F[var_ind/2] = 
            //             maxsat_formula->occurance_list[var_ind] >= maxsat_formula->occurance_list[var_ind - 1] ? 
            //             maxsat_formula->occurance_list[var_ind] : maxsat_formula->occurance_list[var_ind - 1];
            //         }
            //         else {
            //             maxsat_formula->occurance_F[var_ind/2] = f;
            //         }
            //     }
            // }
            // the median heuristic ends here
            // if (maxsat_formula->hard_clause_identifier <= ceil(maxsat_formula->occurance_list[var(maxsat_formula->getSoftClause(i).clause[j])])) {
            //     assert(false);
            // }
        }
    }
    bucket_start = 0;
    myfile.open(stream_maxsat_file);
    poolfile.open(pool_stream_maxsat_file);

    myfile << "p wcnf " + to_string(maxsat_formula->nVars()) + " " + to_string(maxsat_formula->nSoft()) + " " + to_string(maxsat_formula->hard_clause_identifier) << endl;
    if (use_pool) poolfile << "p wcnf " + to_string(maxsat_formula->nVars()) + " " + to_string(BUCKET_SIZE) + " " + to_string(maxsat_formula->hard_clause_identifier) << endl;
    for (auto start_index = bucket_start; start_index < maxsat_formula->nSoft();
            start_index++) {
        myfile << maxsat_formula->getSoftClause(start_index).weight << " ";
        if (use_pool) poolfile << maxsat_formula->getSoftClause(start_index).weight << " ";
        for (int j = 0;
            j < maxsat_formula->getSoftClause(start_index).clause.size(); j++) {
                if (sign(maxsat_formula->getSoftClause(start_index).clause[j])) {
            myfile << "-";
                if (use_pool) poolfile << "-";
                }
                myfile << var(maxsat_formula->getSoftClause(start_index).clause[j]) + 1
                    << " ";
                if (use_pool) poolfile << var(maxsat_formula->getSoftClause(start_index).clause[j]) + 1
                    << " ";
        }
        myfile << "0" << endl;
        if (use_pool) poolfile << "0" << endl;
    }
    // if (decision_heu) {
    //     bias_thre = bias_threshold(maxsat_formula);
    //     gamma = 0;
    //     if (maxsat_formula->bias <= bias_thre) {
    //         gamma = (double) ((maxsat_formula->bias) / (2 * bias_thre));
    //         gamma += 0.5;
    //     }
    // }
    for (int variable = 1; variable <= maxsat_formula->nVars(); variable++) {
        positive_phase = ceil(maxsat_formula->occurance_list[2 * (variable - 1)] - maxsat_formula->temp_occurance_list[2 * (variable - 1)]);
        negative_phase = ceil(maxsat_formula->occurance_list[2 * (variable - 1) + 1] - maxsat_formula->temp_occurance_list[2 * (variable - 1) + 1]);
        if (maxsat_formula->temp_occurance_list[2 * (variable - 1)] <= 1e-5 && maxsat_formula->temp_occurance_list[2 * (variable - 1) + 1] <= 1e-5) {
            // if (maxsat_formula->assignment[variable] == l_True) {
            //     // if (positive_phase > negative_phase) {
            //         // myfile << static_cast<uint64_t>(maxsat_formula->hard_clause_identifier) << " " << variable << " " << 0 << endl;
            //         number_of_hard_clause++;
            //         // myfile << static_cast<uint64_t>(negative_phase) << " " << -variable << " " << 0 << endl;
            //     // }
            // }
            // else if (maxsat_formula->assignment[variable] == l_False) {
            //     // if (positive_phase < negative_phase) {
            //     //     myfile << static_cast<uint64_t>(positive_phase) << " " << variable << " " << 0 << endl;
            //         // myfile << static_cast<uint64_t>(maxsat_formula->hard_clause_identifier) << " " << -variable << " " << 0 << endl;
            //         number_of_hard_clause++;
            //     // }
            // }
        }
        else if (positive_phase > 0 || negative_phase > 0) {
            // if (maxsat_formula->assignment[variable] == l_True) {
            //     if (positive_phase > negative_phase) {
            //         myfile << static_cast<uint64_t>(positive_phase) << " " << variable << " " << 0 << endl;
            //         myfile << static_cast<uint64_t>(negative_phase) << " " << -variable << " " << 0 << endl;
            //     }
            // }
            // else if (maxsat_formula->assignment[variable] == l_False) {
            //     if (positive_phase < negative_phase) {
            //         myfile << static_cast<uint64_t>(positive_phase) << " " << variable << " " << 0 << endl;
            //         myfile << static_cast<uint64_t>(negative_phase) << " " << -variable << " " << 0 << endl;
            //     }
            // }
            // disable assignment heuristic
            if (!use_pool && decision_heu) {
                if (true) {
                    if (maxsat_formula->assignment[variable] == l_True) {
                        if (maxsat_formula->occurance_list[2*(variable - 1)] > F * maxsat_formula->occurance_list[2*(variable - 1) + 1]) {
                            myfile << static_cast<uint64_t>(positive_phase-negative_phase) << " " << variable << " " << 0 << endl;
                        // myfile << static_cast<uint64_t>(negative_phase) << " " << -variable << " " << 0 << endl;
                        }
                    } else if (maxsat_formula->assignment[variable] == l_False) {
                        if (F * maxsat_formula->occurance_list[2*(variable - 1)] < maxsat_formula->occurance_list[2*(variable - 1) + 1]) {
                            // myfile << static_cast<uint64_t>(positive_phase) << " " << variable << " " << 0 << endl;
                            myfile << static_cast<uint64_t>(negative_phase-positive_phase) << " " << -variable << " " << 0 << endl;
                        }
                    }
                }
                // else {
                //     if (maxsat_formula->assignment[variable] == l_True) {
                //         if (maxsat_formula->var_bias[variable - 1] > 0) {
                //             myfile << static_cast<uint64_t>(gamma * positive_phase + 1) << " " << variable << " " << 0 << endl;
                //             // myfile << static_cast<uint64_t>((1 - gamma) * negative_phase) << " " << -variable << " " << 0 << endl;
                //         }
                //     } else if (maxsat_formula->assignment[variable] == l_False) {
                //         if (maxsat_formula->occurance_list[2*(variable - 1)] > maxsat_formula->occurance_list[2*(variable - 1) + 1]) {
                //             // myfile << static_cast<uint64_t>((1 - gamma) * positive_phase) << " " << variable << " " << 0 << endl;
                //             myfile << static_cast<uint64_t>(gamma * negative_phase + 1) << " " << -variable << " " << 0 << endl;
                //         }
                //     }
                // }
            }
            // here the heuristic stops
        }
        // reset the state of temp occurence list
        maxsat_formula->temp_occurance_list[2 * (variable - 1)] = 0;
        maxsat_formula->temp_occurance_list[2 * (variable - 1) + 1] = 0;
    }
    maxsat_formula->temp_occurance_list.clear(true);
    myfile.close();
    stringStream.str("");
    current_time = std::chrono::high_resolution_clock::now();
    remaining_buckets = (nbuckets - maxsat_formula->nSoft() / BUCKET_SIZE + 1);
    remaining_time =  current_time - start_time;
    remaining_time_second = ceil((TIMEOUT - remaining_time.count()) / (remaining_buckets + remaining_buckets));
    timeout = min(SMALL_TIMEOUT, remaining_time_second);
    timeout = (timeout < 10) ? 10 : timeout;
    incompatible.clear(true);
    agreed.clear(true);
    int available_memory = total_memory;
    if (use_fixed_memory)
    {
        int used_memory = currentUsedSizeinVM(maxsat_formula);
        used_memory += sizeof(maxsat_formula->in_bucket[0]) * maxsat_formula->in_bucket.size();
        used_memory = used_memory / (1024 * 1024);
        available_memory = (available_memory > used_memory) ? (available_memory - used_memory) : available_memory;
    }
    cout << "The available memory (1st maxsat call): " << available_memory << endl;
    cout << "Calling maxsat query from clause = " << maxsat_formula->clause_seen_so_far << " to clause = " << maxsat_formula->clause_seen_so_far +
                            maxsat_formula->nSoft() << endl;
    stringStream.str("");
    cout << "The number of clauses in the stream: (1st maxsat)";
    stringStream << "wc -l " << stream_maxsat_file;
    system(stringStream.str().c_str());
    // cout << "The number of hard clauses is: " << number_of_hard_clause << endl;

    stringStream.str("");
    stringStream << "./open-wbo_static -print-model -cpu-lim=" << timeout << " -mem-lim=" << available_memory << " " << stream_maxsat_file + " > " + "result_" + stream_maxsat_file;
    // calling the smapled maxsat query
    cout << stringStream.str() << endl;
    system(stringStream.str().c_str());

    cout << "The memory used already:" << endl;
    stringStream.str("");
    stringStream << "grep \"The already used memory is\" " << "result_" + stream_maxsat_file;
    system(stringStream.str().c_str());

    stringStream.str("");
    std::string open_wbo_maxsat_file = "result_open_wbo_" + stream_maxsat_file;
    stringStream << "mv " << open_wbo_maxsat_file << " result_" + stream_maxsat_file;
    system(stringStream.str().c_str());
    open_wbo_maxsat_file.clear();

    // reading the recent maxsat call
    cout << "Reading the result ..." << endl;
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
                    if (!maxsat_formula->in_bucket[abs(lit)]) {
                        continue;   // no need to update the assignment because it is not in bucket
                    }
                    if (lit < 0) {
                        if (use_pool && maxsat_formula->assignment[abs(lit)] == l_True) {
                            incompatible.push(lit);
                        }
                        else {
                            maxsat_formula->assignment[abs(lit)] = l_False;
                            if (use_pool)
                                agreed.push(lit);
                        }
                    }
                    else {
                        if (use_pool && maxsat_formula->assignment[lit] == l_False) {
                            incompatible.push(lit);
                        }
                        else {
                            maxsat_formula->assignment[lit] = l_True;
                            if (use_pool)
                                agreed.push(lit);
                        }
                    }
                }
            }
        }
    }
    maxsat_formula->in_bucket.clear();
    maxsat_formula->in_bucket.shrink_to_fit();
    resultfile2.close();
    cout << "Already read the result ..." << endl;
    if (!no_assign) {
        cout << " I found no assignment" << endl;
        // exit(1);
        number_of_no_assignment++;
        cout << "c The number of no assignment is: " << number_of_no_assignment << endl;
        if (stop_after_memout) {
            exit(1);
        }
    }
    if (!use_pool) {
        maxsat_formula->clause_seen_so_far += maxsat_formula->nSoft();
        assignfile.open("result_" + stream_maxsat_file);
        assignfile << "v ";
        for (int variable = 1; variable <= maxsat_formula->nVars(); variable++) {
            if (maxsat_formula->assignment[variable] == l_True) {
                assignfile << variable << " ";
            } else if (maxsat_formula->assignment[variable] == l_False) {
                assignfile << -variable << " ";
            }
        }
        assignfile.close();
        cout << " Not using the pool " << endl;
        maxsat_formula->memory_consumed_by_bucket = 0; // the bucket is processed and the bucket size is initialized to zero 
        maxsat_formula->bucket_index++;
        return;
    }
    // if (incompatible.size() > 0) {
    //     // now invoking maxsat query again
    //     // myfile.open(stream_maxsat_file, std::ios_base::app);
    //     cout << "maxsat_formula->nPool(): " << maxsat_formula->nPool() << endl;
        for (int cla_index = 0; cla_index < maxsat_formula->nPool();
                cla_index++) {
            poolfile << maxsat_formula->getPoolClause(cla_index).weight << " ";
            for (int j = 0; j < maxsat_formula->getPoolClause(cla_index).clause.size(); j++) {
                if (sign(maxsat_formula->getPoolClause(cla_index).clause[j])) {
                    poolfile << "-";
                }
                poolfile << var(maxsat_formula->getPoolClause(cla_index).clause[j]) + 1 << " ";
            }
            poolfile << "0" << endl;
        }
    // }
    int c = 0;
    if (median_heu && use_filtering_condition) {
        
        vector<double> temp_f;
        double f = 0;
        for (auto lit_index = 0; lit_index < agreed.size(); lit_index++)
        {
            var_ind = 2 * (abs(agreed[lit_index]) - 1);
            if (agreed[lit_index] > 0)
            {
                f = 0;
                if (maxsat_formula->occurance_list[var_ind] >= maxsat_formula->occurance_list[var_ind + 1])
                {
                    f = maxsat_formula->occurance_list[var_ind] / maxsat_formula->occurance_list[var_ind + 1];
                    if (isinf(f))
                    {
                        f = maxsat_formula->occurance_list[var_ind];
                    }
                    temp_f.push_back(f);
                }
            }
            else if (agreed[lit_index] < 0)
            {
                f = 0;
                if (maxsat_formula->occurance_list[var_ind + 1] >= maxsat_formula->occurance_list[var_ind])
                {
                    f = maxsat_formula->occurance_list[var_ind + 1] / maxsat_formula->occurance_list[var_ind];
                    if (isinf(f))
                    {
                        f = maxsat_formula->occurance_list[var_ind + 1];
                    }
                    temp_f.push_back(f);
                }
            }
        }
        sort(temp_f.begin(), temp_f.end());
        F = temp_f[temp_f.size() * npercentile];
        temp_f.clear();
        temp_f.shrink_to_fit();
    }
    if (agreed.size() > 0 && use_hard)
    {
        // adding agreed literals as hard clauses
        for (auto lit_index = 0; lit_index < agreed.size(); lit_index++)
        {
            int var_ind = 2 * (abs(agreed[lit_index]) - 1);
            if (agreed[lit_index] < 0 && maxsat_formula->occurance_list[var_ind + 1] >= F * maxsat_formula->occurance_list[var_ind])
            {
                // debugfile << agreed[lit_index] << " ";
                // debugfile << maxsat_formula->occurance_list[var_ind + 1] << " ";
                // debugfile << maxsat_formula->occurance_list[var_ind] << endl;
                poolfile << maxsat_formula->hard_clause_identifier << " " << agreed[lit_index] << " 0" << endl;
                c++;
            }
            else if (agreed[lit_index] > 0 && maxsat_formula->occurance_list[var_ind] >= F * maxsat_formula->occurance_list[var_ind + 1])
            {
                // debugfile << agreed[lit_index] << " ";
                // debugfile << maxsat_formula->occurance_list[var_ind] << " ";
                // debugfile << maxsat_formula->occurance_list[var_ind + 1] << endl;
                poolfile << maxsat_formula->hard_clause_identifier << " " << agreed[lit_index] << " 0" << endl;
                c++;
            }
        }
    }
    cout << "Total " << incompatible.size() + (agreed.size() - c) << " (" << c << ") literals are incompatible (compatibles)" << endl;
    bool call_second_maxsat = (incompatible.size() + (agreed.size() - c)) > 0;
    incompatible.clear(true);
    agreed.clear(true);
    poolfile.close();
    // file renaming
    stringStream.str("");
    stringStream << "mv " + pool_stream_maxsat_file + " " + stream_maxsat_file;
    cout << stringStream.str() << endl;
    system(stringStream.str().c_str());
        // calling the new maxsat query
    if (call_second_maxsat && use_pool) {
        current_time = std::chrono::high_resolution_clock::now();
        remaining_buckets = (nbuckets - maxsat_formula->nSoft() / BUCKET_SIZE + 1);
        remaining_time =  current_time - start_time;
        remaining_time_second = ceil((TIMEOUT - remaining_time.count()) / (remaining_buckets + remaining_buckets - 1));
        timeout = min(SMALL_TIMEOUT, remaining_time_second);
        timeout = (timeout < 10) ? 10 : timeout;

        // cout << "The timeout for second maxsat: " << timeout << endl;
        stringStream.str("");

        // mem limit for second maxsat call are same
        // available_memory = total_memory;
        // if (use_fixed_memory)
        // {
        //     int used_memory = currentUsedSizeinVM() / 1024;
        //     available_memory = (available_memory > used_memory) ? (available_memory - used_memory) : available_memory;
        // }
        cout << "The available memory (2nd maxsat call): " << available_memory << endl;
        stringStream << "./open-wbo_static -print-model -cpu-lim=" << timeout << " -mem-lim=" << available_memory << " " << stream_maxsat_file + " > " + "result_" + stream_maxsat_file;
        // calling the new maxsat query
        cout << stringStream.str() << endl;


        // auto start = std::chrono::high_resolution_clock::now();
        system(stringStream.str().c_str());
        // auto end = std::chrono::high_resolution_clock::now();

        // auto int_s = std::chrono::duration_cast<std::chrono::seconds>(end - start);

        // std::cout << "maxsat() elapsed time is " << int_s.count() << " seconds )" << std::endl;
        // cout << "Get result !!!" << endl;
        stringStream.str("");
        open_wbo_maxsat_file = "result_open_wbo_" + stream_maxsat_file;
        stringStream << "mv " << open_wbo_maxsat_file << " result_" + stream_maxsat_file;
        system(stringStream.str().c_str());
        open_wbo_maxsat_file.clear();

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
                        if (!maxsat_formula->seen[abs(lit)]) {
                            continue; // no need to update the assignment because it is not seen
                        }
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
        // cout << "Read result !!!" << endl;
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

    if (maxsat_formula->last_index_in_pool + bound > POOL_SIZE) {
        // this clauses will be replaced
        // int clause_need_replace = ((double) (BUCKET_SIZE) / (BUCKET_SIZE + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
        int clause_need_replace = ceil(((double) (mpz_get_d(maxsat_formula->bucket_clause_weight)) / (mpz_get_d(maxsat_formula->clause_weight_sum))) 
            * maxsat_formula->nPool());

        int clause_already_added = 0;
        int space_in_pool = POOL_SIZE - maxsat_formula->last_index_in_pool;
        for (int cla_index = 0; cla_index < space_in_pool; cla_index++)
        {
            maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(cla_index).weight,
                                             maxsat_formula->getSoftClause(cla_index).clause, maxsat_formula->last_index_in_pool);
            maxsat_formula->last_index_in_pool++;
            if (maxsat_formula->getSoftClause(cla_index).clause.size() <= maxsat_formula->beta)
            {
                mpz_sub_ui(maxsat_formula->bucket_clause_weight, maxsat_formula->bucket_clause_weight, maxsat_formula->getSoftClause(cla_index).weight); // update the weight sum
            }
            mpz_sub_ui(maxsat_formula->clause_weight_sum, maxsat_formula->clause_weight_sum, maxsat_formula->getSoftClause(cla_index).weight);
            //
            clause_already_added++;
        }
        if (clause_already_added > 0)
        {
            maxsat_formula->weight_sampler.erase(maxsat_formula->weight_sampler.begin(), maxsat_formula->weight_sampler.begin() + clause_already_added);
        }
        cout << "clause_already_added: " << clause_already_added << endl;
        int remaining_clause = BUCKET_SIZE;
        // cout << "clause_need_replace: " << clause_need_replace << endl;
        clause_need_replace = min(clause_need_replace, bound) - clause_already_added;
        // cout << "clause_need_replace: " << clause_need_replace << endl;

        // if (false && i + 1 == bound) {
        //     remaining_clause = i + 1;
        //     clause_need_replace = ((double) (remaining_clause) / (remaining_clause + maxsat_formula->clause_seen_so_far)) * maxsat_formula->nPool();
        // }
        unordered_set<uint32_t> replaced_clause_pool = maxsat_formula->pick_k_clauses_from_pool(clause_need_replace);
        unordered_set<uint32_t> replaced_clause_bucket = maxsat_formula->pick_k_clauses(clause_need_replace, true);
        cout << " From bucket index "
                     << maxsat_formula->clause_seen_so_far << " to "
                     << maxsat_formula->clause_seen_so_far +
                            maxsat_formula->nSoft()
                     << " => ";
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
        for (auto cla_index = 0; cla_index < maxsat_formula->nSoft(); cla_index++) {
            maxsat_formula->updatePoolClause(maxsat_formula->getSoftClause(cla_index).weight, 
                maxsat_formula->getSoftClause(cla_index).clause, maxsat_formula->last_index_in_pool);
            maxsat_formula->last_index_in_pool++; // increment the next position in pool (ESSENTIAL)
        }
    }
    maxsat_formula->clause_seen_so_far += maxsat_formula->nSoft();
    mpz_init(maxsat_formula->bucket_clause_weight);
    maxsat_formula->memory_consumed_by_bucket = 0; // the bucket is processed and the bucket size is initialized to zero 
    maxsat_formula->bucket_index++;
    maxsat_formula->weight_sampler.clear();
    maxsat_formula->weight_sampler.shrink_to_fit();
    if (verbose)
        maxsat_formula->status_pool();

    
}
#endif