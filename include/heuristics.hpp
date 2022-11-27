#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "movements.hpp"

#include <omp.h>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <string>
#include <cstring>

struct Data {
    Data() : iter(-1), maxIter(-1), P(nullptr) {}
    Data(int i, int m, float* p) : iter(i), maxIter(m), P(p) {}
    int iter, maxIter;
    float* P;
};

// Greedy construction of a feasible solution
std::tuple<char*, int, char*> GreedyConstruction(
        int m,
        int n,
        const int* C,
        const char* A,
        const float* U,
        float* phi = nullptr,
        Data selection = Data());

// Greedy improvement of a feasible solution through (deep) local search
void GreedyImprovement(
        int m,
        int n,
        const int* C,
        const char* A,
        char* x,
        int z,
        bool deep = true,
        char* column = nullptr);

// ACO for the Set Packing Problem
std::tuple<int, int, int> ACO(
        const int m,
        const int n,
        const int* C,
        const char* A,
        const float* U,
        std::vector<int>& zInits,
        std::vector<int>& zAmels,
        std::vector<int>& zBests,
        float* phi,
        int maxAnts = 15,
        int maxIter = 100,
        bool deep = true,
        bool parallel = true);

#endif /* end of include guard: HEURISTICS_H */
