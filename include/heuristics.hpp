#ifndef HEURISTICS_H
#define HEURISTICS_H

#include "movements.hpp"

#include <string>
#include <cstring>

// Greedy construction of a feasible solution
std::tuple<char*, int, char*> GreedyConstruction(
        int m,
        int n,
        const int* C,
        const char* A,
        const float* U);

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
void ACO(
        const int m,
        const int n,
        const int* C,
        const char* A,
        const float* U,
        std::vector<int>& zInits,
        std::vector<int>& zAmels,
        std::vector<int>& zBests,
        int nbIter = 100,
        bool deep = true,
        bool parallel = true);

#endif /* end of include guard: HEURISTICS_H */
