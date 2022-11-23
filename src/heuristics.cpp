#include "heuristics.hpp"
#include "librarySPP.hpp"

#include <cmath>
#include <omp.h>

std::tuple<char*, int, char*> GreedyConstruction(
        int m,
        int n,
        const int* C,
        const char* A,
        const float* U) {
    bool valid = true;
    int i(0), j(0), s(0);
    char *x = new char[n], *column = new char[m];
    for(i = 0; i < n; i++) x[i] = 0;

    // Indices of utilities in utilities decreasing order
    size_t *u_order = argsort(n, U); // DON'T FORGET TO DELETE
    // We set the variable with the greatest utility to 1
    x[u_order[0]] = 1;
    // Selecting that variable means that we must select the
    // corresponding column in the matrix A and check if the
    // constraints are still verified
    for(j = 0; j < m; j++) {
        column[j] = A[INDEX(u_order[0], j)];
        s += column[j];
    }

    // Repeat the same process with each utility until constraints
    // are eventually violated
    i = 1;
    while(s != m && i < n) {
        for(j = 0, valid = true; j < m && valid; j++)
            valid = !(column[j] & A[INDEX(u_order[i], j)]);
        for(j = 0, s = 0; j < m && valid; s += column[j], j++)
            column[j] += A[INDEX(u_order[i], j)];
        x[u_order[i++]] = valid;
    }

    delete[] u_order;
    return std::make_tuple(x, dot(n, x, C), column);
}

void GreedyImprovement(
        int m,
        int n,
        const int* C,
        const char* A,
        char* x,
        int* z,
        bool deep,
        char* column) {
    int i(2);
    bool (*f[3])(int, int, const int*, const char*, char*, int*, bool, char*) = {
            zero_oneExchange,
            one_oneExchange,
            two_oneExchange
        };

    // We modify x and z directly (no copy)
    while(i >= 0){
        if(!f[i](m, n, C, A, x, z, deep, column)) i--;
    }
}

void ACO(
        const int m,
        const int n,
        const int* C,
        const char* A,
        const float* U,
        std::vector<int>& zInits,
        std::vector<int>& zAmels,
        std::vector<int>& zBests,
        int nbIter,
        bool deep,
        bool parallel) {
    return;
}
