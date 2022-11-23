#ifndef LIBRARYSPP_H
#define LIBRARYSPP_H

#include <tuple>
#include <string>
#include <chrono>
#include <random>
#include <climits>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <filesystem>
#include <cassert>
#include <exception>
#include <stdexcept>

#include <glpk.h>

// Macro wrapping assert to print a message
#define m_assert(expr, msg) assert(( (void)(msg), (expr) ))

// Macros to time an expression
#define __CHRONO_HRC__ std::chrono::high_resolution_clock
#define __DURATION__ std::chrono::duration<float>
// CAUTION : Call INIT_TIMER() just once
#define INIT_TIMER() __CHRONO_HRC__::time_point __m_time_var_a__, __m_time_var_b__;
#define __m_time__() __CHRONO_HRC__::now()
#define TIMED(t, expr) \
    __m_time_var_a__ = __m_time__(); \
    expr; \
    __m_time_var_b__ = __m_time__(); \
    t = __DURATION__(__m_time_var_b__ - __m_time_var_a__).count();

// Macro that returns a 1D index from 2D coordinates (row and col)
#define _INDEX(row, col, nCols) (col * nCols + row)
// Shorter version. CAUTION: only use if  n  is in the scope
#define INDEX(row, col) _INDEX(row, col, n)

// Macros for color printing
#define _CLR     "\u001B[0m"
#define _CLRd    "\033[1;31m"
#define _CLG     "\033[1;32m"
#define _CLB     "\u001B[36m"
#define _CLP     "\033[1;35m"

// Collect the unhidden filenames available in a given folder
std::vector<std::string> getfname(std::string pathtofolder);

// Reads  fname  and returns :
//  *  m  the number of constraints
//  *  n  the number of variables
//  *  C  the vector of coefficients from the objective function
//  *  A  the binary matrix of constraints (as a 1D array)
//  *  U  a vector of utilities computed for each variables
std::tuple<int, int, int*, char*, float*> loadSPP(std::string fname);

// Models the SPP and run GLPK on instance  instance :
void modelSPP(
        std::string fname,
        std::string path = "",
        float* tt = nullptr,
        bool verbose = true);

// Takes  C  ,  A  and  x  and returns :
//  * true if  x  is feasible
//  * false otherwise
bool isFeasible(
        int m,
        int n,
        const int* C,
        const char* A,
        const char* x,
        const char* extColumn = nullptr,
        bool verbose = true);

// Free  C  ,  A  and  U
void freeSPP(int* C, char* A, float* U);

// Computes indirect sort of an array (decreasing order)
template<typename T>
std::vector<int> argsort(int size, const T* arr) {
    // initialize original index locations
    std::vector<int> idx(size);
    std::iota(idx.begin(), idx.begin()+size, 0);

    // sort indexes based on comparing values in arr
    std::stable_sort(idx.begin(), idx.begin()+size,
        [&arr](size_t e1, size_t e2) {
            return arr[e1] > arr[e2];
        }
    );

    return idx;
}

// Computes the dot product of two arrays
template<typename T1, typename T2>
T2 dot(int size, const T1* v1, const T2* v2) {
    T2 product(0);
    for(--size; size >= 0; size--) // Avoids creating extra variable
        product += v1[size] * v2[size];
    return product;
}

void m_print(std::ostream& out);

template<typename T, typename... Args>
void m_print(std::ostream& out, T t, Args... args) {
    out << std::fixed << t; out.flush();
    if(sizeof...(args) != 0) m_print(out, args...);
    out.flush();
}

/* Functions from answer at https://stackoverflow.com/a/16421677 */
template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

#endif /* end of include guard: LIBRARYSPP_H */
