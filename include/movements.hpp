#ifndef MOVEMENTS_H
#define MOVEMENTS_H

#include "librarySPP.hpp"

#include <deque>
#include <forward_list>

template<typename T>
std::forward_list<int> findItems(int size, T* const arr, T target) {
    std::forward_list<int> indices;

    for(--size; size >= 0; size--) // Avoids creating extra variable
        if(arr[size] == target)
            indices.push_front(size);

    return indices;
}

template<typename R>
std::tuple<R, R> find01(int size, char* const arr) {
    R idx0, idx1;

    for(--size; size >= 0; size--) {
        if(arr[size] == 0) idx0.push_front(size);
        if(arr[size] == 1) idx1.push_front(size);
    }

    return std::make_tuple(idx0, idx1);
}

// Implements 01-exchange.
// Returns true if an improved solution is found.
bool zero_oneExchange(
        int m,
        int n,
        const int* C,
        const char* A,
        char* x,
        int* z,
        bool deep,
        char* column = nullptr);

// Implements 11-exchange
// Returns true if an improved solution is found.
bool one_oneExchange(
        int m,
        int n,
        const int* C,
        const char* A,
        char* x,
        int* z,
        bool deep,
        char* column = nullptr);

// Helper function for two_oneExchange()
// compute all non-symmetrical pairs of indices (indices of the variables
// set to 1) recursively so we don't have to store all the combinations
//
// CAUTION : since we use recursion we should be cautious of the program
// stack size. Only things like the return address, the function arguments
// and the variable used in the function are stored in the stack. Since we
// don't initialize too much stuff in the function, the memory used comes (for
// the most part) from storing the return address and function arguments on
// stack. The size of the arguments should be around 70 bytes per call. We
// stop at  recursion_depth == 2  and we make at most N-2 recursive calls
// per level with N the number of ones in x.
// In conclusion : we expect to use around (2N-4) * 70 Bytes of memory.
// With a stack of size 1MB (125000 Bytes) we can use this function up to
// around N = 894 (894 variables set to 1 in x).
bool combinations(
        // heuristic variables
        int m,
        int n,
        const int* C,
        const char* A,
        char* x,
        int* z,
        bool deep,
        char* column,
        bool* stop,
        int* best_z,
        std::tuple<int, int, int>& best_move,
        // indices list (for variables set to 0)
        const std::forward_list<int>& idx0,
        // variables useful for recursion
        int* pair_of_1,
        std::forward_list<int>::iterator start,
        std::forward_list<int>::iterator end,
        int depth);

// Implements 21-exchange
// Returns true if an improved solution is found.
bool two_oneExchange(
        int m,
        int n,
        const int* C,
        const char* A,
        char* x,
        int* z,
        bool deep,
        char* column = nullptr);

#endif /* end of include guard: MOVEMENTS_H */
