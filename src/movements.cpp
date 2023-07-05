#include "movements.hpp"

bool zero_oneExchange(
    int m,
    int n,
    const int* C,
    const char* A,
    char* x,
    int* z,
    bool deep,
    char* column) {
  int c(-1), tmp_z(-1), best_z(*z), best_move(-1);
  std::forward_list<int> idx0 = findItems<char>(n, x, 0);

  for(int i : idx0) {
    x[i] = 1; tmp_z = *z + C[i];
    if(tmp_z > best_z) {
      for(c = 0; c < m && column; c++) column[c] += A[INDEX(i, c)];
      if(isFeasible(m, n, C, A, x, column, false)) {
        if(deep) {
          best_z = tmp_z; best_move = i;
        }
        else return (*z = tmp_z);
      }
      for(c = 0; c < m && column; c++) column[c] -= A[INDEX(i, c)];
    }
    x[i] = 0;
  }

  if(best_move != -1) {
    x[best_move] = 1, *z = best_z;
    for(c = 0; c < m && column; c++)
      column[c] += A[INDEX(best_move, c)];
  }

  return best_move != -1;
}

bool one_oneExchange(
    int m,
    int n,
    const int* C,
    const char* A,
    char* x,
    int* z,
    bool deep,
    char* column) {
  int c(-1), tmp_z(-1), best_z(*z);
  std::tuple<int, int> best_move(-1, -1);
  std::forward_list<int> idx0, idx1;
  std::tie(idx0, idx1) = find01<std::forward_list<int>>(n, x);

  for(int i : idx1) {
    x[i] = 0;
    for(int j : idx0) {
      if(i != j) {
        x[j] = 1, tmp_z = *z - C[i] + C[j];
        if(tmp_z > best_z) {
          for(c = 0; c < m && column; c++)
            column[c] += A[INDEX(j, c)] - A[INDEX(i, c)];
          if(isFeasible(m, n, C, A, x, column, false)) {
            if(deep)
              best_z = tmp_z, best_move = std::make_tuple(i, j);
            else return (*z = tmp_z);
          }
          for(c = 0; c < m && column; c++)
            column[c] += A[INDEX(i, c)] - A[INDEX(j, c)];
        }
        x[j] = 0;
      }
    }
    x[i] = 1;
  }

  int i, j;
  std::tie(i, j) = best_move;
  if((deep = i != -1 && j != -1)) {
    x[i] = 0, x[j] = 1, *z = best_z;
    for(c = 0; c < m && column; c++)
      column[c] += A[INDEX(j, c)] - A[INDEX(i, c)];
  }

  return deep;
}

void combinations(
    // heuristic variables
    int m,
    int n,
    const int *C,
    const char *A,
    char *x,
    int *z,
    bool deep,
    char* column,
    // indices list (for variables set to 0)
    // and additional variables
    bool *stop,
    int *best_z,
    std::tuple<int, int, int>& best_move,
    const std::deque<int> &idx0,
    // variables useful to build up the combinations
    int *pair_of_1,
    std::deque<int>::iterator start,
    std::deque<int>::iterator end,
    int depth) {
  if(*stop) return ;

  // Current combination is ready to be used
  if(depth == 2) {
    int c(-1);
    x[pair_of_1[0]] = 0, x[pair_of_1[1]] = 0;

    for(int k : idx0) {
      if(pair_of_1[0] != k && pair_of_1[1] != k) {
        // Repurpose depth
        x[k] = 1, depth = *z - C[pair_of_1[0]] - C[pair_of_1[1]] + C[k];
        if(depth > *best_z) {
          for(c = 0; c < m && column; c++)
            column[c] +=  A[INDEX(k, c)]
                  - A[INDEX(pair_of_1[0], c)]
                  - A[INDEX(pair_of_1[1], c)];
          if(isFeasible(m, n, C, A, x, column, false)) {
            if(deep) {
              *best_z = depth;
              best_move = std::make_tuple(pair_of_1[0],
                            pair_of_1[1],
                            k);
            } else {
              *stop = true, *z = depth;
              return ;
            }
          }
          for(c = 0; c < m && column; c++)
            column[c] +=  A[INDEX(pair_of_1[0], c)]
                  + A[INDEX(pair_of_1[1], c)]
                  - A[INDEX(k, c)];
        }
        x[k] = 0;
      }
    }

    x[pair_of_1[0]] = 1, x[pair_of_1[1]] = 1;
  } else {
    // Replace index with all possible elements.
    for(std::deque<int>::iterator it = start;
      it != end && end - it+1 >= 2-depth;
      ++it) {
      pair_of_1[depth] = *it;
      combinations(m, n, C, A, x, z, deep, column, stop,
          best_z, best_move, idx0, pair_of_1, it+1, end, depth+1);
    }
  }
}

bool two_oneExchange(
    int m,
    int n,
    const int* C,
    const char* A,
    char* x,
    int* z,
    bool deep,
    char* column) {
  bool stop(false);
  int best_z(*z);
  int pair_of_1[2] = { -1, -1 };
  std::deque<int> idx0, idx1;
  std::tuple<int, int, int> best_move(-1, -1, -1);
  std::tie(idx0, idx1) = find01<std::deque<int>>(n, x);

  combinations(m, n, C, A, x, z, deep, column, &stop, &best_z,
      best_move, idx0, pair_of_1, idx1.begin(), idx1.end(), 0);

  int i(-1), j(-1), k(-1); std::tie(i, j, k) = best_move;
  if((deep = i != -1 && j != -1 && k != -1)) {
    x[i] = 0, x[j] = 0, x[k] = 1, *z = best_z;
    for(int c = 0; c < m && column; c++)
      column[c] +=  A[INDEX(k, c)]
            - A[INDEX(i, c)]
            - A[INDEX(j, c)];
  }

  return deep;
}
