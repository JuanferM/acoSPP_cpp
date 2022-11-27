#include "heuristics.hpp"
#include "librarySPP.hpp"

struct Data;

std::tuple<char*, int, char*> GreedyConstruction(
        int m,
        int n,
        const int* C,
        const char* A,
        const float* U,
        float* phi,
        Data selection) {
    bool valid = true;
    int i(0), j(0), k(0), s(0);
    float sum_phi;
    std::vector<float> previous_probs;
    char *x = new char[n], *column = new char[m];
    for(i = 0; i < n; i++) x[i] = 0;

    // Indices of utilities in utilities decreasing order
    std::vector<int> u_order = phi ? argsort(n, phi) : argsort(n, U);
    // We set the variable with the greatest utility to 1
    x[u_order[0]] = 1;
    // Selecting that variable means that we must select the
    // corresponding column in the matrix A and check if the
    // constraints are still verified
    for(j = 0; j < m; j++) {
        column[j] = A[INDEX(u_order[0], j)];
        s += column[j];
    }

    // Si on a bien l'adresse de P, on modifie sa valeur (mode sélection)
    if(selection.P) {
        *selection.P = log10(selection.iter) / log10(selection.maxIter);
        sum_phi = -1, previous_probs = std::vector<float>(n);
    }

    // Repeat the same process with each utility until constraints
    // are eventually violated
    i = 1;
    while(s != m && i < n) {
        if(selection.P && ((double)rand() / RAND_MAX) > *selection.P) {
            // Init roulette wheel probabilities (once)
            bool init = sum_phi == -1;
            if(init) sum_phi = 0.0f;
            for(k = 0; k < n && init; k++) {
                previous_probs[k] = sum_phi + phi[k];
                sum_phi += phi[k];
            }

            k = 0, init = true;
            float r = (double)rand()/RAND_MAX;
            while(k < n && init) {
                if(r < previous_probs[k]/sum_phi) {
                    // Cette variable est choisie selon la roulette
                    // mais est-ce que son introduction engendre un conflit?
                    for(j = 0, valid = true; j < m; j++)
                        valid = !(column[j] & A[INDEX(k, j)]);
                    init = !valid; // Pas de conflit? Ok. Sinon, on vérifie
                                   // avec les autres
                } else k += 1;
            }

            // S'il n'y a pas de variable candidate on ne fait rien. Sinon,
            if(k < n) {
                x[k] = 1;
                for(j = 0, s = 0; j < m && valid; s += column[j], j++)
                    column[j] += A[INDEX(k, j)];
            }
        } else {
            for(j = 0, valid = true; j < m && valid; j++)
                valid = !(column[j] & A[INDEX(u_order[i], j)]);
            for(j = 0, s = 0; j < m && valid; s += column[j], j++)
                column[j] += A[INDEX(u_order[i], j)];
            x[u_order[i++]] = valid;
        }
    }

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

// TODO add missing parameters (n, z*, z**, lastIterRestart, maxIter, iter)
void managePheromones(
        float* phi,
        char* xbest,
        char* xbest_iter) {

}

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
        int maxAnts,
        int maxIter,
        bool deep,
        bool parallel) {
    float P(0); // ant's curiosity level
    int zinit(-1), zls(-1), zbest(-1);
    char *xbest(nullptr), *column(nullptr), *xbest_iter(nullptr);
    // Elaborate greedy solution
    std::tie(xbest, zinit, column) = GreedyConstruction(m, n, C, A, U);
    zbest = zinit;
    // Local search
    GreedyImprovement(m, n, C, A, xbest, &zbest, deep, column);
    if(column) delete[] column, column = nullptr;
    zls = zbest;

    for(int iter = 0, zbest_iter = -1; iter < maxIter; iter++) {
        xbest_iter = new char[n];

        #pragma omp parallel for if(parallel)
        for(int ant = 0; ant < maxAnts; ant++) {
            char *x(nullptr), *column_ant(nullptr);
            bool condition = ((double)rand() / RAND_MAX) > P;
            Data selection = Data(iter, maxIter, &P);

            std::tie(x, zInits[iter * maxAnts + ant], column_ant) =
                condition ? GreedyConstruction(m, n, C, A, U, phi)
                          : GreedyConstruction(m, n, C, A, U, phi, selection);
            zAmels[iter * maxAnts + ant] = zInits[iter * maxAnts + ant];
            GreedyImprovement(m, n, C, A, x, &zAmels[iter * maxAnts + ant],
                    deep, column_ant);

            #pragma omp single
            if(zAmels[iter * maxAnts + ant] > zbest_iter) {
                memcpy(xbest_iter, x, sizeof(char) * n);
                zbest_iter = zAmels[iter * maxAnts + ant];
                // TODO add vector zBestsIter for plots
                if(zbest_iter > zbest) {
                    if(xbest) delete[] xbest, xbest = nullptr;
                    memcpy(xbest, xbest_iter, sizeof(char) * n);
                }
            }

            if(x) delete[] x, x = nullptr;
        }

        managePheromones(phi, xbest, xbest_iter);
        if(xbest_iter) delete[] xbest_iter, xbest_iter = nullptr;
    }

    /* MOST IMPORTANT SECTION */
    // Ne pas supprimer xbest et la renvoyer si l'on souhaite garder
    // la meilleure solution
    if(xbest) delete[] xbest, xbest = nullptr;
    return std::make_tuple(zinit, zls, zbest);
}
