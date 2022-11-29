#include "heuristics.hpp"
#include "librarySPP.hpp"
#include <cmath>
#include <cstdlib>

struct Data;

std::tuple<char*, int, char*> GreedyConstruction(
        int m,
        int n,
        const int* C,
        const char* A,
        const float* U,
        float* phi,
        Data selection) {
    bool valid(true), valid2(true);
    int i(0), j(0), k(0), s(0);
    float sum_phi;
    std::vector<float> probs;
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
        *selection.P = !selection.iter ? selection.iter
                     : log10(selection.iter) / log10(selection.maxIter);
        // Init roulette wheel probabilities
        sum_phi = 0, probs = std::vector<float>(n);
        for(k = 0; k < n; k++) {
            probs[k] = sum_phi + phi[k];
            sum_phi += phi[k];
        }
    }

    // Repeat the same process with each utility until constraints
    // are eventually violated
    i = 1;
    while(s != m && i < n) {
        if(selection.P && ((float)rand() / (float)RAND_MAX) > *selection.P) {
            float r = (float)rand() / (float)RAND_MAX;

            // Selon la roulette, on vérifie si une variable candidate existe
            for(k = 0, valid = false; k < n && !(valid && valid2); k++) {
                for(j = 0, valid = true; j < m && valid && valid2; j++) {
                    valid = !(column[j] & A[INDEX(k, j)]);
                    valid2 = (r < probs[k]/sum_phi);
                }

                // si la variable candidate existe, la prendre en compte
                for(j = 0, s = 0; j < m && valid && valid2; s += column[j], j++)
                    column[j] += A[INDEX(k, j)];
                if(valid2) x[k] = valid, probs[k] = -1;
            }

            i++;
        } else {
            for(j = 0, valid = true; j < m && valid; j++)
                valid = !(column[j] & A[INDEX(u_order[i], j)]);
            for(j = 0, s = 0; j < m && valid; s += column[j], j++)
                column[j] += A[INDEX(u_order[i], j)];
            if(selection.P) probs[u_order[i]] = -1;
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

void managePheromones(
        const int n,
        const float rhoE,
        const float rhoD,
        const float phiNul,
        const int iter,
        const int maxIter,
        const int itStag,
        float* phi,
        float** phi_bef,
        float** phi_aft,
        char* xbest_iter,
        int* nbRestart,
        bool capturePhi) {
    int i(0);
    bool existPhiNul(false);

    for(i = 0; i < n; i++) {
        // Pheromone evaporation
        phi[i] = phi[i] * rhoE;
        // Pheromone deposition
        if(xbest_iter[i]) phi[i] = phi[i] + rhoD;
        if(phi[i] <= phiNul) existPhiNul = true;
    }

    // Territory disturbance
    // std::cout << itStag << " " << existPhiNul << std::endl;
    if(itStag == 0 && existPhiNul) {
        (*nbRestart)++;
        // If phi_bef not initialized, copy phi into phi_bef
        if(!(*phi_bef) && capturePhi) {
            *phi_bef = new float[n];
            std::copy(phi, phi+n, *phi_bef);
        }
        int pn = rand() % (int)ceil(0.1 * n);

        // Disturb the pheromones
        for(i = 0; i < n; i++) {
            phi[i] = phi[i] * 0.95 * (!iter ? iter : log10(iter)/log10(maxIter));
            if(i < pn) {
                float r = (float)rand() / (float)RAND_MAX,
                      d = (1.0 - (float)iter/(float)maxIter)*0.5 - 0.05;
                phi[rand() % n] = 0.05 + r * d;
            }
            // Offset on the pheromones with low level
            if(phi[i] < 0.1) {
                float r = (float)rand() / (float)RAND_MAX,
                      d = (1.0 - (float)iter/(float)maxIter)*0.5 - 0.05;
                phi[i] = phi[i] + 0.05 + r * d;
            }
        }

        // If phi_aft not initialized, copy phi into phi_aft
        if(!(*phi_aft) && capturePhi) {
            *phi_aft = new float[n];
            std::copy(phi, phi+n, *phi_aft);
        }
    }
}

std::tuple<int, int, int, int> ACO(
        const int m,
        const int n,
        const int* C,
        const char* A,
        const float* U,
        std::vector<int>& zInits,
        std::vector<int>& zAmels,
        std::vector<int>& zBests,
        std::vector<float>& probas,
        float* phi,
        float** phi_bef,
        float** phi_aft,
        const float phiNul,
        const float rhoE,
        const float rhoD,
        const int iterStagnant,
        int maxAnts,
        int maxIter,
        bool deep,
        bool parallel,
        bool restartStop,
        int maxRestart,
        bool capturePhi) {
    float P(0); // ant's curiosity level
    bool keep_going(true);
    int zinit, zls, zbest, zbest_iter;
    int nbRestart(0), iter(0), itStag(iterStagnant);
    char *xbest(nullptr), *column(nullptr), *xbest_iter(nullptr);
    // Elaborate greedy solution
    std::tie(xbest, zinit, column) = GreedyConstruction(m, n, C, A, U);
    zbest = zinit;
    // Local search
    GreedyImprovement(m, n, C, A, xbest, &zbest, deep, column);
    if(column) delete[] column, column = nullptr;
    zls = zbest;

    for(iter = 0; iter < maxIter && keep_going; iter++) {
        xbest_iter = new char[n], zbest_iter = -1;

        // std::cout << iter << std::endl;
        #pragma omp parallel for if(parallel)
        for(int ant = 0; ant < maxAnts; ant++) {
            char *x(nullptr), *column_ant(nullptr);
            bool condition = ((float)rand() / (float)RAND_MAX) >= P;
            Data selection = Data(iter, maxIter, &P);

            // std::cout << "A" << std::endl;
            std::tie(x, zInits[iter * maxAnts + ant], column_ant) =
                condition ? GreedyConstruction(m, n, C, A, U, phi, selection)
                          : GreedyConstruction(m, n, C, A, U, phi);
            // std::cout << "B" << std::endl;
            probas[iter * maxAnts + ant] = P;
            zAmels[iter * maxAnts + ant] = zInits[iter * maxAnts + ant];
            // std::cout << "C" << std::endl;
            GreedyImprovement(m, n, C, A, x, &zAmels[iter * maxAnts + ant],
                    deep, column_ant);
            if(column_ant) delete[] column_ant, column_ant = nullptr;

            #pragma omp critical
            if(zAmels[iter * maxAnts + ant] > zbest_iter) {
                // std::cout << "D" << std::endl;
                std::copy(x, x+n, xbest_iter);
                zbest_iter = zAmels[iter * maxAnts + ant];
                if(zbest_iter > zbest) {
                    // std::cout << "E" << std::endl;
                    itStag = iterStagnant;
                    zbest = zbest_iter;
                    std::copy(xbest_iter, xbest_iter+n, xbest);
                }
            }

            if(x) delete[] x, x = nullptr;
        }

        // std::cout << "F" << std::endl;
        managePheromones(n, rhoE, rhoD, phiNul, iter, maxIter, itStag--,
                        phi, phi_bef, phi_aft, xbest_iter, &nbRestart,
                        capturePhi);
        if(itStag <= 0) itStag = 0;
        if(restartStop && nbRestart == maxRestart) keep_going = false;
        if(xbest_iter) delete[] xbest_iter, xbest_iter = nullptr;
    }

    // Compute zBests using zAmels
    zbest_iter = std::max(zls, zAmels[0]), zBests[0] = zbest_iter;
    for(itStag = 1; itStag < maxIter*maxAnts; itStag++) {
        zbest_iter = std::max(zbest_iter, zAmels[itStag]);
        zBests[itStag] = zbest_iter;
    }

    /* MOST IMPORTANT SECTION */
    // Ne pas supprimer xbest et la renvoyer si l'on souhaite garder
    // la meilleure solution
    if(xbest) delete[] xbest, xbest = nullptr;
    return std::make_tuple(zinit, zls, zbest, iter*maxAnts);
}
