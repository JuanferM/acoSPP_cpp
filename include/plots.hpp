#ifndef PLOTS_H
#define PLOTS_H

#include "librarySPP.hpp"

#include <cmath>
#include <string>
// Modified version of https://github.com/alandefreitas/matplotplusplus
#include <matplot/matplot.h>
#include <matplot/util/common.h>

// Plot l'examen d'un run de ACO sur
// une instance
void plotRunACO(
    const std::string instance,
    const std::vector<int>& zInits,
    const std::vector<int>& zAmels,
    const std::vector<int>& zBests,
    const std::vector<float>& probas,
    const int zinit,
    const int zls,
    const int zbest,
    const int done_iter,
    std::string save_path = "",
    bool silent_mode = false);

// Plot des niveaux de phéromones avant
// et après le premier restart
void plotPhiRunACO(
    const std::string instance,
    const int n,
    const float* phi,
    bool before = true,
    std::string save_path = "",
    bool silent_mode = false);

// Plot le bilan de tous les runs de ACO
// sur une instance (plot exactement NUM_DIVISION
// points avec NUM_DIVISION <= NUM_ITER)
void plotAnalyseACO(
    const std::string instance,
    const std::vector<double>& divs,
    const std::vector<int>& zMin,
    const std::vector<double>& zMoy,
    const std::vector<int>& zMax,
    const int allrunzmin,
    const float allrunzmoy,
    const int allrunzmax,
    std::string save_path = "",
    bool silent_mode = false);

// Plot le bilan CPUt pour chaque instance (le
// temps d'exécution moyen d'un run)
void plotCPUt(
    std::vector<std::string>& fnames,
    std::vector<float>& tMoy,
    std::string save_path = "",
    bool silent_mode = false);

#endif /* end of include guard: PLOTS_H */
