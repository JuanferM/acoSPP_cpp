#ifndef GRASPPLOTS_H
#define GRASPPLOTS_H

#include "librarySPP.hpp"

#include <cmath>
// Modified version of https://github.com/alandefreitas/matplotplusplus
#include <matplot/matplot.h>
#include <matplot/util/common.h>

// Plot l'examen d'un run de GRASP sur
// une instance
void plotRunGRASP(
        const std::string instance,
        const std::vector<int>& zInits,
        const std::vector<int>& zAmels,
        const std::vector<int>& zBests,
        std::string save_path = "",
        bool silent_mode = false);

// Plot des probabilités pour chaque valeur
// de alpha
void plotProbaRunGRASP(
        const std::string instance,
        const std::vector<double>& alpha,
        const std::vector<double>& proba,
        std::string save_path = "",
        bool silent_mode = false);

// Plot le bilan de tous les runs de GRASP
// sur une instance (plot exactement NUM_DIVISION
// points avec NUM_DIVISION <= NUM_ITER)
void plotAnalyseGRASP(
        const std::string instance,
        const std::vector<double>& divs,
        const std::vector<int>& zMin,
        const std::vector<double>& zMoy,
        const std::vector<int>& zMax,
        std::string save_path = "",
        bool silent_mode = false);

// Plot le bilan CPUt pour chaque instance (le
// temps d'exécution moyen d'un run)
void plotCPUt(
        std::vector<std::string>& fnames,
        std::vector<float>& tMoy,
        std::string save_path = "",
        bool silent_mode = false);

#endif /* end of include guard: GRASPPLOTS_H */
