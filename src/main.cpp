#include "plots.hpp"
#include "heuristics.hpp"

#include <omp.h>

// Paramètres GLPK
#define USE_GLPK        false
#define VERBOSE_GLPK    false

// Paramètres OpenMP
#define PARALLEL        true
#define MAX_THREADS     10

// Paramètres ACO
#define NUM_RUN         1
#define NUM_ITER        200

#define NUM_DIVISION    20
#define DEEPSEARCH      true

// Paramètres plot
#define INTERACTIVE     false
#define SILENT_MODE     false
#define PATH_PLOT       ""

int main() {
    // This program will create different sequence of
    // random numbers on every program run
    // Use current time as seed for random generator
    std::string pwd(std::filesystem::current_path());
    std::string path(pwd + "/instances/");
    std::cout.precision(3);

    m_print(std::cout, _CLRd, "Etudiants : MERCIER et PICHON\n", _CLR);
    #if !USE_GLPK
        #if NUM_ITER < 2 // We need at least two iterations or else the plots
            #undef NUM_ITER //break
            #define NUM_ITER 2
        #endif
        #if NUM_RUN < 1
            #undef NUM_RUN
            #define NUM_RUN 1
        #endif
        #if MAX_THREADS < 1
            #undef MAX_THREADS
            #define MAX_THREADS 10
        #endif

        INIT_TIMER(); srand(time(NULL));
        if(PARALLEL) omp_set_num_threads(MAX_THREADS);

        const int _NBD_ = NUM_DIVISION > NUM_ITER ? NUM_ITER : NUM_DIVISION;
        m_print(std::cout, _CLP, "\nnombre de runs\t\t: ", NUM_RUN);
        m_print(std::cout, "\nnombre d'itérations\t: ", NUM_ITER);
        m_print(std::cout, "\nparallélisation\t\t: ", (PARALLEL ? "oui" : "non"));
        if(PARALLEL)
            m_print(std::cout, "\nnombre de threads\t: ", MAX_THREADS);
        m_print(std::cout, "\ndescente profonde\t: ", (DEEPSEARCH ? "oui" : "non"));
        m_print(std::cout, "\nplot des runs en \t: ", _NBD_, " points");
        if(std::string("").compare(PATH_PLOT))
            m_print(std::cout, "\nrépertoire plots \t: ", PATH_PLOT);
        m_print(std::cout, "\nmode silencieux\t\t: ", (SILENT_MODE ? "oui" : "non"));
        m_print(std::cout, "\nmode intéractif\t\t: ", (INTERACTIVE ? "oui" : "non"), "\n\n", _CLR);

        float t(0), *U(nullptr);
        char *A(nullptr);
        int ins(0), run(0), div(0), m(-1), n(-1), *C(nullptr);
        std::vector<int> zInits(NUM_ITER, 0),
                         zAmels(NUM_ITER, 0),
                         zBests(NUM_ITER, 0);
        std::vector<float> tMoy;
        auto divs = matplot::transform(
            matplot::linspace(1, NUM_ITER, _NBD_),
            [](double x) {return (int)x;});
        if(NUM_ITER-1 <= 1) divs[0] = 1;
    #else
        float tt(0.f);;
        m_print(std::cout, _CLP, "\nRÉSOLUTION AVEC GLPK\n\n", _CLR);
    #endif

    std::vector<std::string> fnames = getfname(path);
    for(auto instance : fnames) {
        #if USE_GLPK
            modelSPP(instance, path, &tt, VERBOSE_GLPK);
        #else
            std::vector<int>    zMin(_NBD_, INT_MAX),
                                zMax(_NBD_, INT_MIN);
            std::vector<double> zMoy(_NBD_, 0);
            if(tMoy.size() == 0) {
                for(ins = 0; ins < (int)fnames.size(); ins++)
                    tMoy.push_back(0);
                ins = 0;
            }

            // Load one numerical instance
            std::tie(m, n, C, A, U) = loadSPP(path + instance);
            m_print(std::cout, _CLB, "\nInstance : ", instance, "\n", _CLR);

            m_print(std::cout, "Run exécutés :");
            for(run = 0; run < NUM_RUN; run++) {
                // Run ACO NUM_RUN times
                TIMED(t,
                    ACO(m, n, C, A, U, zInits, zAmels, zBests,
                        NUM_ITER, DEEPSEARCH, PARALLEL);
                );
                tMoy[ins] = (!run) ? t : tMoy[ins]+t;
                // Compute zMax, zMin and zMoy NUM_DIVISION time
                for(div = 0; div < _NBD_; div++) {
                    zMin[div] = std::min(zBests[divs[div]-1], zMin[div]);
                    zMax[div] = std::max(zBests[divs[div]-1], zMax[div]);
                    zMoy[div] += zBests[divs[div]-1];
                }
                m_print(std::cout, " ", run+1);
            }

            // Finish computing average z values
            for(div = 0; div < _NBD_; div++) zMoy[div] /= (double)NUM_RUN;

            // Plots
            m_print(std::cout, "\nPlot du dernier run...\n");
            plotRunACO(instance, zInits, zAmels, zBests, PATH_PLOT, SILENT_MODE);
            // TODO plot des niveaux de phéromones
            // m_print(std::cout, "Plot des probabilités des α pour le dernier run...\n");
            // plotProbaRunGRASP(instance, alpha, proba, PATH_PLOT, SILENT_MODE);
            m_print(std::cout, "Bilan de l'ensemble des runs...\n");
            plotAnalyseACO(instance, divs, zMin, zMoy, zMax, PATH_PLOT, SILENT_MODE);

            /* MOST IMPORTANT SECTIONS */
            freeSPP(C, A, U);
            ins++;
        #endif
    }

    #if USE_GLPK
        glp_free_env();
    #else
        // Finish computing average CPUt values
        for(ins = 0; ins < (int)fnames.size(); ins++)
            tMoy[ins] /= NUM_RUN;

        // Plots
        m_print(std::cout, "\n\nBilan CPUt moyen (par run) pour chaque instance...\n");
        plotCPUt(fnames, tMoy, PATH_PLOT, SILENT_MODE);

        if(INTERACTIVE) {
            m_print(std::cout, _CLG, "\nMODE INTÉRACTIF: Appuyez sur ENTRER pour terminer...\n", _CLR);
            std::cin.get();
        }
    #endif

    return 0;
}
