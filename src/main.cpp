#include "plots.hpp"
#include "heuristics.hpp"

#include <omp.h>

// Paramètres GLPK
#define USE_GLPK        false
#define VERBOSE_GLPK    false

// Paramètres OpenMP
#define PARALLEL        false
#define MAX_THREADS     10

// Paramètres ACO
#define NUM_RUN         1
#define MAX_ANTS        15
#define MAX_ITER        200
#define PHI_INIT        1.0
#define RHO_E           0.8
#define PHI_NUL         0.001
#define ITER_STAGNANT   8
#define DEEPSEARCH      false
#define RESTARTSTOP     true
#define MAX_RESTART     1

// Paramètres plots
#define CAPTURE_PHI     true
#define INTERACTIVE     false
#define SILENT_MODE     false
#define PATH_PLOT       ""
#define NUM_DIVISION    20

int main() {
    // This program will create different sequence of
    // random numbers on every program run
    // Use current time as seed for random generator
    std::string pwd(std::filesystem::current_path());
    std::string path(pwd + "/instances/");
    std::cout.precision(3);

    m_print(std::cout, _CLRd, "Etudiants : MERCIER et PICHON\n", _CLR);
    #if !USE_GLPK
        #if MAX_ITER < 2 // We need at least two iterations or else the plots
            #undef MAX_ITER //break
            #define MAX_ITER 2
        #endif
        #if NUM_RUN < 1
            #undef NUM_RUN
            #define NUM_RUN 1
        #endif
        #if MAX_THREADS < 1
            #undef MAX_THREADS
            #define MAX_THREADS 10
        #endif
        #if ITER_STAGNANT > MAX_ITER
            #undef ITER_STAGNANT
            #define ITER_STAGNANT int(MAX_ITER/10)
        #endif
        #if MAX_ANTS < 1
            #undef MAX_ANTS
            #define MAX_ANTS 1
        #endif

        INIT_TIMER(); srand(time(NULL));
        if(PARALLEL) omp_set_num_threads(MAX_THREADS);

        const int _NBD_ = NUM_DIVISION > MAX_ITER ? MAX_ITER : NUM_DIVISION;
        const float RHO_D = PHI_INIT * (1.0 - RHO_E);
        m_print(std::cout, _CLP, "\nnombre de runs\t\t: ", NUM_RUN);
        m_print(std::cout, "\nnombre de fourmis\t: ", MAX_ANTS);
        m_print(std::cout, "\nnombre d'itérations\t: ", MAX_ITER);
        m_print(std::cout, "\nvaleur initiale ϕ\t: ", PHI_INIT);
        m_print(std::cout, "\ntaux d'évaporation\t: ", RHO_E);
        m_print(std::cout, "\ntaux de dépôt\t\t: ", RHO_D);
        m_print(std::cout, "\nseuil phéromone nul\t: ", PHI_NUL);
        m_print(std::cout, "\nmax. cycles stagnants\t: ", ITER_STAGNANT);
        m_print(std::cout, "\nparallélisation\t\t: ", (PARALLEL ? "oui" : "non"));
        if(PARALLEL)
            m_print(std::cout, "\nnombre de threads\t: ", MAX_THREADS);
        m_print(std::cout, "\ndescente profonde\t: ", (DEEPSEARCH ? "oui" : "non"));
        m_print(std::cout, "\nplot des runs en \t: ", _NBD_, " points");
        if(std::string("").compare(PATH_PLOT))
            m_print(std::cout, "\nrépertoire plots \t: ", PATH_PLOT);
        m_print(std::cout, "\ncapturer ϕ\t\t: ", (CAPTURE_PHI ? "oui" : "non"));
        m_print(std::cout, "\nmode silencieux\t\t: ", (SILENT_MODE ? "oui" : "non"));
        m_print(std::cout, "\nmode intéractif\t\t: ", (INTERACTIVE ? "oui" : "non"), "\n\n", _CLR);

        char *A(nullptr);
        float t(0), *U(nullptr), *phi(nullptr), *phi_bef(nullptr), *phi_aft(nullptr);
        int ins(0), run(0), div(0), zinit(-1), zls(-1), zbest(-1),
            done_iter(0), m(-1), n(-1), *C(nullptr);
        // Vecteur de phéromones
        std::vector<int> zInits(MAX_ITER * MAX_ANTS, 0),
                         zAmels(MAX_ITER * MAX_ANTS, 0),
                         zBests(MAX_ITER * MAX_ANTS, 0);
        std::vector<float> probas(MAX_ITER * MAX_ANTS, 0.0);
        std::vector<float> tMoy;
        std::vector<double> divs;
    #else
        float tt(0.f);;
        m_print(std::cout, _CLP, "\nRÉSOLUTION AVEC GLPK\n\n", _CLR);
    #endif

    std::vector<std::string> fnames = getfname(path);
    for(auto instance : fnames) {
        #if USE_GLPK
            modelSPP(instance, path, &tt, VERBOSE_GLPK);
        #else
            float allrunzmoy(0);
            int allrunzmin(INT_MAX), allrunzmax(INT_MIN);
            std::vector<int>    zMin(_NBD_, INT_MAX),
                                zMax(_NBD_, INT_MIN);
            std::vector<double> zMoy(_NBD_, 0);
            if(tMoy.size() == 0) {
                for(ins = 0; ins < (int)fnames.size(); ins++)
                    tMoy.push_back(0);
                ins = 0;
            }

            // Load one numerical instance (also init phi)
            std::tie(m, n, C, A, U, phi) = loadSPP(path + instance, PHI_INIT);
            m_print(std::cout, _CLB, "\nInstance : ", instance, "\n", _CLR);

            m_print(std::cout, "Run exécutés :");
            for(run = 0; run < NUM_RUN; run++) {
                // Run ACO NUM_RUN times
                TIMED(t,
                    std::tie(zinit, zls, zbest, done_iter) = ACO(m, n, C, A, U,
                        zInits, zAmels, zBests, probas, phi, &phi_bef, &phi_aft,
                        PHI_NUL, RHO_E, RHO_D, ITER_STAGNANT, MAX_ANTS, MAX_ITER,
                        DEEPSEARCH, PARALLEL, RESTARTSTOP, MAX_RESTART, CAPTURE_PHI);
                );
                tMoy[ins] = (!run) ? t : tMoy[ins]+t;
                // Compute zMax, zMin and zMoy NUM_DIVISION time
                divs = matplot::transform(
                    matplot::linspace(1, done_iter/(float)MAX_ANTS, _NBD_),
                    [](double x) {return (int)x;});
                if(MAX_ITER-1 <= 1) divs[0] = 1;
                int idx = 0;
                for(div = 0; div < _NBD_; div++) {
                    idx = (divs[div]-1)*MAX_ANTS;
                    zMin[div] = std::min(zBests[idx], zMin[div]);
                    zMax[div] = std::max(zBests[idx], zMax[div]);
                    zMoy[div] += zBests[idx];
                }
                // Compute allrunzmin, allrunzmoy and allrunzmax
                allrunzmin = std::min(allrunzmin, zbest);
                allrunzmax = std::max(allrunzmax, zbest);
                allrunzmoy += zbest;

                m_print(std::cout, " ", run+1);
            }

            // Finish computing average z values
            allrunzmoy /= (double)NUM_RUN;
            for(div = 0; div < _NBD_; div++) zMoy[div] /= (double)NUM_RUN;

            // Plots
            m_print(std::cout, "\nPlot du dernier run...\n");
            plotRunACO(instance, zInits, zAmels, zBests, probas, zinit,
                        zls, zbest, done_iter, PATH_PLOT, SILENT_MODE);
            if(CAPTURE_PHI && phi_bef && phi_aft) {
                m_print(std::cout, "Plot des niveaux de phéromones avant le premier restart...\n");
                plotPhiRunACO(instance, n, phi_bef, true, PATH_PLOT, SILENT_MODE);
                m_print(std::cout, "Plot des niveaux de phéromones après le premier restart...\n");
                plotPhiRunACO(instance, n, phi_aft, false, PATH_PLOT, SILENT_MODE);
                if(phi_bef) delete[] phi_bef, phi_bef = nullptr;
                if(phi_aft) delete[] phi_aft, phi_aft = nullptr;
            }
            m_print(std::cout, "Bilan de l'ensemble des runs...\n");
            plotAnalyseACO(instance, divs, zMin, zMoy, zMax, allrunzmin,
                    allrunzmoy, allrunzmax, PATH_PLOT, SILENT_MODE);

            /* MOST IMPORTANT SECTIONS */
            freeSPP(C, A, U, phi);
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
