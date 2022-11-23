#include "plots.hpp"
#include "matplot/util/common.h"

void plotRunACO(
        const std::string instance,
        const std::vector<int>& zInits,
        const std::vector<int>& zAmels,
        const std::vector<int>& zBests,
        std::string save_path,
        bool silent_mode) {
    int i(0), n = zInits.size(), ins_i(-1);
    auto X = matplot::linspace(1, n, n);
    std::string ins(instance);
    for(i = 0; i < (int)ins.size(); i++) {
        if(ins[i] == '_')
            ins_i = i, ins.replace(i, 1,  "\\\\\\_"), i+=4;
        if(ins[i] == '.')
            ins = ins.substr(0, i);
    }
    std::string tit("GRASP-SPP | z_{Init} z_{LS} z_{Best} | " + ins);

    double lb = *std::min_element(std::begin(zInits), std::end(zInits)),
           ub = *std::max_element(std::begin(zBests), std::end(zBests));

    if(lb == ub) {
        lb -= 5;
        ub += 5;
    }

    auto fig = matplot::figure(true);
    fig->name("Examen d'un run");
    fig->size(576, 576);
    fig->title(tit);
    matplot::xlabel("Itérations");
    matplot::ylabel("valeurs de z(x)");
    matplot::xticks({1.0, ceil(n/4.0), ceil(n/2.0), ceil((3*n)/4.0), (double)n});
    matplot::axis({0, n+1.0, lb-(int(lb/100)+1)*2, ub+(int(ub/100)+1)*2});
    matplot::plot(X, zBests)
        ->line_width(2)
        .line_width(2)
        .color("green")
        .display_name("meilleures solutions");
    matplot::hold(true); // Allow multiple plot() calls
    for(i = 1; i <= n; i++) {
        matplot::line(i, zInits[i-1], i, zAmels[i-1])
            ->line_width(0.5)
            .color("blue")
            .display_name(""); // Don't show in legend
    }
    matplot::plot(X, zAmels, "o")
        ->marker_size(4)
        .color("green")
        .marker("^")
        .marker_face(true) // Filled
        .display_name("toutes solutions améliorées");
    matplot::plot(X, zInits, "o")
        ->marker_size(2)
        .color("red")
        .marker(".")
        .display_name("toutes solutions construites");
    matplot::legend()
        ->location(matplot::legend::general_alignment::bottomright);
    if(!silent_mode) fig->draw();
    if(ins_i != -1) ins.replace(ins_i, 4, "_");
    if(save_path.compare(""))
        matplot::save(save_path + "run_" + ins + ".png");
}

// TODO remplacer par plot des niveaux de phéromones
// void plotProbaRunACO(
//         const std::string instance,
//         const std::vector<double>& alpha,
//         const std::vector<double>& proba,
//         std::string save_path,
//         bool silent_mode) {
//     int i(0), ins_i(-1);
//     std::string ins(instance);
//     for(i = 0; i < (int)ins.size(); i++) {
//         if(ins[i] == '_')
//             ins_i = i, ins.replace(i, 1,  "\\\\\\_"), i+=4;
//         if(ins[i] == '.')
//             ins = ins.substr(0, i);
//     }
//     std::string tit("GRASP-SPP | proba_{α} | " + ins);
//
//     auto fig = matplot::figure(true);
//     fig->name("Probabilités p des α pour un run");
//     fig->size(576, 576);
//     fig->title(tit);
//     matplot::xlabel("α");
//     matplot::ylabel("P(α)");
//     matplot::ylim({0, 1});
//     matplot::xticks(alpha);
//     matplot::bar(alpha, proba);
//     if(!silent_mode) fig->draw();
//     if(ins_i != -1) ins.replace(ins_i, 4, "_");
//     if(save_path.compare(""))
//         matplot::save(save_path + "proba_" + ins + ".png");
// }

void plotAnalyseACO(
        const std::string instance,
        const std::vector<double>& divs,
        const std::vector<int>& zMin,
        const std::vector<double>& zMoy,
        const std::vector<int>& zMax,
        std::string save_path,
        bool silent_mode) {
    int n = divs.size(), ins_i(-1);
    std::string ins(instance);
    for(int i = 0; i < (int)ins.size(); i++) {
        if(ins[i] == '_')
            ins_i = i, ins.replace(i, 1,  "\\\\\\_"), i+=4;
        if(ins[i] == '.')
            ins = ins.substr(0, i);
    }
    std::string tit("GRASP-SPP | z_{min} z_{moy} z_{max} | " + ins);
    auto yerr1 = matplot::transform(matplot::linspace(0, n-1, n),
            [zMin, zMoy](double x) {
                return zMoy[int(x)]-zMin[(int)x];
            }),
         yerr2 = matplot::transform(matplot::linspace(0, n-1, n),
            [zMoy, zMax](double x) {
                return zMax[(int)x]-zMoy[(int)x];
            }),
         xerr = matplot::linspace(0, 0, n);

    double ub = *std::max_element(zMax.begin(), zMax.end())
                + (*std::max_element(yerr2.begin(), yerr2.end()))/2;

    auto fig = matplot::figure(true);
    fig->name("Bilan tous runs");
    fig->size(576, 576);
    fig->title(tit);
    matplot::xlabel("Itérations");
    matplot::ylabel("valeurs de z(x)");
    matplot::axis({divs[0]-1, divs[n-1]+1, 0, ub+(int(ub/100)+1)*2});
    matplot::xticks(divs);
    matplot::errorbar(divs, zMoy, yerr1, yerr2, xerr, xerr)
        ->line_width(1)
        .color("black")
        .marker("+")
        .display_name("zMin zMax");
    matplot::hold(true); // Allow multiple plot() calls
    matplot::plot(divs, zMoy, "-")
        ->marker_size(4)
        .color("green")
        .marker("o")
        .marker_face(true)
        .display_name("zMoy");
    matplot::legend()
        ->location(matplot::legend::general_alignment::bottomright);
    if(!silent_mode) fig->draw();
    if(ins_i != -1) ins.replace(ins_i, 4, "_");
    if(save_path.compare(""))
        matplot::save(save_path + "analyse_" + ins + ".png");
}

void plotCPUt(
        std::vector<std::string>& fnames,
        std::vector<float>& tMoy,
        std::string save_path,
        bool silent_mode) {
    int n;
    for(n = 0; n < (int)fnames.size(); n++) {
        for(int i = 0; i < (int)fnames[n].size(); i++) {
            if(fnames[n][i] == '_')
                fnames[n].replace(i, 1,  "\\\\\\_"), i+=4;
            if(fnames[n][i] == '.')
                fnames[n] = fnames[n].substr(0, i);
        }
    }
    bool sp(n == 1); // Singlepoint ? (little trick if there is only one
                     // instance)
    auto x = matplot::linspace(1, sp ? n+2 : n, sp ? n+2 : n);
    if(sp) {
        tMoy.push_back(tMoy[0]), tMoy.push_back(0), tMoy[0] = 0;
        fnames.push_back(fnames[0]), fnames.push_back(""), fnames[0] = "";
    }

    auto fig = matplot::figure(true);
    fig->name("Bilan CPUt tous runs");
    fig->size(576, 676);
    fig->title("GRASP-SPP | tMoy");
    matplot::ylabel("CPUt moyen (en s)");
    matplot::xticks(x);
    matplot::xticklabels(fnames);
    matplot::xtickangle(60);
    matplot::plot(x, tMoy, sp ? "o" : "--")
    ->line_width(0.5)
    .marker_size(4)
    .color("blue")
    .marker("o")
    .marker_face(true)
    .display_name("tMoy");
    matplot::legend()
        ->location(matplot::legend::general_alignment::bottomright);
    if(!silent_mode) fig->draw();
    if(sp) {
        tMoy[0] = tMoy[1], tMoy.pop_back(), tMoy.pop_back();
        fnames[0] = fnames[1], fnames.pop_back(), fnames.pop_back();
    }
    if(save_path.compare(""))
        matplot::save(save_path + "CPUt.png");
}
