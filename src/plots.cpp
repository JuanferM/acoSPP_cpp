#include "plots.hpp"
#include "matplot/freestanding/axes_functions.h"
#include "matplot/freestanding/axes_lim.h"
#include "matplot/util/common.h"
#include <string>

void plotRunACO(
        const std::string instance,
        const std::vector<int>& ozInits,
        const std::vector<int>& ozAmels,
        const std::vector<int>& ozBests,
        const std::vector<float>& oprobas,
        const int zinit,
        const int zls,
        const int zbest,
        const int done_iter,
        std::string save_path,
        bool silent_mode) {
    int i(0), n = done_iter, ins_i(-1);
    auto X = matplot::linspace(1, n, n);
    std::vector<int> zInits = std::vector<int>(ozInits.begin(),
                                ozInits.begin()+n),
                     zAmels = std::vector<int>(ozAmels.begin(),
                                ozAmels.begin()+n),
                     zBests = std::vector<int>(ozBests.begin(),
                                ozBests.begin()+n);
    std::vector<float> probas = std::vector<float>(oprobas.begin(),
                                oprobas.begin()+n);
    std::string ins(instance),
                sp(" | "),
                gr("Greedy : " + std::to_string(zinit)),
                ls("LS : " + std::to_string(zls)),
                aco("ACO : " + std::to_string(zbest));
    for(i = 0; i < (int)ins.size(); i++) {
        if(ins[i] == '_')
            ins_i = i, ins.replace(i, 1,  "\\\\\\_"), i+=4;
        if(ins[i] == '.')
            ins = ins.substr(0, i);
    }
    std::string tit("SPP : " + ins + sp + gr + sp + ls + sp + aco);

    double ub = *std::max_element(std::begin(zBests), std::end(zBests));

    auto fig = matplot::figure(true);
    fig->name("Examen d'un run");
    fig->title(tit);
    fig->size(576, 476);
    fig->title_font_size_multiplier(1);
    matplot::xlabel("# itérations x # fourmis");
    matplot::ylabel("valeurs de z(x)");
    matplot::y2label("Probabilité de sélection curieuse/fonceuse");
    matplot::xticks({1.0, ceil(n/4.0), ceil(n/2.0), ceil((3*n)/4.0), (double)n});
    matplot::axis({0, n+1.0, 0, ub+(int(ub/100)+1)*2});
    matplot::y2lim({0, 1});
    matplot::plot(X, zBests)
        ->line_width(2)
        .line_width(2)
        .color("red")
        .display_name("meilleures solutions");
    matplot::hold(true); // Allow multiple plot() calls
    for(i = 1; i <= n; i++) {
        matplot::line(i, zInits[i-1], i, zAmels[i-1])
            ->line_width(0.5)
            .color("black")
            .line_style("--")
            .display_name(""); // Don't show in legend
    }
    matplot::plot(X, std::vector<float>(done_iter, 0), "--")
        ->marker_size(0)
        .line_width(0.5)
        .color("black")
        .display_name("toutes solutions");
    matplot::plot(X, probas, "-")
        ->use_y2(true)
        .line_width(1)
        .color("blue")
        .display_name("Proba. exploitation");
    matplot::legend()
        ->location(matplot::legend::general_alignment::bottomright);
    matplot::gca()->y_axis().color("black");
    matplot::gca()->y2_axis().color("black");
    if(!silent_mode) fig->draw();
    if(ins_i != -1) ins.replace(ins_i, 4, "_");
    if(save_path.compare(""))
        matplot::save(save_path + "run_" + ins + ".png");
}

void plotPhiRunACO(
        const std::string instance,
        const int n,
        const float* phi,
        bool before,
        std::string save_path,
        bool silent_mode) {
    int i(0), ins_i(-1);
    std::string ins(instance),
                st((before) ? "avant" : "après");
    for(i = 0; i < (int)ins.size(); i++) {
        if(ins[i] == '_')
            ins_i = i, ins.replace(i, 1,  "\\\\\\_"), i+=4;
        if(ins[i] == '.')
            ins = ins.substr(0, i);
    }
    std::string tit("ACO-SPP : " + ins + " | phi(x)");

    auto fig = matplot::figure(true);
    fig->name("Etat du vecteur de phéromone " +st+ " application de la procédure de perturbation");
    fig->size(576, 476);
    fig->title(tit);
    fig->title_font_size_multiplier(1);
    matplot::xlabel("x");
    matplot::ylabel("phi(x)");
    matplot::ylim({0, 1});
    for(i = 1; i <= n; i++) {
        matplot::line(i, 0, i, phi[i-1])
            ->line_width(0.5)
            .color("black")
            .line_style("-")
            .display_name(""); // Don't show in legend
    }
    if(!silent_mode) fig->draw();
    if(ins_i != -1) ins.replace(ins_i, 4, "_");
    if(save_path.compare(""))
        matplot::save(save_path + "phi_" +st+"_"+ ins + ".png");
}

void plotAnalyseACO(
        const std::string instance,
        const std::vector<double>& divs,
        const std::vector<int>& zMin,
        const std::vector<double>& zMoy,
        const std::vector<int>& zMax,
        const int allrunzmin,
        const float allrunzmoy,
        const int allrunzmax,
        std::string save_path,
        bool silent_mode) {
    int n = divs.size(), ins_i(-1);
    std::ostringstream mo; mo.precision(2);
    mo << std::fixed << "z_{moy} : " << allrunzmoy;
    std::string ins(instance),
                sp(" | "),
                mi("z_{min} : " + std::to_string(allrunzmin)),
                ma("z_{max} : " + std::to_string(allrunzmax));
    for(int i = 0; i < (int)ins.size(); i++) {
        if(ins[i] == '_')
            ins_i = i, ins.replace(i, 1,  "\\\\\\_"), i+=4;
        if(ins[i] == '.')
            ins = ins.substr(0, i);
    }
    std::string tit("ACO-SPP : " + ins + sp + mi + sp + mo.str() + sp + ma);
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
    fig->title(tit);
    fig->size(576, 476);
    fig->title_font_size_multiplier(1);
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
        .color("red")
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
    int n = 0;
    std::vector<std::string> tMoytxt = std::vector<std::string>(fnames.size());
    for(n = 0; n < (int)fnames.size(); n++) {
        std::ostringstream v; v.precision(3);
        v << std::fixed << tMoy[n];
        tMoytxt[n] = v.str();
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
        tMoytxt.push_back(""); tMoytxt.push_back("");
        std::string tmp = tMoytxt[0];
        tMoytxt[0] = "", tMoytxt[1] = tmp;
    }

    auto fig = matplot::figure(true);
    fig->name("Bilan CPUt tous runs");
    fig->size(576, 676);
    fig->title("ACO-SPP | tMoy");
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
    matplot::text(x, tMoy, tMoytxt)->color("red");
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
