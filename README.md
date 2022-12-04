# acoSPP_cpp
Ant Colony Optimization (ACO) for the Set Packing Problem (SPP) in C++

*Commande pour compiler*:

    make


*Commande pour lancer le programme*:

    ./bin/DM3


*Commande pour effacer les fichiers générer lors de la compilation*:

    make clean


**ATTENTION:** le répertoire `instances` doit se trouver dans le même répertoire
d'où le programme est lancé. Si le programme est lancé depuis `bin` alors `instances`
doit être dans `bin`. Pour utiliser le makefile les répertoires `include`, `src` et
`lib` doivent être présents. Toutes les commandes présentées dans ce README sont
supposées être lancé depuis le répértoire parent de `src` (sauf indication contraire).

### Dépendances et standard C++ utilisé
Le standard C++ utilisé dans ce projet est le C++20, assurez-vous d'avoir une version
de g++ supportant ce standard.

Il est possible de résoudre les instances avec GLPK, assurez-vous d'avoir GLPK
installé sur votre machine (l'installation de la bibliothèque **libglpk-dev** est
fortement recommandée et généralement suffisante pour utiliser GLPK dans notre cas).

Pour paralléliser ACO, nous utilisons OpenMP, assurez-vous d'avoir OpenMP installé
sur votre machine (l'installation de **libomp-dev** est fortement recommandée et
généralement suffisante pour utiliser OpenMP dans notre cas).

Pour générer les plots, nous utilisons une version modifiée de la bibliothèque Matplot++
(https://alandefreitas.github.io/matplotplusplus/) qui a besoin des dépendances
suivantes :
* **gnuplot** (version 5.2.6+)
* **libx11-dev** (recommandé si vous utilisez le protocole de
  fenêtrage X11)

La bibliothèque Matplot++ utilise aussi d'autres dépendances mais la bibliothèque dynamique
_libmatplot.so_ (répertoire `lib`) comporte tous le nécessaire pour compiler Matplot++.
Cette bibliothèque fonctionne pour l'architecture **x86**. Un
script (`configure.sh`, utilisant la commande _cp_) a été créé pour
remplacer _libmatplot.so_ par une version supportant l'architecture **aarch** (ARM) si besoin.
Pour utiliser l'architecture x86 lancez depuis votre terminal :

    ./configure.sh x86

Pour utiliser l'architecture aarch :

    ./configure.sh aarch


### GLPK
Pour lancer GLPK sur les instances veuillez redéfinir USE_GLPK (constante préprocesseur dans le
fichier `src/main.cpp`) à true :
```c
#define USE_GLPK true
```
Sinon pour utiliser notre solution :
```c
#define USE_GLPK false
```

De plus, pour activer tous les messages de GLPK veuillez redéfinir VERBOSE_GLPK (constante préprocesseur dans
le fichier `src/main.cpp`) à true :
```c
#define VERBOSE_GLPK true
```
Sinon pour ne recevoir que les messages normaux :
```c
#define VERBOSE_GLPK false
```

### OpenMP
Pour activer la parallélisation veuillez redéfinir PARALLEL (constante
prépocesseur dans le fichier `src/main.cpp`) à true :
```c
#define PARALLEL true
```
Sinon pour exécuter les itérations séquentiellement :
```c
#define PARALLEL false
```

Pour préciser le nombre de threads maximum qu'OpenMP peut utiliser pour paralléliser ACO veuillez
redéfinir MAX_THREADS (constante préprocesseur dans le fichier `src/main.cpp`) :
```c
#define MAX_THREADS <x>
```
Où _x_ est un nombre entier strictement positif (toute valeur incorrecte sera remplacée par 10).


### Paramètres du ACO
#### Nombre de fourmies
Pour modifier le nombre de fourmies maximum veuillez redéfinir MAX_ANTS (constante
préprocesseur dans le fichier `src/main.cpp`) :
```c
#define MAX_ANTS <x>
```
Où _x_ est un nombre entier.

#### Valeur initiale des phéromones
Pour modifier la valeur initiale des phéromones veuillez redéfinir PHI_INIT (constante
préproceseseur dans le fichier `src/main.cpp`) :
```c
#define PHI_INIT <x>
```
Où _x_ est un nombre entre 0 et 1.

#### Seuil en deça duquel la valeur d'un phéromone est nulle
Pour modifier le seuil en deça duquel la valeur d'un phéromone est considérée
comme nulle veuillez redéfinir PHI_NUL (constante préprocesseur dans le fichier
`src/main.cpp`) :
```c
#define PHI_NUL <x>
```
Où _x_ est un nombre entre 0 et 1.

#### Coefficient d'évaporation
Pour modifier le coefficient d'évaporation des phéromones veuillez redéfinir
RHO_E (constante préprocesseur dans le fichier `src/main.cpp`) :
```c
#define RHO_E <x>
```
Où _x_ est un nombre entre 0 et 1.

#### Coefficient de dépôt
Le coefficient de dépôt des phéromones est déterminé selon la valeur de RHO_E
et selon la formule ```RHO_D = PHI_INIT * (1.0 - RHO_E)```.

#### Nombre d'itérations maximum
Pour modifier le nombre d'itérations maximum veuillez redéfinir MAX_ITER (
constante préprocesseur dans le fichier `src/main.cpp`) :
```c
#define MAX_ITER <x>
```
Où _x_ est un nombre entier strictement positif au moins égal à 2 (toute
valeur incorrecte sera remplacée par 2).

#### Nombre d'itérations avant de considéré qu'ACO stagne
Pour modifier la valeur du nombre d'itérations utilisée pour déterminer
si l'heuristique stagne, veuillez redéfinir ITER_STAGNANT (constante
préprocesseur dans le fichier `src/main.cpp`) :
```c
#define ITER_STAGNANT <x>
```
Où _x_ est un nombre entier strictement positif inférieur à MAX_ITER.

#### Condition d'arrêt dynamique
Pour activer la condition d'arrêt dynamique (arrêt après un certain nombre
de restart) veuillez redéfinir RESTARTSTOP (constante préprocesseur dans le
fichier `src/main.cpp`) :
```c
#define RESTARTSTOP true
```
Pour la désactiver :
```c
#define RESTARTSTOP false
```

#### Nombre de restart maximum
Lorsque la condition d'arrêt dynamique est activée, veuillez redéfinir MAX_RESTART
(constante préprocesseur dans le fichier `src/main.cpp`) pour modifier le nombre
de restart maximum :
```c
#define MAX_RESTART <x>
```
Où _x_ est un nombre entier.

#### Nombre d'exécutions
Pour modifier le nombre d'exécutions (run) de ACO veuillez
redéfinir NUM_RUN (constante préprocesseur dans le fichier
`src/main.cpp`) :
```c
#define NUM_RUN <x>
```
Où _x_ est un nombre entier strictement positif au moins égal à 1 (toute valeur
incorrecte sera remplacée par 1).

#### (Plus profonde) descente
Pour effectuer des améliorations par recherche locale de type plus profonde descente
veuillez redéfinir DEEPSEARCH (constante préprocesseur dans le
fichier `src/main.cpp`) à true :
```c
#define DEEPSEARCH true
```
Sinon pour des descentes "normales" :
```c
#define DEEPSEARCH false
```

### Paramètres plots
#### Nombre de points pour l'affichage des plots
Pour modifier le nombre de points utilisé pour afficher les runs de ACO
veuillez redéfinir NUM_DIVISION (constante préprocesseur dans le
fichier `src/main.cpp`) :
```c
#define NUM_DIVISION <x>
```
Où x est un nombre entier strictement positif entre 1 et NUM_ITER (toute valeur
incorrecte sera remplacée par NUM_ITER).

#### Capturer les niveaux de phéromones avant et après perturbation
Pour capturer les  niveaux de phéromones avant et après le premier restart (si
une perturbation est effectuée), veuillez redéfinir CAPTURE_PHI (constante
préprocesseur dans le fichier `src/main.cpp`) :
```c
#define CAPTURE_PHI true
```
Pour désactiver cette option :
```c
#define CAPTURE_PHI false
```

#### Désactivation des plots
L'affichage des plots peut être désactivé. Veuillez redéfinir SILENT_MODE
(constante préprocesseur dans le fichier `src/main.cpp`) :
```c
#define SILENT_MODE true
```
Il est recommandé d'enregistrer les plots dans un répertoire pour garder
une trace de l'exécution lorsque ce mode est activé (voir ci-dessous). Sinon,
pour afficher les plots :
```c
#define SILENT_MODE false
```

#### Mode intéractif
Les plots affichés sont tous intéractifs. L'actionnement de la molette de la
souris permet de zoomer sur les plots. Cependant, une fois que le programme
est terminé les plots deviennent statiques et il est impossible d'intéragir
avec ces-derniers. Un mode intéractif a été implémenté et permet d'intéragir avec les plots
tant que l'utilisateur n'appuie pas sur une touche depuis le terminal. Pour
activer le mode intéractif, veuillez redéfinir INTERACTIVE (constante
préprocesseur dans le fichier `src/main.cpp`) :
```c
#define INTERACTIVE true
```
Pour le désactiver :
```c
#define INTERACTIVE false
```

#### Enregistrement des plots
Pour enregistrer les plots dans un répertoire veuillez redéfinir PATH_PLOT (constante
préprocesseur dans le fichier `src/main.cpp`) :
```c
#define PATH_PLOT <path>
```
Où _path_ est une chaîne de caractères spécifiant le chemin vers le répertoire (Exemple :
"exp/" ou "."). Pour ne pas enregistrer les plots :
```c
#define PATH_PLOT ""
```


### Listes des instances qui ont été utilisées pour l'expérimentation
* didactic.dat
* pb_100rnd0100.dat
* pb_200rnd0100.dat
* pb_200rnd0300.dat
* pb_200rnd0500.dat
* pb_200rnd0900.dat
* pb_200rnd1500.dat
* pb_500rnd0100.dat
* pb_500rnd0700.dat
* pb_500rnd0900.dat
* pb_500rnd1100.dat
* pb_1000rnd0100.dat
* pb_1000rnd0500.dat
* pb_1000rnd0700.dat
* pb_2000rnd0700.dat
