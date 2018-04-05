# Schéma de Ninomiya-Victoir pour le CIR et application au modèle de Heston

Etudier le schéma d’ordre faible élevé introduit dans [1] et étendu dans [2] pour les modèles affines.
* Appliquer cette méthode pour le pricing d’options asiatiques dans le modèle de Heston. Comparer
avec différentes méthodes de réduction de variance.
* Comment mettre en œuvre cette méthode dans le cadre Quasi-Monte Carlo ? Tester numériquement.


[1] S. Ninomiya and N. Victoir. Weak approximation of stochastic differential equations and application to derivative pricing. Appl. Math. Finance, 15(1-2) :107–121, 2008.
[2] A. Alfonsi. High order discretization schemes for the cir process : Application to affine term structure and heston models. Math. Comput., 79(269) :209–237, 2010.

https://github.com/Gregory-Calvez/proba_num

## Prérequis

Ce projet utilise principalement la STL du C++11.
Certains éléments doivent être installés pour reproduire tous les graphes et résultats du rapport.
* GNUplot. Les graphes sont générés à partir de GNUplot. http://www.gnuplot.info/
* GSL. GNU Scientific Library. https://www.gnu.org/software/gsl/
* NVIDIA Toolkit. https://developer.nvidia.com/cuda-downloads?target_os=Linux&target_arch=x86_64&target_distro=Ubuntu&target_version=1604&target_type=runfilelocal. Les calculs sur GPU du rapport ont été faits avec la version 8.0 du Toolkit NVIIDA pour des raisons de compatibilités avec le GPU. Tout devrait fonctionner avec le Toolkit le plus récent.

## Structure du projet et compilation

### Structure du projet

Le projet se sépare en deux parties distinctes:
* une partie reposant sur la programmation orientée objet en C++. Les sources sont disponibles dans le dossier src/
* une partie reposant sur du calcul parallélisé en CUDA/C. Tout est regroupé dans le dossier cuda/

Un makefile, rédigé à la main, comprend toutes les commandes nécessaires pour reproduire les résultats du rapport.

Des fichiers annexes sont présents, pour les tracés des graphes notamment. Enfin, nous avons conservé dans le dossier du projet des outputs de la console (output/) et des graphes générés (img/) que nous avons réutilisés dans le rapport.

Une documentation de la partie C++ réalisée avec Doxygen est disponible dans le dossier doc/.

### Compilation

Dans le fichier src/main.cpp, nous avons écrit toutes les fonctions permettant de tester la partie C++ du projet.
Des lignes de la fonction main sont commentées, il suffit de dé-commenter la partie que vous voulez tester et de compiler et lancer l'exécutable. Pour ce faire :
'''
make -B run
'''

Si un graphe doit être généré pendant l'exécution et que vous voulez le récupérer, utilsez la commande :
'''
make -B plot
'''

La partie CUDA peut être compilée et lancée à partir du même makefile.
Par exemple, pour reproduire les derniers graphes du rapport, utilisez la commande :
'''
make -B cuda
'''

Attention, l'exécution de certaines partie du code peut être très longue (20 min pour les derniers graphes du rapport en CUDA, ou certains Monte Carlo).


### Test sur la partie C++

Nous décrivons ici les tests que vous pouvez effectuer pour tester notre code. Touts les fonctions sont déjà écrites, il suffit de décommenter certaines partie du main.

* Plotting some examples of trajectories. Vous avez la possibilités de tracer une dizaine de trajectoires simulées par les schémas de discrétisation de Ninomiya-Victoir et de Glasserman. Il faut plotter les graphes. Chaque ligne doit être décommentée séparément.
* Testing the control variates. Si vous décommentez les 4 lignes qui suivent ce commentaire, vous obtiendrez en console les estimations et les intervalles de confiance du prix d'une option avec et sans réduction de variance, pour un schéma d'ordre faible et le schéma de Glasserman, pour des options asiatiques et des options européennes.
* Test for cuda. Cette fonction et cet appel nous a permis de comparer les temps de calcul entre le CUDA et le C++. Il est probablement inutile de le relancer.
* Test for Sobol. Cet appel permet d'optenir un QMC pour estimer le prix d'une option. Les résultats numériques ne sont pas satisfaisants. Cf rapport.
* Test for Alfonsi's graph. Cette portion du main permet de reproduire les graphes p.26-27 de l'article de M. Alfonsi. Les temps de calculs sont absolument prohibitif mais vous pouvez le tester avec un cap sur le nombre d'itération faible. Il y a ici un graph à tracer.

### Test sur la partie CUDA/C

Ici on décrit les tests présents dans le main du fichier cuda/cuda_kernels.cu.

Pour les lancers, le toolkit NVIDIA est nécessaire, la commande du makefile est :
'''
make -B cuda
'''
Cette fois, les graphes sont générées par cette seule commande.


* Test for Alfonsi's graph. Idem que dans la partie C++. Les temps de calculs deviennent raisonnables.
* Example Asian option. Un calcul MC d'une option asiatique en calcul parallélisés. Les résultats sont donnés en console. Les graphes ne représentent rien.

Une remarque sur le nombre de simulations : trois variables globales, définies au début du fichier cuda_kernels.cu permettent de définir le nombre de blocks, le nombre de threads par block et le nombre de simulations par thread pour les calculs de Monte Carlo. Nous les avons choisi pour accélérer le temps de calcul sur un GPU donné, ces paramètres peuvent être peut-être modifiés

## Auteurs
Elie Bohbot
Grégory Calvez
