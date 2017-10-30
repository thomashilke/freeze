# Introduction
Ce dossier contient le code 'freeze' de solidification du bain autour
d'une particule d'alumine, ainsi que les dependances necessaires a la
compilation. 


# Dependances externes
Les packages suivant doivent etre installe et accessibles:
- petsc version 3.7.6 (fichiers de developpement)
  freeze supporte probablement les version ulterieurs de petsc, mais
  ces cas de figures n'ont pas ete testes.
  La version 3.6.3 n'est pas supportee.
  La version 3.7.6 est supportee
  La version 3.7.7 est supportee
  Voir la section plus bas pour installer petsc!
- liblapacke-dev version 3.5.0 (fichiers de developpement)
- clang++ version 3.8.0
- g++ version 4.8.4


# Compilation 
Il faut commence par configurer la variable `PETSC_DIR` dans le fichier
`config.mk`. (voir la section sur l'installation de petsc plus bas) La
variable `PETSC_DIR` doit pointer vers le dossier d'installation de
petsc. Elle a la meme semantique que dans le package petsc lui-meme.

Il suffit ensuite de lancer gnu make:
```shell
  user@hostname freeze> make
```
L'executable et les fichier de configuration se trouves dans le
dossier `freeze/bin/`.

La compilation des differents packages installe des header (`*.hpp`) et
des librairies (`lib*.a`) dans les dossiers `~/.local/lib/` et
`~/.local/include/`. 


# Freeze: execution du code de simulation de formation de gel
Le dossier `freeze/bin/` contient une serie de fichiers de parametres
qui correspondent a differents problemes, dont la description est la
suivante:

- `aluminium-slab.conf`:  gel d'aluminium metallique sur une plaque
  infinie
  
- `cryolite-particle.conf`:  gel de bain sur une particule de cryolite

- `cryolite-particle-remelt.conf`:  gel de bain sur une particule de
  cryolite avec refonte
  
- `cryolite-particle-remelt-study.conf`:  gel de bain sur une particule
  de cryolite avec refonte, resultats consignes dans la these,
  
- `cryolite-slab.conf`:  gel de cryolite sur une plaque infinie

- `cryolite-slab-remelt.conf`:  gel de cryolite sur une plaque infinie
  avec refonte
  
- `double-growth.conf`:  convergence vers une solution a l'equilibre
  non triviale en coordonnees cartesiennes
  
- `double-spherical-growth.conf`:  convergence vers une solution a l'equilibre
  non triviale en coordonnees spheriques
  
- `stephan-exact-neumann-solution.conf`:  convergence L2 de la solution
  exacte de Neumann.

Les resultats consignes dans la these peuvent etre obtenus avec la
command suivante:
```shell
  user@hostname freeze/freeze/bin> ./freeze  cryolite-particle-remelt-study.conf
```
pour l'evolution du rayon de gel en fonction du temps, et 
```shell
  user@hostname freeze/freeze/bin> ./neumann stephan-exact-neumann-solution.conf
```
pour la convergence de l'erreur du schema de Chernoff.


# Description des packages
- freeze: Implementation du schema de Chernoff en une dimension,
  coordonnee cartesienne ou spherique,
  Depend de: `parameters`, `potential`
  
- parameter: Interpretation de fichier de parametres, et gestion de
  collections de parametres,
  Depend de: `spikes`
  
- alucell-db-tools: acces en lecture/ecriture aux fichiers dbfile
  genere par le code d'Alucell,
  Depend de:
  
- lexer: Extraction de lexemes a partir d'une source de caracteres et
  de la definition d'un ensemble de lexemes.
  Depend de:
  
- potential: Template Finite Element Library (tfel); Librairie
  elements finits N dimensions base essentiellement sur des
  techniques de metaprogrammation en C++, oriente vers un haut niveau
  d'abstraction, et une syntaxe aussi proche que possible de la
  formulation mathematique d'un probleme elements finits.
  Depend de: `alucell-db-tools`, `spikes`
  
- spikes: collection de single header libraries.
  Depend de:


# Installation de petsc

## Compilation
La methode standard est d'installer les distributions de petsc dans le
dossier `~/.local/opt/petsc/petsc-*.*.*/`, ou `*.*.*` est a remplacer
par la version, par exemple 3.7.7. A supposer que le tarball se
nomme `~/petsc-3.7.7.tar.gz`, on procede de la maniere suivant pour
l'installation (on suppose l'utilisation d'un interpreteur bash): 
```shell
  user@hostname ~> tar -xvf petsc-3.7.7.tar.gz; cd petsc-3.7.7/
  user@hostname ~> PETSC_DIR=$(HOME)/petsc-3.7.7 ./configure  --prefix=$(HOME)/.local/opt/petsc/petsc-3.7.7/
  user@hostname ~> PETSC_DIR=$(HOME)/petsc-3.7.7 make
  user@hostname ~> PETSC_DIR=$(HOME)/petsc-3.7.7 make install
```
On peut eventuellement exporter la variable `PETSC_DIR` dans le fichier
.bashrc: 
```shell
  export PETSC_DIR=/home/.../.local/opt/petsc/petsc-3.7.7/
```
pour que celle-ci soit definie par defaut dans chaque instance de
l'interpreteur bash.

Un script `install-petsc.sh` est fourni pour automatiser la procedure.


## Configuration du cache de librairies dynamiques
Comme le dossier d'installation des librairies dynamiques de petsc
n'est pas standard, il faut probablement configurer le cache de
librairies dynamiques de linux. Pour ce faire, on a deux solutions:

1. Ajouter le dossier /home/.../.local/opt/petsc/petsc-3.7.7/lib a la 
   variable d'environement LD_LIBRARY_PATH:
   ```shell
     user@hostname ~> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/.../.local/opt/petsc/petsc-3.7.7/lib
   ```
   Cette option a l'avantage de ne pas necessiter de droits d'admin.

2. Reconstruire le cache de ld. Il faut editer le fichier /etc/ld.so.conf,
   et ajouter la ligne suivante:
   ```shell
     /home/.../.local/opt/petsc/petsc-3.7.7/lib
   ```
   Puis reconstruire le cache en executant la commande suivante avec les
   droits d'admin:
   ```shell
     user@hostname ~# ldconfig 
   ```
   On peut verifier que la librairie dynamique figure bien dans le cache
   en executant:
   ```shell
     user@hostname ~# ldconfig -p
   ```
