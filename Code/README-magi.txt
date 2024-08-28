https://www.azlyrics.com/lyrics/donnasummer/ontheradiolongversion.html

https://www.quechoisir.org/outil-speedtest-n64483/ 
____ le 26/1/24 12h32
Réception : 24,31 Mb/s
Envoi : 3,06 Mb/s
Ping moyen : 38 ms (Médiane : 25 ms)
____ le 2/2/24 11h46
Réception : 16,12 Mb/s
Envoi : 3,36 Mb/s
Ping moyen : 39 ms (Médiane : 29 ms)

https://www.youtube.com/watch?v=UEkZ44oj1Ws
https://www.radiofrance.fr/franceinter/podcasts/le-moment-meurice/les-fans-de-didier-raoult-2475015
https://www.radiofrance.fr/franceinter/podcasts/par-jupiter
________________________________________________________________
https://www.facebook.com/profile.php?id=100093615966605
https://deepspacecowboyz.bandcamp.com/
https://www.youtube.com/watch?v=Du3Ax9Rkyyw

https://www.2anes.com/alaffiche/a-partir-du-25-janvier-prochain/
https://dl.ipgp.fr/xciuw
NUM CARTE 6310354338292577
PIN 2085

#!/bin/bash
#SBATCH --job-name=1gE763
#SBATCH --output=1-gamma_E-7-6-3.out
#SBATCH --error=1-gamma_E-7-6-3.err
#SBATCH --partition=MISC-56c
#SBATCH --ntasks=1
#SBATCH --mem=8G
PROB="1-gamma_E-7-6-3"
module load cplex/20.1
cplex -c read ./lp/$PROB.lp optimize write ./sol/$PROB.sol 

#SBATCH --cpus-per-task=8
export OMP_NUM_THREADS=8
#       modified:   ../Code/README-magi.txt
#       modified:   ../Code/decompo-premier.c
#       modified:   ../Code/problemes.c
#       modified:   ../Code/rhogamma.c
#       modified:   ../Code/tableaux.h
#       modified:   designs-annexe.pdf
#       modified:   gammas.pdf
#       modified:   gammas.tex
#       modified:   point-Gamma+Gamma_E.txt
_____________________________________________________________________ mei#V4eM
https://www.youtube.com/watch?v=oK5spIIhhk0
wiki: 
	http://www-magi.univ-paris13.fr/doku.php
_________________________________________________
ci-dessous exemple de script à lancer avec 
	sbatch
autres commandes utiles:
	squeue
	scancel (numéro de job  donné par squeue)
remarques sur les parametres:
	tjrs 1 noeud pour Cplex
	on peut limiter le temps sur MISC-56c (infini par défaut) à
		3h: partition=MISC-56c-SHORT
		ou 3j: partition=MISC-56c-VERYSHORT
NB port: le préciser pour ssh (-p2822) mais scp (-P2822) aussi
TODO :: UTILISER SMP

#!/bin/bash
#SBATCH --job-name=1-gamma_E-8-2-2
#SBATCH --output=1-gamma_E-8-2-2.out
#SBATCH --error=1-gamma_E-8-2-2.err
#SBATCH --partition=MISC-56c
#SBATCH --nodes=1
module load cplex/20.1
srun cplex -c read 1-gamma_E-8-2-2.lp optimize

# job en attente ::
24547  MISC-56c    1g654 sophie.t PD       0:00      1 (Nodes required for job are DOWN, DRAINED or reserved for jobs 
in higher priority partitions)
squeue -j 24547
