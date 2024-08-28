/*	gamma_pl

	TODO :
	- chaînes traitement main : 
		- cont -> int -> pv opt int ; 
		- bin -> int

gcc -c -o reduc.o -Wall reduc.c -I/Users/Soso/Applications/IBM/ILOG/CPLEX_Studio_Community127/cplex/include
gcc reduc.o -L/Applications/IBM/ILOG/CPLEX_Studio_Community127/cplex/lib/x86-64_osx/static_pic -lcplex -lm -lpthread

# /opt/lipn/software/cplex121/examples/x86_debian4.0_4.1
# /opt/lipn/software/cplex102/examples/x86_rhel4.0_3.4
# /opt/lipn/software/cplex101/examples/x86_rhel4.0_3.4
# /opt/lipn/software/cplex/examples/x86_debian4.0_4.1
# /opt/cplex90/examples/i86_linux2_glibc2.3_gcc3.2/

CPLEXDIR = /opt/lipn/software/cplex121
CFLAGS = -I$(CPLEXDIR)/include -m32 -fPIC
LDFLAGS = -L$(CPLEXDIR)/lib/x86_debian4.0_4.1/static_pic -lcplex -m32 -lm -pthread

CPLEXDIR = /opt/lipn/software/cplex
LDFLAGS = -L$(CPLEXDIR)/lib/x86_debian4.0_4.1/static_pic -lcplex -m32 -lm -pthread

CPLEXDIR = /opt/cplex90
CFLAGS = -I$(CPLEXDIR)/include -fPIC
LDFLAGS = -L$(CPLEXDIR)/lib/i86_linux2_glibc2.3_gcc3.2/static_pic -lcplex -lm -pthread

gcc -c -Wall --pedantic -I -I /opt/lipn/software/cplex102/include/ 

LOCAL_INCDIR= ~/usr_local/include/
CFLAGS= -Wall --pedantic -I $(LOCAL_INCDIR)
#CPLEXCFLAGS= -I $(LOCAL_INCDIR)
CPLEXCFLAGS= -I /opt/lipn/software/cplex102/include/

LDFLAGS= -lm -lgen -L $(LOCAL_LIBDIR)
LDFLAGS_GLPK= -lglpk
LDFLAGS_CPLEX= -L /opt/lipn/software/cplex102/lib/x86_rhel4.0_3.4/static_pic/ -lcplex -lpthread

-lm -L/opt/lipn/software/cplex102/lib/x86_rhel4.0_3.4/static_pic/ -lcplex -lpthread
-lm -L/opt/lipn/software/cplex121/lib/x86_debian4.0_4.1/static_pic/ -lcplex -lpthread
*/

#include <stdio.h>	// pour traces
#include <stdlib.h>	// pour EXIT_SUCCESS
#include <math.h> 	// pour pow()
#include "vn.h" 		// representation des vecteurs et des sous-ensembles par des nombres
#include "gamma.h" 	// (env, pl, q, p, k, E_q)

/*	Quelques jeux connus :

	gamma_E(5, 2, 2) :
		- opt_E(5, 2, 2) = 0.2
		- solution entiere de valeur 2/10 =0.2
		- optimum binaire de valeur 1/7

	gamma_E(6, 3, 3) :
		- opt_E(6, 3, 3) = 0.054787821063
		- solution entiere de valeur 153700/2805368 =38425/701342
		- optimum binaire ?
*/
int main() {
	/* __ Variables */
	int status = 0, q =8, p =2, k =2, E_q =1;
	double val_cible =1.0/8;
	long int num_cible =-1, den_cible =-1, R_cible =-1; 
	gamma_pl g;
	resolution resol = RESOL_INT;
	char nom_pl[100];

	#if (DEBUG)
	printf("%s IN\n", __func__);
	#endif

	/* __ PL */

	/* Initialisation */
	switch(resol) {
		case RESOL_CONT:
			status = init_cont(& g, q, p, k, E_q);
			break;

		case RESOL_INT:
			status = init_int(& g, q, p, k, E_q, val_cible);
			break;

		case RESOL_BIN:
			status = init_bin(& g, q, p, k, E_q);
			break;

		case RESOL_INT_PV_OPT:
			status = init_int_pv_opt(& g, q, p, k, E_q, num_cible, den_cible);
			break;

		case RESOL_BIN_PV_OPT:
			status = init_bin_pv_opt(& g, q, p, k, E_q, R_cible);
			break;

		default:
			printf("resolution (%d) non reconnue\n", resol);
			status =STAT_ERROR_PARAM;
	}

	/* Definition du PL (=> on l'ecrit) */
	if(status == 0) {
		status = setPL(& g);
		if (status)  {
			printf("Pb definition PL.\n");
		}
		else {
			sprintf(nom_pl, "%d-gamma%s-%d-%d-%d.lp", resol, (E_q == 1 ? "_E" : ""), q, p, k);
			status =CPXwriteprob (g.env, g.pl, nom_pl, "LP");
		}
	}

	/* __ Resolution du PL (=> on affiche la solution si une solution a ete obtenue) */
	if(status == 0) {
		status = resoudre(& g);
		if(status == 0) {
			affSol(& g);
		}
	}

	/* __ Fin */
	printf ("Status avant liberation memoire : %d\n", status);

	status = liberer(& g);
	printf ("Status apres liberation memoire : %d\n", status);

	return EXIT_SUCCESS;
}

