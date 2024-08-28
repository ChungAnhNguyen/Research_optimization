/*	rho_pl
	rho_E(7, 3, 4) : 0.020576 (inverse 48.600000)	=> sol int : 5/243
*/

#include <stdio.h>	// pour traces
#include <stdlib.h>	// pour EXIT_SUCCESS
#include <math.h> 	// pour pow()
#include "vn.h" 		// representation des vecteurs et des sous-ensembles par des nombres
#include "rho.h" 	// (env, pl, q, p, k, E_q)

/*	Quelques jeux connus :
		rho_E(8, 2, 2) = 1/36, solution entiere donne 4/144
*/
int main() {
	/* __ Variables */
	int status = 0, nu =7, q =3, k =4, E_q =1;
	double val_cible =10.0/486;
	long int num_cible =1, den_cible =4, R_cible =12; 
	rho_pl g;
	resolution resol = RESOL_BIN;
	char nom_pl[100];

	#if (DEBUG)
	printf("%s IN\n", __func__);
	#endif

	/* __ PL */

	/* Initialisation */
	switch(resol) {
		case RESOL_CONT:
			status = init_cont(& g, q, nu, k, E_q);
			break;

		case RESOL_INT:
			status = init_int(& g, q, nu, k, E_q, val_cible);
			break;

		case RESOL_BIN:
			status = init_bin(& g, q, nu, k, E_q);
			break;

		case RESOL_INT_PV_OPT:
			status = init_int_pv_opt(& g, q, nu, k, E_q, num_cible, den_cible);
			break;

		case RESOL_BIN_PV_OPT:
			status = init_bin_pv_opt(& g, q, nu, k, E_q, R_cible);
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
			sprintf(nom_pl, "%d-rho%s-%d-%d-%d.lp", resol, (E_q == 1 ? "_E" : ""), nu, q, k);
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

