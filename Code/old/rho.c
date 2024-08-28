/*	rho_pl
*/

#include "rho.h" 	// (q, nu, k, E_q)
#include "vn.h" 		// representation des vecteurs et des sous-ensembles par des nombres

#if (DEBUG)
#include <stdio.h>	// pour traces
#endif

#include <stdlib.h>	// pour malloc()
#include <string.h>	// pour sprintf()
#include <math.h> 	// pour pow()

/* __________________________________________________________________________ Declaration fonctions statiques
*/

/* Sous-routine universelle : intervalle de travail pour v : {0 ,.., q^q -1} si E_q = 0, {0 ,.., q^(q -1) -1} sinon
*/
static int get_v_max(rho_pl* g);

/* Sous-routine universelle : intervalle de travail pour J : {0 ,.., 2^nu -1}
*/
static int get_J_max(rho_pl* g);

/* Sous-routine universelle : intervalle de travail pour a rhs contrainte bkI : {0 ,.., q^k -1} si E_q = 0, {0 ,.., q^(k -1) -1} sinon
*/
static int get_a_max(rho_pl* g);

/* Sous-routine des fonctions d'initialisation :
	Pre-condition : 		champs resol, val_cible, num_cible, den_cible, R_cible tous initialises
	Post-condition : 	champs atomiques q, nu, k, E_q tous initialises
	Valeur retournee : 	0 si ok, -1 si pb parametres, -2 si pb malloc, code cplex > 0 sinon
	Post-condition en cas de succes :
		- environnement Cplex +PL crees
		- tableau ctr_bkI[][] construit et initialise
	Post-condition en cas d'echec : memoire eventuellement allouee en cours de fonction totalement liberee
*/
static int init(rho_pl* g, int nu, int q, int p, int k, int E_q);

/* Sous-routine des fonctions d'initialisation : verification de la coherence des donnees
	valeur retournee : 0 si tout ok et -1 sinon
*/
static int verification(rho_pl* g);

/* Sous-routine de la fonction de definition du PL
	Contraintes :
	- R = sum_v Q(v)	 								si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	- Q mesure										si RESOL_CONT : 
	- sum_v Q(v) = R_cible	 						si RESOL_BIN_PV_OPT
	- Q(0, 0 ,.., 0) >= (ou =) 1 					si RESOL_INT ou RESOL_BIN ou RESOL_INT_PV_OPT 
	- Q(0, 0 ,.., 0)/R >= val_cible					si RESOL_INT
	- Q(0, 0 ,.., 0)/R >= num_cible/den_cible		si RESOL_INT_PV_OPT
	- Q (ou Q_E) in I_q^k							toujours
	Pre-conditions : g intialisee, parametres coherents
	Post-condition : contraintes crees (sans les coefficients des variables), tableau ctr_bkI rempli
*/
static int rho_ctr(rho_pl* g);

/* Sous-routine de la fonction de definition du PL
	Variables:
		- Q(v), v in Z_nu^q	: toujours, a valeur dans [0, 1], N, {0, 1} selon g->resol
		- R					: toujours, a valeur dans [0, 1] ou N selon g->resol
	Pre-conditions : g intialisee, contraintes toutes generees (sans les coefficients des variables), tableau ctr_bkI rempli, parametres coherents
	Post-condition : PL totalement defini (sinon type et sens optimisation)
*/
static int rho_var(rho_pl* g);

/* Sous-routine de la fonction de definition du PL
	Sens optimisation :
		- Max 		si RESOL_CONT ou RESOL_BIN_PV_OPT
		- Min		si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	Type probleme :
		- LP			si RESOL_CONT
		- MILPL		dans tous les autres cas
	Pre-conditions : PL totalement defini, sauf type et sens optimisation
	Post-condition : type et sens optimisation definis
*/
static int finaliser(rho_pl* g);

/* Sous-routine de la fonction de resolution : affichage synthetique du resultat de la resolution
	Pre-condition : le PL a ete resolu
	NB statuts resolution :
		https://www.tu-chemnitz.de/mathematik/discrete/manuals/cplex/doc/refman/html/appendixB.html
		1 :	CPX_STAT_OPTIMAL
		2 :	CPX_STAT_UNBOUNDED
		3 :	CPX_STAT_INFEASIBLE 
		101 : CPXMIP_OPTIMAL 
		103 : CPXMIP_INFEASIBLE 
		118 : CPXMIP_UNBOUNDED 
*/
static int restitution(rho_pl* g);

/* __________________________________________________________________________ Initialisation / destruction PL
*/

/* ____ Fonctions d'initialisation :
	Post-condition : 	champs atomiques nu, q, k, E_q, resol, val_cible, num_cible, den_cible, R_cible tous initialises
	Valeur retournee : 	0 si ok, -1 si pb parametres, -2 si pb malloc, code cplex > 0 sinon
	Post-condition en cas de succes :
		- environnement Cplex +PL crees
		- tableau ctr_bkI[][] construit et initialise
	Post-condition en cas d'echec : memoire eventuellement allouee en cours de fonction totalement liberee
*/

/* Initialisation resolution continue */
int init_cont(rho_pl* g, int q, int nu, int k, int E_q) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d) IN\n", __func__, g, q, nu, k, E_q);
	#endif

	g->type = TAB_GAMMA;
	g->resol = RESOL_CONT;
	g->val_cible = NO_VAL;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = NO_VAL;

	init(g, nu, q, NO_VAL, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation resolution entiere */
int init_int(rho_pl* g, int q, int nu, int k, int E_q, double val_cible) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %lf) IN\n", __func__, g, q, nu, k, E_q, val_cible);
	#endif

	g->type = TAB_GAMMA;
	g->resol = RESOL_INT;
	g->val_cible = val_cible;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = NO_VAL;

	init(g, nu, q, NO_VAL, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation PV optimalite solution entiere */
int init_int_pv_opt(rho_pl* g, int q, int nu, int k, int E_q, long int num_cible, long int den_cible) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %ld, %ld) IN\n", __func__, g, q, nu, k, E_q, num_cible, den_cible);
	#endif

	g->type = TAB_GAMMA;
	g->resol = RESOL_INT_PV_OPT;
	g->val_cible = NO_VAL;
	g->num_cible = num_cible;
	g->den_cible = den_cible;
	g->R_cible = NO_VAL;

	init(g, nu, q, NO_VAL, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation resolution binaire */
int init_bin(rho_pl* g, int q, int nu, int k, int E_q) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d) IN\n", __func__, g, q, nu, k, E_q);
	#endif

	g->type = TAB_GAMMA;
	g->resol = RESOL_BIN;
	g->val_cible = NO_VAL;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = NO_VAL;

	init(g, nu, q, NO_VAL, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation resolution binaire PV ratio 1/R_cible est opt s.c. R =R_cible */
int init_bin_pv_opt(rho_pl* g, int q, int nu, int k, int E_q, long int R_cible) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %ld) IN\n", __func__, g, q, nu, k, E_q, R_cible);
	#endif

	g->type = TAB_GAMMA;
	g->resol = RESOL_BIN_PV_OPT;
	g->val_cible = NO_VAL;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = R_cible;

	init(g, nu, q, NO_VAL, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Sous-routine des fonctions d'initialisation :
	Pre-condition : 		champs resol, val_cible, num_cible, den_cible, R_cible tous initialises
	Post-condition : 	champs atomiques q, nu, k, E_q tous initialises
	Valeur retournee : 	0 si ok, -1 si pb parametres, -2 si pb malloc, code cplex > 0 sinon
	Post-condition en cas de succes :
		- environnement Cplex +PL crees
		- tableau ctr_bkI[][] construit et initialise
	Post-condition en cas d'echec : memoire eventuellement allouee en cours de fonction totalement liberee
*/
int init(rho_pl* g, int nu, int q, int p, int k, int E_q) {
	int status = 0, J, a, J_max, a_max;
	char nom_pb[100];
	char errmsg[1024];

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %d) IN\n", __func__, g, nu, q, p, k, E_q);
	#endif

	/* ____ Parametres probleme & intialisations */
	g->nu =nu;
	g->q =q;
	g->p =p;
	g->k =k;
	g->E_q =E_q;

	g->ctr_bkI =NULL;
	g->env =NULL;
	g->pl =NULL;

	status =verification(g);

	/* ____ Cplex */
	if (status == 0) {
		/* environnement */
		g->env = CPXopenCPLEX (& status);
		if (g->env == NULL)  {
			CPXgeterrorstring (g->env, status, errmsg);
			printf("Pb CPXopenCPLEX, code erreur = %d, message = %s\n", status, errmsg);
		}
		else {
			/* probleme */
			sprintf(nom_pb, "%d_rho_%s_%d_%d_%d", g->resol, (g->E_q == 1 ? "_E" : ""), g->nu, g->q, g->k);
			g->pl = CPXcreateprob (g->env, &status, nom_pb);
			if (g->pl == NULL)  {
				printf("Pb CPXcreateprob, code erreur = %d\n", status);
			}
		}
	}

	/* ____ Indices contraintes BkI dans le PL */
	if (status == 0) {
		J_max =(int)get_J_max(g);
		a_max =get_a_max(g);

		/* Contraintes BkI */
		g->ctr_bkI =malloc((J_max +1) * sizeof(int*));
		if (g->ctr_bkI == NULL)  {
			printf("Pb malloc g->ctr_bkI.\n");
			status =STAT_ERROR_MALLOC;
		}
		else {
			for(J =0 ; J <= J_max ; J ++) {
				g->ctr_bkI[J] =malloc((a_max +1) * sizeof(int));
				if (g->ctr_bkI[J] == NULL)  {
					printf("Pb malloc g->ctr_bkI[%d].\n", J);
					status =STAT_ERROR_MALLOC;
				}
				else {
					for(a =0 ; a <= a_max ; a ++) {
						g->ctr_bkI[J][a] =NO_VAL;
					}
				}
			}
		}
	}

	/* ____ Le cas echeant : on fait le menage */
	if(status != 0) {
		liberer(g);
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Sous-routine des fonctions d'initialisation : verification de la coherence des donnees
	valeur retournee : 0 si tout ok et -1 sinon
*/
static int verification(rho_pl* g) {
	int status =STAT_ERROR_PARAM;

	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
	#endif

	if(g->k > g->nu || g->q < 2 || g->k < 2) {
		printf("%s Pb donnees : parametres nu (%d) / q (%d) / k (%d) incoherents\n", __func__, g->nu, g->q, g->k);
	}
	else if(g->E_q != 0 && g->E_q != 1) {
		printf("%s Pb donnees : parametre E_q (%d) incoherent\n", __func__, g->E_q);
	}
	else if(g->resol == RESOL_BIN_PV_OPT && g->R_cible <= 1) {
		printf("%s Pb donnees : resolution RESOL_BIN_PV_OPT mais R_cible (%ld) incoherent\n", __func__, g->R_cible);
	}
	else if(g->resol == RESOL_INT && g->val_cible <= 0) {
		printf("%s Pb donnees : resolution RESOL_INT mais val_cible (%lf) incoherent\n", __func__, g->val_cible);
	}
	else if(g->resol == RESOL_INT_PV_OPT && (g->num_cible <= 0 || g->den_cible <= 0 || g->num_cible > g->den_cible) ) {
		printf("%s Pb donnees : resolution RESOL_INT_PV_OPT mais den_cible (%ld)/num_cible (%ld) incoherents\n", __func__, g->num_cible, g->den_cible);
	}
	else {
		status = 0;
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* ____ nettoyage
	Pre-conditions : champs g->q, g->k, g->_q correctement renseignes
	Post-conditions : memoire liberee
*/
int liberer(rho_pl* g) {
	int status = 0, J, J_max =(int)get_J_max(g);
	char errmsg[1024];

	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
	#endif

	/* PL Cplex */
	if (g->pl != NULL) {
		status = CPXfreeprob(g->env, & g->pl);
		if (status) {
			printf("Pb CPXfreeprob, code erreur = %d.\n", status);
		}
	}

	/* Environnement Cplex */
	if (g->env != NULL) {
		status = CPXcloseCPLEX(& g->env);
		if (status) {
			CPXgeterrorstring (g->env, status, errmsg);
			printf("Pb CPXcloseCPLEX, code erreur = %d, message = %s\n", status, errmsg);
		}
	}

	/* Contraintes BkI */
	if (g->ctr_bkI != NULL) { 
		for(J =0 ; J <= J_max ; J ++) {
			if (g->ctr_bkI[J] != NULL)  {
				free(g->ctr_bkI[J]);
			}
		}

		free(g->ctr_bkI);
		g->ctr_bkI = NULL;
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* __________________________________________________________________________ Definition du PL
*/

/* ____ definition du PL
	Variables :
		- Q(v), v in Z_q^nu	: toujours, a valeur dans [0, 1], N, {0, 1} selon g->resol
		- R in N				: si g->resol est RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	Contraintes :
		- Q mesure										si RESOL_CONT : 
		- R = sum_v Q(v)	 								si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
		- Q(0, 0 ,.., 0) >= (ou =) 1 					si RESOL_INT ou RESOL_BIN ou RESOL_INT_PV_OPT 
		- Q(0, 0 ,.., 0)/R >= val_cible					si RESOL_INT
		- Q(0, 0 ,.., 0)/R >= num_cible/den_cible		si RESOL_INT_PV_OPT
		- sum_v Q(v) = R_cible	 						si RESOL_BIN_PV_OPT
		- Q (ou Q_E) in I_q^k							toujours
	Objectif :
		- Max Q(0, 0 ,.., 0)								si RESOL_CONT ou RESOL_BIN_PV_OPT
		- Min R											si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	Pre-conditions : g intialisee, parametres coherents
	Post-condition : PL totalement defini (dont type et sens optimisation), tableau ctr_bkI rempli
*/
int setPL(rho_pl* g) {
	int status = 0;	

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* Contraintes (sans les coefficients des variables) */
	status = rho_ctr(g);
	if (status)  {
		printf("Pb generation contraintes PL.\n");
	}
	else {
		/* variables (avec bornes, coefficient objectif, coefficient dans les contraintes) */
		status = rho_var(g);
		if (status)  {
			printf("Pb generation variables PL.\n");
		}
		else {
			status = finaliser(g);
			if (status)  {
				printf("Pb definition type de probleme.\n");
			}
		}
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d :: %d variables et %d contraintes creees\n", __func__, status, CPXgetnumcols(g->env, g->pl), CPXgetnumrows(g->env, g->pl));
    #endif

	return status;
}

/* Sous-routine de la fonction de definition du PL
	Contraintes :
	- R = sum_v Q(v)	 								si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	- Q mesure										si RESOL_CONT : 
	- sum_v Q(v) = R_cible	 						si RESOL_BIN_PV_OPT
	- Q(0, 0 ,.., 0) >= (ou =) 1 					si RESOL_INT ou RESOL_BIN ou RESOL_INT_PV_OPT 
	- Q(0, 0 ,.., 0)/R >= val_cible					si RESOL_INT
	- Q(0, 0 ,.., 0)/R >= num_cible/den_cible		si RESOL_INT_PV_OPT
	- Q (ou Q_E) in I_q^k							toujours
	Pre-conditions : g intialisee, parametres coherents
	Post-condition : contraintes crees (sans les coefficients des variables), tableau ctr_bkI rempli
*/
static int rho_ctr(rho_pl* g) {
	int status = 0, nb_ctr, J, J_max =(int)get_J_max(g), a, a_max =get_a_max(g);
	double rhs;
	char sens;
	char* nom_ctr = malloc(100 * sizeof(char));

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ R = sum_v Q(v) ssi sum_v Q(v) -R = 0 (contrainte d'indice 0) :
	*/
	sens = 'E';
	rhs = 0;
	sprintf(nom_ctr, "R_def");
	status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);

	if(status == 0) {
		/* __ R (ssi sum_v Q(v)) = une valeur donnee (contrainte d'indice 1) :
			R = 1 (Q mesure)		si RESOL_CONT
			R = R_cible 			si RESOL_BIN_PV_OPT
		*/
		if(g->resol == RESOL_CONT || g->resol == RESOL_BIN_PV_OPT) {
			sens ='E';
			rhs = (g->resol == RESOL_CONT ? 1 : g->R_cible);
			sprintf(nom_ctr, "R_val");
			status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
		}

		/* __ Q(0, 0 ,.., 0) >= 1 (contrainte d'indice 1) :
			Q(0, 0 ,.., 0)  = 1	si RESOL_BIN
			Q(0, 0 ,.., 0) >= 1	si RESOL_INT ou RESOL_INT_PV_OPT
		*/
		else if(g->resol == RESOL_BIN || g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT) {
			sens = (g->resol == RESOL_BIN ? 'E' : 'G');
			rhs = 1;
			sprintf(nom_ctr, "Q_star_G_1");
			status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
		}
	}

	/* __ Q(0, 0 ,.., 0)/R au moins (>= ou >) une certaine valeur (contrainte d'indice 2) :
		Q(0, 0 ,.., 0) -val_cible x R 				>= 0		si RESOL_INT
		num_cible x Q(0, 0 ,.., 0) -den_cible x R 	>= 1		si RESOL_INT_PV_OPT  
	*/
	if( status == 0 && (g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT) ) {
		sens ='G';
		rhs = (g->resol == RESOL_INT ? 0 : 1);
		sprintf(nom_ctr, "%s", (g->resol == RESOL_INT ? "val_cible" : "val_cible_G"));
		status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
	}

	/* __ Q in I_q^k ou Q_E in I_q^k :
		E_q =0 : pour tout sous-ensemble J de taille k de [q] et tout a in Z_q^k, 
					sum_{v in Z_q^q: v_J =a} Q(v) = R/q^k
				<=>	sum_{v in Z_q^q: v_J =a} Q(v) -R/q^k =0
		E_q =1 : pour tout sous-ensemble J de taille k de [q] et tout a in {0} x Z_q^(k -1), 
					sum_{v in U: v_J in {a, a +(1 ,.., 1) ,.., a +(q -1 ,.., q -1)}} Q(v) = R/q^(k -1)
				<=> sum_{v in U: v_J in {a, a +(1 ,.., 1) ,.., a +(q -1 ,.., q -1)}} Q(v) -R/q^(k -1) =0
	*/
	nb_ctr = CPXgetnumrows(g->env, g->pl);
	sens = 'E';
	rhs = 0;

	/* parcours des sous-ensembles J de cardinalite k de {0 ,.., q -1} */
	for (J = 0 ; (J <= J_max) && (status == 0) ; J ++) {
		if(get_card(J, g->nu) == g->k) {
			/* parcours des vecteurs a de Z_q^k si E_q = 0, de {0} x Z_q^(k -1) sinon */
			for (a = 0 ; (a <= a_max) && (status == 0) ; a ++) {
				/* nouvelle contrainte (contrainte numero nb_ctr +1 ssi contrainte d'indice nb_ctr)  */
				sprintf(nom_ctr, "BkI_%d_%d", J, a);
				status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
				if (status == 0) {
					g->ctr_bkI[J][a] =nb_ctr;
					nb_ctr ++;
				}
			}
		}
	}

	/* nettoyage */
	free(nom_ctr);

   	#if (DEBUG)
    printf("%s OUT : status = %d, %d contraintes creees\n", __func__, status, CPXgetnumrows(g->env, g->pl));
    #endif

	return status;
}

/* Sous-routine de la fonction de definition du PL
	Variables:
		- Q(v), v in Z_nu^q	: toujours, a valeur dans [0, 1], N, {0, 1} selon g->resol
		- R					: toujours, a valeur dans [0, 1] ou N selon g->resol
	Pre-conditions : g intialisee, contraintes toutes generees (sans les coefficients des variables), tableau ctr_bkI rempli, parametres coherents
	Post-condition : PL totalement defini (sinon type et sens optimisation)
*/
static int rho_var(rho_pl* g) {
	int status =0, v, J, J_max =(int)get_J_max(g), vJ, ind_var_Q =0, v_max =get_v_max(g); /* v <= q^nu -1 si E_q = 0, q^(nu -1) -1 sinon */
	int coef_bkI =get_a_max(g) +1; /* q^k si E_q =0 et q^(k -1) sinon */
	double coef_obj, lb, ub;
	char type;
	char* nom_var = malloc(100 * sizeof(char));

	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
	#endif

	/* __ R : generee toujours. Participe a :
		- obj 													si RESOL_BIN, RESOL_INT ou RESOL_INT_PV_OPT
		- sum_v Q(v) -R = 0 										toujours
		- R = une valeur donnée									si RESOL_CONT ou RESOL_INT_PV_OPT
		- Q(0 ,.., 0)/R au moins (>= ou >) une certaine valeur	si RESOL_INT ou RESOL_INT_PV_OPT 
		- sum_(v : v_J =a) Q(v) -R/q^k = 0 (ou version E_q)		toujours
	*/

	/* generation */
	sprintf(nom_var, "R");
	lb =0;

	if(g->resol == RESOL_CONT) {
		type = CPX_CONTINUOUS;
		ub =1;
		coef_obj =0;
	}
	else {
		type = CPX_INTEGER;
		ub =CPX_INFBOUND;
		coef_obj =(g->resol == RESOL_CONT || g->resol == RESOL_BIN_PV_OPT ? 0 : 1);
	}

	status = CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom_var);

	/* participation de R aux contraintes (sauf BkI) */
	if(status == 0) {
		/* contrainte (d'indice 0) sum_v P(v) -R = 0 */
		status = CPXchgcoef(g->env, g->pl, IND_ROW_DEF_R, IND_COL_R, -1);

		if (status == 0) {
			/* contrainte (d'indice 1) R =1 (Q mesure) ou R =R_cible (d'indice 1) si RESOL_CONT ou RESOL_BIN_PV_OPT */
			if(g->resol == RESOL_CONT || g->resol == RESOL_BIN_PV_OPT) {
				status = CPXchgcoef(g->env, g->pl, IND_ROW_R_IS_EQUAL_TO, IND_COL_R, 1);
			}
			/* contrainte (d'indice 2) Q(0 ,.., 0)/R au moins (>= ou >) une certaine valeur :
				Q(0, 0 ,.., 0) -val_cible x R 				>= 0		si RESOL_INT
				num_cible x Q(0, 0 ,.., 0) -den_cible x R 	>= 1		si RESOL_INT_PV_OPT  
			*/
			else if (g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT) {
				status = CPXchgcoef(g->env, g->pl, IND_ROW_LB, IND_COL_R, -(g->resol == RESOL_INT ? g->val_cible : g->den_cible) );
			}
		}
	}

	/* __ Q(v) : generees toujours. Participe a :
		- obj 													si v = (0 ,..., 0) et RESOL_CONT ou RESOL_BIN_PV_OPT
		- Q(0 ,.., 0) >= (ou =) 1								si v = (0 ,..., 0) et RESOL_BIN, RESOL_INT ou RESOL_INT_PV_OPT 
		- Q(0 ,.., 0)/R au moins (>= ou >) une certaine valeur	si v = (0 ,..., 0) et RESOL_INT ou RESOL_INT_PV_OPT 
		- sum_v Q(v) -R = 0 										toujours
		- sum_(v : v_J =a) Q(v) -R/q^k = 0 (ou version E_q)		toujours
	*/
	if(status == 0) {
		if(g->resol == RESOL_CONT) {
			type = CPX_CONTINUOUS;
			lb =0;
			ub =1;
		}
		else if(g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT || g->resol == RESOL_BIN_PV_OPT) {
			type = CPX_INTEGER;
			lb =0;
			ub =CPX_INFBOUND;
		}
		else if(g->resol == RESOL_BIN) {
			type = CPX_BINARY;
			lb =0;
			ub =1;
		}
	}

	/* parcours des vecteurs v de Z_q^nu si E_q = 0, de {0} x Z_q^(nu -1) si E_q = 1 */
	for (v = 0 ; (v <= v_max) && (status == 0) ; v ++) {
		ind_var_Q ++;

		/* generation : 
			- Q(v) genere pour tout v in 0..v_max 
			- seul Q(0 ,.., 0) peut participer à la fonction ojectif, et il ne le fait que pour RESOL_CONT et RESOL_BIN_PV_OPT
		*/
		sprintf(nom_var, "Q%d", v);
		coef_obj = ( (v == 0) && (g->resol == RESOL_CONT || g->resol == RESOL_BIN_PV_OPT) );
		status = CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom_var);

		/* participation specifique de Q(0 ,.., 0) aux contraintes */
		if(status == 0 && v == 0) {
			/* contrainte (d'indice 1) Q(0 ,.., 0) >= (ou =) 1 : RESOL_BIN, RESOL_INT ou RESOL_INT_PV_OPT
			*/
			if(g->resol == RESOL_BIN || g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT) {
				status = CPXchgcoef(g->env, g->pl, IND_ROW_Q_G_1, ind_var_Q, 1);
			}

			/* contrainte (d'indice 2) Q(0 ,.., 0)/R au moins (>= ou >) une certaine valeur :
				Q(0 ,.., 0) -val_cible x R 				>= 0		si RESOL_INT
				num_cible x Q(0 ,.., 0) -den_cible x R 	>= 1		si RESOL_INT_PV_OPT  
			*/
			if(status == 0 && (g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT)) {
				status = CPXchgcoef(g->env, g->pl, IND_ROW_LB, ind_var_Q, (g->resol == RESOL_INT ? 1 : g->num_cible));
			}
		}

		/* participation de Q(v) a la contrainte (d'indice 0) sum_v Q(v) -R =0 */
		if(status == 0) {
			status = CPXchgcoef(g->env, g->pl, IND_ROW_DEF_R, ind_var_Q, 1);
		}

		/* participation de Q(v) aux contraintes BkI d'indices (J =(j_1 ,.., j_k), a =(a_1 ,.., a_t)) vérifiant :
			- v_J = a 														si E_q = 0
			- v_J in {a, a +(1 ,.., 1) ,.., a +(q -1 ,.., q -1)} ou a_1 =0 	si E_q = 1
		*/

		/* parcours des sous-ensembles J de taille k de Z_q */
		for (J = 0 ; (J <= J_max) && (status == 0) ; J ++) {
			if(get_card(J, g->nu) == g->k) {
				/* recuperation de v_J si E_q = 0, de v_J -(v_{j_1} ,.., v_{j_1}) si E_q = 1 */
				vJ =get_v_J(v, J, g->q, g->nu, g->k, g->E_q);

				/* +coef_bkI pour Q et -1 pour R, ou coef_bkI =q^k si E_q =0 et q^(k -1) sinon */
				status = CPXchgcoef(g->env, g->pl, g->ctr_bkI[J][vJ], ind_var_Q, coef_bkI);		/* Q(v) */

				/* coefficent de R : TODO optimiser (ici on le fait plusieurs fois ...) */
				if(status == 0) {
					status = CPXchgcoef(g->env, g->pl, g->ctr_bkI[J][vJ], IND_COL_R, -1);		/* R */
				}
			}
		}
	}

	/* nettoyage */
	free(nom_var);

	#if (DEBUG)
    printf("%s OUT : status = %d, %d variables creees\n", __func__, status, CPXgetnumcols(g->env, g->pl));
	#endif

	return status;
}

/* Sous-routine de la fonction de definition du PL
	Sens optimisation :
		- Max 		si RESOL_CONT ou RESOL_BIN_PV_OPT
		- Min		si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	Type probleme :
		- LP			si RESOL_CONT
		- MILPL		dans tous les autres cas
	Pre-conditions : PL totalement defini, sauf type et sens optimisation
	Post-condition : type et sens optimisation definis
*/
static int finaliser(rho_pl* g) {
	int status =0;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Sens optimisation :
		Max Q(0, 0 ,.., 0)	si RESOL_CONT ou RESOL_BIN_PV_OPT
		Min R				si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	*/
	CPXchgobjsen(g->env, g->pl, (g->resol == RESOL_CONT || g->resol == RESOL_BIN_PV_OPT ? CPX_MAX : CPX_MIN ) );

	/* __ Resolution
	*/

	/* Continue */
	if(g->resol == RESOL_CONT) {
		status = CPXchgprobtype(g->env, g->pl, CPXPROB_LP);
	}
	/* Entiere */
	else {
		status = CPXchgprobtype(g->env, g->pl, CPXPROB_MILP);
	}

	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
	#endif

	return status;
}

/* __________________________________________________________________________ Resolution du PL
*/

/* ____ resolution PL :
	Pre-condition : PL totalement defini (sinon type et sens optimisation), tableau ctr_bkI rempli
	Post-conditions : 
		- sens obj et type de probleme definis
		- PL resolu
	Affichages : restitution valeur optimale

	TODO Cplex : comprendre necessite de specifier pb continu qd on la deja fait a la construction du PL, et que l'on a ensuite insere des variables continues
*/
int resoudre(rho_pl* g) {
	int status =0;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Sens optimisation :
		Max Q(0, 0 ,.., 0)	si RESOL_CONT ou RESOL_BIN_PV_OPT
		Min R				si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	*/
	CPXchgobjsen(g->env, g->pl, (g->resol == RESOL_CONT || g->resol == RESOL_BIN_PV_OPT ? CPX_MAX : CPX_MIN ) );

	/* __ Resolution
	*/

	/* Continue */
	if(g->resol == RESOL_CONT) {
		status = CPXchgprobtype(g->env, g->pl, CPXPROB_LP);
		if (status)  {
			printf("\t******** Pb type de probleme.\n");
		}
		else {
			status = CPXlpopt(g->env, g->pl);
			if (status)  {
				printf("\t******** Pb optimisation.\n");
			}
		}
	}
	/* Entiere */
	else {
		status = CPXchgprobtype(g->env, g->pl, CPXPROB_MILP);
		if (status)  {
			printf("\t******** Pb type de probleme.\n");
		}
		else {
			status = CPXmipopt(g->env, g->pl);
			if (status)  {
				printf("\t******** Pb optimisation.\n");
			}
		}
	}

	/* __ Restitution
	*/
	if (status == 0)  {
		status = restitution(g);
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Sous-routine de la fonction de resolution : affichage synthetique du resultat de la resolution
	Pre-condition : le PL a ete resolu
	NB statuts resolution :
		https://www.tu-chemnitz.de/mathematik/discrete/manuals/cplex/doc/refman/html/appendixB.html
		1 :	CPX_STAT_OPTIMAL
		2 :	CPX_STAT_UNBOUNDED
		3 :	CPX_STAT_INFEASIBLE 
		101 : CPXMIP_OPTIMAL 
		103 : CPXMIP_INFEASIBLE 
		118 : CPXMIP_UNBOUNDED 
*/
static int restitution(rho_pl* g) {
	int status =0, stat_sol;
	double val_opt, val_Q_star;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Probleme resolu */
	printf("\t******** rho_%s(%d, %d, %d) ", (g->E_q == 1 ? "_E" : ""), g->nu, g->q, g->k);
	switch(g->resol) {
		case RESOL_CONT:
			printf("resolution CONTINUE (%d)\n", g->resol);
			break;

		case RESOL_INT:
			printf("resolution ENTIERE (%d) valeur cible >= %lf\n", g->resol, g->val_cible);
			break;

		case RESOL_BIN:
			printf("resolution BINAIRE (%d)\n", g->resol);
			break;

		case RESOL_INT_PV_OPT:
			printf("resolution ENTIERE (%d) valeur cible > %ld/%ld = %lf\n", g->resol, g->num_cible, g->den_cible, (double)g->num_cible/g->den_cible);
			break;

		case RESOL_BIN_PV_OPT:
			printf("resolution ENTIERE (%d) nb lignes cible = %ld\n", g->resol, g->R_cible);
			break;

		default:
			printf("resolution (%d) non reconnue\n", g->resol);
			break;
	}

	/* __ Solution obtenue */
	stat_sol = CPXgetstat(g->env, g->pl);
	if(stat_sol == CPX_STAT_UNBOUNDED || stat_sol == CPXMIP_UNBOUNDED) {
		printf("\t\t**** probleme de valeur infinie.\n");
	}
	else if(stat_sol == CPX_STAT_INFEASIBLE || stat_sol == CPXMIP_INFEASIBLE) {
		printf("\t\t**** probleme non realisable.\n");
	}
	else if(stat_sol != CPX_STAT_OPTIMAL && stat_sol != CPXMIP_OPTIMAL) {
		printf("\t\t**** statut de resolution %d non traite (voir codes Cplex).\n", stat_sol);
	}
	else {
		printf("\t\t**** solution optimale trouvee ");

		status = CPXgetobjval(g->env, g->pl, &val_opt);
		if(status)  {
			printf("MAIS Pb recuperation val opt\n");
		}
		else if(g->resol == RESOL_CONT) {
			printf("de valeur %lf (inverse %lf) ", val_opt, 1/val_opt);
		}
		else if(g->resol == RESOL_BIN_PV_OPT) {
			printf("de valeur %lf/%ld = %lf (inverse %lf)", val_opt, g->R_cible, val_opt/g->R_cible, g->R_cible/val_opt);
		}
		else if(g->resol == RESOL_BIN) {
			printf("de valeur 1/%lf = %lf", val_opt, 1.0/val_opt);
		}
		else if(g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT) {
			printf("utilisant %lf lignes ", val_opt);

			status = CPXgetx(g->env, g->pl, &val_Q_star, IND_COL_Q_STAR, IND_COL_Q_STAR);
			if(status != 0) {
				printf("MAIS Pb recuperation val Q(0, 0 ,.., 0).\n");
			}
			else {
				printf(" de valeur %lf/%lf = %lf (inverse %lf)\n", val_Q_star, val_opt, val_Q_star/val_opt, val_opt/val_Q_star);
			}
		}
		else {
			printf("MAIS resolution (%d) non reconnue.\n", g->resol);
		}
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* __________________________________________________________________________ Affichage solution
*/

/* affichage de la solution
*/
int affSol(rho_pl* g) {
	int v, status =0, ind_var;
	int v_max =get_v_max(g);		// intervalle pour v : {0 ,.., q^nu -1} si E_q = 0, {0 ,.., q^(nu -1) -1} sinon
	double valQ;
	char nom_var[100];

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Parcours des vecteurs v in (Z_q)^q représentés par des nombres i in Z_(q^q) dont v est la décomposition en base q */
	for (v = 0 ; (v <= v_max) && (status == 0) ; v ++) {
		sprintf(nom_var, "Q%d", v);
		status = CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
		if(status == 0) {
			status = CPXgetx(g->env, g->pl, &valQ, ind_var, ind_var);
		}

		if (status == 0 && valQ > 0) {
			/* Affichage si P(v) > 0 ou Q(v) > 0 */
			printf("%d en base %d:\t", v, g->q);
			affiche_v(v, g->q, g->nu);
			printf("\tval: %lf\t%s\n", valQ, (v == 0 ? "(*)" : ""));
		}
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* __________________________________________________________________________ Divers (sous-routines)
*/

/* Intervalle de travail pour v : {0 ,.., q^nu -1} si E_q = 0, {0 ,.., q^(nu -1) -1} sinon
*/
static int get_v_max(rho_pl* g) {
	return pow(g->q, g->nu -g->E_q) -1;
}

/* Intervalle de travail pour J : {0 ,.., 2^nu -1}
*/
static int get_J_max(rho_pl* g) {
	return pow(g->q, g->nu -g->E_q) -1;
}

/* Intervalle de travail pour a rhs contrainte bkI : {0 ,.., q^k -1} si E_q = 0, {0 ,.., q^(k -1) -1} sinon
*/
static int get_a_max(rho_pl* g) {
	return pow(g->q, g->k -g->E_q) -1;
}

