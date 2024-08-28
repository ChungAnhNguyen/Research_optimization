/*	commun-pl
*/

#if (DEBUG)
#include <stdio.h>	// pour traces
#endif

#include <string.h>	// pour sprintf()
#include <math.h> 	// pour pow()

#define IND_COL_R 0
#define IND_COL_Q_STAR 1

#define IND_ROW_DEF_R 0
#define IND_ROW_R_IS_EQUAL_TO 1
#define IND_ROW_Q_G_1 1
#define IND_ROW_LB 2

enum e_resolution {
	RESOL_CONT =0,
	RESOL_INT,
	RESOL_INT_PV_OPT,
	RESOL_BIN,
	RESOL_BIN_PV_OPT
};

enum e_tableaux {
	TAB_RHO,
	TAB_GAMMA
};

struct s_tableaux{
	int nu;					// =q pour gamma
	int q;
	int k;
	int p;					// non exploite pour rho
	int E_q;
	enum e_resolution resol;	// choix de resolution : continu, entier, ...
	double val_cible; 		// valeur cible : quand recherche solution entiere minimisant R quand valeur continue connue
	long int num_cible; 	// num_cible et den_cible : quand on veut demontrer solution entiere Q/R =den_cible / num_cible est opt (du fait arrondis valeurs irationnelles)
	long int den_cible; 	
	long int R_cible;	// R_cible : quand on veut demontrer qu'avec le nombre de lignes min pour des tableaux simples, on ne peut faire mieux que 1/R_cible 	
	int **ctr_bkI;		// associe a (J, a) l'indice de la contrainte associee dans le PL => pour simplifier l'access, tableau de taille 2^nu (tous les sous-ens de Z_q)
	CPXENVptr env;
	CPXLPptr pl;
	enum e_tableaux type;
};

typedef struct s_tableaux tableaux;

#define STAT_ERROR_MALLOC -2
#define STAT_ERROR_PARAM -1
#define NO_VAL -1	/* utilisee pour les champs val_cible et ctr_bkI[J][a]  */

/* __________________________________________________________________________ Declaration fonctions statiques
*/

/* Sous-routine des fonctions d'initialisation :
	Pre-condition : 		champs resol, val_cible, num_cible, den_cible, R_cible tous initialises
	Post-condition : 	champs atomiques q, p, k, E_q tous initialises
	Valeur retournee : 	0 si ok, -1 si pb parametres, -2 si pb malloc, code cplex > 0 sinon
	Post-condition en cas de succes :
		- environnement Cplex +PL crees
		- tableau ctr_bkI[][] construit et initialise
	Post-condition en cas d'echec : memoire eventuellement allouee en cours de fonction totalement liberee
*/
static int init(tableaux* t, int nu, int q, int p, int k, int E_q);

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
static int finaliser(gamma_pl* g);

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
static int restitution(gamma_pl* g);

/* __________________________________________________________________________ Initialisation / destruction PL
*/

static void init_param(int n, int q, int p, int k, int E_q) {
	g->q =q;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = NO_VAL;
}

/* ____ Fonctions d'initialisation :
	Post-condition : 	champs atomiques q, p, k, E_q, resol, val_cible, num_cible, den_cible, R_cible tous initialises
	Valeur retournee : 	0 si ok, -1 si pb parametres, -2 si pb malloc, code cplex > 0 sinon
	Post-condition en cas de succes :
		- environnement Cplex +PL crees
		- tableau ctr_bkI[][] construit et initialise
	Post-condition en cas d'echec : memoire eventuellement allouee en cours de fonction totalement liberee
*/

/* Initialisation resolution continue */
int init_cont(tableaux* g, int n, int q, int p, int k, int E_q) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %d) IN\n", __func__, g, nu, q, p, k, E_q);
	#endif

	g->resol = RESOL_CONT;
	g->val_cible = NO_VAL;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = NO_VAL;

	init(g, q, p, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation resolution entiere */
int init_int(gamma_pl* g, int n, int q, int p, int k, int E_q, double val_cible) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %lf) IN\n", __func__, g, q, p, k, E_q, val_cible);
	#endif

	g->resol = RESOL_INT;
	g->val_cible = val_cible;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = NO_VAL;

	init(g, q, p, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation PV optimalite solution entiere */
int init_int_pv_opt(gamma_pl* g, int n, int q, int p, int k, int E_q, long int num_cible, long int den_cible) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %ld, %ld) IN\n", __func__, g, q, p, k, E_q, num_cible, den_cible);
	#endif

	g->resol = RESOL_INT_PV_OPT;
	g->val_cible = NO_VAL;
	g->num_cible = num_cible;
	g->den_cible = den_cible;
	g->R_cible = NO_VAL;

	init(g, q, p, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation resolution binaire */
int init_bin(gamma_pl* g, int n, int q, int p, int k, int E_q) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d) IN\n", __func__, g, q, p, k, E_q);
	#endif

	g->resol = RESOL_BIN;
	g->val_cible = NO_VAL;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = NO_VAL;

	init(g, q, p, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Initialisation resolution binaire PV ratio 1/R_cible est opt s.c. R =R_cible */
int init_bin_pv_opt(gamma_pl* g, int n, int q, int p, int k, int E_q, long int R_cible) {
	int status =0;

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d, %ld) IN\n", __func__, g, q, p, k, E_q, R_cible);
	#endif

	g->resol = RESOL_BIN_PV_OPT;
	g->val_cible = NO_VAL;
	g->num_cible = NO_VAL;
	g->den_cible = NO_VAL;
	g->R_cible = R_cible;

	init(g, q, p, k, E_q);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Sous-routine des fonctions d'initialisation :
	Pre-condition : 		champs resol, val_cible, num_cible, den_cible, R_cible tous initialises
	Post-condition : 	champs atomiques q, p, k, E_q tous initialises
	Valeur retournee : 	0 si ok, -1 si pb parametres, -2 si pb malloc, code cplex > 0 sinon
	Post-condition en cas de succes :
		- environnement Cplex +PL crees
		- tableau ctr_bkI[][] construit et initialise
	Post-condition en cas d'echec : memoire eventuellement allouee en cours de fonction totalement liberee
*/
static int init(gamma_pl* g, int q, int p, int k, int E_q) {
	int status = 0, J, nb_J =(int)pow(2, q), a, a_max;
	char nom_pb[100];
	char errmsg[1024];

	#if (DEBUG)
    printf("%s(%p, %d, %d, %d, %d) IN\n", __func__, g, q, p, k, E_q);
	#endif

	/* ____ Parametres probleme & intialisations */
	g->q =q;
	g->p =p;
	g->k =k;
	g->E_q =E_q;

	g->ind_Q_star =NO_VAL;
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
			sprintf(nom_pb, "%d_gamma_%s_%d_%d_%d", g->resol, (g->E_q == 1 ? "_E" : ""), g->q, g->p, g->k);
			g->pl = CPXcreateprob (g->env, &status, nom_pb);
			if (g->pl == NULL)  {
				printf("Pb CPXcreateprob, code erreur = %d\n", status);
			}
		}
	}

	/* ____ Indices contraintes BkI dans le PL */
	if (status == 0) {
		/* Contraintes BkI */
		g->ctr_bkI =malloc(nb_J * sizeof(int*));
		if (g->ctr_bkI == NULL)  {
			printf("Pb malloc g->ctr_bkI.\n");
			status =STAT_ERROR_MALLOC;
		}
		else {
			a_max =get_a_max(g);
			for(J =0 ; J < nb_J ; J ++) {
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

/* ____ nettoyage
	Pre-conditions : champs g->q, g->k, g->_q correctement renseignes
	Post-conditions : memoire liberee
*/
int liberer(gamma_pl* g) {
	int status = 0, J, nb_J =(int)pow(2, g->q);
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
		for(J =0 ; J < nb_J ; J ++) {
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
int setPL(gamma_pl* g) {
	int status = 0;	

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* Contraintes (sans les coefficients des variables) */
	status = ctr(g);
	if (status)  {
		printf("Pb generation contraintes PL.\n");
	}
	else {
		/* variables (avec bornes, coefficient objectif, coefficient dans les contraintes) */
		if (g->type == TAB_GAMMA)
			status = gamma_var(g);
		else if (g->type == TAB_RHO)
			status = rho_var(g);
		
		if (status) 
			printf("Pb generation variables PL.\n");
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
	Sens optimisation :
		- Max 		si RESOL_CONT ou RESOL_BIN_PV_OPT
		- Min		si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	Type probleme :
		- LP			si RESOL_CONT
		- MILPL		dans tous les autres cas
	Pre-conditions : PL totalement defini, sauf type et sens optimisation
	Post-condition : type et sens optimisation definis
*/
static int finaliser(gamma_pl* g) {
	int status =0;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Sens optimisation :
		Max Q(0, 1 ,.., q -1)	si RESOL_CONT ou RESOL_BIN_PV_OPT
		Min R					si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	*/
	CPXchgobjsen(g->env, g->pl, (g->resol == RESOL_CONT || g->resol == RESOL_BIN_PV_OPT ? CPX_MAX : CPX_MIN ) );

	/* __ Resolution
	*/

	/* Continue */
	if(g->resol == RESOL_CONT)
		status = CPXchgprobtype(g->env, g->pl, CPXPROB_LP);
	/* Entiere */
	else
		status = CPXchgprobtype(g->env, g->pl, CPXPROB_MILP);

	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
	#endif

	return status;
}

/* __________________________________________________________________________ Resolution du PL
*/

/* ____ resolution PL :
	Pre-condition : PL totalement defini (sinon type et sens optimisation), champ ind_Q_star instancie, tableau ctr_bkI rempli
	Post-conditions : 
		- sens obj et type de probleme definis
		- PL resolu
	Affichages : restitution valeur optimale

	TODO Cplex : comprendre necessite de specifier pb continu qd on la deja fait a la construction du PL, et que l'on a ensuite insere des variables continues
*/
int resoudre(gamma_pl* g) {
	int status =0;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Sens optimisation :
		Max Q(0, 1 ,.., q -1)	si RESOL_CONT ou RESOL_BIN_PV_OPT
		Min R					si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
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
static int restitution(gamma_pl* g) {
	int status =0, stat_sol;
	double val_opt, val_Q_star;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Probleme resolu */
	printf("\t******** gamma_%s(%d, %d, %d) ", (g->E_q == 1 ? "_E" : ""), g->q, g->p, g->k);
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

			status = CPXgetx(g->env, g->pl, &val_Q_star, g->ind_Q_star, g->ind_Q_star);
			if(status != 0) {
				printf("MAIS Pb recuperation val Q(0, 1 ,.., q -1).\n");
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

