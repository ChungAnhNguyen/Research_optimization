/* NB modif juillet 21:
int allouer_tableau(tableau * T, unsigned n);	-> appelée par: enumerer_nombre_incidence
void liberer_tableau(tableau * T);				-> appelée par: ctr_bkI, var_gen
void generer_nombre_incidence()					-> appelée par: generer_nombre_incidence, enumerer_nombre_incidence
void enumerer_nombre_incidence()				-> appelée par: ctr_bkI, var_gen
*/

#include "tableaux.h"	/* tableaux */
#include "vn.h" 		// representation des vecteurs et des sous-ensembles par des nombres

#include <stdlib.h>		// pour malloc(), itoa()
#include <string.h>		// pour sprintf()
#include <math.h> 		// pour pow()

#if (DEBUG)
#include <stdio.h>		// pour traces
#endif

/* ____________________________________________________________________________________________ Fonctions statiques
*/

/* _____________ Vérification des paramètres d'optimisation
	valeur retournee: 0 si tout ok et -1 (STAT_ERROR_PARAM) sinon
*/
static int check_param_optim(tableaux* g);

/* _____________ sous-routine formation nom de fichier
*/
static void get_nom_fic(tableaux* g, char nom_pl[100], const char* ext);

/* ____________________________________________________________________________________________ Gestion mémoire
*/
int liberer(tableaux* g) {
	int status =0;
	char errmsg[1024];

	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
	#endif

	/* PL Cplex */
	if (g->pl != NULL) {
		status =CPXfreeprob(g->env, & g->pl);
		if (status) {
			printf("Pb CPXfreeprob, code erreur =%d.\n", status);
		}
	}

	/* Environnement Cplex */
	if (g->env != NULL) {
		status =CPXcloseCPLEX(& g->env);
		if (status) {
			CPXgeterrorstring (g->env, status, errmsg);
			printf("Pb CPXcloseCPLEX, code erreur =%d, message =%s\n", status, errmsg);
		}
	}

	/* Contraintes BkI 
		(uniquement pour rho et gamma) 
	*/
	if (g->type == TAB_RHO || g->type == TAB_GAMMA)
		libererBkI(g);

   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* ____________________________________________________________________________________________ Initialisation / vérification des paramètres
*/

/* Initialisation & vérification des champs de la structure (dimensions) :: rho */
int rho_init_param(tableaux* g, int nu, int q, int k, int E_q) {
	int status =STAT_ERROR_PARAM;
	
	g->type =TAB_RHO;

	g->nu =nu;	/* dimension */
	g->q =q;	/* alphabet */
	g->k =k;	/* force */
	g->E_q =E_q;	/* Rho ou Rho_E */
	
	/* test des paramètres */
	if(g->k >g->nu || g->q <1 || g->k <1)	/* on veut nu >= k et min(q, k) >= 1 */
		printf("%s Pb données rho: paramètres nu (%d) / q (%d) / k (%d) incohérents\n", __func__, g->nu, g->q, g->k);
	else if(g->E_q != 0 && g->E_q != 1)		/* on veut E_q in {0, 1} */
		printf("%s Pb données: paramètre E_q (%d) incoherent\n", __func__, g->E_q);
	else
		status =0;

	return status;
}

/* Initialisation & vérification des champs de la structure (dimensions) :: gamma */
int gamma_init_param(tableaux* g, int q, int p, int k, int E_q, int p_phi) {
	int status =STAT_ERROR_PARAM;

	g->type =TAB_GAMMA;

	g->nu =q;	/* dimension */
	g->q =q;	/* alphabet */
	g->p =p;	/* restriction */
	g->k =k;	/* force */
	g->E_q =E_q;	/* Gamma ou Gamma_E */

	if(p_phi != NO_VAL)		/* possible restriction supplémentaire sur Phi */
		g->p_phi =p_phi;
	else
		g->p_phi =q;

	/* test des paramètres */
	if(g->p >g->q || g->k >g->p || g->k <1 || g->p_phi >g->q || g->p_phi <0)	/* on veut q >=p >=k >=1 et q >=p_phi >=0 */
		printf("%s Pb données gamma: paramètres q (%d) / p (%d) / k (%d) / p_phi (%d) incohérents\n", __func__, g->q, g->p, g->k, g->p_phi);
	else if(g->E_q != 0 && g->E_q !=1)
		printf("%s Pb données: paramètre E_q (%d) incohérent\n", __func__, g->E_q);
	else
		status =0;

	return status;
}

/* Initialisation & vérification des champs de la structure (dimensions) :: delta (réguliers ou non) */
int delta_init_param(tableaux* g, int nu, int p, int k, int p_phi) {
	int status =STAT_ERROR_PARAM;
	
	g->type =TAB_DELTA;
	
	g->nu =nu;		/* dimension */
	g->p =p;		/* restriction */
	g->k =k;		/* force */
	g->q =2;		/* alphabet (champ non exploité) */
	g->E_q =0;	/* initialisation technique: pour appel fait à la fonction "get_v_J" dans la fonction "var_gen" */

	if(p_phi != NO_VAL)		/* possible restriction supplémentaire sur Phi */
		g->p_phi =p_phi;
	else
		g->p_phi =nu -1;

	/* test des paramètres */
	if(g->p >g->nu || g->k >g->p || g->k <1 || g->p_phi >=g->nu || g->p_phi <0)	/* test nu >=p >=k >=1 et 0 <= p_phi <nu */
		printf("%s Pb données delta: paramètres nu (%d) / p (%d) / k (%d) / p_phi (%d) incohérents\n", __func__, g->nu, g->p, g->k, g->p_phi);
	else
		status =0;

	#if(TRACE_DELTA)
	printf("%s :: p_phi =%d, g->p_phi =%d\n", __func__, p_phi, g->p_phi);
	#endif
	
	return status;
}

/* Initialisation & vérification des champs de la structure (dimensions) :: delta réguliers */
int reg_init_param(tableaux* g, int nu, int p, int k, int p_phi) {
	int status =delta_init_param(g, nu, p, k, p_phi);	
	g->type =TAB_REG;
	
	return status;
}

/* _____________ Initialisation des paramètres d'optimisation
	valeur retournee: -1 (STAT_ERROR_PARAM) si mode de résolution non reconnu et 0 sinon
*/
int init_param_optim(tableaux* g, enum e_resolution resol, long int num_cible, long int den_cible, long double val_cible, long int R_cible) {
	int status =0;

	/* initialisation des champs  */
	g->resol =NO_VAL;
	g->val_cible =NO_VAL;
	g->num_cible =NO_VAL;
	g->den_cible =NO_VAL;
	g->R_cible =NO_VAL;

	/* __ attribution des valeurs */
	if(resol == RESOL_CONT)
		g->resol =RESOL_CONT;
	else if (resol == RESOL_INT) {
		g->resol =RESOL_INT;
		g->val_cible =val_cible;
		if(R_cible >0)
			g->R_cible =R_cible;
	} else if (resol == RESOL_BIN && g->type != TAB_REG) {
		g->resol =RESOL_BIN;
		if(R_cible >0)
			g->R_cible =R_cible;
	} else if (resol == RESOL_INT_PV_OPT) {
		g->resol =RESOL_INT_PV_OPT;
		g->num_cible =num_cible;
		g->den_cible =den_cible;
	} else if (resol == RESOL_MAX_FREQ && g->type != TAB_REG) {
		g->resol =RESOL_MAX_FREQ;
		g->R_cible =R_cible;
	} else if (resol == RESOL_R && g->type != TAB_REG) {
		g->resol =RESOL_R;
		if(R_cible >0)
			g->R_cible =R_cible;
	} else
		status =STAT_ERROR_PARAM;

	/* __ vérifications (+ éventuel ajout lb connue pour gamma) */
	if(status == 0) {
		/* si TAB_GAMMA et résolution RESOL_CONT *DO_LB*: on renseigne num_cible et den_cible pour exploiter une borne inférieure connue */
		if(g->type == TAB_GAMMA && g->resol == RESOL_CONT && DO_LB) {
			/* bornes inférieures 
				gamma_E(q, q, q) =gamma(q, q, q) =1 
				gamma_E(q, p, 2) >=gamma(q, p, 2) =delta(q, p, 2) =floor(p/2)ceil(p/2)/[(q -floor(p/2))(q -ceil(p/2))]
				gamma_E(q, p, k) >=gamma(q, p, k) >=gamma(q -p +k, k, k) =delta(q -p +k, k, k) =2/(T(q -p +k, k) +1) si q >k
			*/
			if(g->q == g->k)
				g->num_cible =g->den_cible =1;
			else if (g->k == 2) {
				if(g->p % 2 == 0) {
					g->num_cible =g->p/2 * g->p/2;
					g->den_cible =(g->nu -g->p/2) * (g->nu -g->p/2);
				} else {
					g->num_cible =(g->p -1)/2 * (g->p +1)/2;
					g->den_cible =(g->nu -(g->p -1)/2) * (g->nu -(g->p +1)/2);
				} 
			} else {
				g->num_cible =2;
				g->den_cible =Tqk(g->q -g->p +g->k, g->k) +1;
			}
		}

		status =check_param_optim(g);
	}

	return status;
}

/* _____________ Vérification des paramètres d'optimisation
	valeur retournee: 0 si tout ok et -1 (STAT_ERROR_PARAM) sinon
*/
static int check_param_optim(tableaux* g) {
	int status =STAT_ERROR_PARAM;

	if(g->resol == RESOL_MAX_FREQ && g->R_cible <1)
		printf("%s Pb données: résolution RESOL_MAX_FREQ mais R_cible (%ld) incohérent\n", __func__, g->R_cible);
	else if(g->resol == RESOL_INT && g->val_cible <=0)
		printf("%s Pb données: résolution RESOL_INT mais val_cible (%Lf) incohérent\n", __func__, g->val_cible);
	else if( g->resol == RESOL_INT_PV_OPT && (g->num_cible <=0 || g->den_cible <=0 || g->num_cible >g->den_cible) )
		printf("%s Pb données: résolution RESOL_INT_PV_OPT mais den_cible (%ld)/num_cible (%ld) incohérents\n", __func__, g->num_cible, g->den_cible);
	else if(g->resol == RESOL_INT_PV_OPT && g->R_cible != NO_VAL)
		printf("%s Pb données: option #lignes min ou cible incompatible avec résolution RESOL_INT_PV_OPT\n", __func__);
	/* arrivé ici: paramètres obligatoires sont (chaque fois le cas échéant) 
						>=1 pour R_cible et >0 pour val_cible, num_cible, den_cible  */
	else if(g->type == TAB_RHO && g->resol == RESOL_INT && g->R_cible>0) {
		if (g->E_q == 0 && mod(g->R_cible, (int)pow(g->q, g->k)) != 0)
			printf("%s Pb données: résolution RESOL_INT rho mais R_cible (%ld) n'est pas un multiple de q^k (%d)\n", __func__, g->R_cible, (int)pow(g->q, g->k));
		else if (g->E_q == 1 && mod(g->R_cible, (int)pow(g->q, g->k -1)) != 0)
			printf("%s Pb données: résolution RESOL_INT rho_E mais R_cible (%ld) n'est pas un multiple de q^(k -1) (%d)\n", __func__, g->R_cible, (int)pow(g->q, g->k -1));
		else /* arrivé ici: TAB_RHO RESOL_INT avec paramètre R_cible >0 vérifiant R =\lambda \times q^k */
			status =0;
	} else if((g->type == TAB_GAMMA || g->type == TAB_DELTA) && g->resol == RESOL_CONT && DO_LB) {
		if(g->num_cible <=0 || g->den_cible <=0 || g->num_cible >g->den_cible)
			printf("%s Pb données: résolution RESOL_CONT gamma mais den_cible (%ld)/num_cible (%ld) incohérents\n", __func__, g->num_cible, g->den_cible);
		else /* arrivé ici: TAB_GAMMA ou TAB_DELTA RESOL_CONT avec paramètres den_cible >= num_cible >0 */
			status =0;
	} else /* arrivé ici: 
			TAB_RHO RESOL_INT sans R_cible, ou TAB_RHO autre résolution, ou (TAB_GAMMA ou TAB_DELTA) toute résolution sauf RESOL_CONT, ou TAB_REG toute résolution
			paramètres obligatoires sont (chaque fois le cas échéant) >=1 pour R_cible et >0 pour val_cible, num_cible, den_cible */
		status =0;

	return status;
}

/* ___________________________________________________________________________________________________________________ Definition du PL
*/

/* Construction environnement Cplex & PL */
int init_cplex(tableaux* g) {
	int status =0;
	char nom_pb[100];
	char errmsg[1024];

	/* ____ Cplex */
	if (status == 0) {
		/* environnement */
		g->env =CPXopenCPLEX (& status);
		if (g->env == NULL)  {
			CPXgeterrorstring (g->env, status, errmsg);
			printf("Pb CPXopenCPLEX, code erreur =%d, message =%s\n", status, errmsg);
		}
		else {
			/* problème */
			sprintf(nom_pb, "DEF");
			if (g->type == TAB_RHO)
				sprintf(nom_pb, "%d_rho_%s_%d_%d_%d", g->resol, (g->E_q == 1 ? "_E": ""), g->nu, g->q, g->k);
			else if (g->type == TAB_GAMMA)
				sprintf(nom_pb, "%d_gamma_%s_%d_%d_%d", g->resol, (g->E_q == 1 ? "_E": ""), g->q, g->p, g->k);
			else if (g->type == TAB_REG)
				sprintf(nom_pb, "%d_reg_%d_%d_%d", g->resol, g->nu, g->p, g->k);
			else /* alors type == TAB_DELTA */
				sprintf(nom_pb, "%d_delta_%d_%d_%d", g->resol, g->nu, g->p, g->k);

			g->pl =CPXcreateprob (g->env, &status, nom_pb);
			if (g->pl == NULL)  {
				printf("Pb CPXcreateprob, code erreur =%d\n", status);
			}
		}
	}

	return status;
}

/* Sens et type optimisation PL */
int finaliser(tableaux* g) {
	int status =0;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Sens optimisation :
		Max Q(0, 1 ,.., q -1)	si RESOL_CONT ou RESOL_MAX_FREQ
		Min R					si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	*/
	CPXchgobjsen(g->env, g->pl, (g->resol == RESOL_CONT || g->resol == RESOL_MAX_FREQ ? CPX_MAX: CPX_MIN ) );

	/* __ Resolution
	*/

	/* Continue */
	if(g->resol == RESOL_CONT)
		status =CPXchgprobtype(g->env, g->pl, CPXPROB_LP);
	/* Entiere */
	else
		status =CPXchgprobtype(g->env, g->pl, CPXPROB_MILP);

	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
	#endif

	return status;
}

/* Borne inferieure nombre de lignes */
int set_R_min(tableaux* g, long int R_min) {
	char borne ='L';
	double lb =(double)R_min; 
	int status;

   	#if (DEBUG)
    printf("%s(%p, %ld) IN\n", __func__, g, R_min);
    #endif
    
	status =CPXchgbds(g->env, g->pl, 1, &(g->ind_R), &borne, &lb);

   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* Borne supérieure nombre de lignes */
int set_R_max(tableaux* g, long int R_max) {
	char borne ='U';
	double lb =(double)R_max; 
	int status;

   	#if (DEBUG)
    printf("%s(%p, %ld) IN\n", __func__, g, R_max);
    #endif
    
	status =CPXchgbds(g->env, g->pl, 1, &(g->ind_R), &borne, &lb);

   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* ___________________________________________________________________________________________________________________ Resolution du PL
*/

/* Resolution PL :
	Pre-condition: PL totalement defini (sinon type et sens optimisation), champ ind_Q_star instancie, tableau ctr_bkI rempli
	Post-conditions: 
		- sens obj et type de problème definis
		- PL resolu
	Affichages: restitution valeur optimale

	TODO Cplex: comprendre necessite de specifier pb continu qd on la deja fait a la construction du PL, et que l'on a ensuite insere des variables continues
*/
int resoudre(tableaux* g) {
	int status =0;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Sens optimisation :
		Max Q(v*)	si RESOL_CONT ou RESOL_MAX_FREQ
		Min R		si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	*/
	CPXchgobjsen(g->env, g->pl, (g->resol == RESOL_CONT || g->resol == RESOL_MAX_FREQ ? CPX_MAX: CPX_MIN ) );

	/* __ Resolution
	*/

	/* Continue */
	if(g->resol == RESOL_CONT) {
		status =CPXchgprobtype(g->env, g->pl, CPXPROB_LP);
		if (status)
			printf("\t******** Pb type de problème.\n");
		else {
			status =CPXlpopt(g->env, g->pl);
			if (status)
				printf("\t******** Pb optimisation.\n");
		}
	}
	/* Entiere */
	else {
		status =CPXchgprobtype(g->env, g->pl, CPXPROB_MILP);
		if (status)
			printf("\t******** Pb type de problème.\n");
		else {
			status =CPXmipopt(g->env, g->pl);
			if (status) 
				printf("\t******** Pb optimisation.\n");
		}
	}

	/* __ Restitution
	*/
	if (status == 0) 
		status =restitution(g);

   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* ____________________________________________________________________________________________ Sauvegardes / affichages
*/

/* _____________ sous-routine formation nom de fichier
*/
static void get_nom_fic(tableaux* g, char nom_pl[100], const char* ext) {
	char tmp[100];

	/* résolution */
	sprintf(nom_pl, "%d-", g->resol);	

	/* type de problème */
	if(g->type == TAB_RHO)
		strcat(nom_pl, "rho");
	else if (g->type == TAB_GAMMA)
		strcat(nom_pl, "gamma");
	else if (g->type == TAB_REG)
		strcat(nom_pl, "reg");
	else /* g->type == TAB_DELTA */
		strcat(nom_pl, "delta");

	/* le cas échéant, E_q (alors nécessairement TAB_RHO ou TAB_GAMMA) */
	if( (g->type == TAB_RHO || g->type == TAB_GAMMA) && g->E_q == 1)
		strcat(nom_pl, "_E");
		
	/* dimensions du problème */
	if(g->type == TAB_RHO)
		sprintf(tmp, "-%d-%d-%d", g->nu, g->q, g->k);
	else if (g->type == TAB_GAMMA) {
		if(g->p_phi <g->q -1)
			sprintf(tmp, "-%d-%d-%d-pmax=%d", g->q, g->p, g->k, g->p_phi);
		else
			sprintf(tmp, "-%d-%d-%d", g->q, g->p, g->k);
	} else { /* g->type == TAB_DELTA ou g->type == TAB_REG */
		if(g->p_phi <g->nu -1)
			sprintf(tmp, "-%d-%d-%d-pmax=%d", g->nu, g->p, g->k, g->p_phi);
		else
			sprintf(tmp, "-%d-%d-%d", g->nu, g->p, g->k);
	}
	strcat(nom_pl, tmp);

	/* éventuel paramètre complémentaire "val_cible" */
	if( (g->R_cible >0) && ( (g->resol == RESOL_R) || (g->resol == RESOL_BIN) || (g->resol == RESOL_MAX_FREQ) || (g->resol == RESOL_INT) ) ) {
		sprintf(tmp, "-R=%ld", g->R_cible);
		strcat(nom_pl, tmp);
	}

	/* extension du fichier */
	sprintf(tmp, ".%s", ext);
	strcat(nom_pl, tmp);
}

/* _____________ PL avant résolution
*/

/* Ecrit le PL au format LP */
int writePL(tableaux* g) {
	int status =0;
	char nom_pl[100];

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	get_nom_fic(g, nom_pl, "lp");

	status =CPXwriteprob(g->env, g->pl, nom_pl, "LP");

   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* _____________ Statut résolution après résolution
*/

/* Affichage synthetique du resultat de la résolution */
int restitution(tableaux* g) {
	int status =0, stat_sol;
	double val_opt, val_R, val_R_star;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ problème resolu */
	if (g->type == TAB_RHO)
		printf("******** rho%s(%d, %d, %d) ", (g->E_q == 1 ? "_E": ""), g->nu, g->q, g->k);
	else if (g->type == TAB_GAMMA)
		printf("******** gamma%s(%d, %d, %d) ", (g->E_q == 1 ? "_E": ""), g->q, g->p, g->k);
	else if (g->type == TAB_REG)
		printf("******** reg(%d, %d, %d) ", g->nu, g->p, g->k);
	else /* g->type == TAB_DELTA */
		printf("******** delta(%d, %d, %d) ", g->nu, g->p, g->k);

	/* __ Résolution */
	if(g->resol == RESOL_CONT)
		printf("résolution CONTINUE (%d)", g->resol);
	else if(g->resol == RESOL_INT) 
		printf("résolution ENTIERE (%d) valeur cible >=%Lf", g->resol, g->val_cible);
	else if(g->resol == RESOL_INT_PV_OPT)
		printf("résolution ENTIERE (%d) valeur cible >%ld/%ld =%Lf", g->resol, g->num_cible, g->den_cible, (long double)g->num_cible/g->den_cible);
	else if(g->resol == RESOL_BIN)
		printf("résolution BINAIRE (%d)", g->resol);
	else if(g->resol == RESOL_MAX_FREQ)
		printf("résolution ENTIERE (%d) nb lignes cible =%ld", g->resol, g->R_cible);
	else if(g->resol == RESOL_R)
		printf("résolution R (%d)", g->resol);
	else
		printf("résolution (%d) non reconnue\n", g->resol);

	/* __ Options supplémentaires */
	if( (g->resol == RESOL_INT || g->resol == RESOL_BIN || g->resol == RESOL_R) && g->R_cible >0)
		printf(" R_cible =%ld", g->R_cible);
	if(g->type == TAB_GAMMA && g->resol == RESOL_CONT && g->num_cible >0 && g->den_cible >0)
		printf(" valeur cible >=%ld/%ld =%Lf", g->num_cible, g->den_cible, (long double)g->num_cible/g->den_cible);
	if( (g->type == TAB_GAMMA || g->type == TAB_DELTA || g->type == TAB_REG) && g->p_phi <(g->type == TAB_GAMMA ? g->q :g->nu -1))
		printf(" p_phi =%d", g->p_phi);
	printf("\n");

	/* __ Solution obtenue */
	stat_sol =CPXgetstat(g->env, g->pl);
	if(stat_sol == CPX_STAT_UNBOUNDED || stat_sol == CPXMIP_UNBOUNDED)
		printf("\t**** problème de valeur infinie.\n");
	else if(stat_sol == CPX_STAT_INFEASIBLE || stat_sol == CPXMIP_INFEASIBLE)
		printf("\t**** problème non réalisable.\n");
	else if(stat_sol != CPX_STAT_OPTIMAL && stat_sol != CPXMIP_OPTIMAL)
		printf("\t**** statut de résolution %d non traite (voir codes Cplex).\n", stat_sol);
	else {
		printf("\t**** solution optimale trouvée ");
		status =CPXgetobjval(g->env, g->pl, &val_opt);
		if(status)
			printf("MAIS Pb récupération val opt\n");
		else{	/* récupération de R^* */
			status =CPXgetx(g->env, g->pl, &val_R_star, g->ind_R_star, g->ind_R_star);
			if(status)
				printf("MAIS Pb récupération val R^*\n");
			else{	/* récupération de R */
				status =CPXgetx(g->env, g->pl, &val_R, g->ind_R, g->ind_R);
				if(status)
					printf("MAIS Pb récupération val R\n");
				else if(g->resol == RESOL_CONT)
					printf("de valeur %lf (inverse %lf)\n", val_opt, 1.0/val_opt);
				else if(g->resol == RESOL_BIN && g->type != TAB_REG)
					printf("de valeur %lf =1/%lf\n", 1.0/val_opt, val_R);
				else if(g->type == TAB_REG)	/* solutions entières Delta régulières */
					printf("de valeur %lf (inverse %lf) =%lf/%lf\n", 2*val_R_star/val_R, val_R/(2*val_R_star), val_R_star, val_R/2);
				else /* solutions entières Rho ou Gamma ou Delta */
					printf("de valeur %lf (inverse %lf) =%lf/%lf\n", val_R_star/val_R, val_R/val_R_star, val_R_star, val_R);
			}
		}
	}

   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* _____________ Solutions après résolution
*/

/* sauvegarde la solution dans un fichier */
int writeSolPL(tableaux* g) {
	int status =0;
	char nom_pl[100];

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif
	
	/* __ on forme le nom du fichier */
	get_nom_fic(g, nom_pl, "sol");

	/* __ on écrit le fichier */
	status =CPXsolwrite(g->env, g->pl, nom_pl);
	
   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* affichage des solutions */
int affSol(tableaux* g) {
	int stat_sol;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	stat_sol =CPXgetstat(g->env, g->pl);
	if(stat_sol != CPX_STAT_OPTIMAL && stat_sol != CPXMIP_OPTIMAL) {
 	  	printf("\tpas de solution courante optimale, on n'affiche pas l'éventuelle solution en cours\n");
		return 0;
	}

	if(g->type == TAB_RHO)
		return rho_affSol(g);

	if(g->type == TAB_GAMMA || g->type == TAB_DELTA)
		return gamma_delta_affSol(g);
	
	if(g->type == TAB_REG)
		return reg_affSol(g);

	return STAT_ERROR_PARAM;

   	#if (DEBUG)
    printf("%s OUT\n", __func__);
    #endif
}
