/* tableaux :: implémentation des fonctions spécifiques aux solutions Delta régulières

	Fonctions :
		int reg_setPL(tableaux* g);
		int reg_affSol(tableaux* g);
*/

#include "tableaux.h" 	// (nu, q, p, k, E_q)
#include "vn.h" 		// representation des vecteurs et des sous-ensembles par des nombres

#include <stdlib.h>		// pour malloc(), itoa()
#include <string.h>		// pour sprintf()
#include <math.h> 		// pour pow()

#if (DEBUG)
#include <stdio.h>		// pour traces
#endif

/* ____________________________________________________________________________________________ Declaration des fonctions statiques
*/

/* _____________ Sous-routines définition du PL */
static int reg_generer_var(tableaux* g);
static int reg_generer_ctr(tableaux* g);
static int reg_populer_ctr(tableaux* g);
static int reg_finaliser(tableaux* g);

/* ____________________________________________________________________________________________ Définition du PL (hors sous-routines)
*/

int reg_setPL(tableaux* g) {
	int status;	

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	status =init_cplex(g);
	if(status)
		printf("Pb initialisation PL.\n");
 	else {
		status =reg_generer_var(g);
		if(status)
			printf("Pb reg_generer_var\n");
		else {
			status =reg_generer_ctr(g);
			if(status)
				printf("Pb reg_generer_ctr\n");
			else {
				status =reg_populer_ctr(g);
				if(status)
					printf("Pb reg_populer_ctr\n");
				else {
					status =reg_finaliser(g);
					if(status)
						printf("Pb reg_finaliser_pl\n");
				}
			}
		}
	}

   	#if (DEBUG)
    printf("%s OUT : status = %d :: %d variables et %d contraintes creees\n", __func__, status, CPXgetnumcols(g->env, g->pl), CPXgetnumrows(g->env, g->pl));
    #endif

	return status;
}

/* ____________________________________________________________________________________________ Sous-routines définition du PL
*/

/* "Finalise" le problème en fonction du type de résolution */
static int reg_finaliser(tableaux* g) {
	int status =0;	

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ sens objectif & type de problème */
	if(g->resol == RESOL_CONT){
		CPXchgobjsen(g->env, g->pl, CPX_MAX);
		status =CPXchgprobtype(g->env, g->pl, CPXPROB_LP);
	} else { /* g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT */
		CPXchgobjsen(g->env, g->pl, CPX_MIN);
		status =CPXchgprobtype(g->env, g->pl, CPXPROB_MILP);
	}

   	#if (DEBUG)
    printf("%s OUT: status =%d\n", __func__, status);
    #endif

	return status;
}

/* Génération des variables
	(NB contraintes R =1 et z_i <=0, i in {p +1 .. nu} données par les domaines de définition des variables)
			
	Domaine des variables:
	- R: 						{1} 	si RESOL_CONT,		N\{0, 1} si RESOL_INT ou RESOL_INT_PV_OPT
	- R_nu, i in {0 .. nu}: 	[0, 1]	si RESOL_CONT,		N\{0} si RESOL_INT ou RESOL_INT_PV_OPT
	- R_i, i in {0 .. nu}: 		[0, 1]	si RESOL_CONT,		N si RESOL_INT ou RESOL_INT_PV_OPT
	- z_nu, i in {p .. nu}: 	[-1, 0]	si RESOL_CONT,		-N\{0} si RESOL_INT ou RESOL_INT_PV_OPT
	- z_i, i in {0 .. p}: 		[-1, 1]	si RESOL_CONT,		Z si RESOL_INT ou RESOL_INT_PV_OPT
	- z_i, i in {p .. nu}: 		[-1, 0]	si RESOL_CONT,		-N si RESOL_INT ou RESOL_INT_PV_OPT

	Coefficients dans fonction objectif:
	2 pour R_vu	(et tous les autres coefficients nuls) si RESOL_CONT
	1 pour R (et tous les autres coefficients nuls) si RESOL_INT ou RESOL_INT_PV_OPT
*/
static int reg_generer_var(tableaux* g) {
	int status =0;	
	int j;
	double lb, ub, coef_obj;
	char type =CPX_CONTINUOUS;	/* toutes les variables sont continues */
	char* nom_var =malloc(100 * sizeof(char));

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ generation de R:
		Variable réelle de valeur 1 si RESOL_CONT, à valeur dans N\{0, 1} si RESOL_INT ou RESOL_INT_PV_OPT
		De coefficient 1 dans la fonction objectif si RESOL_CONT (0 sinon)
	*/
	if(g->resol == RESOL_CONT) {
		type =CPX_CONTINUOUS;
		lb =ub =1;
		coef_obj =0;
	} else { /* cas (g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT) */
		type =CPX_INTEGER;
		lb =2;
		ub =CPX_INFBOUND;
		coef_obj =1;
	}
	/* on génère la variable R et retient son indice */
	sprintf(nom_var, "R");
	g->ind_R =CPXgetnumcols(g->env, g->pl);
	status =CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom_var);
	if(status)
		goto SORTIE_FONCTION;

	/* __ generation de R_j, j in {0 .. nu} :
		Variables réelles sur [0, 1] si RESOL_CONT, entières >=0 si RESOL_INT ou RESOL_INT_PV_OPT
		Ne participent pas à la fonction objectif, 
			à l'exception de R_nu qui y participe avec un coefficient 2 *si* RESOL_CONT
	 	Domaine:						RESOL_CONT		RESOL_INT ou RESOL_INT_PV_OPT
			si j =nu:					0 <=R_j <=1		1 <=R_j <=infty 
			si nu >j >max(p, p_phi):		R_j  =0			R_j  =0
			si j <= max(p, p_phi):		0 <=R_j <=1		0 <=R_j <=infty
	*/
	type =(g->resol == CPX_CONTINUOUS ? 0 : CPX_INTEGER);
	for(j =g->nu ; j >=0 && (status == 0); j --) {
		coef_obj =( (j == g->nu && g->resol == RESOL_CONT) ? 2 : 0);
		lb =( (j == g->nu && g->resol != RESOL_CONT) ? 1 : 0);
		if(j <g->nu && j >g->p && j >g->p_phi)
			ub =0;
		else
			ub =(g->resol == RESOL_CONT ? 1 : CPX_INFBOUND);
		sprintf(nom_var, "R_%d", j);
		if(j == g->nu)	/* on retient l'indice de la variable R_nu */
			g->ind_R_star=CPXgetnumcols(g->env, g->pl);
		status =CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom_var);
	}	
	if(status)
		goto SORTIE_FONCTION;

	/* __ generation de z_j, j in {0 .. nu}:
		Variables réelles sur [-1, 1] si RESOL_CONT, entières si RESOL_INT ou RESOL_INT_PV_OPT
		Ne participent jamais à la fonction objectif
	 	Domaine:						RESOL_CONT		RESOL_INT ou RESOL_INT_PV_OPT
			si j =nu:					-1 <=z_j <=0	-infty <=z_j <=-1 
			si nu >j >max(p, p_phi):		 z_j  =0			 z_j  =0
			si p_phi <j <=p:			 0 <=z_j <=1		 0 <=z_j <=infty
			si p <j <=p_phi:			-1 <=z_j <=0	-infty <=z_j <=0
			si nu >min(p, p_phi) >=j:	-1 <=z_j <=1	-infty <=z_j <=infty
	*/
	type =(g->resol == RESOL_CONT ? CPX_CONTINUOUS : CPX_INTEGER);
	coef_obj =0;
	for(j =g->nu ; j >=0 && (status == 0); j --) {
		if(j == g->nu){
			lb =(g->resol == RESOL_CONT ? -1 : -CPX_INFBOUND);
			ub =(g->resol == RESOL_CONT ? 0 : -1);
		} else if (j >g->p_phi && j >g->p)
			ub =lb =0;
		else if (j >g->p_phi && j <=g->p){
			lb =0;
			ub =(g->resol == RESOL_CONT ? 1 : CPX_INFBOUND);
		} else if (j >g->p && j <=g->p_phi){
			lb =(g->resol == RESOL_CONT ? -1 : -CPX_INFBOUND);
			ub =0;
		} else {	/* j <=min (p, p_phi) */
			lb =(g->resol == RESOL_CONT ? -1 : -CPX_INFBOUND);
			ub =(g->resol == RESOL_CONT ? 1 : CPX_INFBOUND);
		}
		sprintf(nom_var, "z_%d", j);
		status =CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom_var);
	}
	if(status)
		goto SORTIE_FONCTION;

SORTIE_FONCTION:
	/* nettoyage */
	free(nom_var);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Génération des contraintes 
	(NB les contraintes 
				R =1 si RESOL_CONT,
			et	z_i <=0, i in {p +1 .. nu} pour tous les types de résolution
		sont représentées par les domaines de définition des variables)

	NB les contraintes ici sont indépendantes du type de résolution considéré
*/
static int reg_generer_ctr(tableaux* g) {
	int status =0, j;	
	double rhs;	/* toutes les contraintes ont 0 pour rhs, *sauf pê val_cible* si RESOL_INT_PV_OPT */
	char sens;
	char* nom_ctr =malloc(100 * sizeof(char));

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* -- Definition de la variable dérivée R
		comme sum_{i =0}^nu R_i:	sum_{i =0}^nu R_i -R	=0 
	*/
	rhs =0;
	sens ='E';
	sprintf(nom_ctr, "R_def");
	status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
	if(status)
		goto SORTIE_FONCTION;

	/* __ Contrainte additionnelle liant R et R_nu si RESOL_INT ou RESOL_INT_PV_OPT
		2R_nu 				-val_cible xR 	>=0		si RESOL_INT
		2den_cible xR_nu 	-num_cible xR 	>=1		si RESOL_INT_PV_OPT
	*/
	if(g->resol != RESOL_CONT) {
		sens ='G';
		rhs =(g->resol == RESOL_INT_PV_OPT ? 1 : 0);
		sprintf(nom_ctr, "%s", "val_cible");
		status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
	}

	/* __ Definition des variables dérivées z_i, i in {nu} U {0 .. k -1} selon les z_j, j in {k .. nu -1} :
			z_nu: 					sum_{j =k}^{nu -1} C_{nu -k}^{j -k} z_j +z_nu =0
			z_i, i in {0 .. k -1}: 	(-1)^{k -i} sum_{j =k}^{nu -1} C_{nu -1 -j}^{j -i} C_{j -i -1}^{k -i -1} z_j -z_i =0
	*/
	rhs =0;
	sens ='E';
	sprintf(nom_ctr, "z_%d_def", g->nu);
	status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
	for(j =0 ; j < g->k && (status == 0); j ++) {
		sprintf(nom_ctr, "z_%d_def", j);
		status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
	}
	if(status)
		goto SORTIE_FONCTION;

	/* __ Definition des variables dérivées 
		R_j comme -C_nu^j z_j pour j in {p +1 .. nu}:	z_j +R_j =0
	*/
	for(j =g->p +1 ; j <= g->nu && (status == 0); j ++) {
		sprintf(nom_ctr, "R_%d_def", j);
		status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
	}
	if(status)
		goto SORTIE_FONCTION;

	/* __ Bornes inférieures des variables dérivées 
		R_j, j in {0 .. p}:	R_j [+-] C_nu^j z_j >=0, j in {0 .. p}
	*/
	rhs =0;	
	sens ='G';
	for(j =0 ; j <= g->p && (status == 0); j ++) {
		sprintf(nom_ctr, "R_%d_dm", j);
		status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
		if(status == 0) {
			sprintf(nom_ctr, "R_%d_dp", j);
			status = CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
		}
	}

SORTIE_FONCTION:
	/* nettoyage */
	free(nom_ctr);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* Population de la matrice des contraintes
	(NB les contraintes 
				R =1 si RESOL_CONT,
			et	z_i <=0, i in {p +1 .. nu} pour tous les types de résolution
		sont représentées par les domaines de définition des variables)

	NB les contraintes ici sont indépendantes du type de résolution considéré
*/
static int reg_populer_ctr(tableaux* g) {
	int status = 0;	
	char* nom_ctr = malloc(100 * sizeof(char));
	char* nom_var = malloc(100 * sizeof(char));
	int j, i, ind_var, ind_ctr;
	int puiss;
	double coef;

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Definition de R comme sum_{i =0}^nu R_i:	sum_{i =0}^nu R_i -R	=0 */

	/* on récupère l'indice de la contrainte dans le PL */
	sprintf(nom_ctr, "R_def");
	status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_ctr);
	if(status == 0) {
		/* traitement de la variable R (coefficient -1) */
		status =CPXchgcoef(g->env, g->pl, ind_ctr, g->ind_R, -1);
		/* traitement des variables R_0 .. R_nu (coefficient 1) */
		for(j =0 ; j <=g->nu && (status == 0); j ++) {
			sprintf(nom_var, "R_%d", j);
			status = CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
			if(status == 0)
				status = CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, 1);
		}
	}
	if(status)
		goto SORTIE_FONCTION;

	/* __ Definition de z_nu:	sum_{j =k}^{nu -1} C_{nu -k}^{j -k} z_j +z_nu =0 */

	/* on récupère l'indice de la contrainte dans le PL */
	sprintf(nom_ctr, "z_%d_def", g->nu);
	status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_ctr);
	if(status == 0) {
		/* traitement de la variable z_nu (coefficient 1) */
		sprintf(nom_var, "z_%d", g->nu);
		status = CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
		if(status == 0) {
			status= CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, 1);
			/* traitement des variables z_k .. z_{nu -1} (coefficient C_{nu -k}^{j -k}) */
			for(j =g->k ; j <g->nu && (status == 0); j ++) {
				coef =(double)binom(g->nu -g->k, j -g->k);
				sprintf(nom_var, "z_%d", j);
				status= CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
				if(status == 0)
					status= CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, coef);
			}
		}
	}
	if(status)
		goto SORTIE_FONCTION;

	/* __ Definition des z_i, i in {0 .. k -1} :
		(-1)^{k -i} sum_{j =k}^{nu -1} C_{nu -1 -i}^{j -i} C_{j -i -1}^{k -i -1} z_j -z_i =0 */

	/* traitement des contraintes z_0_def .. z_{k -1}_def */
	for(i =0, puiss =(g->k % 2 == 0 ? 1 : -1) ; i <g->k && (status == 0); i ++, puiss =-puiss) {
		sprintf(nom_ctr, "z_%d_def", i);
		status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_ctr);
		if(status == 0) {
			/* traitement de la variable z_i (coefficient -1) */
			sprintf(nom_var, "z_%d", i);
			status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
			if(status == 0) {
				status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, -1);
				/* traitement des variables z_k .. z_{nu -1} (coefficient (-1)^{k -i} C_{nu -1 -i}^{j -i}C_{j -1 -i}^{k -1 -i}) */
				for(j =g->k ; j <g->nu && (status == 0); j ++) {
					sprintf(nom_var, "z_%d", j);
					status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
					if(status == 0)
						status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, puiss * binom(g->nu -1 -i, j -i) * binom(j -1 -i, g->k -1 -i));
				}
			}
		}
	}
	if(status)
		goto SORTIE_FONCTION;

	/* __ Definition de R_j comme -C_nu^j z_j pour j in {p +1 .. nu}:	C_nu^j z_j +R_j =0 */

	/* traitement des contraintes R_{p +1}_def .. R_nu_def */
	for(j =g->p +1 ; j <= g->nu && (status == 0); j ++) {
		sprintf(nom_ctr, "R_%d_def", j);
		status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_ctr);
		if(status == 0) {
			/* traitement de la variable R_j (coefficient 1) */
			sprintf(nom_var, "R_%d", j);
			status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
			if(status == 0) {
				status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, 1);
				if(status == 0) {
					/* traitement de la variable z_j (coefficient C_nu^j) */
					sprintf(nom_var, "z_%d", j);
					status = CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
					if(status == 0)
						status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, binom(g->nu, j));
				}
			}
		}
	}
	if(status)
		goto SORTIE_FONCTION;

	/* __ Bornes inférieures des variables dérivées R_j, j in {0 .. p} : R_j -C_nu^j z_j >=0, j in {0 .. p} */
	
	/* traitement des contraintes R_0_dp .. R_p_dp */
	for(j =0 ; j <= g->p && (status == 0); j ++) {
		sprintf(nom_ctr, "R_%d_dp", j);
		status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_ctr);
		if(status == 0) {
			/* traitement de la variable R_j (coefficient 1) */
			sprintf(nom_var, "R_%d", j);
			status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
			if(status == 0) {
				status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, 1);
				if(status == 0) {
					/* traitement de la variable z_j (coefficient -C_nu^j) */
					sprintf(nom_var, "z_%d", j);
					status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
					if(status == 0)
						status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, -binom(g->nu, j));
				}
			}
		}
	}
	if(status)
		goto SORTIE_FONCTION;

	/* traitement des contraintes R_0_dm .. R_p_dm */
	for(j =0 ; j <= g->p && (status == 0); j ++) {
		sprintf(nom_ctr, "R_%d_dm", j);
		status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_ctr);
		if(status == 0) {
			/* traitement de la variable R_j (coefficient 1) */
			sprintf(nom_var, "R_%d", j);
			status = CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
			if(status == 0) {
				status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, 1);
				if(status == 0) {
					/* traitement de la variable z_j (coefficient C_nu^j) */
					sprintf(nom_var, "z_%d", j);
					status = CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
					if(status == 0)
						status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var, binom(g->nu, j));
				}
			}
		}
	}
	if(status)
		goto SORTIE_FONCTION;
	
	/* __ Contrainte additionnelle liant R et R_nu si RESOL_INT ou RESOL_INT_PV_OPT
		2R_nu 				-val_cible xR 	>=0		si RESOL_INT
		2den_cible xR_nu 	-num_cible xR 	>=1		si RESOL_INT_PV_OPT
	*/
	if(status == 0 && g->resol != RESOL_CONT) {
		strcpy(nom_ctr, "val_cible");	/* on récupère l'indice de la contrainte val_cible */
		status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_ctr);
		if(status == 0){
			/* R: -val_cible si RESOL_INT et -num_cible sinon */
			status =CPXchgcoef(g->env, g->pl, ind_ctr, g->ind_R, -1*(g->resol == RESOL_INT ? g->val_cible : g->num_cible) );
			if(status == 0)
				/* R_nu: 2 si RESOL_INT et 2*den_cible si RESOL_INT_PV_OPT */
				status =CPXchgcoef(g->env, g->pl, ind_ctr, g->ind_R_star, 2*(g->resol == RESOL_INT ? 1 : g->den_cible) );
		}
	}

SORTIE_FONCTION:
	/* nettoyage */
	free(nom_var);
	free(nom_ctr);

   	#if (DEBUG)
    printf("%s OUT : status = %d\n", __func__, status);
    #endif

	return status;
}

/* ____________________________________________________________________________________________ Sauvegardes / affichages
*/
int reg_affSol(tableaux* g){
	int status =0, ind_var, j;
	double val_z, val_R, val_conj =0.0;
	char nom_var[100];

	/* __ Valeur conjecturée */
	if (g->k == 2) {
		if(g->p % 2 == 0)
			val_conj =(g->p/2.0 * g->p/2)/((g->nu -g->p/2) * (g->nu -g->p/2));
		else
			val_conj =((g->p +1)/2.0 * (g->p -1)/2)/((g->nu -(g->p +1)/2) * (g->nu -(g->p -1)/2));
	} else if (g->p == g->k)
		val_conj =( g->nu == g->k ? 1.0 : 2.0/(Tqk(g->nu, g->k) +1) );
	
	if (val_conj && g->p_phi == g->nu -1)
		printf("Borne inf. %lf\n", val_conj);

	/* __ Parcours des composantes z_nu .. z_0 et R_nu .. R_0 */
	for (j =g->nu ; j >=0 && (status == 0) ; j --) {
		if(j <g->nu && j >g->p && j >g->p_phi)	/* variables toujours nulles */
			continue;
		/* z_j */
		sprintf(nom_var, "z_%d", j);
		status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
		if(status == 0) {
			status =CPXgetx(g->env, g->pl, &val_z, ind_var, ind_var);
			if(status == 0) {
				/* R_j */
				sprintf(nom_var, "R_%d", j);
				status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
				if(status == 0) {
					status =CPXgetx(g->env, g->pl, &val_R, ind_var, ind_var);
					if(status == 0) {
						/* affichage z_j et R_j */
						if(val_z != 0 || val_R != 0) {
							/* reconnaissance de variables particulières (k, p) */
							if(j == g->k && j == g->p)
								printf("\t%s", "_kp_ ");
							else if(j == g->k)
								printf("\t%s", "_k__ ");
							else if(j == g->p)
								printf("\t%s", "_p__ ");
							else if(j == g->nu)
								printf("\t%s", "     ");
							else
								printf("\t%s", "     ");
							/* classe de la variable (K, H, L pour j <k, P, N, O pour j >= k) */
							if(j == g->nu)
								printf("%s pour j =%d:\t", "***", j);
							else if(j >= g->k) {
								if(val_z >0)
									printf("%s pour j =%d:\t", "_P_", j);
								else if(val_z <0)
									printf("%s pour j =%d:\t", "-N-", j);
								else /* val_z == 0 */
									printf("%s pour j =%d:\t", "-O-", j);
							}
							else {/* j < k */
								if(val_z == 0)
									printf("%s pour j =%d:\t", "_L_", j);
								else if( (j%2 == g->k%2 && val_z >0) || (j%2 != g->k%2 && val_z <0) )
									printf("%s pour j =%d:\t", "_K_", j);
								else	/* (j%2 == g->k%2 && val_z <0) || (j%2 != g->k%2 && val_z >0) */
									printf("%s pour j =%d:\t", "_H_", j);
							}
							/* valeur des variables */
							if(val_z >=0)
								printf("z =%+lf\t R =%+lf\n", val_z, val_R);
							else
								printf("z =%+lf\t-R =%+lf\n", val_z, -val_R);
						}
					}
				}
			}
		}
	}

	return status;
}
