/* tableaux :: implémentation des fonctions spécifiques aux solutions Rho, Gamma et Delta

	Fonctions :
		int setPL(tableaux* g):
		void libererBkI(tableaux* g);
		int rho_affSol(tableaux* g);
		int gamma_delta_affSol(tableaux* g);
*/

#include "tableaux.h"	/* tableaux */
#include "vn.h" 		// representation des vecteurs et des sous-ensembles par des nombres

#include <stdlib.h>		// pour malloc(), itoa()
#include <string.h>		// pour sprintf()
#include <math.h> 		// pour pow()

#if (DEBUG)
#include <stdio.h>		// pour traces
#endif

/* ________________________________________________________________________ Declaration des fonctions statiques
*/

/* _____________ Accesseurs en lecture
*/

/* Renvoie v* =(0 .. 0) si TAB_RHO, (0, 1 .. q -1) si TAB_GAMMA, (1 .. 1) si TAB_DELTA
*/
static int get_v_star(tableaux* g);

/* Renvoie le plus grand indice de variable, sachant que l'intervalle pour v est: 
	- si TAB_RHO et TAB_GAMMA :: {0 ,.., q^q -1} si E_q =0, {0 ,.., q^(q -1) -1} sinon
		(suppose pour TAB_GAMMA: champ de dimension nu bien initialisé au champ q de taille de l'alphabet)
	- si DELTA :: {0 ,.., 2^nu -1}
*/
static int get_v_max(tableaux* g);

/* Intervalle de travail pour J: {0 ,.., 2^nu -1}
	(suppose pour TAB_GAMMA: champ de dimension vu bien initialisé au champ q de taille de l'alphabet)
*/
static int get_J_max(tableaux* g);

/* Intervalle de travail pour a rhs contrainte bkI:
	- si TAB_RHO et TAB_GAMMA :: {0 ,.., q^k -1} si E_q =0, {0 ,.., q^(k -1) -1} sinon
	- si DELTA :: {0 ,.., 2^k -1}
*/
static int get_a_max(tableaux* g);

/* _____________ Sous-routines définition du PL (hors contraintes bkI)
*/

/* Contraintes base (generation sans les coefficients des variables) :
	- contrainte 0:	 sum_v Q(v) 		    -R	 =0	si TAB_RHO
					 sum_v P(v) 		    -R	 =0	si TAB_GAMMA ou TAB_DELTA 
	- contrainte 1:	Q(v*) 		 -val_cible xR	>=0	si RESOL_INT
				den_cible xQ(v*) -num_cible xR 	>=1	si RESOL_INT_PV_OPT
				den_cible xQ(v*) -num_cible xR 	>=0	si (TAB_GAMMA ou TAB_DELTA) et RESOL_CONT (lower bound 2/(T(q -p +k, k) +1))

	NB TAB_GAMMA : pour BkI comme pour definition de R, sum_v P(v) (+tôt que sum_v Q(v)) car + econome
	NB on utilise les lb sur les variables R et Q(v*) pour representer les contraintes 
		R {=1, =R_cible} (selon le contexte)
		Q(v*) >=1 (le cas échéant)
*/
static int ctr_base(tableaux* g);

/* Variables R et Q(v*) (dont leur coefficient dans certaines des contraintes de base)	
	R est a valeur dans [0, 1] si RESOL_CONT et N sinon. 
	R est impliquée dans les contraintes (hors contraintes BkI) :
	- représentées par lb sur R:	 			 R 	 =1		si RESOL_CONT
									 	  		 R 	 =R_cible	si MAX_FREQ
	- contrainte 0:		 sum_v Q(v) 		    -R	 =0		si TAB_RHO
						 sum_v P(v) 		    -R	 =0		si TAB_GAMMA ou TAB_DELTA
	- contrainte 1:			   Q(v*) -val_cible xR	>=0		si RESOL_INT
					den_cible xQ(v*) -num_cible xR 	>=1		si RESOL_INT_PV_OPT
					den_cible xQ(v*) -num_cible xR 	>=0		si (TAB_GAMMA ou TAB_DELTA) et RESOL_CONT (LB 2/(T(q -p +k, k) +1))
	- contrainte optionnelle :			  	 	 R	 =R_cible	si RESOL_INT, RESOL_INT_PV_OPT, RESOL_R, RESOL_BIN et TAB_RHO

	Q(v*) est spécifiquement impliquée dans les contraintes (soit, hors contraintes définition R et BkI) :
	- représentées par lb sur Q(v*) :
							   Q(v*)				>=1 	si RESOL_BIN, RESOL_INT, RESOL_INT_PV_OPT, RESOL_R
	- contrainte 1:	Q(v*) 			 -val_cible xR	>=0		si RESOL_INT
					den_cible xQ(v*) -num_cible xR 	>=1		si RESOL_INT_PV_OPT
					den_cible xQ(v*) -num_cible xR 	>=0		si (TAB_GAMMA ou TAB_DELTA) et RESOL_CONT (LB 2/(T(q -p +k, k) +1))
*/
static int var_base(tableaux* g);

/* Variables Q(v), v in Z_q^nu (toujours) et (si TAB_GAMMA) P(v), v in U
	Generation des variables:
	- Q(v) pour v in Z_q^nu\{v*}		toujours
	- P(v) pour v in U					si TAB_GAMMA ou TAB_DELTA
	Definition des coefficients dans les contraintes de definition de R et BkI des variables
	- Q(v), v in Z_q^nu					toujours
	- P(v), v in U						si TAB_GAMMA ou TAB_DELTA
	Contraintes Q(v*) -Q(v) >=0 		si TAB_RHO et pas RESOL_BIN
		soit dans les contraintes
		0. sum_v Q(v) -R	=0			si TAB_RHO
		1. sum_v P(v) -R	=0			si TAB_GAMMA ou TAB_DELTA
		2. Q(v*) -Q(v) >=0				si TAB_RHO et pas RESOL_BIN
		3. sum_{v_J =a} Q(v) -R/q^k	=0	si TAB_RHO et E_q =0 (contraintes similaires si E_q =1) 
		4. sum_{v_J =a} (P(v) -Q(v))=0	si TAB_GAMMA et E_q =0 (contraintes similaires si E_q =1) ou TAB_DELTA
*/
static int var_gen(tableaux* g);

/* _____________ Sous-routines contraintes bkI
*/

/* Allocation memoire & initialsation tableau ctr_bkI
	En cas d'erreur : toute la memoire (dont l'environnement Cplex et le PL) est lberee
*/
static int init_bkI(tableaux* g);

/* Structure permettant de contenir les variables avant écriture dans les pl
	Le choix du type 'unsigned long long' permet de repousser au maximum les overflows
*/
typedef struct tableau_st {
	unsigned long long * elements;
	unsigned longueur;
}tableau;

/* Contraintes BkI (generation sans les coefficients des variables, sinon R si TAB_RHO):
		sum_{v_J =a} Q(v) -R/q^k	=0		si TAB_RHO et E_q =0 (contraintes similaires si E_q =1) 
		sum_{v_J =a} (P(v) -Q(v)	=0		si TAB_GAMMA et E_q =0 (contraintes similaires si E_q =1) ou TAB_DELTA
	Post-conditions : tableau ctr_bkI instancié et si TAB_RHO, coefficient de R dans contraintes BkI mis à jour
	Dans ce cas, on écrit comme suit les contraintes :
		q^k 		x sum_{v_J =a} Q(v) 										-R	=0	si TAB_RHO et E_q =0
		q^(k -1) 	x sum_{v_J in {a, a +(1 ,..,1) ,.., (q -1 ,.., q -1)}} Q(v) -R	=0	si TAB_RHO et E_q =1
*/
static int ctr_bkI(tableaux* g);

/* Les deux fonctions suivantes calculent tous les J utilisés dans la fonction ctr_bkI et les stocke dans le tableau *Tab_ni.
	Cette fonction a une complexité polynomiale
	Préconditions : Tab_ni instancié
					n, k >= 1
					(tout cela si seulement enumerer_nombre_incidence est utilisée)
	Post-conditions : Tab_ni->elements alloué dynamiquement
*/
static void generer_nombre_incidence(tableau * Tab_ni, long long unsigned nombre_incidence, unsigned n, unsigned k, unsigned indice_a_choisir, unsigned * indice_tab);	// -> appelée par: generer_nombre_incidence, enumerer_nombre_incidence
static void enumerer_nombre_incidence(tableau * Tab_ni, unsigned n, unsigned k);	// -> appelée par: ctr_bkI, var_gen

/* Alloue dynamiquement le champs elements du tableau *T
	Préconditions : T instancié
	Post-conditions : (T->longueur =n et T->elements alloué) si malloc != NULL
*/
static int allouer_tableau(tableau * T, unsigned n);	// -> appelée par: enumerer_nombre_incidence

/* libère le champs elements du tableau *T
	Préconditions : elements alloué
	Post-conditions : elements libérés et longueur =0
*/
static void liberer_tableau(tableau * T);				// -> appelée par: ctr_bkI, var_gen

/* ________________________________________________________________________ Accesseurs en lecture
*/

/* Renvoie v* =(0 .. 0) si TAB_RHO, (0, 1 .. q -1) si TAB_GAMMA, (1 .. 1) si TAB_DELTA
*/
static int get_v_star(tableaux* g) {
	if(g->type == TAB_RHO)
		return 0;

	if(g->type == TAB_GAMMA)
		return get_v_012etc(g->q);

	/* g->type == TAB_DELTA */
	return get_v_max(g); 
}

/* Renvoie le plus grand indice de variable, sachant que l'intervalle pour v est: 
	- si TAB_RHO et TAB_GAMMA :: {0 ,.., q^q -1} si E_q =0, {0 ,.., q^(q -1) -1} sinon
		(suppose pour TAB_GAMMA: champ de dimension nu bien initialisé au champ q de taille de l'alphabet)
	- si DELTA :: {0 ,.., 2^nu -1}
*/
static int get_v_max(tableaux* g) {
	if(g->type == TAB_DELTA)
		return get_J_max(g);

	/* g->type == TAB_RHO ou g->type == TAB_GAMMA */
	return pow(g->q, g->nu -g->E_q) -1;
}

/* Intervalle de travail pour J: {0 ,.., 2^nu -1}
	(suppose pour TAB_GAMMA: champ de dimension vu bien initialisé au champ q de taille de l'alphabet)
*/
static int get_J_max(tableaux* g) {
	return pow(2, g->nu) -1;
}

/* Intervalle de travail pour a rhs contrainte bkI:
	- si TAB_RHO et TAB_GAMMA :: {0 ,.., q^k -1} si E_q =0, {0 ,.., q^(k -1) -1} sinon
	- si DELTA :: {0 ,.., 2^k -1}
*/
static int get_a_max(tableaux* g) {
	if(g->type == TAB_DELTA)
		return pow(2, g->k) -1;

	/* g->type == TAB_RHO ou g->type == TAB_GAMMA */
	return pow(g->q, g->k -g->E_q) -1;
}

/* ________________________________________________________________________ Définition du PL (hors sous-routines)
*/

/* definition du PL */
int setPL(tableaux* g) {
	int status =0;	

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* Construction environnement Cplex & PL */
	status =init_cplex(g);
	if(status)
		printf("Pb initialisation PL.\n");
 	else {
		/* Allocation mémoire  contraintes bki */
		status =init_bkI(g);		
		if(status)
			printf("Pb allocation mémoire contraintes bkI.\n");
		else {
			/* Contraintes de base (sans les coefficients des variables) */
			status =ctr_base(g);
			if(status)
				printf("Pb generation contraintes base PL.\n");
			else {
				/* Variables R, Q(v*) (avec certains coefficients des variables) */
				status =var_base(g);
				if(status)
					printf("Pb generation variables base PL.\n");
				else {
					/* Contraintes BkI (avec les coefficients de R) */
					status =ctr_bkI(g);
					if(status)
						printf("Pb generation contraintes BkI PL.\n");
					else {
						/* Variables autres que R, Q(v*) et coefficients qui restaient a definir dans les contraintes */
						status =var_gen(g);
						if (status) 
							printf("Pb generation variables generales PL.\n");
						else {
							/* Finalisation PL (sens optimisation, type pb) */
							status =finaliser(g);
							if (status) 
								printf("Pb definition type de probleme.\n");
						}
					}
				}
			}
		}
	}
	
   	#if (DEBUG)
    printf("%s OUT: status =%d :: %d variables et %d contraintes creees\n", __func__, status, CPXgetnumcols(g->env, g->pl), CPXgetnumrows(g->env, g->pl));
    #endif

	return status;
}

/* ________________________________________________________________________ Sous-routines définition du PL (hors contraintes bkI)
*/

/* Contraintes base (generation sans les coefficients des variables) */
static int ctr_base(tableaux* g) {
	int status =0;
	double rhs;
	char sens;
	char* nom_ctr =malloc(100 * sizeof(char));

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Definition de R comme sum_v Q(v) ou sum_v P(v): toujours
		sum_v Q(v) -R	=0		si TAB_RHO
		sum_v P(v) -R	=0		si TAB_GAMMA ou TAB_DELTA
	*/
	sens ='E';
	rhs =0;
	sprintf(nom_ctr, "R_def");
	status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);

	/* __ LB Q(v*)/R au moins (>= ou >) une certaine valeur: si RESOL_INT ou RESOL_INT_PV_OPT
			ou (TAB_GAMMA ou TAB_DELTA) et RESOL_CONT (lower bound 2/(T(q -p +k, k) +1))
		Q(v*) 				-val_cible xR 	>=0		si RESOL_INT
		den_cible xQ(v*)	-num_cible xR 	>=0		si (TAB_GAMMA ou TAB_DELTA) et RESOL_CONT *et DO_LB*
		den_cible  Q(v*) 	-num_cible xR 	>=1		si RESOL_INT_PV_OPT
	*/
	if(status == 0 && ( g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT || 
			(g->resol == RESOL_CONT && (g->type == TAB_GAMMA || g->type == TAB_DELTA) && DO_LB) ) ) {
		sens ='G';
		rhs =(g->resol == RESOL_INT_PV_OPT ? 1 : 0);
		sprintf(nom_ctr, "%s", "val_cible");
		status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
	}

	/* nettoyage */
	free(nom_ctr);

   	#if (DEBUG)
    printf("%s OUT: status =%d, %d contraintes creees\n", __func__, status, CPXgetnumrows(g->env, g->pl));
    #endif

	return status;
}

/* Variables R et Q(v*) (dont leur coefficient dans certaines des contraintes de base) */
static int var_base(tableaux* g) {
	int status =0, v_star =get_v_star(g), ind_ctr; 
	double coef_obj, lb, ub;
	char type;
	char* nom =malloc(100 *sizeof(char));

	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
	#endif

	/* __ generation de R :
		Domaine:
			in {1} 			si RESOL_CONT
			in {R_cible} 	si RESOL_MAX_FREQ 
			in N\{0}		si RESOL_BIN, RESOL_INT, RESOL_INT_PV_OPT
			in {R_cible} 	si R_cible >0 et RESOL_INT, RESOL_R, RESOL_BIN
		R participe a obj (alors coef 1) si RESOL_BIN, RESOL_INT, RESOL_INT_PV_OPT, RESOL_R */
	sprintf(nom, "R");

	/* obj */
	coef_obj =(g->resol == RESOL_CONT || g->resol == RESOL_MAX_FREQ ? 0 : 1);

	/* type */
	type =(g->resol == RESOL_CONT ? CPX_CONTINUOUS : CPX_INTEGER);

	/* bornes */
	if(g->resol == RESOL_CONT)	/* Q (et le cas échéant, P) mesure ssi R =1 */ 
		lb =ub =1;
	else if(g->R_cible >0)		/* R_cible >0: valeur exacte si RESOL_MAX_FREQ (obl) ou RESOL_BIN, RESOL_INT ou RESOL_R (opt) */
		lb =ub =g->R_cible;
	else { 						/* R_cible == NO_VAL: RESOL_INT_PV_OPT notamment (jamais de R_cible) */
		lb =1;
		ub =CPX_INFBOUND;
	}

	/* on génère R et l'on retient son indice dans ind_R */
	g->ind_R =CPXgetnumcols(g->env, g->pl);
	status =CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom);
	if(status != 0)
		goto FIN_FONCTION;

	/* __ generation de Q(v*) :
		in [0, 1] 		si RESOL_CONT
		in {1} 			si RESOL_BIN
		in N^*			si RESOL_INT, RESOL_INT_PV_OPT, RESOL_MAX_FREQ
		>=1				si RESOL_INT, RESOL_BIN, RESOL_INT ou RESOL_INT_PV_OPT
		=> representation par lb des contraintes
			Q(v*) >= 1	si 
		Q(v*) participe a obj (alors coef 1) si RESOL_CONT ou RESOL_MAX_FREQ */
	sprintf(nom, "Q");
	itoa(v_star, g->q, g->nu, nom +1);

	if(g->resol == RESOL_CONT) {
		type =CPX_CONTINUOUS;
		lb =0;
		ub =1;
		coef_obj =1;
	} else if(g->resol == RESOL_MAX_FREQ) {
		type =CPX_INTEGER;
		lb =1;
		ub =CPX_INFBOUND;
		coef_obj =1;
	} else if(g->resol == RESOL_BIN) {
		type =CPX_BINARY;
		lb =ub =1;
		coef_obj =0;
	} else /* RESOL_INT, RESOL_INT_PV_OPT, RESOL_R */ {
		type =CPX_INTEGER;
		lb =1;
		ub =CPX_INFBOUND;
		coef_obj =0;
	}

	/* on génère Q(v*) et l'on retient son indice dans ind_R_star */
	g->ind_R_star =CPXgetnumcols(g->env, g->pl);
	status =CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom);
	if(status != 0)
		goto FIN_FONCTION;

	/* __ participation de R à la contrainte R_def de definition de R:
			sum_v Q(v) 	-R	=0	 si TAB_RHO
			sum_v P(v) 	-R	=0	 si TAB_GAMMA ou TAB_DELTA */
	strcpy(nom, "R_def");	/* il faut récupérer l'indice de la contrainte R_def */
	status =CPXgetrowindex(g->env, g->pl, nom, &ind_ctr);
	if(status == 0) {
		status =CPXchgcoef(g->env, g->pl, ind_ctr, g->ind_R, -1);

		/* __ participation de R et Q(v*) à la contrainte "val_cible" Q(v*)/R au moins (>= ou >) une certaine valeur :
				Q(v*) 				-val_cible xR 	>=0		si RESOL_INT
				den_cible xQ(v*)	-num_cible xR 	>=0		si (TAB_GAMMA ou TAB_DELTA) et RESOL_CONT *et DO_LB*
				den_cible xQ(v*) 	-num_cible xR 	>=1		si RESOL_INT_PV_OPT
		*/
		if( status == 0 && ( g->resol == RESOL_INT || g->resol == RESOL_INT_PV_OPT ||
				(g->resol == RESOL_CONT && (g->type == TAB_GAMMA || g->type == TAB_DELTA) && DO_LB) ) ) {
			strcpy(nom, "val_cible");	/* il faut récupérer l'indice de la contrainte val_cible */
			status =CPXgetrowindex(g->env, g->pl, nom, &ind_ctr);
			if(status == 0){
				/* R: -val_cible si RESOL_INT et -num_cible sinon */
				status =CPXchgcoef(g->env, g->pl, ind_ctr, g->ind_R, -1*(g->resol == RESOL_INT ? g->val_cible : g->num_cible) );
				/* Q(v*): 1 si RESOL_INT et den_cible sinon */
				status =CPXchgcoef(g->env, g->pl, ind_ctr, g->ind_R_star, (g->resol == RESOL_INT ? 1 : g->den_cible) );
			}
		}
	}

FIN_FONCTION:
	/* nettoyage */
	free(nom);

	#if (DEBUG)
    printf("%s OUT: status =%d, %d variables creees\n", __func__, status, CPXgetnumcols(g->env, g->pl));
	#endif

	return status;
}

/* Variables Q(v), v in Z_q^nu (toujours) et P(v), v in U (si TAB_GAMMA)
	Generation des variables:
		- Q(v) pour v in Z_q^nu\{v*}	toujours, sauf si TAB_GAMMA et #val utilisées par v >p_phi
		- P(v) pour v in Z_q^nu			si TAB_GAMMA et #val utilisées par v <=p
	Definition des coefficients dans les contraintes de definition de R et BkI des variables
		- Q(v), v in Z_q^nu				toujours
		- P(v), v in U					si TAB_GAMMA
	Contraintes Q(v*) -Q(v) >=0 		si TAB_RHO et pas RESOL_BIN
		soit dans les contraintes
		0. sum_v Q(v) -R	=0			si TAB_RHO
		1. sum_v P(v) -R	=0			si TAB_GAMMA
		2. Q(v*) -Q(v) >=0				si TAB_RHO et pas RESOL_BIN
		3. sum_{v_J =a} Q(v) -R/q^k	= 0	si TAB_RHO et E_q =0 (contraintes similaires si E_q =1) 
		4. sum_{v_J =a} (P(v) -Q(v)	= 0	si TAB_GAMMA et E_q =0 (contraintes similaires si E_q =1)
*/
static int var_gen(tableaux* g) {
	int status =0, v, J, vJ, ind_var_P, ind_var_Q, ind_ctr, ind_R_def, div;
	int v_star =get_v_star(g), v_max =get_v_max(g), nb_var =CPXgetnumcols(g->env, g->pl);
	double coef_obj =0, lb =0, ub, coef_Q_ctr;
	char type, sens ='G';
	char *nom_var =malloc(100 * sizeof(char)), *nom_ctr =malloc(100 * sizeof(char));
	double rhs =0;

	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
	#endif

	/* on récupère l'indice de la contrainte R_def de définition de R (dont on aura besoin plus tard) */ 	
	strcpy(nom_ctr, "R_def");
	status =CPXgetrowindex(g->env, g->pl, nom_ctr, &ind_R_def);
	if(status) {
		printf("Pb récupération indice contrainte \"R_def\".\n");
		goto FIN_FONCTION;
	}

	/* __ Propriétés des variables générées 
				Q(v), v <>v* (et respectant p_phi si TAB_GAMMA ou TAB_DELTA) 
			et	P(v), v in Z_q^q (et respectant p si TAB_GAMMA ou TAB_DELTA):
		- elles ne participent jamais à la fonction objectif
		- elles ont toutes 0 comme borne inférieure
		- si RESOL_CONT: 				continues in [0, 1]
		- si RESOL_BIN:  				binaires in {0, 1}
		- dans tous les autres cas:		entières in N */
	if(g->resol == RESOL_CONT) {
		type =CPX_CONTINUOUS;
		ub =1;
	} else if(g->resol == RESOL_BIN) {
		type =CPX_BINARY;
		ub =1;
	} else {
		type =CPX_INTEGER;
		ub =CPX_INFBOUND;
	}

	/* coefficient de Q(v) dans contraintes BkI */
	if(g->type == TAB_RHO)	
		coef_Q_ctr =(g->E_q == 0 ? pow(g->q, g->k) : pow(g->q, g->k -1) );
	else /* g->type == TAB_GAMMA ou g->type == TAB_DELTA */
		coef_Q_ctr =-1;

	/* ___ Parcours des vecteurs v */
	for (v =0 ; (v <=v_max) && (status == 0) ; v ++) {
		ind_var_Q =ind_var_P =NO_VAL;
		if(g->type == TAB_GAMMA)		/* nombre de valeurs distinctes prises par le vecteur de Z_q^q associé à v */
			div =get_nb_val(v, g->q, g->q);
		else if(g->type == TAB_DELTA)	/* nombre de composantes non nulles du vecteur de {0, 1}^nu associé à v */
			div =get_card(v, g->nu);

		#if(TRACE_DELTA)
		printf("%s :: v =%d (v_star = %d), poids(v) =%d\n", __func__, v, v_star, div);
		#endif

		/* __ génération de Q(v), P(v) et leur participation à la contrainte définissant R */

		if(v == v_star)	/* variable Q(v^*) déjà générée */
			ind_var_Q =g->ind_R_star;
		/* on ne génère ni Q(v), ni P(v) si 
				TAB_GAMMA et v <>(0 .. q -1) utilise un nombre <nu, >max(p, p_phi) valeurs distinctes, 
			ou 	TAB_DELTA et v de poids <nu, >max(p, p_phi) */
		else if (g->type == TAB_GAMMA && div >g->p_phi && div >g->p && div <g->q)
			continue;
		else if (g->type == TAB_DELTA && div >g->p_phi && div >g->p && div <g->nu)
			continue;
		/* on génère Q(v) si
				TAB_RHO 
			ou 	TAB_GAMMA et v <>(0 .. q -1) utilise un nombre <=p_phi de valeurs distinctes
			ou 	TAB_DELTA et v <>(1 .. 1) a un nombre <=p_phi de composantes non nulles */
		else if(g->type == TAB_RHO || (g->type == TAB_GAMMA && div <=g->p_phi) || (g->type == TAB_DELTA && div <=g->p_phi)) {
			sprintf(nom_var, "Q");
			itoa(v, g->q, g->nu, nom_var +1);
			status =CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom_var);
			ind_var_Q =nb_var;
			nb_var ++;

			#if(TRACE_DELTA)
			printf("\t%s :: variable Q générée\n", __func__);
			#endif
		}

		/* participation de Q(v) à la contrainte définissant R: si TAB_RHO pour tout v, coefficient alors 1 */
		if(status == 0 && g->type == TAB_RHO) {
			status =CPXchgcoef(g->env, g->pl, ind_R_def, ind_var_Q, 1);

			/* participation de Q(v) a la contrainte Q(v*) -Q(v) >=0: si TAB_RHO, v <>v* et resolution non binaire */
			if ((status == 0) && (g->resol != RESOL_BIN) && (v != v_star)) {
				sprintf(nom_ctr, "QstargeqQ%d", v);
				status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom_ctr);
				if(status == 0) {
					ind_ctr =CPXgetnumrows(g->env, g->pl) -1;
					status =CPXchgcoef(g->env, g->pl, ind_ctr, g->ind_R_star, 1);
					if(status == 0)
						status =CPXchgcoef(g->env, g->pl, ind_ctr, ind_var_Q, -1);
				}
			}
		}

		/* generation de P(v) et sa participation alors à la contrainte définissant R:
				si 	TAB_GAMMA et les composantes de v prennent au plus p valeurs distinctes, 
				ou 	TAB_DELTA et v de poids <=p */
		ind_var_P =-1;
		if(status == 0 && ( (g->type == TAB_GAMMA && div <=g->p) || (g->type == TAB_DELTA && div <=g->p) ) ) {
			sprintf(nom_var, "P");
			itoa(v, g->q, g->nu, nom_var +1);
			status =CPXnewcols(g->env, g->pl, 1, &coef_obj, &lb, &ub, &type, &nom_var);
			ind_var_P =nb_var;
			nb_var ++;

			#if(TRACE_DELTA)
			printf("\t%s :: variable P générée\n", __func__);
			#endif

			/* P(v) de coefficient 1 dans la contrainte définissant R */
			if(status == 0)
				status =CPXchgcoef(g->env, g->pl, ind_R_def, ind_var_P, 1);
		}

		/* __ participation de Q(v) (toujours) et P(v) (si TAB_GAMMA ou TAB_DELTA et v in U) 
				aux contraintes BkI d'indices (J =(j_1 ,.., j_k), a =(a_1 ,.., a_t)) vérifiant :
				v_J =a 															si E_q =0 ou TAB_DELTA
				v_J in {a, a +(1 ,.., 1) ,.., a +(q -1 ,.., q -1)} ou a_1 =0 	si E_q =1 et TAB_GAMMA
			- pour rho: 	coefficient q^k ou q^(k -1) pour Q(v) selon E_q
			- pour gamma: 	coefficient 1 pour P(v) si v in U et -1 pour Q(v), quel que soit E_q
			- pour delta: 	coefficient 1 pour P(v) si v in U et -1 pour Q(v)
		*/
		if(status == 0 && (ind_var_P != NO_VAL || ind_var_Q != NO_VAL)) {
			tableau Tab_ni;

			#if(TRACE_DELTA)
			printf("\t%s :: coefficients dans contraintes BkI\n", __func__);
			#endif
	
			/* énumération des sous-ensembles J de Z_q */
			enumerer_nombre_incidence(&Tab_ni, g->nu, g->k);
			for (unsigned i =0; i <Tab_ni.longueur; ++i) {
				J =Tab_ni.elements[i]; 
				/* recuperation de v_J si E_q =0 ou TAB_DELTA, de v_J -(v_{j_1} ,.., v_{j_1}) si E_q =1 et TAB_GAMMA */
				vJ =get_v_J(v, J, g->q, g->nu, g->k, g->E_q);

				/* coefficient -1, q^k ou q^(k -1) pour Q */
				if(ind_var_Q != NO_VAL)
					status =CPXchgcoef(g->env, g->pl, g->ctr_bkI[J][vJ], ind_var_Q, coef_Q_ctr);

				/* coefficient +1 (le cas echeant) pour P */
				if(status == 0 && g->type != TAB_RHO && ind_var_P != NO_VAL)
					status =CPXchgcoef(g->env, g->pl, g->ctr_bkI[J][vJ], ind_var_P, 1);	
			}
			liberer_tableau(&Tab_ni);
		}
	}

FIN_FONCTION:
	/* nettoyage */
	free(nom_var);
	free(nom_ctr);

	#if (DEBUG)
    printf("%s OUT: status =%d, le PL contient %d variables au total\n", __func__, status, CPXgetnumcols(g->env, g->pl));
	#endif

	return status;
}

/* ________________________________________________________________________ Sous-routines contraintes bkI
*/

/* Libération mémoire
*/
void libererBkI(tableaux* g) {
	long int J, J_max =get_J_max(g);

	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
	#endif

	/* Contraintes BkI */
	if (g->ctr_bkI != NULL) { 
		for(J =0 ; J <= J_max ; J ++) {
			if (g->ctr_bkI[J] != NULL)  {
				free(g->ctr_bkI[J]);
			}
		}

		free(g->ctr_bkI);
		g->ctr_bkI =NULL;
	}

   	#if (DEBUG)
    printf("%s OUT \n", __func__);
    #endif
}

/* Allocation memoire & initialsation tableau ctr_bkI
	En cas d'erreur : toute la memoire (dont l'environnement Cplex et le PL) est libérée
*/
static int init_bkI(tableaux* g) {
	int status =0;
	int J, J_max =get_J_max(g), a, a_max =get_a_max(g);

	/* ____ Indices contraintes BkI dans le PL */
	g->ctr_bkI =malloc((J_max +1) * sizeof(int*));
	if (g->ctr_bkI == NULL)  {
		printf("Pb malloc g->ctr_bkI.\n");
		status =STAT_ERROR_MALLOC;
	}
	else {
		for(J =0 ; J <=J_max ; J ++) {
			g->ctr_bkI[J] =malloc((a_max +1) * sizeof(int));
			if (g->ctr_bkI[J] == NULL)  {
				printf("Pb malloc g->ctr_bkI[%d].\n", J);
				status =STAT_ERROR_MALLOC;
			}
			else {
				for(a =0 ; a <=a_max ; a ++) {
					g->ctr_bkI[J][a] =NO_VAL;
				}
			}
		}
	}

	/* ____ Le cas echeant : on fait le menage */
	if(status != 0) {
		liberer(g);
	}

	return status;
}

/* Generation des contraintes BkI (sans les coefficients des variables, sinon R si TAB_RHO) */
static int ctr_bkI(tableaux* g) {
	int status =0, nb_ctr, J, a, a_max =get_a_max(g);
	double rhs;
	char sens;
	char* nom =malloc(100 * sizeof(char));

   	#if (DEBUG)
    printf("%s(%p) IN\n", __func__, g);
    #endif

	/* __ Q ou Q_E in I_q^k / P -Q ou (P -Q)_E in I_q^k :
		E_q =0: 1 contrainte pour tout sous-ensemble J de taille k de [q] et tout a in Z_q^k
		E_q =1: 1 contrainte pour tout sous-ensemble J de taille k de [q] et tout a in {0} x Z_q^(k -1)
	*/
	nb_ctr =CPXgetnumrows(g->env, g->pl);
	sens ='E';
	rhs =0;
	
	tableau Tab_ni;

	/* énumération des sous-ensembles J de Z_q */
	if(status == 0) {
		enumerer_nombre_incidence(&Tab_ni, g->nu, g->k);
		for (unsigned i =0; i < Tab_ni.longueur; ++i) {
			J =Tab_ni.elements[i];
			/* parcours des vecteurs a de Z_q^k si E_q =0, de {0} x Z_q^(k -1) sinon */
			for (a =0 ; (a <= a_max) && (status == 0) ; a ++) {
				/* nouvelle contrainte (contrainte numero nb_ctr +1 ssi contrainte d'indice nb_ctr) */
				sprintf(nom, "BkI_%d_%d", J, a);
				status =CPXnewrows(g->env, g->pl, 1, &rhs, &sens, NULL, &nom);
				if (status == 0) {
					/* mise a jour ctr_bkI[J][a] */
					g->ctr_bkI[J][a] =nb_ctr;

					/* si TAB_RHO: on ajoute le coefficient de R dans la contrainte */
					if(g->type == TAB_RHO)
						status =CPXchgcoef(g->env, g->pl, g->ctr_bkI[J][a], g->ind_R, -1);

					/* contrainte suivante */
					nb_ctr ++;
				}
			}
		}
	}

	/* nettoyage */	
	liberer_tableau(&Tab_ni);
	free(nom);

   	#if (DEBUG)
    printf("%s OUT: status =%d, le PL contient %d contraintes\n", __func__, status, CPXgetnumrows(g->env, g->pl));
    #endif

	return status;
}

/* Les deux fonctions suivantes calculent tous les J utilisés dans la fonction ctr_bkI et les stocke dans le tableau *Tab_ni */
static void generer_nombre_incidence(tableau * Tab_ni, long long unsigned nombre_incidence, unsigned n, unsigned k, unsigned indice_a_choisir, unsigned * indice_tab) {
    if (k == 0) {
        Tab_ni->elements[*indice_tab] =nombre_incidence;
        (*indice_tab)++;
        return;
    }
    if (n -indice_a_choisir == k) {
		nombre_incidence += pow(2, indice_a_choisir)*(pow(2, k) -1);
        Tab_ni->elements[*indice_tab] =nombre_incidence;
        (*indice_tab)++;
        return;
    }
    unsigned long long V0, V1;
	
	V0 =nombre_incidence;
	V1 =nombre_incidence + pow(2, indice_a_choisir);
	
    generer_nombre_incidence(Tab_ni, V0, n, k, indice_a_choisir +1, indice_tab);
    generer_nombre_incidence(Tab_ni, V1, n, k -1, indice_a_choisir +1, indice_tab);
}

static void enumerer_nombre_incidence(tableau * Tab_ni, unsigned n, unsigned k) {
	if (allouer_tableau(Tab_ni, binom(n, k)) == -1) {
		fprintf(stderr, "Le tableau n'a pas pu être alloué.\n");
		exit(-1);
	}
	unsigned long long nombre_incidence =0;
	unsigned indice_tab =0;
	generer_nombre_incidence(Tab_ni, nombre_incidence, n, k, 0, &indice_tab);
}

/* Alloue dynamiquement le champs elements du tableau *T */
static int allouer_tableau(tableau * T, unsigned n) {
    // Le tableau est supposé déjà alloué: on n'en a pas besoin de manière dynamique
    if ((T->elements =(unsigned long long *)(malloc(n * sizeof(unsigned long long)))) == NULL)
        return -1;  // Erreur d'allocation
        
    T->longueur =n;

    return 0;   // Tout se passe bien
}

/* libère le champs elements du tableau *T */
static void liberer_tableau(tableau * T) {
    T->longueur =0;
    free(T->elements);
}

/* ________________________________________________________________________ Sauvegardes / affichages
*/

/* _____________ Solutions après résolution (selon le type de problème)
*/

/* Affichage rho */
int rho_affSol(tableaux* g) {
	int status =0, ind_var, v, v_max =get_v_max(g);
	double valQ;
	char nom_var[100];

	/* __ Parcours des vecteurs v in (Z_q)^nu ou {0} x Z_q^(nu -1) selon E_q */
	for (v =0 ; (v <= v_max) && (status == 0) ; v ++) {
		/* recuperation indice +valeur Q(v) */
		sprintf(nom_var, "Q");
		itoa(v, g->q, g->nu, nom_var +1);
		status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
		if(status == 0)
			status =CPXgetx(g->env, g->pl, &valQ, ind_var, ind_var);

		/* affichage si Q(v) >0 */
		if (status == 0 && valQ >0) {
			affiche_v(v, g->q, g->nu);
			printf("\t%lf", valQ);
			printf("%s\n", (ind_var == g->ind_R_star ? " (*)" : ""));
		}
	}

	return status;
}

/* Affichage gamma ou delta (non supposé régulier) */
int gamma_delta_affSol(tableaux* g) {
	int status =0, ind_var, v, v_max =get_v_max(g), v_star =get_v_star(g), div;
	double valP, valQ;
	char nom_var[100];

	/* __ Parcours des vecteurs v in (Z_q)^q ou {0} x Z_q^(q -1) selon E_q */
	for (v =0 ; (v <=v_max) && (status == 0) ; v ++) {
		/* une ligne v n'est considérée que si elle est 0 .. q -1 (v autorisée dans Phi), 
				ou son nombre de composantes est <=p (v autorisée dans Psi) ou <=p_phi (v autorisée dans Phi) */
		if(g->type == TAB_GAMMA)
			div =get_nb_val(v, g->q, g->q);	/* nb valeurs distinctes proses par les composantes de v in Z_q^q */
		else /* g->type == TAB_DELTA */
			div =get_card(v, g->nu);		/* poids de v in {0, 1}^nu */

		if (v != v_star && div >g->p_phi && div >g->p)
			continue;

		valP =valQ =0;
	
		/* recuperation indice +valeur P(v) si les composantes de v prennent au plus p valeurs */
		if (div <=g->p) {
			sprintf(nom_var, "P");
			itoa(v, g->q, g->nu, nom_var +1);
			status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
			if(status == 0)
				status =CPXgetx(g->env, g->pl, &valP, ind_var, ind_var);
			else
				printf("\tproblème récupération index variable %s\n", nom_var);
		}

		/* recuperation indice +valeur Q(v) si les composantes de v prennent soit q, soit au plus p_phi valeurs */
		if(v == v_star)
			status =CPXgetx(g->env, g->pl, &valQ, g->ind_R_star, g->ind_R_star);
		else if (div <=g->p_phi) {
			sprintf(nom_var, "Q");
			itoa(v, g->q, g->nu, nom_var +1);
			status =CPXgetcolindex(g->env, g->pl, nom_var, &ind_var);
			if(status == 0)
				status =CPXgetx(g->env, g->pl, &valQ, ind_var, ind_var);
			else
				printf("\tproblème récupération index variable %s\n", nom_var);
		}

		/* Affichage si P(v) >0 ou Q(v) >0 */
		if (status == 0 && (valP +valQ >0)) {
			affiche_v(v, g->q, g->nu);
			printf("\tPsi: %lf\tPhi: %lf", valP, valQ);
			printf("%s\n", (v == get_v_star(g) ? " (*)" : ""));
		}
	}

	return status;
}
