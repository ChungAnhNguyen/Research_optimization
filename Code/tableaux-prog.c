#include <stdio.h>		// pour affichages
#include <stdlib.h>		// pour EXIT_SUCCESS
#include "tableaux.h"	/* tableaux */
#include "vn.h" 		// Tqk

// Modification 19/07/2021 : prendre en argument de la fonction main tous les paramètres pour les résolutions
/* Fonction principale qui prend en argument les entrées utilisateur dans l'ordre suivant :
		- type: 
			TAB_RHO normal 		0
			TAB_RHO E_q 		1
			TAB_GAMMA normal 	2
			TAB_GAMMA E_q	 	3
			TAB_REG 			4
			TAB_DELTA 			5
		- resolution: 
			RESOL_CONT 			0
			RESOL_INT			1
			RESOL_INT_PV_OPT	2 
			RESOL_BIN			3
			RESOL_MAX_FREQ		4 
			RESOL_R				5
	puis
		- nu si TAB_RHO		q si TAB_GAMMA 		nu si TAB_DELTA ou TAB_REG
		- q  si TAB_RHO		p si TAB_GAMMA		p  si TAB_DELTA ou TAB_REG
		- k
puis (paramètres obligatoires optimisation):
	Si RESOL_INT
		- val_cible: valeur ciblée optimale *sous la forme val ou num/den*
	Si RESOL_INT_PV_OPT
		- num_cible: numérateur cible
		- den_cible: dénominateur cible
	Si RESOL_MAX_FREQ
		- R_cible: nombre de lignes ciblé
puis (paramètres optionnels optimisation):
	Si RESOL_INT ou RESOL_BIN ou RESOL_R
		- R_cible: nombre de lignes ciblé
	Si TAB_GAMMA ou TAB_DELTA ou TAB_REG
		- p_phi (*requiert* pour RESOL_INT ou RESOL_BIN ou RESOL_R de *spécifier R_cible* (à 0 si non exploité))
*/

/* Définition des constantes sur les arguments */
#define ARGCMIN_CONT_BIN_MINR 6
#define ARGCMIN_INT_MAXFREQ 7
#define ARGCMIN_PV_OPT 8
// Indices des arguments
#define INDICE_TYPE 1
#define INDICE_RESOL 2
#define INDICE_DIM 3
#define INDICE_Q_RHO 4
#define INDICE_P 4
#define INDICE_K 5
#define INDICE_VAL_CIBLE 6
#define INDICE_NUM_CIBLE 6
#define INDICE_DEN_CIBLE 7
#define INDICE_R_CIBLE 6
#define IND_ARGOPT_R 0
#define IND_ARGOPT_PMAX_CONT_PVOPT_FREQ 0
#define IND_ARGOPT_PMAX_INT_BIN_MINR 1

/* lecture des données du problème ; en cas de succès, valeurs:
	type, nu, k 	toujours renseignées
	E_q, q			renseignée si type TAB_RHO ou TAB_GAMMA
	p				renseignée si type TAB_GAMMA ou TAB_DELTA ou TAB_REG
*/
int get_probleme(int argc, char * argv[], enum e_type* type, int* E_q, int* nu, int* q, int* p, int* k);

/* lecture du type et des paramètres (obligatoires) d'optimisation ; en cas de succès, valeurs:
	resol					toujours renseignée
	val_cible 				renseignée si RESOL_INT
	num_cible, den_cible	renseignées si RESOL_INT_PV_OPT
	R_cible					renseignée si RESOL_MAX_FREQ
*/
int get_optim(int argc, char * argv[], enum e_resolution* resol, long double* val_cible, long int* num_cible, long int* den_cible, long int* R_cible);

/* __ lecture des paramètres optionnels d'optimisation */

/* Nombre de lignes (selon le contexte, min ou max ou à atteindre):
	- possible 1er argument optionnel pour RESOL_INT, RESOL_BIN, RESOL_R
	(de valeur 0 si une valeur n'est mise que pour pouvoir spécifier d'autres paramètres optionnels à la suite)
*/
void get_Rcible(int argc, char * argv[], enum e_resolution resol, long int* R_cible);	/* pour type de résolution <> CONT, PV_OPT */

/* Pour TAB_GAMMA, TAB_DELTA et TAB_REG :: éventuelle limitation p_phi sur les mots de Phi 
	(outre 0 1 .. (q -1) si GAMMA, 1 1  .. 1 si DELTA) de Phi:
	- possible 2è argument optionnel pour RESOL_INT, RESOL_BIN, RESOL_R
	- possible 1er argument optionnel pour RESOL_CONT, RESOL_INT_PV_OPT, RESOL_BIN, RESOL_MAX_FREQ
*/
void get_p_phi(int argc, char * argv[], enum e_resolution resol, int* p_phi);	/* suppose type de problème GAMMA ou DELTA */

int main(int argc, char * argv[]) {
	for (unsigned i =0; i <argc; printf("%s ", argv[i++]));
	printf("\n");

	/* __ Variables */
	int status =0;
	tableaux g;
	g.ctr_bkI =NULL;	/* en cas de pb *avant* de faire évzntuellement appel à set_PL: on n'essaiera pas de libérer la mémoire */
	enum e_type type;
	int E_q =NO_VAL, nu, q, k, p =0, p_phi =NO_VAL;
	enum e_resolution resol;
	long int num_cible, den_cible, R_cible =0;
	long double val_cible;
	
	if (argc <ARGCMIN_CONT_BIN_MINR) {
		fprintf(stderr, "Il n'y a pas assez d'arguments renseignés.\n");
		return 0;
	}

	/* __ Interprétation des arguments du main */
	
	/* lecture des données du problème ; en cas de succès, valeurs:
		type, nu, k 	toujours renseignées
		E_q, q			renseignée si type TAB_RHO ou TAB_GAMMA
		p				renseignée si type TAB_GAMMA ou TAB_DELTA ou TAB_REG
	*/
	status =get_probleme(argc, argv, &type, &E_q, &nu, &q, &p, &k);
	if(status)
		return 0;

	/* résolution (type et paramètres obligatoires) 
	*/
	status =get_optim(argc, argv, &resol, &val_cible, &num_cible, &den_cible, &R_cible);
	if(status)
		return 0;

	/* paramètres optionnels */
	if(resol != RESOL_CONT && resol != RESOL_INT_PV_OPT)
		get_Rcible(argc, argv, resol, &R_cible);
	
	if(type == TAB_GAMMA || type == TAB_DELTA || type == TAB_REG)
		get_p_phi(argc, argv, resol, &p_phi);

	#if (DEBUG)
	printf("\n");
	#endif

	/* __ Initialisation de la structure */
	/* les dimensions du problème */
	if(type == TAB_RHO)	/* champs nu, q, k, E_q de g */
		status =rho_init_param(&g, nu, q, k, E_q);
	else if(type == TAB_GAMMA)	/* champs q, p, k, E_q, p_max (optionnel) et aussi nu de g */
		status =gamma_init_param(&g, q, p, k, E_q, p_phi);
	else if (type == TAB_REG) /* champs nu, p, k et aussi q de g */
		status =reg_init_param(&g, nu, p, k, p_phi);
	else /* type == TABLE_DELTA */ /* champs nu, p, k et aussi q de g */
		status =delta_init_param(&g, nu, p, k, p_phi);
	if(status != 0) {
		fprintf(stderr, "Dimensions incohérentes.\n");
		return 0;
	}
	
	/* sa resolution (champs num_cible, den_cible, num_cible, den_cible, R_min de g) 
		on en profite si gamma (ou gamme_E) et cont on aide (pê) la résolution en ajoutant une borne inf connue
	*/
	status =init_param_optim(&g, resol, num_cible, den_cible, val_cible, R_cible);
	if (status != 0){
		fprintf(stderr, "Resolution (%d) non reconnue ou non disponible ou de paramètres fournis incohérents.\n", resol);
		return 0;
	}

	/* __ Definition du PL (=> on l'ecrit) */
	if(type == TAB_REG)
		status =reg_setPL(& g);
	else /* tous les autres types */
		status =setPL(& g);
	if(status) {
		printf("Pb definition PL.\n");
		goto TERMINATE;
	}

	status =writePL(& g);
	if(status) {
		printf("Pb ecriture PL.\n");
		goto TERMINATE;
	}

	/* __ Resolution du PL (=> on affiche la solution si une solution a ete obtenue) */
	status =resoudre(& g);
	if(status) {
		printf("Pb resolution PL.\n");
		goto TERMINATE;
	}
	
	status =writeSolPL(& g);
	if(status) {
		printf("Pb ecriture solution PL.\n");
		goto TERMINATE;
	}

	status =affSol(& g);
	if(status) {
		printf("Pb affichage solution PL.\n");
		goto TERMINATE;
	}

	/* __ Fin */
TERMINATE:
	status =liberer(& g);
	if(status)
		printf("Pb libération mémoire.\n");

	return EXIT_SUCCESS;
}

/* lecture des données du problème ; en cas de succès, valeurs:
	type, nu, k 	toujours renseignées
	E_q, q			renseignée si type TAB_RHO ou TAB_GAMMA
	p				renseignée si type TAB_GAMMA ou TAB_DELTA ou AB_REG
*/
int get_probleme(int argc, char * argv[], enum e_type* type, int* E_q, int* nu, int* q, int* p, int* k) {
	int status =0;

	/* __ type, nu, k, E_q */
	*type =atoi(argv[INDICE_TYPE]);
	
	if(*type >= 4) /* 4 => TAB_REG (2) ; 5 => TAB_DELTA (3) */
		*type =*type -2;	
	else { /* 0, 1 => TAB_RHO (0) ; 2, 3 => TAB_GAMMA (1) ; 1, 3 => E_q == 1 */
		*E_q =*type % 2;
		*type =*type/2;
	}
	*nu =atoi(argv[INDICE_DIM]);
	*k =atoi(argv[INDICE_K]);
	
	/* __ q (toujours), p, p_phi (si TAB_GAMMA ou TAB_DELTA ou TAB_REG) */
	if (*type == TAB_RHO) {	
		*q =atoi(argv[INDICE_Q_RHO]);
		#if (DEBUG)
		printf("rho: E_q =%d nu =%d q =%d k =%d\t", *E_q, *nu, *q, *k);
		#endif
	} else if (*type == TAB_GAMMA || *type == TAB_DELTA || *type == TAB_REG) {
		*p =atoi(argv[INDICE_P]);
		/* traitement selon TAB_GAMMA d'une part, TAB_DELTA ou TAB_REG de l'autre */
		if (*type == TAB_GAMMA) {
			*q =*nu;
			#if (DEBUG)
			printf("gamma: E_q =%d q =%d p =%d k =%d\t", *E_q, *q, *p, *k);
			#endif
		} else if (*type == TAB_DELTA || *type == TAB_REG) {
			#if (DEBUG)
			printf("delta: nu =%d p =%d k =%d\t", *nu, *p, *k);
			#endif
		}
	} else { /* type de design non reconnu */
		fprintf(stderr, "Type de problème %d non reconnu.\n", *type);
		status =STAT_ERROR_PARAM;
	}
	
	return status;
}

/* lecture du type et des paramètres (obligatoires) d'optimisation ; en cas de succès, valeurs:
	resol					toujours renseignée
	val_cible 				renseignée si RESOL_INT
	num_cible, den_cible	renseignées si RESOL_INT_PV_OPT
	R_cible					renseignée si RESOL_MAX_FREQ
*/
int get_optim(int argc, char * argv[], enum e_resolution* resol, long double* val_cible, long int* num_cible, long int* den_cible, long int* R_cible) {
	long double val_cible_den =0;
	int status =0;
	
	*resol =atoi(argv[INDICE_RESOL]);

	/* __ type de résolution & paramètres obligatoire */
	if (*resol == RESOL_CONT || *resol == RESOL_BIN || *resol == RESOL_R) {	/* résolutions sans paramètre obligatoire */
		#if (DEBUG)
			if(*resol == RESOL_CONT)
				printf("Résolution continue\t");
			if(*resol == RESOL_BIN)
				printf("Résolution binaire\t");
			if(*resol == RESOL_R)
				printf("Résolution nombre de lignes min\t");
		#endif
	}
	else if (*resol == RESOL_INT) {	/* résolution entière: il faut spécifier une val_cible */
		if (argc < ARGCMIN_INT_MAXFREQ) {
			fprintf(stderr, "Val cible non précisée dans la résolution entière.\n");
			status =STAT_ERROR_PARAM;
		} else {
			sscanf(argv[INDICE_VAL_CIBLE], "%Lf/%Lf", val_cible, &val_cible_den);
			if(val_cible_den)
				*val_cible =*val_cible/val_cible_den;
			#if (DEBUG)
			printf("Résolution entière val_cible =%Lf\t", *val_cible);
			#endif
		}
	} else if (*resol == RESOL_INT_PV_OPT) {	/* résolution preuve l'optimalité: il faut spécifier une paire (num_cible, den_cible) */
		if (argc < ARGCMIN_PV_OPT) {
			fprintf(stderr, "Vous n'avez pas précisé le numérateur cible et/ou le dénominateur cible dans le cadre de la preuve de l'optimalité du problème à valeurs entières.\n");
			status =STAT_ERROR_PARAM;
		} else {
			*num_cible =atol(argv[INDICE_NUM_CIBLE]);
			*den_cible =atol(argv[INDICE_DEN_CIBLE]);
			#if (DEBUG)
			printf("Résolution preuve d'optimalité num_cible =%ld den_cible =%ld\t", *num_cible, *den_cible);
			#endif
		}
	} 	else if (*resol == RESOL_MAX_FREQ) {	/* résolution fréquence max: il faut spécifier un nombre R_cible de lignes */
		if (argc <ARGCMIN_INT_MAXFREQ) {
			fprintf(stderr, "Vous n'avez pas précisé le nombre de lignes cible dans le cadre de la maximisation de la fréquence.\n");
			status =STAT_ERROR_PARAM;
		}
		else {
			*R_cible =atol(argv[INDICE_R_CIBLE]);
			#if (DEBUG)
			printf("Résolution fréquence max R_cible =%ld\t", *R_cible);
			#endif
		}
	} else {	/* type de résolution incorrect */
		fprintf(stderr, "Type de résolution %d non reconnu.\n", *resol);
		status =STAT_ERROR_PARAM;
	}

	return status;
}

/* __ lecture des paramètres optionnels d'optimisation */

/* Nombre de lignes (selon le contexte, min ou max ou à atteindre):
	- possible 1er argument optionnel pour RESOL_INT, RESOL_BIN, RESOL_R
	(de valeur 0 si une valeur n'est mise que pour pouvoir spécifier d'autres paramètres optionnels à la suite)
*/
void get_Rcible(int argc, char * argv[], enum e_resolution resol, long int* R_cible) {
	if( (resol == RESOL_INT) && (argc > ARGCMIN_INT_MAXFREQ +IND_ARGOPT_R) ) {
		*R_cible =atol(argv[ARGCMIN_INT_MAXFREQ +IND_ARGOPT_R]);
		#if (DEBUG)
		printf("Paramètre optionnel R_cible =%ld\t", *R_cible);
		#endif
	} else if( ((resol == RESOL_BIN) || (resol == RESOL_R)) && (argc > ARGCMIN_CONT_BIN_MINR +IND_ARGOPT_R) ) {
		*R_cible =atol(argv[ARGCMIN_CONT_BIN_MINR +IND_ARGOPT_R]);
		#if (DEBUG)
		printf("Paramètre optionnel R_cible =%ld\t", *R_cible);
		#endif
	}
} 

/* Pour TAB_GAMMA, TAB_DELTA et TAB_REG :: éventuelle limitation p_phi sur les mots de Phi 
	(outre 0 1 .. (q -1) si GAMMA, 1 1  .. 1 si DELTA):
	- possible 2è argument optionnel pour RESOL_INT, RESOL_BIN, RESOL_R
	- possible 1er argument optionnel pour RESOL_CONT, RESOL_INT_PV_OPT, RESOL_BIN, RESOL_MAX_FREQ
*/
void get_p_phi(int argc, char * argv[], enum e_resolution resol, int* p_phi) {
	if( (resol == RESOL_CONT) && (argc > ARGCMIN_CONT_BIN_MINR +IND_ARGOPT_PMAX_CONT_PVOPT_FREQ) ) {	/* résolution CONT: p_phi 1er argument optionnel */
		*p_phi =atol(argv[ARGCMIN_CONT_BIN_MINR +IND_ARGOPT_PMAX_CONT_PVOPT_FREQ]);
		#if (DEBUG)
		printf("Paramètre optionnel p_phi =%d", *p_phi);
		#endif
	} else if( (resol == RESOL_INT_PV_OPT) && (argc > ARGCMIN_PV_OPT +IND_ARGOPT_PMAX_CONT_PVOPT_FREQ) ) {	/* résolution PV_OPT: p_phi 1er argument optionnel */
		*p_phi =atol(argv[ARGCMIN_PV_OPT +IND_ARGOPT_PMAX_CONT_PVOPT_FREQ]);
		#if (DEBUG)
		printf("Paramètre optionnel p_phi =%d", *p_phi);
		#endif
	} else if( (resol == RESOL_MAX_FREQ) && (argc > ARGCMIN_INT_MAXFREQ +IND_ARGOPT_PMAX_CONT_PVOPT_FREQ) ) {	/* résolution MAX_FREQ: p_phi 1er argument optionnel */
		*p_phi =atol(argv[ARGCMIN_INT_MAXFREQ +IND_ARGOPT_PMAX_CONT_PVOPT_FREQ]);
		#if (DEBUG)
		printf("Paramètre optionnel p_phi =%d", *p_phi);
		#endif
	} else if( (resol == RESOL_INT) && (argc > ARGCMIN_INT_MAXFREQ +IND_ARGOPT_PMAX_INT_BIN_MINR) ) {	/* résolution INT: p_phi 2è argument optionnel */
		*p_phi =atol(argv[ARGCMIN_INT_MAXFREQ +IND_ARGOPT_PMAX_INT_BIN_MINR]);
		#if (DEBUG)
		printf("Paramètre optionnel p_phi =%d", *p_phi);
		#endif
	} else if( ((resol == RESOL_BIN) || (resol == RESOL_R)) && (argc > ARGCMIN_CONT_BIN_MINR +IND_ARGOPT_PMAX_INT_BIN_MINR) ) {	/* résolutions BIN et MIN_R: p_phi 2è argument optionnel */
		*p_phi =atol(argv[ARGCMIN_CONT_BIN_MINR +IND_ARGOPT_PMAX_INT_BIN_MINR]);
		#if (DEBUG)
		printf("Paramètre optionnel p_phi =%d", *p_phi);
		#endif
	}
} 
