/*	rho_pl :

	Données :	
		- nu >= 2, q >= 2, k in {2 ,.., nu} entiers,
		- E_q in {0, 1}
		- resol

	Variables : 
		- Q : Z_q^nu -> [0, 1], N ou {0, 1} selon le cas considéré ; Z_q^nu restreint aux vecteurs de {0} x Z_q^(nu -1) si E_q =1
		- R : a valeur dans [0, 1] ou N selon le cas considéré

	Les problemes considérés pour E_q = 0 (le cas E_q = 1 est similaire) :

	- valeur optimale :
		Max Q(0, 0 ,.., 0)
		s.c. 	sum_v Q(v) = 1
				Q in I_q^k
				Q(v) in [0, 1] 

	- solution binaire (tableau simple) optimale :
		Min R
		s.c.		sum_v Q(v) -R		= 0
				Q(0, 0 ,.., 0)		= 1
				Q in I_q^k
				Q(v) in {0, 1} 

	- solution entiere de valeur au moins une valeur donnée (on ne met pas '=' à cause des arrondis nombres irrationnels) :
		Min R
		s.c.		sum_v Q(v) -R 			 		 = 0
				Q(0, 0 ,.., 0)					>= 1
				Q(0, 0 ,.., 0) -val_cible x R 	>= 0
				Q in I_q^k

	- solution entiere de valeur strictement supérieure à une valeur donnée :
		Min R
		s.c.		sum_v Q(v) -R 						  		 = 0
				Q(0, 0 ,.., 0)							 	>= 1
				num_cible x Q(0, 0 ,.., 0) -den_cible x R 	>= 1
				Q in I_q^k

	- meilleure solution entiere pour un nombre de lignes fixé :
		Max Q(0, 0 ,.., 0)
		s.c.		sum_v Q(v) = R_cible
				Q in I_q^k
*/

#define DEBUG 1

#include <ilcplex/cplex.h>

#define IND_COL_R 0
#define IND_COL_Q_STAR 1

#define IND_ROW_DEF_R 0
#define IND_ROW_R_IS_EQUAL_TO 1
#define IND_ROW_Q_G_1 1
#define IND_ROW_LB 2

enum e_type {
	TAB_RHO,
	TAB_GAMMA
};

enum e_resolution {
	RESOL_CONT =0,
	RESOL_INT,
	RESOL_INT_PV_OPT,
	RESOL_BIN,
	RESOL_BIN_PV_OPT
};

typedef enum e_resolution resolution ;

struct s_rho_pl{
	int nu;
	int q;
	int p;
	int k;
	int E_q;
	enum e_type type; 	// rho ou gamma
	enum e_resolution resol;	// choix de resolution : continu, entier, ...
	double val_cible; 	// valeur cible : quand recherche solution entiere minimisant R quand valeur continue connue
	long int num_cible; 	// num_cible et den_cible : quand on veut demontrer solution entiere Q/R =den_cible / num_cible est opt (du fait arrondis valeurs irationnelles)
	long int den_cible; 	
	long int R_cible;	// R_cible : quand on veut demontrer qu'avec le nombre de lignes min pour des tableaux simples, on ne peut faire mieux que 1/R_cible 	
	int **ctr_bkI;		// associe a (J, a) l'indice de la contrainte associee dans le PL => pour simplifier l'access, tableau de taille 2^q (tous les sous-ens de Z_q)
	CPXENVptr env;
	CPXLPptr pl;
};

typedef struct s_rho_pl rho_pl;

#define STAT_ERROR_MALLOC -2
#define STAT_ERROR_PARAM -1
#define NO_VAL -1	/* utilisee pour les champs val_cible et ctr_bkI[J][a]  */

/* ____ Fonctions d'initialisation :
	Post-condition : 	champs atomiques nu, q, k, E_q, resol, val_cible, num_cible, den_cible, R_cible tous initialises
	Valeur retournee : 	0 si ok, -1 si pb parametres, -2 si pb malloc, code cplex > 0 sinon
	Post-condition en cas de succes :
		- environnement Cplex +PL crees
		- tableau ctr_bkI[][] construit et initialise
	Post-condition en cas d'echec : memoire eventuellement allouee en cours de fonction totalement liberee
*/
int init_cont(rho_pl* g, int q, int nu, int k, int E_q);
int init_int(rho_pl* g, int q, int nu, int k, int E_q, double val_cible);
int init_int_pv_opt(rho_pl* g, int q, int nu, int k, int E_q, long int num_cible, long int den_cible);
int init_bin(rho_pl* g, int q, int nu, int k, int E_q);
int init_bin_pv_opt(rho_pl* g, int q, int nu, int k, int E_q, long int R_cible);

/* ____ nettoyage
	Pre-conditions : champs g->q, g->k, g->_q correctement renseignes
	Post-conditions : memoire liberee
*/
int liberer(rho_pl* g);

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
int setPL(rho_pl* g);

/* ____ resolution PL :
	- sens obj +type probleme (continu)
	- resolution
	- restitution valeur optimale
*/
int resoudre(rho_pl* g);

/* ____ affichage de la solution
*/
int affSol(rho_pl* g);

