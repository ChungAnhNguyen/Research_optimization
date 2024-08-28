/*	tableaux :: calcul 
		d'OA(., nu, q, k) (dont DS), 
		de Gamma(q, p, k) (dont Gamma_E) et 
		de Delta(nu, p, k) (réguliers pour ces derniers) 
	optimaux pour différents critères d'optimisation

	On note v* =(0 ,.., 0) pour rho et (0, 1 ,.., q -1) pour gamma.

	_____________________________________________________________________________ Rho: v* =(0 ,..., 0)
	Données:	
	- nu >=2, q >=2, k in {2 ,.., nu} entiers,
	- E_q in {0, 1}
	- resol

	Variables: 
		- Q: 	Z_q^nu -> [0, 1], N ou {0, 1} selon le cas considéré
				représentant le nombre d'occurrences ou la fréquence des vecteurs de Z_q^nu dans le tableau
				Si E_q =1, on se restreint sans perte de généralité à {0} x Z_q^(nu -1)
		- R: 	variable atomique à valeur dans [0, 1] ou N selon le cas considéré
				représentant le nombre de lignes ou la fréquence cumulée dans le tableau

	Les problèmes considérés pour E_q =0 (le cas E_q =1 est similaire):

	Tous considèrent le bloc (B) de contraintes:
		sum_v Q(v) 	 -R =0		// contrainte liant les variables Q à la variable R
		Q(v*) -Q(v) >= 0 pour tout v in Z_q^nu	
								// fixation arbitraire d'un vecteur le plus fréquent (pour éviter les symétries inutiles) 
		sum_{v in Z_q^nu: v_J =a} Q(v)	-R/q^k	=0 pour tout sous-ensemble J de cardinalité k de [nu] et tout a in Z_q^k
								// Q est "balanced k-wise independent"

	0. valeur optimale:
		Max Q(v*)
		s.c. 	(B)
				R 								 =1
				Q(v) 							in [0, 1] 
	(objet: exhiber rho(nu, q, k))

	1. solution entière de valeur au moins une valeur donnée (on ne met pas '=' à cause des arrondis nombres irrationnels):
		Min R
		s.c. 	(B)
				Q(v*)							>=1
				Q(v*) -val_cible xR 			>=0
				Q(v) 							in N
				+contrainte optionnelle: R =R_cible
	(objet: exhiber un OA réalisant rho(nu, q, k), et tant qu'à faire on en exhibe un qui est minimal en nombre de lignes)

	2. solution entière de valeur strictement supérieure à une valeur donnée:
		Min R
		s.c. 	(B)
				num_cible xQ(v*) -den_cible xR 	>=1
				Q(v) 							in N
	(objet: démonter qu'une valeur num_cible/den_cible donnée est bien la valeur de rho(nu, q, k))

	3. solution binaire (tableau simple) optimale:
		Min R	/!\ pour rho_E: on autorise un unique représentant par classe de vecteurs {u, u +1 ,.., u +q -1}  /!\
		s.c. 	(B)
				Q(v*)			 				>=1
				Q(v) 							in {0, 1} 

				+contrainte optionnelle: R =R_cible
	(objet: exhiber un OA minimal en nombre de lignes parmi ceux qui n'admettent pas de répétition de lignes)

	4. meilleure solution entière pour un nombre de lignes fixé:
		Max Q(v*)
		s.c. 	(B)
				R 								 =R_cible
				Q(v) 							in N
	(objet: exhiber, parmi les OA qui sont minimaux en nombre de lignes, un OA dans lequel la plus grande fréquence est maximale)

	5. plus petit nombre de lignes:
		Min R
		s.c.	(B)
				Q(v*)							>=1
				Q(v) 							in N
				+contrainte optionnelle: R =R_cible
	(objet: exhiber un OA minimal en nombre de lignes)

	_____________________________________________________________________________ Gamma / Delta: v* =(0, 1 .. q -1) / (1 .. 1)
	Données:
	- q >p >=k 3 entiers (et on note nu =q) / nu >p >=k (et on note q =2)
	et on note U ={v in Z_q^q | |{v_1 ,.., v_q}| <=p} / les mots de poids au plus p de {0, 1}^nu

	Variables:
	- P: U 		-> [0, 1], N ou {0, 1} selon le cas considéré; U restreint aux vecteurs de {0} xZ_q^(q -1) si E_q =1 et TAB_GAMMA
	- Q: Z_q^nu -> [0, 1], N ou {0, 1} selon le cas considéré; {0} xZ_q^(q -1) si E_q =1 et TAB_GAMMA
				
	Contraintes:
	- sum_{v in Z_q^nu} P(v) =sum_{v in Z_q^nu} Q(v) =1
	- P(v) >0 => v in U
	- pour tout J =(j_1 ,..., j_k) et tout a in Z_q^k, sum_{v in Z_q^nu: v_J =a} P(v) -sum_{v in Z_q^nu: v_J =a} Q(v) =0

	Objectif:
	- max Q(v*)

	PL:
	- variables P(v): v dont les composantes prennent au plus p valeurs distinctes, et v_1 =0 si E_q =1
	- variables Q(v): v_1 =0 si E_q =1
	- numérotation des variables:
		- cas entier ou binaire: variable R numérotée 0
		- dans tous les cas: variables Q(v), P(v) numérotées ensuite
	- numérotation des contraintes:
		- cas continu: contrainte P mesure numérotée 0 (NB contrainte Q mesure est redondante)
		- cas entier ou binaire: contraintes R =sum_v P(v) et Q(v^*) >= 1 numérotées 0 et 1
		- cas entier: contrainte Q(v^*)/R >= val_cible numérotée 2
		- dans tous les cas: contraintes bki numerotees ensuite

	Les problèmes considérés pour E_q =0 (le cas E_q =1 est similaire):

	Tous considèrent le bloc (B) de contraintes:
		sum_(v in U) P(v) -R =0	// contrainte liant les variables P (et donc Q) à la variable R
		sum_{v in U: v_J =a} P(v) -sum_{v in Z_q^nu: v_J =a} Q(v) =0
									// P -Q est "balanced k-wise independent"

	0 valeur optimale:
		Max Q(v*)
		s.c. 	(B)
				R 									 =1
				P(v), Q(v) 							in [0, 1] 
			+on espère aider avec 
				(T(q -p +k, k) +1) xQ(v*) -2 xR 	>= 0

	1 solution entière de valeur au moins une valeur donnée (on ne met pas '=' à cause des arrondis nombres irrationnels):
		Min R
		s.c. 	(B)
				Q(v*)								>=1
				Q(v*) -val_cible xR 				>=0
				P(v), Q(v) 							in N 
			+contrainte optionnelle: R =R_cible

	2 solution entière de valeur strictement supérieure à une valeur donnée:
		Min R
		s.c. 	(B)
				Q(v*)							 	>=1
				num_cible xQ(v*) -den_cible xR 		>=1
				P(v), Q(v) 							in N 
			+contrainte optionnelle: R =R_cible

	3 solution binaire (tableaux simples) optimale:
		Min R		/!\ pour gamma_E: on autorise un unique représentant par classe de vecteurs {u, u +1 ,.., u +q -1}  /!\
		s.c. 	(B)
				Q(v*)								>=1
				P(v), Q(v) 							in {0, 1} 
			+contrainte optionnelle: R =R_cible

	4 meilleure solution entière pour un nombre de lignes fixé:
		Max Q(v*)
		s.c. 	(B)
				R 									 =R_cible
				P(v), Q(v) 							in N 
				
	5. plus petit nombre de lignes:
		Min R
		s.c.	(B)
				Q(v*)								>=1
				P(v), Q(v) 							in N

	NB les variables P(v), Q(v) correspondant à un vecteur v dont le nombre de composantes est >max(p, p_phi) ne sont pas générées

	_____________________________________________________________________________ Delta solution régulières (TODO)


	_____________________________________________________________________________ Bornes inf.
	/!\ BORNES INF 
			PARTIELLEMENT INTÉGRÉES. 
			ACTIVABLES OU DÉSACTIVABLES PAR DO_LB 
	/!\

	Pour rho, on dispose des bornes (exactes ou) inférieures:
		TODO 
	on sait aussi:
	- R est un multiple de q^k si E_q =0 et de q^(k -1) si E_q =1 							(non exploite pour l'instant)
	- on dispose des bornes triviales 	F(nu, q, k) <= q^nu et E(nu, q, k) <= q^(nu -1) 	(non exploitees pour l'instant)

	Pour gamma et delta, on dispose des bornes (exactes ou) inférieures:	(EXPLOITÉES POUR GAMMA SEULEMENT)
		gamma_E(q, q, q) =gamma(q, q, q) =1 
		gamma_E(q, p, 2) >=gamma(q, p, 2) =delta(q, p, 2) =floor(p/2)ceil(p/2)/[(q -floor(p/2))(q -ceil(p/2))]
		gamma_E(q, p, k) >=gamma(q, p, k) >=gamma(q -p +k, k, k) =delta(q -p +k, k, k) =2/(T(q -p +k, k) +1) si q >k
			où T(q -p +k, k, k) =sum_{r =0}^k binom(q -p +k -1 -r, k -r) x binom(q -p +k, r) 

	_____________________________________________________________________________ ATTENTION
	________________ p_phi pour gamma et delta
	- p_phi à interpréter différemment dans les delta et dans les gamma 
		(ne serait-ce que parce que dans les Delta, il n'existe qu'un vecteur de poids p tandis que dans les gamma, il existe plusieurs vecteurs dont les composantes sont toutes différentes)
	- au sein de gamma même, interprétation différente selon le type de résolution
		-> si résolution continue, *la borne inf* contraint *aussi* la solution

	________________ champs multiforme (nu, q)
	- on utilise nu pour indiquer la dimension des rho *et des delta* (vs. q pour les gamma)
	- delta modélisé de façon différente de rho et gamma (solutions régulières => on s'intéresse aux familles de lignes d'un poids donné plutôt qu'à aux lignes individuellement)

	________________ paramètres d'optimisation dans le PL 
		(pas nécessairement actifs ni interprétés de la même façon selon le design)
	pour rho et gammma, 6 types de résolution 
		(continue, entière, pv opt, binaire, binaire, max fréq., min R), 
	avec les paramètres obligatoires ::
	- si résolution continue :: pour GAMMA, borne inf automatiquement intégrée sous forme de fraction num_cible/den_cible
	- si résolution entière :: valeur min
	- si résolution pv opt :: fraction min
	- si résolution max fréq. :: #lignes cible
	et les paramètres optionnels ::
	- pour GAMMA et TOUTES les résolutions :: p_phi
	- pour les résolutions INT, BIN et MinR ::  #lignes cible
	pour delta, 3 types de résolution 
		(continue, entière, pv opt):
	avec les paramètres obligatoires ::
	- si résolution entière :: valeur cible
	- si résolution pv opt :: fraction cible
	et les paramètres optionnels ::
	- si résolution continue :: p_phi

	NB #lignes min ou cible géré par les bbor,nes inf et sup sur la variable R (et non comme unne contrainte)

	_____________________________________________________________________________ Organisation du code / implémentation	
	________________ variables "interdites"
	reg.c ::		définition du PL et affichage de la solution pour les delta *réguliers*
	rgd.c ::		didem pour les rho, gamma et delta *généraux*, +libération mémoire contraintes BkI
	tableaux.c ::	tout le reste
	(+ ces fichiers utilisent vn.c)

	________________ variables "interdites"
	dans rho et gamma on ne génère pas les variables superflues (ex. P(v) pour gamma et v interdit dans Psi), tandis que das delta, on met simplement ces variables à 0 (variables z_j et R_j associées aux familles des mots de poids j pour j >max{p, p_phi}) 
*/

#ifndef TABLEAUX
#define TABLEAUX

#define TRACE_DELTA 0
#define DEBUG 0
#define DO_LB 0	/* pour la résolution continue: intégration d'une born inférieure connue (quand elle existe) */

#include <ilcplex/cplex.h>

#define NO_VAL -1	/* utilisee pour les champs: resol, val_cible, num_cible, den_cible, R_cible, p, p_phi, ctr_bkI[J][a] */

enum e_type {
	TAB_RHO,
	TAB_GAMMA,
	TAB_REG,
	TAB_DELTA
};

enum e_resolution {
	RESOL_CONT =0,
	RESOL_INT,
	RESOL_INT_PV_OPT,
	RESOL_BIN,
	RESOL_MAX_FREQ,
	RESOL_R
};

struct s_tableaux{
	int nu;				// le cas échéant, >=q (rho) 	ou >=p (delta)
	int q;				// le cas échéant, >=p (gamma)
	int p;				// le cas échéant, >=k (gamma, delta)
	int k;				// >= 1 (toujours)
	int E_q;
	enum e_type type; 	// rho ou gamma
	enum e_resolution resol;	// choix de resolution: continu, entier, ...
	int p_phi;	// pour gamma et delta (doit être >= 0 et <=q pour gamma, <= nu -1 pour delta)
	long double val_cible; 	// valeur cible: quand recherche solution entiere minimisant R quand valeur continue connue
	long int num_cible; 	// num_cible et den_cible: quand on veut demontrer solution entiere Q/R =den_cible / num_cible est opt (du fait arrondis valeurs irrationnelles)
	long int den_cible; 	
	long int R_cible;	// R_cible: quand on veut demontrer qu'avec le nombre de lignes min pour des tableaux simples, on ne peut faire mieux que 1/R_cible; aussi utilisé comme paramètre optionnel de certaines résolutions pour forcer un nombre exact de lignes
	int **ctr_bkI;		// associe a (J, a) l'indice de la contrainte associee dans le PL => pour simplifier l'acces, tableau de taille 2^q (tous les sous-ens de Z_q)
	CPXENVptr env;
	CPXLPptr pl;
	int ind_R_star;	// indice dans le PL de la variable R*
	int ind_R;	// indice dans le PL de la variable R
};

typedef struct s_tableaux tableaux;

#define STAT_ERROR_MALLOC -2
#define STAT_ERROR_PARAM -1

/* ____________________________________________________________________________________________ Initialisation / vérification des paramètres
*/

/* _____________ Initialisation et vérification des champs de la structure (dimensions) 
	Post-condition: sont initialisés les champs atomiques 
	- type, nu, q, k, E_q																			si rho
	- type, nu (=q, exploité), q, p, k, E_q, p_phi (=q si valeur NO_VAL transmise en paramètre)		si gamma (p_phi parametre optionnel)
	- type, nu, q (=2, non exploité), p, k, p_phi (=nu -1 si valeur NO_VAL transmise en paramètre)	si delta (p_phi parametre optionnel)
	
	NB delta: paramètre d'optimisation R_cible non exploité encore (dans le PL)
*/
int rho_init_param(tableaux* g, int nu, int q, int k, int E_q);
int gamma_init_param(tableaux* g, int q, int p, int k, int E_q, int p_phi);	// p_phi parametre optionnel
int delta_init_param(tableaux* g, int nu, int p, int k, int p_phi);			// p_phi parametre optionnel
int reg_init_param(tableaux* g, int nu, int p, int k, int p_phi);			// p_phi parametre optionnel

/* _____________ Initialisation et vérification des paramètres d'optimisation
		(champs resol, val_cible, num_cible, den_cible, R_cible selon le mode de resolution)
	paramètres:
	- val_cible: paramètre obligatoire pour RESOL_INT (on veut R^* / R =val_cible)
	- R_cible: paramètre optionnel pour RESOL_INT, RESOL_BIN et RESOL_R et RESOL_BIN, obligatoire pour FREQ (dans tous les cas on veut R =R_cible)
	- num_cible, den_cible: paramètres obligatoires pour PV_OPT (on veut R^* /R >=num_cible/den_cible unfeas)
	
	cas gamma (ou gamme_E): si résolution cont, on ajoute une borne inf connue en définissant num_cible et den_cible

	valeur retournee: -1 (STAT_ERROR_PARAM) si mode de résolution non reconnu et 0 sinon
*/
int init_param_optim(tableaux* g, enum e_resolution resol, long int num_cible, long int den_cible, long double val_cible, long int R_cible);

/* ____________________________________________________________________________________________ Gestion mémoire
*/

/* ____ nettoyage
	Pre-conditions: champs g->q, g->k, g->_q correctement renseignes
	Post-conditions: memoire liberee
*/
int liberer(tableaux* g);

/* rho et gamma :: libère la mémoire prise pour les contraintes Bki 
	(appelée par liberer) */
void libererBkI(tableaux* g);

/* ____________________________________________________________________________________________ Définition PL
*/

/* Construction environnement Cplex & PL
	(sous-routine commune à setPL *et reg_setPL*)
*/
int init_cplex(tableaux* g);

/* Sens et type optimisation PL
	(sous-routine commune à setPL *et reg_setPL*)
	Sens optimisation:
		- Max 		si RESOL_CONT ou RESOL_MAX_FREQ
		- Min		si RESOL_BIN ou RESOL_INT ou RESOL_INT_PV_OPT
	Type probleme:
		- LP		si RESOL_CONT
		- MILPL		dans tous les autres cas
	Pre-conditions: PL totalement defini, sauf type et sens optimisation
	Post-condition: type et sens optimisation definis
*/
int finaliser(tableaux* g);

/* Borne inferieure nombre de lignes */
int set_R_min(tableaux* g, long int R_min);

/* Borne supérieure nombre de lignes */
int set_R_max(tableaux* g, long int R_max);

/* ____ definition du PL :: rho et gamma 
	Variables:
		- R: 					toujours, a valeur dans [0, 1], [1, +infty[ ou N^*
		- Q(v), v in Z_q^nu: 	toujours, a valeur dans [0, 1], N, {0, 1} selon g->resol
		- P(v), v in U: 		si gamma, a valeur dans [0, 1], N, {0, 1} selon g->resol
	Contraintes:
		- R =sum_v P(v)	 						si gamma
		- R =sum_v Q(v)	 						si rho
		- R =1									si RESOL_CONT
		- R =R_cible	 						si RESOL_MAX_FREQ ou (RESOL_INT ou RESOL_R ou RESOL_BIN) et R_cible > 0
		- Q(v*) >= (ou =) 1 					si RESOL_INT ou RESOL_BIN ou RESOL_INT_PV_OPT ou RESOL_R
		- Q(v*) -Q(v) >= 0 						si rho et pas RESOL_BIN
		- Q(v*)/R >= val_cible					si RESOL_INT
		- Q(v*)/R >= num_cible/den_cible +1		si RESOL_INT_PV_OPT
		- Q(v*)/R >= num_cible/den_cible			si gamma et RESOL_CONT
		- Q ou Q_E / P -Q (ou P_E -Q_E) I_q^k	toujours
	Objectif:
		- Max Q(v*)								si RESOL_CONT ou RESOL_MAX_FREQ
		- Min R									si RESOL_BIN ou RESOL_INT, RESOL_INT_PV_OPT, RESOL_R
	NB on utilise les lb sur les variables R et Q(v*) pour representer les contraintes
					R 	 =1			: cas RESOL_CONT
					R 	 =R_cible	: cas MAX_FREQ, éventuellement RESOL_INT ou RESOL_R ou RESOL_BIN
			Q(v*)		>= 1 		: cas RESOL_BIN, RESOL_INT, RESOL_INT_PV_OPT, RESOL_R
	Pre-conditions: g intialisee, parametres coherents
	Post-condition: PL totalement defini (dont type et sens optimisation), tableau ctr_bkI rempli
*/
int setPL(tableaux* g);


/* Definition du PL :: delta (solutions régulières)
	Variables:		/!\ ici R représente le nombre de lignes de Psi *plus* le nombre de lignes de Phi /!\
		- R: 					de valeur 1 			si RESOL_CONT, 	dans N\{0} si RESOL_INT ou RESOL_INT_PV_OPT
		- R_i, i in {0 .. nu}: 	à valeur dans [0, 1]	si RESOL_CONT,	dans N si RESOL_INT ou RESOL_INT_PV_OPT
		- z_i, i in {0 .. nu}: 	à valeur dans [-1, 1]	si RESOL_CONT,	dans Z si RESOL_INT ou RESOL_INT_PV_OPT,  <=0 si i >p
	NB on peut considérer que les variables de décision sont les variables {z_j, j =k ,.., nu -1} 
		(la valeur, ou du moins la valeur optimale des autres variables en découlant)
	Contraintes:
		z_nu =sum_{j =k}^{nu -1} C_{nu -k}^{j -k} z_j
		z_i  =(-1)^{k -i} sum_{j =k}^{nu -1} C_{nu -1 -i}^{j -i} C_{j -1 -i}^{k -1 -i} z_j
		R_j  =-C_nu^j z_j, j in {p +1 .. nu}
		R_i >=C_nu^j z_i >=-R_i, i in {0 .. p}
		R	 =sum_{j =0}^nu R_j
	En plus selon type de résolution:
		2R_nu 				-val_cible xR 	>=0		si RESOL_INT
		2den_cible xR_nu 	-num_cible xR 	>=1		si RESOL_INT_PV_OPT
	Objectif :
		- Max 2R_nu		en posant R =1 pour RESOL_INT
		- Min R			en ajoutant la contrainte 2 R_nu >= R val_cible 				pour RESOL_INT
		- Min R			en ajoutant la contrainte 2 R_nu den_cible >= R num_cible +1 	pour RESOL_INT_PV_OPT
	Pre-conditions: g intialisee, parametres coherents
	Post-condition: PL totalement defini (dont type et sens optimisation)
	
	NB les variables z_j, R_j d'indice j >max(p, p_phi), <nu sont générées malgré tout, de domaine {0} 
*/
int reg_setPL(tableaux* g);

/* ____________________________________________________________________________________________ Résolution PL
*/

/* _____________ Resolution PL:
	- sens obj +type probleme (continu)
	- resolution
	- restitution valeur optimale
*/
int resoudre(tableaux* g);

/* ____________________________________________________________________________________________ Sauvegardes / affichages
*/

/* _____________ PL avant résolution
*/

/* Ecrit le PL au format LP */
int writePL(tableaux* g);

/* _____________ Statut résolution après résolution
*/

/* Affichage synthetique du resultat de la resolution
	Pre-condition: le PL a ete resolu
	NB statuts resolution :
		https://www.tu-chemnitz.de/mathematik/discrete/manuals/cplex/doc/refman/html/appendixB.html
		1 :	CPX_STAT_OPTIMAL
		2 :	CPX_STAT_UNBOUNDED
		3 :	CPX_STAT_INFEASIBLE 
		101: CPXMIP_OPTIMAL 
		103: CPXMIP_INFEASIBLE 
		118: CPXMIP_UNBOUNDED 
*/
int restitution(tableaux* g);

/* _____________ Solutions après résolution
*/

/* Sauvegarde l'éventuelle solution du PL */
int writeSolPL(tableaux* g);

/* Affichage solution 
	(fait appel aux fonctions dédiées rho, gamma, delta) 
*/
int affSol(tableaux* g);

/* Affichage des solutions (dédié selon le type de problème)
	/!\ fonctions définies dans rhogamma.c et delta.c /!\
*/
int rho_affSol(tableaux* g);
int gamma_delta_affSol(tableaux* g);
int reg_affSol(tableaux* g);

#endif
