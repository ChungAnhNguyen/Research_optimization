/*	vn : manipulation des entiers pour représenter des vecteurs de Z_q^nu / des sous-ensembles de [nu]
*/

/***
	Modification du 17/07/2021
		Remplacement de la fonction pow par la fonction puissance qui traite les entiers (cela permet de ne pas avoir d'éventuelles erreurs d'arrondi)
			==> Tests OK sur gamma(3, 2, 2), gamma(3, 4, 2) et gamma_E(4, 2, 2), gamma_E(5, 3, 2), gamma_E(6, 3, 2)
			
	Ajout 7/8/2021 : fonction pgcd
			==> TODO : fonction de test

	Modification 29/8/21 : fonction binom (aout de cas non traités précédemment)
			==> TODO : tester ok (pas de régression)
***/

/* renvoie le plus grand diviseur commun à la valeur absolue de deux entiers a et b
*/
int pgcd(int a, int b);

/* renvoie a modulo q (l'operateur % considère (-2) % 3 == -2 ?!?) 
	Ce nombre est dans {0 ,.., q -1}
*/
int mod(long int a, int q);

/* renvoie fact(n)
	Pre-conditions : n >= 1
*/
int fact(int n);

/* renvoie binom(q, k)
	Pre-conditions : q >= k >= 2
		ST 29/8/21 finalement pas, j'ai besoin des autres cas aussi
*/
long long int binom(int q, int k);

/* renvoie T(q, k) = sum_{r =0}^k binom(q, r) binom(q -1 -r, k -r)
	Pre-conditions : q > k >= 0
*/
long long int Tqk(int q, int k);

/* renvoie V(q, p, k) = sum_{r =0}^k binom(q, r) binom(q -1 -r, p -r) binom(p -1 -r, k -r)
	Pre-conditions : q > p > k >= 0
*/
long long int Vqpk(int q, int p, int k);

/* renvoie le nombre dont la decomposition en base q est (0, 1 ,..., q -1) 
	Ce nombre est dans {0 ,.., q^q -1}
*/
int get_v_012etc(int q);

/* renvoie le coefficient v_r de q^(nu -r) dans la décomposition v_1 q^(nu -1) +...+ v_nu q^0 d'un nombre v in {0 ,..., q^nu -1} en base q 
	Ce coefficient est dans {0 ,.., q -1}
*/
int get_v_r(int v, int r, int q, int nu);

/* renvoie le nombre de valeurs distinctes prises par la décomposition (v_1 ,..., v_nu) d'un nombre v in {0 ,..., q^nu -1} en base q */
int get_nb_val(int v, int q, int nu);

/* renvoie la taille du sous-ensemble de [nu] dont le vecteur d'incidence est la décomposition binaire du nombre J in {0 ,..., 2^nu -1} */
int get_card(int J, int nu);

/* pour J sous-ensemble de cardinalité t de [nu], a in Z_q^t et v in Z_q^nu, cette fonction renvoie 1 si 
		- E_q = 0 et v_J =a ou
		- E_q = 1 et v_J in {a, a +(1 ,.., 1) ,.., a +(q -1 ,.., q -1)}
	et 0 sinon 
	Pre-condition : si E_q == 1, alors a_1 = 0 
*/
int is_vJ_equal_to_a(int v, int J, int a, int q, int nu, int t, int E_q);

/* pour v =(v_1 ,.., v_nu) in Z_q^nu et J =(j_1 ,.., j_t) sous-ensemble de taille t de [nu], renvoie le nombre a =(a_1 ,.., a_t) de Z_q^t vérifiant :
	- v_J =a si E_q =0
	- v_J -(v_{j_1} ,.., v_{j_1}) si E_q = 1
*/
int get_v_J(int v, int J, int q, int nu, int t, int E_q);

/* affiche la decomposition d'un nombre v in {0 ,..., q^nu -1} en base q */
void affiche_v(int v, int q, int nu);

/* ecrit dans buffer la decomposition en base q de v
	Pré-conditions : q <= 10 et buffer est de taille suffisante
*/
void itoa(int v, int q, int nu, char* buffer);

/* affiche taille rho(nu, q, k) */
void display_taille_rho(int nu, int q, int k);

/* affiche taille rho_E(nu, q, k) */
void display_taille_rho_E(int nu, int q, int k);

/* affiche taille gamma(q, p, k) */
void display_taille_gamma(int q, int p, int k);

/* affiche taille gamma_E(q, p, k) */
void display_taille_gamma_E(int q, int p, int k);

/* Cette fonction renvoie le nombre de partitions en k >= 1 sous-ensembles non vides d'un ensemble de taille n >= 1 */
int get_nbpart(int n, int k);


/* Ajout 17/07/2021 */

long long unsigned puissance(unsigned base, unsigned exposant);

/* Fin Ajout 17/07/2021 */
