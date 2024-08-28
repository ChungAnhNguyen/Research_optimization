/*	vn : manipulation des entiers pour représenter des vecteurs de Z_q^nu / des sous-ensembles de [nu]
*/

#include "vn.h"

#include <stdio.h>	// pour printf() (fonctions d'affichage)

static void affiche_nombre(long int x);

/* renvoie le plus grand diviseur commun à la valeur absolue de deux entiers
*/
int pgcd(int a, int b) {
	int r;

	/* on transforme (a, b) en (max{|a|,|b|}, min{|a|,|b|}) */
	if(a <0)
		a =-a;
	if(b <0)
		b =-b;
	if(a < b) {
		r =a;
		a =b;
		b =r;
	}

	/* Calcul de pgcd(a, b) où a >= b >= 0 par l'algorithme d'Euclide :
		pour 2 entiers naturels a >= b, pccg(a, b) est égal à
			a 					si b =0
			pgcd(b, a mod b)	sinon
	*/
	while(b >0) {
		r = a %b;
		a = b;
		b = r;
	}

	return a;
}

/* renvoie a modulo q (l'operateur % considère (-2) % 3 == -2 ?!?) 
	Ce nombre est dans {0 ,.., q -1}
*/
int mod(long int a, int q) {
	while (a < 0) {
		a =a +q;
	}

	return (a % q);
}

/* renvoie fact(n)
	Pre-conditions : n >= 1
*/
int fact(int n) {
	int res =1;

	while(n >= 2) {
		res = res * n;
		n --;
	}

	return res;
}

/* renvoie binom(q, k)
	Pre-conditions : q >= k >= 2
		ST 29/8/21 finalement pas, j'ai besoin des autres cas aussi
*/
long long int binom(int q, int k) {
	int res =1, r;
	
	// DEBUG printf("%s(%d, %d)\n)", __func__, q, k);
	
	/* __ ajouts ST 29/8/21 */
	if(k > q || k < 0)		/* k < q ou k < 0 => 0 */
		return 0;
	if(k ==0 || q == k)		/* 0 <= k <= q et (k == 0 ou k == q) */
		return 1;
	if(k ==1 || q == k -1)	/* 1 <= k <= q -1 et (k == 1 ou k == q -1) */
		return q;
	/* fin ajouts ST 29/8/21 */

	/* __ ajouts ST 30/8/22 */
	if(q -k < k)
		k =q -k;
	/* fin ajouts ST 30/8/22 */

	/* q x..x (q -k +1) */
	for (r =0 ; r <= k -1 ; r ++) {
		res = res *(q -r);
	}

	/* k x..x 2 */
	for (r =k ; r >= 2 ; r --) {
		res = res /r;
	}

	return res;
}

/* renvoie T(q, k) =sum_{r =0}^k binom(q -1 -r, k -r) binom(q, r)
	Pre-conditions : q >k >=0
*/
long long int Tqk(int q, int k) {
	long long int res =0, r;

	for (r =0 ; r <= k ; r ++)
		res = res +binom(q, r) *binom(q -1 -r, k -r);

	return res;
}

/* renvoie V(q, p, k) = sum_{r =0}^k binom(q, r) binom(q -1 -r, p -r) binom(p -1 -r, k -r)
	Pre-conditions : q > p > k >= 0
*/
long long int Vqpk(int q, int p, int k) {
	long long int res =0, r;

	for (r =0 ; r <= k ; r ++)
		res = res +binom(q, r) *binom(q -1 -r, p -r) *binom(p -1 -r, k -r);

	return res;
}

/* renvoie le nombre dont la decomposition en base q est (0, 1 ,..., q -1) 
	Ce nombre est dans {0 ,.., q^q -1}
*/
int get_v_012etc(int q) {
	int v =0, r;

	for (r = 1 ; r <= q ; r ++)
		v = v +(r -1) * (int)puissance(q, q -r);

	return v;
}

/* renvoie le coefficient v_r de q^(nu -r) dans la décomposition v_1 q^(nu -1) +...+ v_nu q^0 d'un nombre v in {0 ,..., q^nu -1} en base q 
	Ce coefficient est dans {0 ,.., q -1}
*/
int get_v_r(int v, int r, int q, int nu) {
	int v_r;

	v_r = v/(int)puissance(q, nu -r);	/* le resultat v_r appartient a {0 ,.., q -1}  */
	v_r = v_r % q;

	return v_r;
}

/* renvoie le nombre de valeurs distinctes prises par la décomposition (v_1 ,..., v_nu) d'un nombre v in {0 ,..., q^nu -1} en base q */
int get_nb_val(int v, int q, int nu) {
	int nb_val = 0, a, r;

	/* parcours des entiers {0 ,..., q -1} */
	for (a = 0 ; a <= q -1 ; a ++) {
		/* parcours de (v_1 ,..., v_nu) */
		for (r = 1 ; r <= nu ; r ++) {
			if (a == get_v_r(v, r, q, nu)) {
				nb_val ++;
				break;
			}
		}
	}

	return nb_val;
}

/* renvoie la taille du sous-ensemble de [nu] dont le vecteur d'incidence est la décomposition binaire du nombre J in {0 ,..., 2^nu -1} */
int get_card(int J, int nu) {
	int card = 0, r;

	for (r = 1 ; r <= nu ; r ++) {
		card = card +get_v_r(J, r, 2, nu);
	}

	return card;
}

/* pour J sous-ensemble de cardinalité t de [nu], a in Z_q^t et v in Z_q^nu, cette fonction renvoie 1 si 
		- E_q = 0 et v_J =a ou
		- E_q = 1 et v_J in {a, a +(1 ,.., 1) ,.., a +(q -1 ,.., q -1)}
	et 0 sinon 
	Pre-condition : si E_q == 1, alors a_1 = 0 
*/
int is_vJ_equal_to_a(int v, int J, int a, int q, int nu, int t, int E_q) {
	int egal = 1, r, s =0, v_j1;

	/* parcours des composantes v_r de v pour r in J =(j_1 ,..., j_t) 
		- E_q = 0 : on test v_J = a
		- E_q = 1 : on test v_J -(v_{j_1} ,..., v_{j_1}) = a -(a_1 ,..., a_1)
	*/

	for (r = 1 ; (r <= nu) && egal ; r ++) {
		if(get_v_r(J, r, 2, nu) == 1)	/* J in {0, 1}^nu => on teste sa r-ème composante */  {
			s ++;

			/* pour E_q = 1 : on recupere v_{j_1} */
			if (s == 1)	/* r est le plus petit indice de J */ {
				v_j1 = get_v_r(v, r, q, nu);		/* v in Z_q^nu et l'on prend sa r-ème composante */
			}

			/* r =j_s : test v_{j_s} = a_s si E_q =0, v_{j_s} -v_{j_1} = a_s si E_q =1 (NB les pre-conditions impliquent alors a_s = a_s -a_1) */
			if (E_q == 0) {
				egal = (get_v_r(v, r, q, nu) == get_v_r(a, s, q, t));
			}
			else {
				egal = (mod(get_v_r(v, r, q, nu) -v_j1, q) == get_v_r(a, s, q, t));
			}
		}
	}

	return egal;
}

/* pour v =(v_1 ,.., v_nu) in Z_q^nu et J =(j_1 ,.., j_t) sous-ensemble de taille t de [nu], renvoie le nombre a =(a_1 ,.., a_t) de Z_q^t vérifiant :
	- v_J =a si E_q =0
	- v_J -(v_{j_1} ,.., v_{j_1}) si E_q = 1
*/
int get_v_J(int v, int J, int q, int nu, int t, int E_q) {
	int vJ =0, r, s =0, v_j1;

	/* parcours des composantes v_r de v pour r in J =(j_1 ,..., j_t) */
	for (r = 1 ; (r <= nu) && (s <= t) ; r ++) {
		if(get_v_r(J, r, 2, nu) == 1)	/* J in {0, 1}^nu => on teste sa r-ème composante */  {
			s ++;	/* apres mise a jour de s : r =j_s */

			/* si E_q = 1 et s = 1 : on recupere v_{j_1} */
			if (s == 1) {
				v_j1 = get_v_r(v, r, q, nu); /* v in Z_q^nu et l'on prend sa r-ème composante */
			}

			/* on ajoute a vJ (v_r x q^(t -s)) si E_q = 0, ((v_r -v_j1) x q^(t -s)) si E_q =1 */
			if (E_q == 0) {
				vJ =vJ +get_v_r(v, r, q, nu) * (int)puissance(q, t -s);
			}
			else {
				vJ =vJ +mod(get_v_r(v, r, q, nu) -v_j1, q) * (int)puissance(q, t -s);
			}
		}
	}

	return vJ;
}

/* affiche la decomposition d'un nombre v in {0 ,..., q^nu -1} en base q */
void affiche_v(int v, int q, int nu) {
	int r;

	for (r = 1 ; r <= nu ; r ++)
		printf("%d ", get_v_r(v, r, q, nu));
}


/* ecrit dans buffer la decomposition en base q de v
	Pré-conditions : q <= 10 et buffer est de taille suffisante
*/
void itoa(int v, int q, int nu, char* buffer) {
	int r;
	// printf("%s(%d, %d, %d, %p) IN\n", __func__, v, q, nu, buffer);

	for (r = 1 ; r <= nu ; r ++)
		buffer[r -1] = '0' +(char)get_v_r(v, r, q, nu);

	buffer[nu] ='\0';

	// printf("%s OUT\n", __func__);
}

/* affiche taille rho(nu, q, k) */
void display_taille_rho(int nu, int q, int k) {
	long int nb_var =0, nb_ctr =0;

	/* Nb var : les vecteurs de Z_q^nu + R */
	nb_var =(long int)puissance(q, nu) +1;

	/* Nb ctr : 
		- les paires (J, a) ou J sous-ensemble de taille k de [nu] et a in Z_q^k 
		- R
		- eventuellement : LB
	*/
	nb_ctr =2 +( binom(nu, k) * (long int)puissance(q, k) );	/* choix de J et de a */

	/* affichage */
	printf("rho(%d, %d, %d) comporte ", nu, q, k);
	affiche_nombre(nb_var);
	printf(" variables et ");
	affiche_nombre(nb_ctr);
	printf(" contraintes.\n");
}

/* affiche taille rho_E(nu, q, k) */
void display_taille_rho_E(int nu, int q, int k) {
	long int nb_var =0, nb_ctr =0;

	/* Nb var : les vecteurs de {0} x Z_q^(nu -1) + R */
	nb_var =(long int)puissance(q, nu -1) +1;

	/* Nb ctr : 
		- les paires (J, a) ou J sous-ensemble de taille k de [nu] et a in {0} x Z_q^(k -1) 
		- R
		- eventuellement : LB
	*/
	nb_ctr =binom(nu, k) * (long int)puissance(q, k -1) +2;

	/* affichage */
	printf("rho_E(%d, %d, %d) comporte ", nu, q, k);
	affiche_nombre(nb_var);
	printf(" variables et ");
	affiche_nombre(nb_ctr);
	printf(" contraintes.\n");
}

/* affiche taille gamma(q, p, k) */
void display_taille_gamma(int q, int p, int k) {
	long int nb_var_Q, nb_var_P =0, nb_ctr =0;
	int r;

	/* Nb var : les vecteurs de Z_q^q + les tels vecteurs prenant au plus p valeurs distinctes + R */
	nb_var_Q =(long int)puissance(q, q);

	/* Pour r de 1 a p : 1. partition de q en r sous-ensembles. 2. choix de r valeurs parmi q. 3. attribution des r valeurs aux sous-ensembles */
	for(r =1 ; r <= p ; r ++)
		nb_var_P =nb_var_P +(get_nbpart(q, r) * binom(q, r) * fact(r));	

	/* Nb ctr : 
		- les paires (J, a) ou J sous-ensemble de taille k de [q] et a in Z_q^k 
		- R
		- eventuellement : LB
	*/
	nb_ctr =binom(q, k) * (long int)puissance(q, k) +2;

	/* affichage */
	printf("gamma(%d, %d, %d) comporte ", q, p, k);
	affiche_nombre(nb_var_Q +nb_var_P +1);
	printf(" variables (dont %ld pour P) et ", nb_var_P);
	affiche_nombre(nb_ctr);
	printf(" contraintes.\n");
}

/* affiche taille gamma_E(q, p, k) */
void display_taille_gamma_E(int q, int p, int k) {
	long int nb_var_Q, nb_var_P =0, nb_ctr =0;
	int r;

	/* Nb var : les vecteurs de {0} x Z_q^(q -1) +les tels vecteurs prenant au plus p valeurs distinctes + R */
	nb_var_Q =(long int)puissance(q, q -1);

	/* Pour r de 1 a p : 
		1. partition de [q -1] en r sous-ensembles. 
		2. choix de r valeurs parmi q, *sinon le choix de r valeurs parmi {1 ,.., q -1} si r =p* 
		3. attribution des r valeurs aux sous-ensembles 
	*/
	for(r =1 ; r <= p -1 ; r ++) 
		nb_var_P =nb_var_P +(get_nbpart(q -1, r) * binom(q, r) * fact(r));

	nb_var_P =nb_var_P +(get_nbpart(q -1, p) * binom(q -1, p -1) * fact(p));

	/* Nb ctr : 
		- les paires (J, a) ou J sous-ensemble de taille k de [q] et a in {0} x Z_q^k 
		- R
		- eventuellement : LB
	*/
	nb_ctr =binom(q, k) * (long int)puissance(q, k -1) +2;

	/* affichage */
	printf("gamma_E(%d, %d, %d) comporte ", q, p, k);
	affiche_nombre(nb_var_Q +nb_var_P +1);
	printf(" variables (dont %ld pour P) et ", nb_var_P);
	affiche_nombre(nb_ctr);
	printf(" contraintes.\n");
}

/* DEBUG mais keep (sympa)
	void afftab(int* t, int n, int k);
	void afftab(int* t, int n, int k) {
		int i, j;

		for(i =1 ; i <= n ; i ++) {
			for(j =1 ; j <= i && j <= k ; j ++) {
				printf("%d ", t[k*(i -1) +(j -1)]);
			}

			printf("\n");
		}
	}
*/

/* Cette fonction renvoie le nombre de partitions en k >= 1 sous-ensembles non vides d'un ensemble de taille n >= 1 */
int get_nbpart(int n, int k) {
	int tab[n * k]; /* attention au declage : tan[O] est l'entree pour n =1 */
	int i, j;

	/* initialisation : 
		nbpart(i, 1) = nbpart(i, i) = 1 pour tout i >= 1
	*/
	for(i =1 ; i <= n ; i ++) {
		for(j =1 ; j <= i && j <= k ; j ++) {
			if(j == 1 || j == i)
				tab[k*(i -1) +(j -1)] =1;
			else
				tab[k*(i -1) +(j -1)] =0;
		}
	}

	// DEBUG mais keep (sympa)	// afftab(tab, n, k);

	/* calcul ligne par ligne de gaucha a droite selon la ligne du dessus et la case de gauche :
		nbpart(i, j) = nbpart(i -1, j -1) +j x nbpart(i -1, j) 
	*/
	for(i =3 ; i <= n ; i ++) {
		for(j =2 ; j <= i -1 && j <= k ; j ++) {
			tab[k*(i -1) +(j -1)] =tab[k*(i -2) +(j -2)] +j*tab[k*(i -2) +(j -1)];

			// DEBUG mais keep (sympa)	// afftab(tab, n, k);
			// getchar();
		}
	}

	/* on retourne le resultat */
	return tab[k*(n -1) +(k -1)];
}

static void affiche_nombre(long int x) {
	long int y =1;
	int prems =1;

	/* y = 1000^d où d tel que 1000^(d +1) > x >= 1000^d */
	while(y*1000 < x)
		y =y *1000;

	while(y >= 1) {
		if(prems == 1) {
			printf("%d", (int)(x/y));
			prems =0;
		}
		else {
			printf("%03d", (int)(x/y));
		}

		if(y > 1) /* derniere it° on ne met pas le . à la fin */
			printf(".");

		x =x%y;
		y =y/1000;
	}
}


long long unsigned puissance(unsigned base, unsigned exposant) {
	if (exposant == 0)
		return 1;
	if (exposant == 1)
		return base;
	long long unsigned a = puissance(base, exposant/2);
	return (exposant % 2) ? a*a*base : a*a;
}

