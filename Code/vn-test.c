/*	vn (manipulation des entiers pour représenter des vecteurs de Z_q^nu / des sous-ensembles de [nu]) : fichier test
*/

#include "vn.h"

#include <stdlib.h>	// pour EXIT_SUCCESS
#include <stdio.h>	// pour printf(), getchar()
#include <math.h>	// pour pow()

/* ________ Fonctions de test (declarations)
*/

/* test : 
	- mod()
	- binom()	
	- Tqk()
*/
void test_math();

/* test : 
	- get_v_r() (par fonction d'affichage affiche_v)
	- itoa()
	- get_nb_val()	
*/
void test_decompo();

/* test : 
	- get_v_012etc()	
*/
void test_v_012etc();

/* test : 
	- get_card()	
*/
void test_card();

/* test : 
	- is_vJ_equal_to_a()
*/
void test_vJ_is_a();

/* test : 
	- get_vJ()
*/
void test_get_vJ();

/* test : 
	- get_nbpart()
*/
void test_get_nbpart();

/* test : 
	- display_taille_rho()
	- display_taille_rho_E()
	- display_taille_gamma()
	- display_taille_gamma_E()
*/
void test_taille_pbs();

/* ________ Programme de test
*/

int main () {
	test_decompo();
	getchar();
	
	test_taille_pbs();
	getchar();

	test_get_nbpart();
	getchar();

	test_math();
	getchar();

	test_v_012etc();
	getchar();

	test_card();
	getchar();

	test_vJ_is_a();
	getchar();

	test_get_vJ();

	return EXIT_SUCCESS;
}

/* ________ Fonctions de test (definitions)
*/

/* test : 
	- mod()
	- binom()	
*/
void test_math() {
	int q =5, a, k =3, p, q_max =10;

	printf("%s: q = %d et k = %d\n", __func__, q, k);

	/* modulo */
	for(a =-3*q -1 ; a <= 3*q +1 ; a ++) {
		printf("\t%d mod %d = %d\n", a, q, mod(a, q));
	}

	/* coefficient binomial */
	for(a =2 ; a <= q ; a ++) {
		printf("\tC_%d^%d = %d\n", q, a, binom(q, a));
	}

	/* Tqk */
	for(k =2 ; k <= q ; k ++) {
		printf("\tT(q, %d) :\n", k);
		for(p =k +1 ; p <= q_max ; p ++) {
			printf("\t\tT(%d, %d) = %d\t(T(%d, %d) +1)/2 = %f\n", p, k, Tqk(p, k), p, k, (Tqk(p, k) +1)/2.0);
		}
	}
}

/* test : 
	- get_v_r() (par fonction d'affichage affiche_v)
	- itoa()
	- get_nb_val()	
*/
void test_decompo() {
	int v, q =3, nu =3;
	char tab[100];

	printf("%s: q = %d et nu = %d\n", __func__, q, nu);

	/* parcours des nombres de Z_(q^nu) : leur decomposition en base q, le nombre de valeurs distinctes prises par les composantes de cette décomposition */
	for (v = 0 ; v <= pow(q, nu) -1 ; v ++) {
		printf("v = %d :\t\t", v);
		affiche_v(v, q, nu);
		itoa(v, q, nu, tab);
		printf(" (verif itoa : %s)\t\tnb val differentes: %d\n", tab, get_nb_val(v, q, nu));
	}
}

/* test : 
	- get_v_012etc()	
*/
void test_v_012etc() {
	int v, q =3;

	printf("%s: q = %d\n", __func__, q);

	/* parcours des nombres de Z_(q^q) : leur decomposition en base q, cette decomposition est-elle (0, 1 ,.., q -1) ? */
	for (v = 0 ; v <= pow(q, q) -1 ; v ++) {
		printf("v = %d :\t\t", v);
		affiche_v(v, q, q);
		printf("\t\test (0, 1 ,.., %d) ? %s\n", q -1, (v == get_v_012etc(q) ? "OUI" : "NON"));
	}
}

/* test : 
	- get_card()	
*/
void test_card() {
	int J, nu =4;

	printf("%s: nu = %d\n", __func__, nu);

	/* parcours des nombres de Z_(2^nu) : leur decomposition binaire, le nombre de bits non nuls de cette decomposition */
	for (J = 0 ; J <= pow(2, nu) -1 ; J ++) {
		printf("J = %d :\t\t", J);
		affiche_v(J, 2, nu);
		printf("\t\tcardinalite: %d\n", get_card(J, nu));
	}
}

/* test : 
	- is_vJ_equal_to_a()
*/
void test_vJ_is_a() {
	int v, J, a, q =3, nu =3, t =2;

	printf("%s: q = %d, nu = %d, t = %d : vecteurs v verifiant v_J = a\n", __func__, q, nu, t);

	/* E_q = 0 : parcours des sous-ensembles J de cardinalite t de [nu] */
	for (J = 0 ; J <= pow(2, nu) -1 ; J ++) {
		if(get_card(J, nu) == t) {
			printf("________ ensemble J (%d)\t", J);
			affiche_v(J, 2, nu);
			printf("\tde cardinalite: %d\n", get_card(J, nu));

			/* parcours des nombres a de (Z_q)^t */
			for (a = 0 ; a <= pow(q, t) -1 ; a ++) {
				printf("\t____ vecteur a (%d)\t", a);
				affiche_v(a, q, t);

				/* parcours des nombres v de (Z_q)^nu */
				for (v = 0 ; v <= pow(q, nu) -1 ; v ++) {
					if(is_vJ_equal_to_a(v, J, a, q, nu, t, 0) == 1) {
						printf("\n\t\t __ vecteur v (%d)\t", v);
						affiche_v(v, q, nu);
					}
				}

				printf("\n");
			}
		}
	}
	getchar();

	printf("%s: q = %d, nu = %d, t = %d : vecteurs v verifiant v_J in {a, a +(1 ,.., 1) ,.., a +(q -1 ,.., q -1)}\n", __func__, q, nu, t);

	/* E_q = 1 : parcours des sous-ensembles J de cardinalite t de [nu] */
	for (J = 0 ; J <= pow(2, nu) -1 ; J ++) {
		if(get_card(J, nu) == t) {
			printf("________ ensemble J (%d)\t", J);
			affiche_v(J, 2, nu);
			printf("\tde cardinalite: %d\n", get_card(J, nu));

			/* parcours des nombres a de (Z_q)^t */
			for (a = 0 ; a <= pow(q, t) -1 ; a ++) {
				printf("\t____ vecteur a (%d)\t", a);
				affiche_v(a, q, t);

				/* parcours des nombres v de (Z_q)^nu */
				for (v = 0 ; v <= pow(q, nu) -1 ; v ++) {
					if(is_vJ_equal_to_a(v, J, a, q, nu, t, 1) == 1) {
						printf("\n\t\t __ vecteur v (%d)\t", v);
						affiche_v(v, q, nu);
					}
				}

				printf("\n");
			}
		}
	}
}

/* test : 
	- get_vJ()
*/
void test_get_vJ() {
	int v, J, vJ, vJ0, q =3, nu =3, t =2;

	printf("%s: q = %d, nu = %d, t = %d : v_J pour v in Z_q^nu et J sous-ensemble de taille t de [nu]\n", __func__, q, nu, t);

	/* parcours des nombres v de (Z_q)^nu */
	for (v = 0 ; v <= pow(q, nu) -1 ; v ++) {
		printf("\t________ vecteur v (num %d) = ", v);
		affiche_v(v, q, nu);
		printf("\n");

		/* parcours des sous-ensembles J de cardinalite t de [nu] */
		for (J = 0 ; J <= pow(2, nu) -1 ; J ++) {
			if(get_card(J, nu) == t) {
				printf("\t\tsous-ensemble J (num %d) = ", J);
				affiche_v(J, 2, nu);

				/* recuperation de v_J et v_J -v_{j_1} */
				vJ = get_v_J(v, J, q, nu, t, 0);
				printf("\tv_J (num %d) = ", vJ);
				affiche_v(vJ, q, t);

				vJ0 = get_v_J(v, J, q, nu, t, 1);
				printf("\tnormalise (num %d) = ", vJ0);
				affiche_v(vJ0, q, t);

				printf("\n");
			}
		}

		getchar();
	}
}


/* test : 
	- get_nbpart()
*/
void test_get_nbpart() {
	int n =10, k, nb;

	printf("%s: n = %d\n", __func__, n);

	for (k = 1 ; k <= n ; k ++) {
		nb = get_nbpart(n, k);
		printf("\til existe %d partitions d'un ensemble de taille %d en %d sous-ensembles non vides\n", nb, n, k);
	}
}

/* test : 
	- display_taille_rho()
	- display_taille_rho_E()
	- display_taille_gamma()
	- display_taille_gamma_E()
*/
void test_taille_pbs() {
	int nu, q, k, p;

	/* test rho */
	q =3;
	k =3;
	printf("%s: test rho et rho_E pour q = %d et k = %d\n", __func__, q, k);

	for (nu = k +1 ; nu <= 10 ; nu ++) {
		display_taille_rho(nu, q, k);
		display_taille_rho_E(nu, q, k);
	}
	getchar();

	/* test gamma */
	k =2;
	p =k;
	printf("%s: test gamma pour k = %d et p = %d\n", __func__, k, p);

	for (q = p +1 ; q <= 10 ; q ++) {
		display_taille_gamma(q, p, k);
		display_taille_gamma_E(q, p, k);
	}

	printf("%s: test gamma_E(8, 7, 2)\n", __func__);
	display_taille_gamma_E(8, 7, 2);
}
