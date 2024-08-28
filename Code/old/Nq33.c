/* Depuis U, V pour N(q, p, k), construire N(q +1, p, k) ? Cas : k = p = 3.

	- N(q, 3, 3) -> N(q +1, 3, 3) : si l'on note N le nombre de lignes de U et V, et en supposant v_1 = (0 ,..., q -1) (on ne considère qu'une occurrence de opt) :
		1. Faire u^q_[N] = u^0_[N] et v^q_[N] = (q, v^0_(2 ,..., N))
		2. Pour tous (j < h) in Z_q, insérer dans U la ligne (q ,..., q, j, q ,..., q, h, q ,..., q, q) et dans V la ligne (q ,..., q, j, q ,..., q, h, q ,..., q, 0)
		3. Pour tous j in Z_q, insérer q -2 fois dans U la ligne (q ,..., q, j, q ,..., q, 0) et dans V la ligne (q ,..., q, j, q ,..., q, 0)
		4. Insérer (q -2)^2 -C_(q -2)^2 = C_(q -1)^2 fois dans U la ligne (q ,..., q, q) et dans V la ligne (q ,..., q, 0)

	# 1. => pour j <> h in Z_q et i in 1, .., N, 						1 			occurrence 	(u_i^j, u_i^h, u_i^q) = (j, h, 0) au lieu de (j, h, q) dans V
	# 2. => pour j <> h in Z_q et i in N +1, ..., N +C_q^2, 				1 			occurrence 	(u_i^j, u_i^h, u_i^q) = (j, h, q) au lieu de (j, h, 0) dans V
																 		q -2 		occurrences 	(u_i^j, u_i^h, u_i^q) = (j, q, q) au lieu de (j, q, 0) dans V
																 		C_(q -2)^2 	occurrences (u_i^j, u_i^h, u_i^q) = (q, q, q) au lieu de (j, q, 0) dans V
	# 3. => pour j <> h in Z_q et i in N +C_q^2 +1 ,..., N +(q -2)q , 	q -2 		occurrences 	(u_i^j, u_i^h, u_i^q) = (j, q, 0) au lieu de (j, q, q) dans V
																 		(q -2)^2 	occurrences (u_i^j, u_i^h, u_i^q) = (q, q, 0) au lieu de (j, q, q) dans V

	=> N(q +1, 3, 3) 	<= N(q, 3, 3) +C_q^2 +q(q -2) +C_(q -1)^2
						 = N(q, 3, 3) +q(q -2) +(1/2)(q -1)(q + q -2)
						 = N(q, 3, 3) +(q^2 -2q) +(q -1)^2
						 = N(q, 3, 3) +2(q -1)^2 -1


	Ex. N(3, 3, 3) -> N(4, 3, 3) :

		Init.	0 1 2		|	0 1 2	

		1.		0 1 2 0		|	0 1 2 3

		2.		0 1 2 0		|	0 1 2 3
				0 1 3 3 		|	0 1 3 0
				0 3 2 3 		|	0 3 2 0
				3 1 2 3 		|	3 1 2 0

		3.		0 1 2 0		|	0 1 2 3
				0 1 3 3 		|	0 1 3 0
				0 3 2 3 		|	0 3 2 0
				3 1 2 3 		|	3 1 2 0
				0 3 3 0 		|	0 3 3 3
				3 1 3 0 		|	3 1 3 3
				3 3 2 0 		|	3 3 2 3

		4.		0 1 2 0		|	0 1 2 3
				0 1 3 3 		|	0 1 3 0
				0 3 2 3 		|	0 3 2 0
				3 1 2 3 		|	3 1 2 0
				0 3 3 0 		|	0 3 3 3
				3 1 3 0 		|	3 1 3 3
				3 3 2 0 		|	3 3 2 3
				3 3 3 3		|	3 3 3 0
*/


#include<stdlib.h>
#include<stdio.h>

#define K 3
#define P 3
#define Q_MAX 7
#define N_MAX (Q_MAX * Q_MAX * Q_MAX)

int Nq33_construire(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q);
void Nq33_q_plus_1(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int* N);
int Nq33_evaluer(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N);
void Nq33_afficher(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N);
static void Nq33_init(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX]);

int main() {
	int U[N_MAX][Q_MAX];
	int V[N_MAX][Q_MAX];
	int N, q = Q_MAX;
//	int est_valide = 0;

	Nq33_init(U, V);

	N = Nq33_construire(U, V, q);
//	Nq33_afficher(U, V, q, N);
//	est_valide = Nq33_evaluer(U, V, q, N);

//	printf("Pour q = %d : N == %d, ok = %d\n", q, N, est_valide);

	return EXIT_SUCCESS;
}

int Nq33_construire(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q) {
	int r, N, N_conj;

	/* initialisation */
	for (r = 0 ; r < P ; r ++) {
		U[0][r] = r;
		V[0][r] = r;
	}

	N = 1;
	printf("\t%s :: pour q = %d, N == %d, N_conj == %d, validité == %d\n", __func__, P, N, 1, Nq33_evaluer(U, V, P, N));

	/* construction récursive */
	for (r = P ; r < q ; r ++) {
		Nq33_q_plus_1(U, V, r, &N);
		N_conj = ((r+1) -2) * (r+1) * (2*(r+1) -5) / 3;
		printf("\t%s :: pour q = %d, N == %d, N_conj == %d, validité == %d\n", __func__, r +1, N, N_conj, Nq33_evaluer(U, V, r +1, N));
	}

	/* on retourne la taille de la paire (U, V) ainsi construite pour N(Q_MAX, 2, 2) */
	return N;
}

/* Passage d'une paire (U, V) vérifiant |U| = |V| = *N valide pour q >= 2 à une paire (U, V) vérifiant |U| = |V| = *N valide pour q +1 (*N est mis à jour) */
void Nq33_q_plus_1(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int* N) {
	int i, j, h, l;

	/* 1. u^q_[N] = u^0_[N] et v^q_[N] = (q, v^0_(2 ,..., N))
	*/
	U[0][q] = U[0][0];
	V[0][q] = q;

	for (i = 1 ; i < *N ; i ++) {
		U[i][q] = U[i][0];
		V[i][q] = V[i][0];
	}

	/* 2. pour tous j,h in Z_q vérifiant j < h, on insère les lignes :
		(q ,..., q, j, q ,..., q, h, q ,..., q, q) dans U
		(q ,..., q, j, q ,..., q, h, q ,..., q, 0) dans V
	*/
	for (j = 0 ; j <= q -2 ; j ++) {
		for (h = j +1 ; h <= q -1 ; h ++) {
			for (l = 0 ; l <= j -1 ; l ++) {
				U[*N][l] = q;
				V[*N][l] = q;
			}

			U[*N][j] = j;
			V[*N][j] = j;

			for (l = j +1 ; l <= h -1 ; l ++) {
				U[*N][l] = q;
				V[*N][l] = q;
			}

			U[*N][h] = h;
			V[*N][h] = h;

			for (l = h +1 ; l <= q -1 ; l ++) {
				U[*N][l] = q;
				V[*N][l] = q;
			}

			U[*N][q] = q;
			V[*N][q] = 0;

			*N = *N +1;
		}
	}

	/* 3. pour tout j in Z_q, on insère (q -2) fois les lignes :
		(q ,..., q, j, q ,..., q, 0) dans U
		(q ,..., q, j, q ,..., q, q) dans V
	*/
	for (j = 0 ; j <= q -1 ; j ++) {
		for (i = *N ; i < *N +q -2 ; i ++) {
			for (h = 0 ; h <= j -1 ; h ++) {
				U[i][h] = q;
				V[i][h] = q;
			}

			U[i][j] = j;
			V[i][j] = j;

			for (h = j +1 ; h < q ; h ++) {
				U[i][h] = q;
				V[i][h] = q;
			}

			U[i][q] = 0;
			V[i][q] = q;
		}

		*N = *N +(q -2);
	}

	/* 4. on insère (q -2)^2 -C_(q -2)^2 = C_(q -1)^2 fois les lignes :
		(q ,..., q, q) dans U
		(q ,..., q, 0) dans V
	*/
	for (i = *N ; i < *N +(q -1)*(q -2)/2 ; i ++) {
		for (j = 0 ; j < q ; j ++) {
			U[i][j] = q;
			V[i][j] = q;
		}

		U[i][q] = q;
		V[i][q] = 0;
	}

	*N = *N +(q -1)*(q -2)/2;
}

int Nq33_evaluer(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N) {
	int i, j, h, l, a, b, c, nb_U, nb_V;
	int ok = 1;

	/* U_[(q -1)^2]^[Z_q] n'a bien que 3 valeurs par ligne ?		Ici : trivial */
	/* U_1^[Z_q]  =  (0 ,..., q -1) ? 							Ici : trivial */
	for (j = 0 ; j <= q -3 ; j++)
		for (h = j +1 ; h <= q -2 && ok ; h++)
			for (l = h +1 ; l <= q -1 && ok ; l++)
				for (a = 0 ; a <= q -1 && ok ; a ++)
					for (b = 0 ; b <= q -1 && ok ; b ++)
						for (c = 0 ; c <= q -1 && ok ; c ++)
						{
							nb_U = 0;
							nb_V = 0;
							for (i = 0 ; i < N ; i ++)
							{
								if (U[i][j] == a && U[i][h] == b && U[i][l] == c)
									nb_U ++;

								if (V[i][j] == a && V[i][h] == b && V[i][l] == c)
									nb_V ++;
							}

							if (nb_U != nb_V)
								ok = 0;
						}
	
	return ok;
}

void Nq33_afficher(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N) {
	int i, j;

	for (i = 0 ; i < N ; i ++) {
		/* u_i */
		for (j = 0 ; j < q ; j ++) {
			printf("%d ", U[i][j]);
		}

		printf("\t\t");

		/* v_i */
		for (j = 0 ; j < q ; j ++) {
			printf("%d ", V[i][j]);
		}

		printf("\n");
	}
}

static void Nq33_init(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX]) {
	int i, j;

	for (i = 0 ; i < N_MAX ; i ++)
		for (j = 0 ; j < Q_MAX ; j ++) {
			U[i][j] = -1;
			V[i][j] = -1;
		}
}

