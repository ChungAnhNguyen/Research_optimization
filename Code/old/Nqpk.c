/* Depuis U, V pour N(q, p, k), construire N(q +1, p, k) ?

	EN COURS
*/

#include<stdlib.h>
#include<stdio.h>

#define K 3
#define P 3
#define Q_MAX 20

int binom(int n, int k);
int bmin(int q, int p, int k);

void nqkp_construire(int U[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int V[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int p, int k);
void nqkp_evaluer(int U[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int V[(Q_MAX -1) * (Q_MAX -1)][Q_MAX]);
void nqkp_afficher(int U[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int V[(Q_MAX -1) * (Q_MAX -1)][Q_MAX]);

int main() {
	int U[(Q_MAX -1) * (Q_MAX -1)][Q_MAX];
	int V[(Q_MAX -1) * (Q_MAX -1)][Q_MAX];

	nqkp_construire(U, V, P, K);
	nqkp_afficher(U, V);
	nqkp_evaluer(U, V);

	return EXIT_SUCCESS;
}

/* C_n^k = n * (n -1) *...* (n -k +1) / (k * (k -1) *...* 2) */
int binom(int n, int k) {
	int num = 1, den = 1, i;

	for (i = n ; i >= n -k +1 ; i --)
		num = num * i;

	for (i = k ; i >= 2 ; i --)
		den = den * i;

	return num/den;
}

/* Pour q >= p >= k, retourner i = argmin{C_q^(k -1) + C_(q -1)^(k -1), C_q^(p -1) + C_(q -1)^(p -1)} */
int bmin(int q, int p, int k) {
	int nbk, nbp;

	nbk = binom(q, k -1) + binom(q -1, k -1);
	nbp = binom(q, p -1) + binom(q -1, p -1);

	return (nbk <= nbp ? k : p);
}

void nqkp_construire(int U[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int V[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int p, int k) {
	int q, i, j, a, b, J, n;

	/* initialisation : N(p, p, k) : U = V = {(0 ,..., p -1)} */
	for (j = 0 ; j < p ; j ++) {
		U[0][j] = j;
		V[0][j] = j;
	}

	n = 1;

	/* construction récursive : passage de q à q +1 pour q de p à Q_MAX -1 */
	for (q = p ; q < Q_MAX ; q ++) {
		/* u^q_[n] = u^0_[n] et v^q_[n] = v^0_[n], sauf sur 1è ligne où v^q_1 = q */
		U[0][q] = U[0][0];
		V[0][q] = q;
		for (i = 1 ; i < n ; i ++) {
			U[i][q] = U[i][0];
			V[i][q] = V[i][0];
		}

		/* on place dans U et V, pour les colonnes 0 ... q, respectivement :
			(q ,..., q, a, q ,..., q, q) et (q ,..., q, a, q ,..., q, 0)		sur les lignes i = (q -2)^2 +a, a in Z_(q -1)
			(q ,..., q, q, q ,..., q, 0) et (q ,..., q, q, q ,..., q, q)		sur les lignes (q -2)^2 +q -1 ,..., (q -2)^2 +q -1 +q -2 */ // HERE
		b = bmin(q, p, k);
		for (J = 0 ; J < 2^q ; J ++) {
			for (j = 0 ; j < q ; j ++)
				vec[j] = (i/q^(q -r)) % q;
		}

		for (a = 0 ; a <= q -2 ; a ++) {
			i = (q -2)*(q -2) +a;
			for (j = 0 ; j <= q -1 ; j ++) {
				if (j == a) {
					U[i][j] = a;
					V[i][j] = a;
				}
				else if (j < q -1) {
					U[i][j] = q -1;
					V[i][j] = q -1;
				}
				else {
					U[i][j] = q -1;
					V[i][j] = 0;
				}
			}
		}
		for (i = (q -2)*(q -2) +q -1 ; i < (q -1)*(q -1) ; i ++) {
			for (j = 0 ; j <= q -2 ; j ++) {
				U[i][j] = q -1;
				V[i][j] = q -1;
			}

			U[i][q -1] = 0;
			V[i][q -1] = q -1;
		}
	}
}

void nqkp_evaluer(int U[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int V[(Q_MAX -1) * (Q_MAX -1)][Q_MAX]) {
	int q, r, s, a, b, i, ok, nb_U, nb_V;

	/* TODO : faire une fonction prenant q en paramètre  */
	for (q = P ; q <= Q_MAX ; q ++) {
		ok = 1;

		/* U_[(q -1)^2]^[Z_q] n'a bien que 2 valeurs par ligne ?		Ici : trivial */
		/* U_1^[Z_q]  =  (0 ,..., q -1) ? 							Ici : trivial */
		/* Pour tous r <> s in [q] et tous a <> b in Z_q, |{i in [(q -1)^2] : (u_i^r, u_i^s) = (a, b)}| = |{i in [(q -1)^2] : (u_i^r, u_i^s) = (a, b)}| ? */
		for (r = 0 ; r <= q -2 ; r++) {
			for (s = r +1 ; s <= q -1 && ok ; s++) {
				for (a = 0 ; a <= q -1 && ok ; a ++) {
					for (b = 0 ; b <= q -1 && ok ; b ++) {
						nb_U = 0;
						nb_V = 0;
						for (i = 0 ; i < (q -1) * (q -1) ; i ++) {
							if (U[i][r] == a && U[i][r] == b)
								nb_U ++;

							if (V[i][r] == a && V[i][r] == b)
								nb_V ++;
						}

						if (nb_U != nb_V)
							ok = 0;
					}
				}
			}
		}
		
		printf("Pour q = %d, ok == %d\n", q, ok);
	}
}

void nqkp_afficher(int U[(Q_MAX -1) * (Q_MAX -1)][Q_MAX], int V[(Q_MAX -1) * (Q_MAX -1)][Q_MAX]) {
	int i, j;

	for (i = 0 ; i < (Q_MAX -1) * (Q_MAX -1) ; i ++) {
		/* u_i */
		for (j = 0 ; j < Q_MAX ; j ++) {
			printf("%d ", U[i][j]);
		}

		printf("\t\t");

		/* v_i */
		for (j = 0 ; j < Q_MAX ; j ++) {
			printf("%d ", V[i][j]);
		}

		printf("\n");
	}
}

