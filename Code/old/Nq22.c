/* Depuis U, V pour N(q, p, k), construire N(q +1, p, k) ? Cas k = p = 2

	N(q, 2, 2) -> N(q +1, 2, 2) : si l'on note N le nombre de lignes de U et V, et en supposant v_1 = (0 ,..., q -1) (on ne considère qu'une occurrence de opt) :
		1. Faire u^q_[N] = u^0_[N] et v^q_[N] = (q, v^0_(2 ,..., N))
			=> pour a, b, c in Z_q: 
				# i : (u_i^a, u_i^q) = (b, c)  == 	# i : (u_i^a, u_i^0) = (b, c) si a <> 0
													# i : u_i^0 = b si a = 0 et b = c
													0	sinon
				# i : (v_i^a, v_i^q) = (b, c)  == 	# i <> 1 : TODO
		2. Pour tout j in Z_q, insérer dans U la ligne (q ,..., q, j, q ,..., q, q) et dans V la ligne (q ,..., q, j, q ,..., q, 0)
		3. Insérer q -1 fois dans U la ligne (q ,..., q, 0) et dans V la ligne (q ,..., q, q)

		=> N(q +1, 2, 2) <= N(q, 2, 2) +(2q -1)
		=> N(q, 2, 2) <= N(2, 2, 2) + sum_{r = 2}^{q -1} (2r -1) = 1 +2*((q -1)q/2 -1) -(q -2) = (q -1)q -q +1 = (q -1)^2 

		NB : x occurrences de opt => même manip. Si l'on note N_x(q, p, k) la taille minimale d'une paire de familles (U, V) pour x occurrences de opt, alors :
		   N_x(q +1, 2, 2) <= N_x(q, 2, 2) +(2q -1)*x
		=> N_x(q, 2, 2) <= N_x(2, 2, 2) + x * sum_{r = 2}^{q -1} (2r -1) = x * ( 1 +sum_{r = 2}^{q -1} (2r -1)) = x * (q -1)^2 
*/

#include<stdlib.h>
#include<stdio.h>

#define K 2
#define P 2
#define Q_MAX 7
#define N_MAX ((Q_MAX -1) * (Q_MAX -1))

int Nq22_construire(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q);
void Nq22_q_plus_1(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int* N);
int Nq22_evaluer(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N);
void Nq22_afficher(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N);
static void Nq22_init(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX]);

int main() {
	int U[N_MAX][Q_MAX];
	int V[N_MAX][Q_MAX];
	int N, q = Q_MAX, est_valide = 0;

	Nq22_init(U, V);

	N = Nq22_construire(U, V, q);
	Nq22_afficher(U, V, q, N);
	est_valide = Nq22_evaluer(U, V, q, N);

	printf("Pour q = %d : N == %d, ok = %d\n", q, N, est_valide);

	return EXIT_SUCCESS;
}

int Nq22_construire(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q) {
	int r, N;

	/* initialisation */
	for (r = 0 ; r < P ; r ++) {
		U[0][r] = r;
		V[0][r] = r;
	}

	N = 1;

	/* construction récursive */
	for (r = P ; r < q ; r ++) {
		Nq22_q_plus_1(U, V, r, &N);
		printf("\t%s :: pour q = %d, N == %d\n", __func__, r +1, N);
	}

	/* on retourne la taille de la paire (U, V) ainsi construite pour N(Q_MAX, 2, 2) */
	return N;
}

/* Passage d'une paire (U, V) vérifiant |U| = |V| = *N valide pour q >= 2 à une paire (U, V) vérifiant |U| = |V| = *N valide pour q +1 (*N est mis à jour) */
void Nq22_q_plus_1(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int* N) {
	int i, j, h;

	/* 1. u^q_[N] = u^0_[N] et v^q_[N] = (q, v^0_(2 ,..., N))
	*/
	U[0][q] = U[0][0];
	V[0][q] = q;

	for (i = 1 ; i < *N ; i ++) {
		U[i][q] = U[i][0];
		V[i][q] = V[i][0];
	}

	/* 2. pour tout j in Z_q, on insère les lignes :
		(q ,..., q, j, q ,..., q, q) dans U
		(q ,..., q, j, q ,..., q, 0) dans V
	*/
	for (j = 0 ; j <= q -1 ; j ++) {
		for (h = 0 ; h <= j -1 ; h ++) {
			U[*N][h] = q;
			V[*N][h] = q;
		}

		U[*N][j] = j;
		V[*N][j] = j;

		for (h = j +1 ; h < q ; h ++) {
			U[*N][h] = q;
			V[*N][h] = q;
		}

		U[*N][q] = q;
		V[*N][q] = 0;

		*N = *N +1;
	}

	/* 3. on insère (q -1) fois les lignes :
		(q ,..., q, 0) dans U
		(q ,..., q, q) dans V
	*/
	for (i = *N ; i < *N +(q -1) ; i ++) {
		for (j = 0 ; j < q ; j ++) {
			U[i][j] = q;
			V[i][j] = q;
		}

		U[i][q] = 0;
		V[i][q] = q;
	}

	*N = *N +(q -1);
}

int Nq22_evaluer(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N) {
	int i, j, h, a, b, nb_U, nb_V;
	int ok = 1;

	/* U_[(q -1)^2]^[Z_q] n'a bien que 2 valeurs par ligne ?		Ici : trivial */
	/* U_1^[Z_q]  =  (0 ,..., q -1) ? 							Ici : trivial */
	/* Pour tous r <> s in [q] et tous a <> b in Z_q, |{i in [(q -1)^2] : (u_i^r, u_i^s) = (a, b)}| = |{i in [(q -1)^2] : (u_i^r, u_i^s) = (a, b)}| ? */
	for (j = 0 ; j <= q -2 ; j++)
		for (h = j +1 ; h <= q -1 && ok ; h++)
			for (a = 0 ; a <= q -1 && ok ; a ++)
				for (b = 0 ; b <= q -1 && ok ; b ++) {
					nb_U = 0;
					nb_V = 0;
					for (i = 0 ; i < N ; i ++) {
						if (U[i][j] == a && U[i][h] == b)
							nb_U ++;

						if (V[i][j] == a && V[i][h] == b)
							nb_V ++;
					}

					if (nb_U != nb_V)
						ok = 0;
				}
	
	return ok;
}

void Nq22_afficher(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX], int q, int N) {
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

static void Nq22_init(int U[N_MAX][Q_MAX], int V[N_MAX][Q_MAX]) {
	int i, j;

	for (i = 0 ; i < N_MAX ; i ++)
		for (j = 0 ; j < Q_MAX ; j ++) {
			U[i][j] = -1;
			V[i][j] = -1;
		}
}
