/* 	gamma_E(8, 3, 2) : 
	- valeur relaxation 0.2301980198 (la précision glpk est de 10 chiffres sig.)
	- 93/404  = (3 * 31)/(4 * 101) est un bon candidat : valeur 0.2301980198019802 (scilab ne va pas au-dessus de 16 chiffres sig.)

	gamma_E(6, 3, 3) :
	- valeur relaxation 0.05478782106 (la précision glpk est de 10 chiffres sig.)
	- 38425/701342  = (25 * 29 *53)/(2 * 173 * 2027) est un bon candidat : valeur 0.0547878210630477 (scilab ne va pas au-dessus de 15 chiffres sig.)
	- il existe une solution entière de valeur 153700/2805368 =38425/701342
	- pour ce que cela vaut, il n'existe pas de solution entière de valeur vérifiant Q(0, 1 ,.., q -1) * 2805368 >= R * 153700 +1

	gamma_E(6, 4, 4) :
	- valeur relaxation 0.03159029059 (la précision glpk est de 10 chiffres sig.)
	- 76895/2434134  = (25 * 29 *53)/(2 * 173 * 2027) n'est pas un bon candidat : valeur 0.031590290427725 (scilab ne va pas au-dessus de 15 chiffres sig.)

	gamma_E(8, 6, 2) : 0.5714285714 -> 4/7
*/
/*	decomposition en nombres premiers

	Pb : valeur fractionnaire vs. la précision glpk est de 10 chiffres significatifs 

	// x = 0.05478782106;	// target i =701342 on trouve i = 151311 avec eps =0.00001, 14328 avec eps =0.00001 etc. 

	// gamma_E(8, 3, 2) : x = 2.301980198e-1, 	i = 404 (eps =0.00001)	-> 93/404
	// gamma_E(6, 3, 3) : x = 5.478782106e-2, 	i = 151311 (eps =0.00001) -> 310799 (eps =0.000004) -> 1091885 (eps =0.0000021) -> 59822/1091885 ?
	//	-> on a exhibé une solution +petite démontrée optimale ; 1 -modulo(i*x, 1) = 0.00000213748
	//	-> je recommence avec 1 -mod > eps ET mod > eps ; i = 79744 -> 390543 (je ne parviens pas à retrouver le dénominateur i = 701342)
	/	Conclusion : la démarche est ko
	// gamma_E(7, 3, 3) : x = 3.428432865e-2, 	i = 4171 ; 61865 ; 276307 ; 383528 ; 874277
	// gamma_E(6, 4, 4) : x = 3.159029059e-2, 	i = 238491 ; 256123 ; 529878 ; 1077388 ; 1624898 -> 51331/1624898 ?
	// gamma_E(7, 5, 4) : x = 5.8898096000e-2,	i = 42599 ; 117542	-> 6923/117542 ?
	x =5.478782106e-2;
	eps =0.00000213749;
	i =390543 ;
	while((1 -modulo(i*x, 1) > eps) & (modulo(i*x, 1) > eps))
		i = i +1;
	end
	printf("x =%.11f et %d x = %.11f\n", x, i, i*x);
	printf("\tmodulo(i*x, 1) = %.11f, 1 -modulo(i*x, 1) = %.11f, eps = %.11f\n", modulo(i*x, 1), 1 -modulo(i*x, 1), eps);
*/

#include <stdlib.h>	// pour EXIT_SUCCESS
#include <stdio.h>	// pour printf()
#include <math.h>	// pour floor(), ceil()

/* ________ Fonctions (declaration)
*/

/* test nombre premier
*/
int est_premier(int p);

/* decomposition en facteurs premiers
*/
void decomposition(int x);

/* recherche multiple
*/
void mult(double x);

/* recherche fractions
*/
void frac(double x, long int d_max);

/* ________ Programme
*/
int main () {
	double x;
	long int n;
	long int d;
/*
	printf("\n__ gamma_E(7, 5, 5) :\n");
	x =1.2817776230e-2;
	frac(x, 10000000);

	return 0;
*/	
	printf("\n__ gamma_E(6, 4, 4) :\n");
	x =3.159029059e-2;	
	n =8355;
	d =264480;
	decomposition(n);
	decomposition(d);
	printf("x = %.9e, %ld/%ld = %.13e arrondi %.9e\n", x, n, d, (double)n/d, (double)n/d);

	printf("\n__ gamma_E(7, 3, 3) :\n");
	x =3.428432865e-2;
	n =3676;
	d =107221;
	decomposition(n);
	decomposition(d);
	printf("x = %.9e, %ld/%ld = %.13e arrondi %.9e\n", x, n, d, (double)n/d, (double)n/d);

	printf("\n__ gamma_E(7, 5, 4) :\n");
	x =5.8898096000e-2;
	n =6923;
	d =117542;
	decomposition(n);
	decomposition(d);
	printf("x = %.9e, %ld/%ld = %.13e arrondi %.9e\n", x, n, d, (double)n/d, (double)n/d);

/*
	printf("\n__ rho(10, 4, 2) :\n");
	x =2.901785714e-2;
	n =13;
	d =448;
	decomposition(n);
	decomposition(d);
	printf("x = %.9e, %ld/%ld = %.13e arrondi %.9e\n", x, n, d, (double)n/d, (double)n/d);

	printf("\n__ rho_E(10, 4, 2) :\n");
	x =8.783783784e-2;
	n =13;
	d =148;
	decomposition(n);
	decomposition(d);
	printf("x = %.9e, %ld/%ld = %.13e arrondi %.9e\n", x, n, d, (double)n/d, (double)n/d);
*/
	return EXIT_SUCCESS;
}

/* ________ Fonctions (definition)
*/

/* test nombre premier
*/
int est_premier(int p) {
	int res = 1, d;

	for (d = 2 ; d * d <= p && res ; d ++) {
		if (p % d == 0) {
			res = 0;
		}
	}

	return res;
}

/* decomposition en facteurs premiers
*/
void decomposition(int x) {
	int p;

	printf("%d = ", x);

	p = 2;
	while (x > 1) {
		if (x % p == 0) {
			printf("%d x ", p);
			x = x /p;
		}
		else {
			do {
				p ++;
			}
			while (! est_premier(p));
		}
	}

	printf("1\n");
}

/* recherche multiple
*/
#define P_MAX 20000
void mult(double x) {
	long int p;

	for (p =2 ; p <= P_MAX ; p ++) {
		if(est_premier(p)) {
			printf("%ld * %.15f = %.15f\n", p, x, p * x);
			getchar();
		}
	}
}

/* recherche fractions
*/
void frac(double x, long int d_max) {
	long int n, d, l =floor(1/x), u =ceil(1/x);

	printf("%s : x = %lf, 1/x = %lf, l = %ld, u = %ld\n", __func__, x, 1/x, l, u);

	for (d =2 ; d <= d_max ; d ++) {
		// printf("\td = %ld : recherche pour n de ceil(%ld/%ld) = %ld a floor(%ld/%ld) = %ld\n", d, d, u, (long int)ceil((double)d/u), d, l, d/l);
		for(n =ceil((double)d/u) ; n <= d/l && n >= 1 && n <= d ; n ++) {
			if(x == (double)n/d) {
				printf("%lf = %ld / %ld\n", x, n, d);
				return;
			}
		}
	}
}

