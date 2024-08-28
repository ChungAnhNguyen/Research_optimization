/*	problèmes : infos sur les pbs résolus .... ou pas */
#include "vn.h"
#include <stdlib.h> // pour EXIT_SUCCESS
#include <stdio.h>	// pour printf(), getchar()
/* ________ Fonctions (declarations)
*/

/* sous-routine de display_bornes renvoyant 
        C_base^b *C_base^{a -h) -C_base^a *C_base^(b -h)
    pré-condition: a > b, base > 0
*/
int f(int base, int a, int b, int h);
/* Valeurs optimales connues pour gamma(q, p, k) */
char* gamma_get_opt(int q, int p, int k);
/* Eval bornes gamma */
void display_bornes();
/* Tqk */
void display_Tqk();
/* gamma(q, 3, 2) */
void display_gammaqp2();
/* bornes inf données par la construction récursive */
void display_BI_Tqk(int p, int k);
/* "généralisation" Tqk */
void display_Tqk_plus();
/* fait / à faire */
void todo_rho();
/* fait / à faire */
void todo_gamma();
/* une solution specifique */
void aff_gamma(int* P, int* Q, int R, int q);
/* ________ Programme de test
*/
int main () {
/*
	int R =16, q =5;
	int P[] ={93, 131, 312, 624, 656, 784, 793, 1033, 1431, 1512, 1814, 1862, 2317, 2331, 2369, 2498};
	int Q[] ={81, 781, 1381, 2306, 783, 1034, 262, 1562, 1812, 1864, 2367, 668, 2343, 194, 623, 2499};
	int R =4, q =6; // gamma(6, 4, 2)
	int P[] ={1919, 2076, 43337, 43602};
	int Q[] ={1865, 2130, 43391, 43548};
	aff_gamma(P, Q, R, q);
	getchar();
	// return EXIT_SUCCESS;
*/
	todo_gamma();
/*
	getchar();
	// return EXIT_SUCCESS;
	todo_rho();
	getchar();
	// return EXIT_SUCCESS;
//	display_Tqk();
//	display_Tqk_plus();
//  display_BI_Tqk(3, 2);
// display_gammaqp2();
    display_bornes();
*/
	return EXIT_SUCCESS;
}
/* ________ Fonctions (definitions)

*/

/* Valeurs optimales connues pour gamma(q, p, k) */
#define NBCAR 100
char* gamma_get_opt(int q, int p, int k) {
	char res[NBCAR] ="valeur optimale non connue";

	/* k = 2 */	
	if (k == 2 && p == 2 && q >= p +1 && q <= 7)
		sprintf(res, "%d/%d", 1, (q -1)*(q -1));
	else if (k == 2 && p == 3 && q >= p +1 && q <= 7)
		sprintf(res, "%d/%lld", 1, binom(q -1, 2));
	else if (k == 2 && p == 4 && q >= p +1 && q <= 7)
		sprintf(res, "%d/%d", (q % 2 ? 4 : 1), (q % 2 ? (q -2)*(q -2) : (q -2)*(q -2)/4));
	else if (k == 2 && p == 5 && q >= p +1 && q <= 7)
		sprintf(res, "%d/%lld", (q % 3 ? 3 : 1), (q % 3 ? binom(q -2, 2) : binom(q -2, 2)/3));
	else if (k == 2 && p == 6 && q == 7)
		sprintf(res, "%s", "9/16");
	/* k > 2 et p =k */	
	else if (k == 3 && p == k && q <= 7)
		sprintf(res, "%d/%lld", 1, (Tqk(q, k) +1)/2);
	else if (k == 4 && p == k && q <= 6)
		sprintf(res, "%d/%lld", 1, (Tqk(q, k) +1)/2);
	else if (k == 5 && p == k && q <= 6)
		sprintf(res, "%d/%lld", 1, (Tqk(q, k) +1)/2);
	/* autres resultats PV opt */	
	else if (k == 3 && p == 4 && q == 5)
		sprintf(res, "%s", "1/5");
	else if (k == 3 && p == 4 && q == 6)
		sprintf(res, "%s", "2/27");
	else if (k == 3 && p == 5 && q == 6)
		sprintf(res, "%s", "1/4");
		
	return res;
}

/* sous-routine de display_bornes renvoyant 
        C_base^b *C_base^{a -h) -C_base^a *C_base^(b -h)
    pré-condition: a > b, base > 0
*/
int f(int base, int a, int b, int h) {
    return binom(base, b) *binom(base, a -h) -binom(base, a) *binom(base, b -h);
}

/* Eval bornes gamma */
void display_bornes() {
	int q, p, k, qmax =12;
	
    for(k =2 ; k <= 2 ; k ++) {
        for(p =k ; p < qmax ; p ++) {
            for(q =p +1 ; q <= qmax ; q ++) {
            /* pour k =2, q, p, i fixés, et j variant de 0 à p -1, affichage du ratio apporté par la solution régulière où 
                    Psi contient $U^p$ et $U^j$, et Phi contient $U^q$ et $U^i$
            */
                printf(" ________ ratio pour (q, p, k) =(%d, %d, %d):\n", q, p, k);
	            long int delta, numopt =0, denopt =1;
                int i =(p % 2 ? (p +1)/2 : p/2), j =0;
                
                /* calcul en (j =0 et i =ceil(p/2))) */
                numopt =binom(q -2, p -2) *binom(q -2, i -1) -binom(q -2, p -1)*binom(q -2, i -2);
                denopt =numopt +binom(q, i) *binom(q -2, p -1);
                delta =pgcd(numopt, denopt);
                printf("- **** (i, j) = (%d, %d):\t", i, j);
                printf("borne %ld/%ld\t", numopt/delta, denopt/delta);
                printf("vs. val opt %s\n", gamma_get_opt(q, p, k));

                /* calcul pour les autres (i, j) */
                for(i =1 ; i < p ; i ++) {
                    for(j =0 ; j < i ; j ++) {
	                    long int num, den, fij, fpi, fpj;
                        int lq, li, lp, lj;

                        if(i == (p +(p % 2))/2 && j == 0)
                            continue;

                        fij =f(q -2, i, j, 1);
                        fpi =f(q -2, p, i, 1);
                        fpj =f(q -2, p, j, 1);
                        delta =pgcd(fij, fpi);
                        delta =pgcd(delta, fpj);
                        lp =fij/delta;
                        li =fpj *lp/fij;
                        lj =fpi *lp/fij;
                        lq =binom(q -2, p -2) *lp +binom(q -2, j -2) *lj -binom(q -2, i -2) *li;
                        num =lq;
                        den =num +binom(q, i) *li;

                        if(num * denopt >= numopt * den) {
                            numopt =num;
                            denopt =den;

                            delta =pgcd(numopt, denopt);
                            printf("- pour (i,j) = (%d, %d)\t", i, j);
                            printf("borne %ld/%ld\t", numopt/delta, denopt/delta);
                            printf("vs. val opt %s\n", gamma_get_opt(q, p, k));
                        }
                    }
                }

/*              val1 =Tqk(q -1, k -1);
                val2 =(binom(q -1, p -1) +Vqpk(q -1, p -1, k -2))/binom(q -k, p -k);
                printf("- passage de q -1 à q:\t\tval1 =%ld\tval2 =%ld\tval1/val2 =%lf\n", val1, val2, (double)val1/val2);  
*/
                /* gamma(q, p, k) : construction récursive */
            /*  num =1;
                for(den =1, int j =p ; j <= q -1 ; j ++)
                    den =den +(binom(j, p -1) +Vqpk(j, p -1, k -2))/binom(j -k +1, p -k);
            */
                /* gamma(q, p, k) : solution déduite (par complétion des lignes) de la construction récursive pour gamma(q -p +k, k, k) */
            /*  num =1;
                den =(Tqk(q -p +k, k) +1)/2;
            */
                /* gamma(q, p, 2) : solution régulière avec U1 dans phi / construction récursive */
            /*  num =(p -1);
                den =(q -1)*(q -p +1);
            */
                /* gamma(q, p, 2) : solution régulière avec U2 dans phi */
            /*  num =2*(p -2);
                den =(q -2)*(q -p +2);
            */
                /* gamma(q, p, 2) : solution régulière avec Ul dans phi */
            /*  int i = (p % 2 ? (p +1)/2 : p/2);
                num =binom(q -2, p -2)*binom(q -2, i -1) -binom(q -2, p -1)*binom(q -2, i -2);
                den =num +binom(q, i)*binom(q -2, p -1);
                int delta =pgcd(num, den);
                printf("pour (q, p, k) =(%d, %d, %d):\t", q, p, k);
                printf("borne %ld/%ld\t", num/delta, den/delta);
                printf("vs. val opt %s\n", gamma_get_opt(q, p, k));
            */
            }
        }
    }
}

/* Tqk */
void display_Tqk() {
	int k, q, val;
	for(k =2 ; k <= 6 ; k ++) {
		printf("____ T(q, k) pour k = %d:\n", k);
		for(q = k +1 ; q <= 10 ; q ++) {
			val =Tqk(q, k);
			printf("\tT(%d, %d) = %d\t=> gamma(%d, %d, %d) >= 1/%f\n", q, k, val, q, k, k, (Tqk(q, k) +1)/2.0);
		}
	}
}

/* gamma(q, 3, 2) */
void display_gammaqp2() {
	int qmax =7, p_phi =qmax, p, q, Qstar, R, delta, num, den, deltaBis;

	for(p =3 ; p <= p_phi ; p++) {
        printf("____ estimation hypothetique pour gamma(q, %d, 2):\n", p);
	    for(q =p ; q <= qmax ; q++) {
            printf("\tdonnees pour q =%d: ", q);
            Qstar = (q -2) *binom(q -2, p -2) -binom(q -2, p -1);
            R =Qstar +binom(q, 2) *binom(q -2, p -1);
            delta =pgcd(Qstar, R);
            printf("\tgamma(%d, %d, 2) >= %d/%d =%d/%d\n", q, p, Qstar, R, Qstar/delta, R/delta);
        }
    }

	/* ____ bornes inférieures connues pour (q, p, 2), q > p >= 2 */
	for(p =2 ; p <= p_phi ; p++) {
	    q =p;
	    Qstar =1;
	    R =1;
        printf("____ bornes inférieures connues pour gamma(q, %d, 2):\n", p);
        printf("\tgamma(%d, %d, 2) == %d/%d\n", q, p, Qstar, R);
        while(q +1 <= qmax) {
            /* calcul de Qstar et R */
            delta =pgcd(Qstar, binom(q -1, p -2));
            R =R * binom(q -1, p -2)/delta +Qstar *( binom(q, p -1) +binom(q -1, p -1) )/delta;
            Qstar = Qstar * binom(q -1, p -2)/delta;
            q ++;
            /* calcul alternatif du ratio obtenu */
            num =p -1;
            den =(q -1)*(q -p +1);
            deltaBis =pgcd(num, den);
            /* restitution */
            delta =pgcd(Qstar, R);
            printf("\tgamma(%d, %d, 2) >= %d/%d =%d/%d == %d/%d\n", q, p, Qstar, R, Qstar/delta, R/delta, num/deltaBis, den/deltaBis);
        }
    }
}

/* bornes inf données par la construction récursive */
void display_BI_Tqk(int p, int k) {
	int q =p +1, qmax =8, Qstar, R;

	/* récupération des données (connues) d'une paire optimale pour (p +1, p, k) */
	if(p == k && k >= 2) {
		Qstar =1;
		R =1;
	    q =p;
	}
	else if(p == 3 && k == 2) {
		Qstar =2;
		R =6;
	}
	else if(p == 4 && k == 3) {
		Qstar =3;
		R =15;
	}
/*	else if(p == 5 && k == 4) {
		Qstar =?;
		R =?;
	}*/
	else if(p == 4 && k == 2) {
		Qstar =8;
		R =18;
	}
	else if(p == 5 && k == 3) {
		Qstar =6;
		R =24;
	}
	else if(p == 7 && k == 4) {
		Qstar =2;
		R =6;
	}
	else if(p == 6 && k == 2) {
		Qstar =9;
		R =16;
	}
	else {
	    printf("paire (p =%d, k =%d) de paramètres non traitée.\n", p, k);
	    return;
	}
	
	/* bornes inférieures déduites pour (q, p, k), q > p */
	printf("____ bornes inférieures connues pour gamma(q, %d, %d), %d >= q >= %d:\n", p, k, qmax, q);
	printf("\tgamma(%d, %d, %d) == %d/%d\n", q, p, k, Qstar, R);
	for( ; q <= qmax ; q ++) {
		R = R +Qstar * Tqk(q, k -1);
	    printf("\tgamma(%d, %d, %d) >= %d/%d\n", q, p, k, Qstar, R);
	}
}

/* rapport obtenu pour les voisinages (nu -> q, d -> p, k -> k) */
void display_Tqk_plus() {
	int k, qbase, q, p, num, den;
	for(k =2 ; k <= 7 ; k ++) {
		for(qbase = k +1 ; qbase <= 8 ; qbase ++) {
			printf("____ Les gamma(%d +kappa, %d +kappa, %d):\n", qbase, k, k);
			for(q =qbase, p =k ; q <= 7 ; q++, p ++) {
				/* calcul pour (q, p, k) =(qbase +kappa, k +kappa, k) :
					num =C_(q -k)-^(p -k)
					den =[C_q^p (T(p, k -1) +1) - C_(q -k)^(p -k) (T(q, k -1) -1)]/2
				*/
				num = binom(q -k, p -k);
				den = binom(q, p) * (Tqk(p, k -1) +1);
				den = den -binom(q -k, p -k) * (Tqk(q, k -1) -1);
				den = den/2;
				printf("\tpour gamma(%d, %d, %d) on trouve : %d/%d\n", q, p, k, num, den);
			}
		}
	}
}

/* fait / à faire */
void todo_rho() {
	printf("________________________________ q =2 (on donne les infos pour rho_E uniquement):\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho_E(17, 2, 2);
	display_taille_rho_E(15, 2, 4);
	display_taille_rho_E(14, 2, 6);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho_E(14, 2, 2);
	display_taille_rho_E(12, 2, 4);
	display_taille_rho_E(12, 2, 6);
	printf("________ Prochaine résolution Min R à réaliser :\n");
	display_taille_rho_E(14, 2, 2);
	display_taille_rho_E(10, 2, 4);
	display_taille_rho_E(11, 2, 6);
	printf("________ Prochaine résolution R bin à réaliser :\n");
	display_taille_rho_E(14, 2, 2);
	display_taille_rho_E(10, 2, 4);
	display_taille_rho_E(11, 2, 6);
	printf("________________________________ q =3 et rho_E:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho_E(14, 3, 2);
	display_taille_rho_E(13, 3, 3);
	display_taille_rho_E(11, 3, 4);
	display_taille_rho_E(11, 3, 5);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho_E(10, 3, 2);		// EN COURS MACPRO -> NON : j'ai dû éteindre la machine ...
	display_taille_rho_E(8, 3, 3);
	display_taille_rho_E(8, 3, 4);
	display_taille_rho_E(8, 3, 5);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho_E(11, 3, 2);
	display_taille_rho_E(9, 3, 3);
	display_taille_rho_E(8, 3, 4);
	display_taille_rho_E(7, 3, 5);
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho_E(11, 3, 2);
	display_taille_rho_E(9, 3, 3);
	display_taille_rho_E(8, 3, 4);
	display_taille_rho_E(7, 3, 5);
	printf("________________________________ q =3 et rho:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho(14, 3, 2);
	display_taille_rho(12, 3, 3);
	display_taille_rho(10, 3, 4);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho(9, 3, 2);
	display_taille_rho(7, 3, 3);
	display_taille_rho(7, 3, 4);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho(8, 3, 2);
	display_taille_rho(7, 3, 3);
	display_taille_rho(7, 3, 4);
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho(9, 3, 2);
	display_taille_rho(7, 3, 3);
	display_taille_rho(7, 3, 4);
	printf("________________________________ q =4 et rho_E:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho_E(12, 4, 2);
	display_taille_rho_E(11, 4, 3);
	display_taille_rho_E(9, 4, 4);
	printf("____ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho_E(9, 4, 2);
	display_taille_rho_E(7, 4, 3);
	display_taille_rho_E(6, 4, 3);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho_E(9, 4, 2);
	display_taille_rho_E(7, 4, 3);
	display_taille_rho_E(6, 4, 3);
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho_E(8, 4, 2);
	display_taille_rho_E(7, 4, 3);
	display_taille_rho_E(6, 4, 3);
	printf("________________________________ q =4 et rho:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho(11, 4, 2);
	display_taille_rho(10, 4, 3);
	display_taille_rho(8, 4, 4);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho(7, 4, 2);
	display_taille_rho(7, 4, 3);
	display_taille_rho(6, 4, 3);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho(7, 4, 2);
	display_taille_rho(7, 4, 3);
	display_taille_rho(6, 4, 3);
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho(7, 4, 2);
	display_taille_rho(7, 4, 3);
	display_taille_rho(6, 4, 3);
	printf("________________________________ q =5 et rho_E:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho_E(10, 5, 2);
	display_taille_rho_E(9, 5, 3);
	display_taille_rho_E(8, 5, 4);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho_E(8, 5, 2);
	display_taille_rho_E(7, 5, 3);
	display_taille_rho_E(6, 5, 4);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho_E(8, 5, 2);
	display_taille_rho_E(7, 5, 3);
	display_taille_rho_E(6, 5, 4);
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho_E(8, 5, 2);
	display_taille_rho_E(7, 5, 3);	// EN COURS MACST
	display_taille_rho_E(6, 5, 4);
	printf("________________________________ q =5 et rho:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho(9, 5, 2);
	display_taille_rho(8, 5, 3);
	display_taille_rho(7, 5, 4);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho(7, 5, 2);
	display_taille_rho(7, 5, 3);
	display_taille_rho(7, 5, 4);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho(7, 5, 2);
	display_taille_rho(7, 5, 3);
	display_taille_rho(7, 5, 4);
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho(7, 5, 2);
	display_taille_rho(7, 5, 3);
	display_taille_rho(7, 5, 4);
	printf("________________________________ q =6 et rho_E:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho_E(10, 6, 2);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho_E(7, 6, 2);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho_E(6, 6, 2);	// EN COURS MACST
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho_E(7, 6, 2);
	printf("________________________________ q =6 et rho:\n");
	printf("________ Prochaines résolutions continues à réaliser (puis à démontrer optimale):\n");
	display_taille_rho(8, 6, 2);
	printf("________ Prochaines résolutions entières à réaliser :\n");
	display_taille_rho(5, 6, 2);
	printf("________ Prochaines résolutions Min R à réaliser :\n");
	display_taille_rho(5, 6, 2);
	printf("________ Prochaines résolutions R bin à réaliser :\n");
	display_taille_rho(5, 6, 2);
}
/* fait / à faire */
void todo_gamma() {
	printf("\n________________________ Cas dont il faut démontrer l'optimalité:\n");
	printf("____ Pour k donné, plus petit q pour lequel gamma(q, k, k) reste à traiter:\n");
	display_taille_gamma(8, 2, 2);
	display_taille_gamma(8, 3, 3);
	display_taille_gamma(6, 4, 4);		// NB glpk : ko après + de 4 jours ...
	display_taille_gamma(6, 5, 5);
	display_taille_gamma(7, 6, 6);
	printf("n____ Pour k et un écart k -p >0 donnés, plus petit q pour lequel gamma(q, p, k) reste à traiter:\n");
	display_taille_gamma(7, 6, 2);		// écart 4
	printf("____ Pour k donné, plus petit q pour lequel gamma_E(q, k, k) reste à traiter:\n");
	display_taille_gamma_E(7, 3, 3);
	display_taille_gamma_E(6, 4, 4);
	printf("n____ Pour k et un écart k -p >0 donnés, plus petit q pour lequel gamma_E(q, p, k) reste à traiter:\n");
	display_taille_gamma_E(8, 3, 2);	// écart 1
	display_taille_gamma_E(7, 5, 4);	// écart 1
	display_taille_gamma_E(8, 4, 2);	// écart 2
	display_taille_gamma_E(8, 5, 2);	// écart 3
	display_taille_gamma_E(7, 6, 3);	// écart 3
	display_taille_gamma_E(8, 6, 2);	// écart 4
	display_taille_gamma_E(8, 7, 2);	// écart 5
	printf("\n________________________ Cas à résoudre en continu:\n");
	printf("____ Pour k et un écart k -p >0 donnés, plus petit q pour lequel gamma(q, p, k) reste à traiter:\n");
	display_taille_gamma(7, 4, 3);		// écart 1
	display_taille_gamma(6, 5, 4);		// écart 1
	display_taille_gamma(7, 6, 5);		// écart 1
	display_taille_gamma(7, 5, 3);		// écart 2
	display_taille_gamma(7, 6, 4);		// écart 2
	display_taille_gamma(7, 6, 3);		// écart 3
	printf("____ Pour k donné, plus petit q pour lequel gamma_E(q, k, k) reste à traiter:\n");
	display_taille_gamma_E(9, 2, 2);	// Cplex fait kill ...
	display_taille_gamma_E(8, 3, 3);
	display_taille_gamma_E(8, 4, 4);
	display_taille_gamma_E(7, 5, 5);	// écart 
	printf("n____ Pour k et un écart k -p >0 donnés, plus petit q pour lequel gamma_E(q, p, k) reste à traiter:\n");
	display_taille_gamma_E(8, 4, 3);	// écart 1
	display_taille_gamma_E(8, 5, 3);	// écart 2
	display_taille_gamma_E(7, 6, 3);	// écart 3
	display_taille_gamma_E(8, 7, 3);	// écart 4
	display_taille_gamma_E(8, 5, 4);	// écart 1
	display_taille_gamma_E(7, 6, 4);	// écart 2
	display_taille_gamma_E(8, 7, 4);	// écart 3
	display_taille_gamma_E(7, 6, 5);
	display_taille_gamma_E(8, 7, 5);
	display_taille_gamma_E(7, 6, 6);
	display_taille_gamma_E(8, 7, 6);
	display_taille_gamma_E(8, 7, 7);
	printf("\n________________________ Cas resolus en continu pour lesquels une solution entiere reste a exhiber:\n");
	printf("____ gamma(q, p, k) ou p > k et q <= 7 (on donne tous les cas) :\n");
	display_taille_gamma(6, 5, 2);	// en cours Dives 3 (NB on peut se restreindre à p_phi =5)
	display_taille_gamma(7, 3, 2);
	display_taille_gamma(7, 4, 2);
	display_taille_gamma(7, 5, 2);
	display_taille_gamma(7, 6, 2);
	display_taille_gamma(6, 4, 3);	// (NB on peut se restreindre à p_phi =4)
	display_taille_gamma(6, 5, 3);	// NB on sait au - qu'il existe une solution vérifiant R =24 et R* = 6 (rapport 1/4)
	display_taille_gamma(6, 5, 4);	// en attente connaissance valeur optimale

	printf("____ gamma_E(q, p, k) ou q <= 8 (on donne tous les cas) :\n");
	display_taille_gamma_E(7, 4, 2);
	display_taille_gamma_E(7, 5, 2);
	display_taille_gamma_E(7, 6, 2);
	display_taille_gamma_E(8, 2, 2);		// Cplex 20 sort en erreur "out of memory"
	display_taille_gamma_E(8, 3, 2);
	display_taille_gamma_E(8, 4, 2);
	display_taille_gamma_E(8, 5, 2);
	display_taille_gamma_E(8, 6, 2);
	display_taille_gamma_E(8, 7, 2);
	display_taille_gamma_E(7, 3, 3);		// TODO
	display_taille_gamma_E(7, 4, 3);		// Min R entiers et binaire : cplex 12 sort en erreur "out of memory status" (int : sans sol tv ; bin : sol R =4096, gap 99,22%)
	display_taille_gamma_E(7, 5, 3);
	display_taille_gamma_E(6, 4, 4);
	display_taille_gamma_E(7, 4, 4);		// résolution continue: valeur 1.3964734051e-0.2, temps 736655.25 s.
	display_taille_gamma_E(7, 5, 4);
// q =8
}

/* une solution specifique */
void aff_gamma(int* P, int* Q, int R, int q) {
	int i;
	for (i = 0 ; i < R ; i ++) {
		affiche_v(P[i], q, q);
		printf("\t\t");
		affiche_v(Q[i], q, q);
		printf("\n");
	}
}
