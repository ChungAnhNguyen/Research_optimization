#!/bin/sh

clear;	echo "___________________________ appels rho 3 2 2"

echo "_________ résolution continue:"
./tableaux-prog 0 0 3 2 			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 0 0 3 2 2			;read;clear	

echo "_________ résolution entière:"
./tableaux-prog 0 1 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 0 1 3 2 2 1/4		;echo " "		#VAL >=1/4
./tableaux-prog 0 1 3 2 2 0.25 8	;echo " "		#VAL >=1/4, R =8
./tableaux-prog 0 1 3 2 2 1/3.0		;echo " "		#VAL >=1/3
./tableaux-prog 0 1 3 2 2 0.25 5	;read;clear		#VAL >=1/4, R =5

echo "_________ résolution PV opt:"
./tableaux-prog 0 2 3 2 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 0 2 3 2 2 1 4		;echo " "		#VAL >1/4
./tableaux-prog 0 2 3 2 2 1 6		;read;clear		#VAL >1/6

echo "_________ résolution binaire:"
./tableaux-prog 0 3 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 0 3 3 2 2			;echo " "
./tableaux-prog 0 3 3 2 2 4			;echo " "		#R =4
./tableaux-prog 0 3 3 2 2 8			;echo " "		#R =8
./tableaux-prog 0 3 3 2 2 5			;read;clear		#R =5

echo "_________ résolution fréquence max:"
./tableaux-prog 0 4 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 0 4 3 2 2 4			;echo " "		#R =4
./tableaux-prog 0 4 3 2 2 8			;echo " "		#R =8		
./tableaux-prog 0 4 3 2 2 5			;read;clear		#R =5

echo "_________ résolution min R:"
./tableaux-prog 0 5 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 0 5 3 2 2			;echo " "		#
./tableaux-prog 0 5 3 2 2 4			;echo " "		#R =4
./tableaux-prog 0 5 3 2 2 8			;echo " "		#R =8
./tableaux-prog 0 5 3 2 2 5			;read			#R =5

clear;	echo "___________________________ appels rho_E 3 2 2"

echo "_________ résolution continue:"
./tableaux-prog 1 0 3 2 			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 1 0 3 2 2			;read;clear	

echo "_________ résolution entière:"
./tableaux-prog 1 1 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 1 1 3 2 2 1/4		;echo " "		#VAL >=1/4
./tableaux-prog 1 1 3 2 2 0.25 8	;echo " "		#VAL >=1/4, R =8
./tableaux-prog 1 1 3 2 2 1/3.0		;echo " "		#VAL >=1/3
./tableaux-prog 1 1 3 2 2 0.25 5	;read;clear		#VAL >=1/4, R =5

echo "_________ résolution PV opt:"
./tableaux-prog 1 2 3 2 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 1 2 3 2 2 1 4		;echo " "		#VAL >1/4
./tableaux-prog 1 2 3 2 2 1 6		;read;clear		#VAL >1/6

echo "_________ résolution binaire:"
./tableaux-prog 1 3 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 1 3 3 2 2			;echo " "
./tableaux-prog 1 3 3 2 2 4			;echo " "		#R =4
./tableaux-prog 1 3 3 2 2 8			;echo " "		#R =8
./tableaux-prog 1 3 3 2 2 5			;read;clear		#R =5

echo "_________ résolution fréquence max:"
./tableaux-prog 1 4 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 1 4 3 2 2 4			;echo " "		#R =4
./tableaux-prog 1 4 3 2 2 8			;echo " "		#R =8		
./tableaux-prog 1 4 3 2 2 5			;read;clear		#R =5

echo "_________ résolution min R:"
./tableaux-prog 1 5 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 1 5 3 2 2			;echo " "		#
./tableaux-prog 1 5 3 2 2 4			;echo " "		#R =4
./tableaux-prog 1 5 3 2 2 8			;echo " "		#R =8
./tableaux-prog 1 5 3 2 2 5			;read			#R =5
