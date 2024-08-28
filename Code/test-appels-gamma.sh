#!/bin/sh

clear;	echo "___________________________ appels gamma 3 2 2"

echo "_________ résolution continue:"
./tableaux-prog 2 0 3 2 			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 2 0 3 2 2			;echo " "	
./tableaux-prog 2 0 3 2 2 2			;echo " "		#pmax =2
./tableaux-prog 2 0 3 2 2 1			;read;clear		#pmax =1	

echo "_________ résolution entière:"
./tableaux-prog 2 1 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 2 1 3 2 2 1/4		;echo " "		#VAL >=1/4
./tableaux-prog 2 1 3 2 2 1/3.0 	;echo " "		#VAL >=1/3
./tableaux-prog 2 1 3 2 2 0.25 5	;echo " "		#VAL >=1/4, R >=5
./tableaux-prog 2 1 3 2 2 0.25 4 2	;echo " "		#VAL >=1/4, R >=4, pmax =2
./tableaux-prog 2 1 3 2 2 0.25 0 1	;read;clear;	#VAL >=1/4, pmax =1

echo "_________ résolution PV opt:"
./tableaux-prog 2 2 3 2 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 2 2 3 2 2 1 6		;echo " "		#VAL >1/6
./tableaux-prog 2 2 3 2 2 1 4		;echo " "		#VAL >1/4
./tableaux-prog 2 2 3 2 2 1 6 2		;echo " "		#VAL >1/6, pmax =2
./tableaux-prog 2 2 3 2 2 1 4 1		;read;clear		#VAL >1/4, pmax =1

echo "_________ résolution binaire:"
./tableaux-prog 2 3 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 2 3 3 2 2			;echo " "
./tableaux-prog 2 3 3 2 2 5			;echo " "		#R >=5
./tableaux-prog 2 3 3 2 2 2	2		;echo " "		#R >=2, pmax =2
./tableaux-prog 2 3 3 2 2 0	1		;read;clear		#pmax =1

echo "_________ résolution fréquence max:"
./tableaux-prog 2 4 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 2 4 3 2 2 5			;echo " "		#R >=5
./tableaux-prog 2 4 3 2 2 2	2		;echo " "		#R >=2, pmax =2
./tableaux-prog 2 4 3 2 2 4	1		;read;clear		#R >=4, pmax =1

echo "_________ résolution min R:"
./tableaux-prog 2 5 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 2 5 3 2 2			;echo " "		#
./tableaux-prog 2 5 3 2 2 5			;echo " "		#R >=5
./tableaux-prog 2 5 3 2 2 4	2		;echo " "		#R >=4, pmax =2
./tableaux-prog 2 5 3 2 2 0	1		;read			#pmax =1

clear; echo "___________________________ appels gamma_E 3 2 2"

echo "_________ résolution continue:"
./tableaux-prog 3 0 3 2 			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 3 0 3 2 2			;echo " "	
./tableaux-prog 3 0 3 2 2 2			;echo " "		#pmax =2
./tableaux-prog 3 0 3 2 2 1			;read;clear		#pmax =1	

echo "_________ résolution entière:"
./tableaux-prog 3 1 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 3 1 3 2 2 1/3		;echo " "		#VAL >=1/3
./tableaux-prog 3 1 3 2 2 1/2.0 	;echo " "		#VAL >=1/2
./tableaux-prog 3 1 3 2 2 0.33 5	;echo " "		#VAL >=1/3, R >=5
./tableaux-prog 3 1 3 2 2 0.33 2 2	;echo " "		#VAL >=1/3, R >=2, pmax =2
./tableaux-prog 3 1 3 2 2 0.33 0 1	;read;clear;	#VAL >=1/3, pmax =1

echo "_________ résolution PV opt:"
./tableaux-prog 3 2 3 2 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 3 2 3 2 2 1 4		;echo " "		#VAL >1/4
./tableaux-prog 3 2 3 2 2 1 3		;echo " "		#VAL >1/3
./tableaux-prog 3 2 3 2 2 1 4 2		;echo " "		#VAL >1/4, pmax =2
./tableaux-prog 3 2 3 2 2 1 4 1		;read;clear		#VAL >1/3, pmax =1

echo "_________ résolution binaire:"
./tableaux-prog 3 3 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 3 3 3 2 2			;echo " "
./tableaux-prog 3 3 3 2 2 4			;echo " "		#R =4
./tableaux-prog 3 3 3 2 2 4 2		;echo " "		#R =2, pmax =2
./tableaux-prog 3 3 3 2 2 0 1		;read;clear		#pmax =1

echo "_________ résolution fréquence max:"
./tableaux-prog 3 4 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 3 4 3 2 2 5			;echo " "		#R >=5
./tableaux-prog 3 4 3 2 2 2 2		;echo " "		#R >=2, pmax =2
./tableaux-prog 3 4 3 2 2 0 1		;read;clear		#pmax =1

echo "_________ résolution min R:"
./tableaux-prog 3 5 3 2				;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 3 5 3 2 2			;echo " "		#
./tableaux-prog 3 5 3 2 2 5			;echo " "		#R >=5
./tableaux-prog 3 5 3 2 2 2	2		;echo " "		#R >=2, pmax =2
./tableaux-prog 3 5 3 2 2 0 1		;read			#pmax =1
