#!/bin/sh

clear;	echo "___________________________ appels delta 3 2 2"

echo "_________ résolution continue:"
./tableaux-prog 4 0 3 2 			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 4 0 3 2 2			;echo " "	
./tableaux-prog 4 0 3 2 2 1			;echo " "		#pmax =1
./tableaux-prog 4 0 3 2 2 0			;read;clear		#pmax =0	

echo "_________ résolution entière:"
./tableaux-prog 4 1 3 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 4 1 3 2 2 1/4		;echo " "		#VAL >=1/4
./tableaux-prog 4 1 3 2 2 1/3.0 	;echo " "		#VAL >=1/3
./tableaux-prog 4 1 3 2 2 0.25 5	;echo " "		#VAL >=1/4, R =5
./tableaux-prog 4 1 3 2 2 0.25 4 1	;echo " "		#VAL >=1/4, R =4, pmax =1
./tableaux-prog 4 1 3 2 2 0.25 0 0	;read;clear;	#VAL >=1/4, pmax =0

echo "_________ résolution PV opt:"
./tableaux-prog 4 2 3 2 2 2			;echo " "		#IL MANQUE DES ARGUMÙENTS
./tableaux-prog 4 2 3 2 2 1 6		;echo " "		#VAL >1/6
./tableaux-prog 4 2 3 2 2 1 4		;echo " "		#VAL >1/4
./tableaux-prog 4 2 3 2 2 1 6 1		;echo " "		#VAL >1/6, pmax =1
./tableaux-prog 4 2 3 2 2 1 4 0		;read;clear		#VAL >1/4, pmax =0
