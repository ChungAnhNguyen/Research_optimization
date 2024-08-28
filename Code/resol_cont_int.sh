#!/bin/sh

printf 'TAB_RHO (0) ou TAB_GAMMA (1) : '
IFS= read -r typ
printf 'E_q (0/1) : '
IFS= read -r E_q
printf 'Assignation de résolution à RESOL_CONT\n'
resolution=0
pmax=0
fic_sortie=val_cible.out
printf 'Création de fichiers temporaires\n'
touch $fic_sortie

case $typ in
0)
	printf 'Entrer nu, q et k : '
	read -r nu q k
	./tableaux-prog $typ $E_q $resolution $nu $q $k
	;;
1)
	printf 'Entrer q, p, k :'
	read -r q p k
	printf 'Voulez-vous préciser pmax ? (Y,n)'
	IFS= read -r accept
	if [ $accept = 'Y' ]; then
		printf 'Entrer pmax'
		IFS= read -r pmax
	fi
	./tableaux-prog $typ $E_q $resolution $q $p $k $pmax
	;;
*)
	printf 'Type tableau incorrect.\n' >&2
	exit 2
esac

IFS= read -r val_cible <$fic_sortie
rm $fic_sortie
resolution=1
R_cible=0

case $typ in
0)
	./tableaux-prog $typ $E_q $resolution $nu $q $k $val_cible $R_cible
	;;
1)
	./tableaux-prog $typ $E_q $resolution $q $p $k $val_cible $R_cible $pmax
	;;
esac
