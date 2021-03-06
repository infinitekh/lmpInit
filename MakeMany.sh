#!/bin/bash -
#===============================================================================
#
#          FILE: MakeMany.sh
#
#         USAGE: ./MakeMany.sh
#
#   DESCRIPTION: 
#
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: KIM Hyeok (kh), ekh0324@gmail.com
#  ORGANIZATION: Konkuk University
#       CREATED: 2019년 03월 05일 17시 25분 59초
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error
let j=1
while read N1 Box  diameter N2 
do
if [ ! -f data.srd.$j ]
then
./lmpInit.out <<ENDL
${diameter}
${N1}
${N2}
0,${Box}
ENDL
mv default_output.out  data.srd.$j 
echo "Diameter= ${diameter}, replica_num = $j, N1=${N1}, N2=${N2}, Box=${Box} cube"
else
	echo "data.srd.$j" is exist.
#	exit 1;
fi
let j++
done <<< "2	16						7.25566335719562						19456
4	16						5.75882382296972						19456
8	16						4.57078149734083						19456
12	16						3.99294542465508						19456
16	16						3.62783167859781						19456
32	16						2.87941191148486						19456
2	32						14.5113267143912						155648
4	32						11.5176476459394						155648
8	32						9.14156299468166						155648
12	32						7.98589084931016						155648
16	32						7.25566335719562						155648
32	32						5.75882382296972						155648
2	16						9.14156299468166						18432
4	16						7.25566335719562						18432
8	16						5.75882382296972						18432
12	16						5.03079599160436						18432
16	16						4.57078149734083						18432
32	16						3.62783167859781						18432
2	32						18.2831259893633						147456
4	32						14.5113267143912						147456
8	32						11.5176476459394						147456
12	32						10.0615919832087						147456
16	32						9.14156299468166						147456
32	32						7.25566335719562						147456 "

tsp -S4
let j=1
while read N1 Box  diameter N2 
do
if [  -f data.srd.$j ]
then
tsp /home/kh/bin/lmp_serial -i rein_make.B -v dia ${diameter} -v lbd 0.1 -v replica $j -v ran1 $RANDOM 
echo "Diameter= ${diameter}, replica_num = $j, N1=${N1}, N2=${N2}, Box=${Box} cube"
fi
let j++
done <<< "2	16						7.25566335719562						19456
4	16						5.75882382296972						19456
8	16						4.57078149734083						19456
12	16						3.99294542465508						19456
16	16						3.62783167859781						19456
32	16						2.87941191148486						19456
2	32						14.5113267143912						155648
4	32						11.5176476459394						155648
8	32						9.14156299468166						155648
12	32						7.98589084931016						155648
16	32						7.25566335719562						155648
32	32						5.75882382296972						155648
2	16						9.14156299468166						18432
4	16						7.25566335719562						18432
8	16						5.75882382296972						18432
12	16						5.03079599160436						18432
16	16						4.57078149734083						18432
32	16						3.62783167859781						18432
2	32						18.2831259893633						147456
4	32						14.5113267143912						147456
8	32						11.5176476459394						147456
12	32						10.0615919832087						147456
16	32						9.14156299468166						147456
32	32						7.25566335719562						147456 "
