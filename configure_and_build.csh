#! /usr/bin/env tcsh
#
# script to configure and  compile the MAGNETOCOSMICS code,  written by L.Desorgher 06/02/2005
# usage: type "./configur_and_build.sh " 


IGRF_TABLE=$PWD/igrfdata/igrf_dgrf_table.txt

lines='#! /bin/bash \n'
lines="$lines""#Setup script for the MAGNETOCOSMICS code\n"
lines="$lines""#Written by L. Desorgher\n"
lines="$lines""export IGRF_TABLE=$IGRF_TABLE\nexport PATH=\$PATH:$PWD/bin\n"


linestcl='#! /bin/tcsh \n'
linestcl="$linestcl""#Setup script for MAGNETOCOSMICS \n"
linestcl="$linestcl""#Written by L. Desorgher\n"
linestcl="$linestcl""setenv IGRF_TABLE $IGRF_TABLE\nsetenv PATH  \${PATH}:$PWD/bin \n "

printf  "$lines" >setupMAGNETOCOSMICS.sh
printf  "$linestcl" >setupMAGNETOCOSMICS.csh

echo "making libFortranMagField.a ..."
cd fortran
make
cd ..

if test "$USE_UNILIB" != ""; then
	echo "making libUnilib30.a ..."
	cd unilib30
	make
	cd ..
fi
echo "making MAGNETOCOSMICS ..."
make