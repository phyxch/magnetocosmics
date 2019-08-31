#! /usr/bin/env bash
#
# script to configure and  compile the static MAGNETOCOSMICS code,  written by L.Desorgher 06/02/2005
# usage: type "./configure_and_build_static.sh " 


IGRF_TABLE=$PWD/igrfdata/igrf_dgrf_table.txt

lines='#! /bin/bash \n'
lines="$lines""#Setup script for the MAGNETOCOSMICS code\n"
lines="$lines""#Written by L. Desorgher\n"
lines="$lines""export IGRF_TABLE=$IGRF_TABLE\n export PATH=\$PATH:$PWD/bin"


linestcl='#! /bin/tcsh \n'
linestcl="$linestcl""#Setup script for MAGNETOCOSMICS that uses ROOT\n"
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
env BUILD_STATIC=1 make	
 


		
			
	
  







