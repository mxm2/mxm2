#!/bin/bash

PSF_LINK="http://departments.icmab.es/leem/siesta/Databases/Pseudopotentials/Pseudos_GGA_Abinit"

for arg in $*; do
	wget "${PSF_LINK}/${arg}_html/${arg}.psf"
done

