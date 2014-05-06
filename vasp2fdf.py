#!/usr/bin/python

import sys
import os
import pymatgen as pm

def get_incar_dict(incar_file):
	ICDict = {}

	with open(incar_file) as f:
		for line in f:
			line = line.rstrip('\r\n')
			ICDict.update(dict([line.split(' = ')]))

	return ICDict

def get_poscar(poscar_file):
	with open(poscar_file) as f:
		poscar = f.readlines()

	pclist = []

	atoms = poscar[6][:-1].split(' ')
	atom_names = poscar[5][:-1].split(' ')
	elems = {}
	n = 1
	for el in atom_names:
		elems[el] = n
		n += 1

	pclist.append(''.join([el + atoms[n - 1] for el, n in elems.iteritems()]))

	pclist.append(str(sum([int(atom) for atom in atoms])))
	pclist.append(str(len(atoms)))

	pclist.append([str(val) + '  ' + str(pm.Element(key).number) + ' ' + key for key, val in elems.iteritems()])

	pclist.append(poscar[1][:-1])
	pclist += [vec[:-1] for vec in poscar[2:5]]

	for coord in poscar[8:]:
		c_list = coord.split(' ')
		c_list[3] = str(elems[c_list[3][:-1]])
		pclist.append(' '.join(c_list))
	print pclist
	return pclist

def get_kpoints(kpoints_file):
	try:
		with open(kpoints_file) as f:
			kpoints = []
			for line in f:
				if line[0] == 'M' or line[0] == 'm' or line[0] == 'G' or line[0] == 'g':
					if line[0] == 'M' or line[0] == 'm':
						dkplus = 0
					else:
						dkplus = -0.5
					try:
						kpoints.append(f.next()[:-1].split(' '))
						kpoints.append([str(int(dk) + dkplus) for dk in f.next()[:-1].split(' ')])
						print kpoints
					except StopIteration:
						kpoints.append(['-0.5', '-0.5', '-0.5'])
					return kpoints
	except:
		return None
	return None

def to_fdf(ICDict, PCoord, KPoints):
	Fdf = []

	Fdf.append("SystemName " + PCoord[0])
	Fdf.append("SystemLabel " + PCoord[0])
	Fdf.append("NumberOfAtoms " + PCoord[1])
	Fdf.append("NumberOfSpecies " + PCoord[2])

	Fdf.append("%block ChemicalSpeciesLabel")
	for l in PCoord[3]:
		Fdf.append(l)
	Fdf.append("%endblock ChemicalSpeciesLabel\n")

	Fdf.append("XC.functional GGA\n")
	Fdf.append("XC.authors PBE\n")

	Fdf.append("LatticeConstant " + PCoord[4] + " Ang")
	Fdf.append("%block LatticeVectors")
	for v in PCoord[5:8]:
		Fdf.append(v)
	Fdf.append("%endblock LatticeVectors\n")

	Fdf.append("AtomicCoordinatesFormat  Ang")

	Fdf.append("%block AtomicCoordinatesAndAtomicSpecies")
	for c in PCoord[8:]:
		Fdf.append(c)
	Fdf.append("%endblock AtomicCoordinatesAndAtomicSpecies\n")

	Fdf.append("ElectronicTemperature %(SIGMA)s eV" % ICDict) 

	if int(ICDict["ISMEAR"]) in range(-5, 0):
		Fdf.append("OccupationFunction FD")
	elif int(ICDict["ISMEAR"]) >= 1:
		Fdf.append("OccupationFunction MP")
		Fdf.append("OccupationMPOrder %(ISMEAR)s" % ICDict)	

	if PCoord:
		Fdf.append("")

	if "MAGMOM" in ICDict.keys():
		MagList = ICDict["MAGMOM"].split(' ')
		FdfSpin = []
		for mag in MagList:
			if '*' in mag:
				m = mag.split('*')
				FdfSpin += int(m[0])*[m[1]]
			else:
				FdfSpin += mag

		if set(FdfSpin) != {'0.6'}:
			Fdf.append("SpinPolarised .true.")
			Fdf.append("%block DM.initspin")
			Fdf += [str(i + 1) + '    ' + l for i, l in enumerate(FdfSpin)]
			Fdf.append("%endblock DM.initspin")

	if KPoints:
		Fdf.append("%block kgrid_Monkhorst_pack")
		Fdf.append(KPoints[0][0] + ' 0 0 ' + KPoints[1][0])
		Fdf.append('0 ' + KPoints[0][1] + ' 0 ' + KPoints[1][1])
		Fdf.append('0 0 ' + KPoints[0][2] + ' ' + KPoints[1][2])
		Fdf.append("%endblock kgrid_Monkhorst_pack\n")

	Fdf.append("MeshCutoff		125.0 Ry")
	Fdf.append("WriteCoorXmol		yes")
	return Fdf

dir_name = sys.argv[1].rstrip('/')
if dir_name == '.':
	file_name = os.getcwd().split('/')[-1]
else:
	file_name = dir_name
kp = get_kpoints(dir_name + '/KPOINTS')
FdfList = to_fdf(get_incar_dict(dir_name + '/INCAR'), get_poscar(dir_name + '/POSCAR'), kp)

with open(file_name + '.fdf', 'w') as f:
	for line in FdfList:
		f.write(line + '\n')

print "Fdf File: " + file_name + ".fdf"
