#!/usr/bin/python

import sys
import pymatgen as pm
import re

m = pm.MPRester("UvSgztFQ7CcO6fb5")

def GetStructure(id):
	st = m.get_structure_by_material_id(id, False)
	
	return st

def GetData(id):
	return m.get_entries(id)

def MakeBlock(block_name, content):
	block = "%block " + block_name + "\n"
	for el in content:
		block += " ".join(el) + "\n"
	block += "%endblock " + block_name + "\n"
	return block

def StructToFDF(st):
	elems = dict(zip(st.symbol_set, range(1, len(st.symbol_set) + 1)))
	anums = st.atomic_numbers
	anums = [v for i,v in enumerate(anums) if v not in anums[:i]]

	ChemSpec = [(str(n), str(anum), elem) for elem, n, anum in zip(elems.iterkeys(), elems.itervalues(), anums)]
	ChemSpeclBlock = MakeBlock("ChemicalSpeciesLabel", ChemSpec)

	LatticeDict = st.to_dict["lattice"]
	LatticeParam = [[str(round(LatticeDict[key], 3)) for key in 'a', 'b', 'c', 'alpha', 'beta', 'gamma']]

	LatticeParamBlock = MakeBlock("LatticeParameters", LatticeParam)

	SitesDict = st.to_dict["sites"]
	# for site in SitesDict:
		# AtomicCoord = [[[str(round(coord, 3)) for coord in site['xyz']]] + [str(elems[site["label"]])]]

	AtomicCoord = [[str(site["xyz"][0]), str(site["xyz"][1]), str(site["xyz"][2]), str(elems[site["label"]])] for site in SitesDict]
	AtomicCoordBlock = MakeBlock("AtomicCoordinatesAndAtomicSpecies", AtomicCoord)

	Fdf = []
	Fdf.append(["SystemName", st.formula])
	Fdf.append(["SystemLabel", sys.argv[1]])
	Fdf.append(["NumberOfAtoms", str(st.num_sites)])
	Fdf.append(["NumberOfSpecies", str(st.ntypesp)])
	Fdf.append(["AtomicCoordinatesFormat", "Ang"])
	Fdf.append(["LatticeConstant", "1 Ang"])
	Fdf.append(["MeshCutoff", "125.0 Ry"])
	Fdf.append(["WriteCoorXmol", "yes"])
	Fdf.append(["\n", ChemSpeclBlock])
	Fdf.append(["\n", LatticeParamBlock])
	Fdf.append(["\n", AtomicCoordBlock])

	return Fdf

def WriteToFDF(Fdf):
	f = open(sys.argv[1] + ".fdf", "w")

	for elem in Fdf:
		f.write(" ".join(elem) + "\n")

	f.write("\n")
	f.close()

WriteToFDF(StructToFDF(GetStructure(sys.argv[1])))