This repo contains the necessary files for reproducing
a comparison between 3DRISM and SuperWater placed water atoms.

layout
------
* pdb\_lists: text files with the PDB names used for comparison
* rism\_and\_workup\_code: code that was used to run 3drism and compare results
* super\_configs: .yaml files for running superwater
* run\_files: bash scripts used to run jobs on HTCondore
* submit\_files: .sub files for HTCondor

Dependancies
------------
* PDBFixer -> add missing residues
* OpenMM 8.0 -> parameterize proteins
* Parmed -> write .top files
* epipy and EPISOL -> run 3DRISM
* numpy -> misc.
* scipy -> misc.
