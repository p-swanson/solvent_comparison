from lzma import open as L_open
from pickle import dump as P_dump
import numpy as np 
from scipy import ndimage
## personal packages
from download_rcsb import *
from epipy import *
from utils import *
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist
###
import warnings

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn
#import pymol 
import time
import numpy as np
import argparse
import MDAnalysis as md
from MDAnalysis.lib.distances import capped_distance
import MDAnalysis.transformations as mdt
import getpass #DOES NOT WORK ON CHTC
import datetime
import os
#################
# Mk.I.I.II
#This is the beta version of run_episol
# takes a pdbfile input
# fixes it, add hydrogens, then generates a topology and
# idc.solute file before running
#################
parser = argparse.ArgumentParser()
parser.add_argument("name",help="AF PDBID",type=str)
parser.add_argument("DENS",help="density gt")#,type=str)
args = parser.parse_args()
print(f"running on: {getpass.getuser()}")
print(f"# PARSED: {args.name}")
#################
# -------------- PARAMS
GT = float(args.DENS) # place at regions with density greater than this value
N_STEPS = 10
TEMP = 300
USE_NO_MASK = True
SIGMA = 0.5 # uncertainty distance in 3drism grid
# -------------- END PARAMS
print("# USING NO_MASK",USE_NO_MASK)
print("# USING TEMP",GT)
print("# USING N_STEPS",GT)
print("# USING RHO GT",GT)
download_rcsb(args.name)
print("# DOWNLOAD COMPLETE")

wot = md.Universe(f"fixed_{args.name}_centroid.pdb") # in this verison water and protein are merged
wot.atoms.write(f"combined_{args.name}.pdb") # for consistency
#################
buffer = 6 +1 #Angstroms. 
# using a (MINIMUM) buffer of 7 A on each dimension so we will have
# minimum image at LEAST 14 A away (with cutoff of 1nm)
# move the RCSB PDB to a new name,, 
get_res = md.Universe(f'{args.name}.pdb')
get_res.atoms.write(f'full_{args.name}.pdb')
# using USE_NO_MASK = True overides this
num_waters_to_place = len(get_res.select_atoms('resname HOH')) 
#num_waters_to_place = len(get_res.select_atoms('resname HOH')) 

#num_waters_to_place = 10*len(get_res.select_atoms('protein and name CA')) 
print("# PLACING",num_waters_to_place,"WATERS")
grid_spacing = 1 #A
#################
start1 = time.time()
# generate topology file and fix PDB file (adding hydrogens and missing residues)
[gridx,gridy,gridz],n_atoms,n_residues = pdb2top(args.name,buffer) # output is fixed_{pdb_name}.gro
start = time.time()
print(f'generated topology for {args.name} took:{start-start1}')
# take fixed coordinate file and center in a box with the minimal distance buffer applied
###############################
print("transformations done")
print("Going to RUn RISM")
run = epipy(f'fixed_{args.name}.gro',f'{args.name}.top',gen_idc=True)
#run.log = f'{args.name}_res_1A.log'
run.report(f'{args.name}_res_1A')

run.cmd_path = '/home/release/'
run.solvent_path = '/home/release/solvent/'
run.delvv = 0.5#''
run.ndiis = 15
run.rc = 0.9
run.err_tol = 1e-08
run.dynamic_delvv = 0.5#''#0.5
print('grids:',run.grid)
run.rism(step=800,resolution=1,args=('guv,excess'))
run.kernel(nt=1,v=2)
print(f'Finished run for {args.name}')
stop = time.time()
print(f"Generating topoology + Calculation took:{stop-start}")
######################################################################
# now select the top N peaks and place waters
prot = md.Universe(f'fixed_{args.name}.gro')
out_waters = run.placement(num_waters_to_place=num_waters_to_place,grid_spacing=grid_spacing,gt=GT,ALL=USE_NO_MASK)

###### GET SASA FOR EACH RESIDUE AND TOTAL SASA
import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
from pymol import cmd, stored

pymol.finish_launching()

cmd.set('dot_solvent', 1)
cmd.set('dot_density', 3)

cmd.load(f'fixed_{args.name}.gro')  # use the name of your pdb file
stored.residues = []
cmd.iterate('name ca', 'stored.residues.append(resi)')

sasa_per_residue = []
sasa_res_names = []
for i in stored.residues:
    sasa_res_names.append(f'{i}')
    sasa_per_residue.append(cmd.get_area('resi %s' % i))

sasa_dict = dict(zip(res_names,sasa_per_residue))
del sasa_res_names
del sasa_per_residue
print('SASA Calculated')
print("TOTAL SASA:",cmd.get_area('all'))
#-------------------------------------------------- BUILD WATER STRUCTURE FILE
n_atoms = len(out_waters)
seg_index = np.zeros(n_atoms)
res_index = seg_index+1
wot = md.Universe.empty(n_atoms=n_atoms,n_residues=n_atoms,
                        atom_resindex=res_index,residue_segindex=seg_index,trajectory=True)

wot.add_TopologyAttr('name',['O']*n_atoms)
wot.add_TopologyAttr('resname',['HOH']*n_atoms)
wot.add_TopologyAttr('segid',['HOH'])
wot.add_TopologyAttr('resid',[i for i in range(n_atoms)])
wot.add_TopologyAttr('ids',[i for i in range(n_atoms)])
# write out
wot.atoms.positions = out_waters
combined = md.Merge(prot.atoms,wot.atoms)
combined.atoms.write("out.gro")
combined.dimensions = prot.dimensions
combined.add_TopologyAttr('tempfactors',np.append(np.zeros(len(prot.atoms)),dens_out))
combined.atoms.write(f"{args.name}_rism_placed_waters_no_mask.pdb")
## ffrom PSE3 EXP:
A,B = -671.598836677146,-2.44628209193456

with open(f'{args.name}_res_1A.log','r') as rr:
    for line in rr:
      t = line.split()
      if t[0] == 'total':
        G_solve = float(t[10])+float(t[8])*A+B
rr.close()
############ WRITE ERROR AND FREE ENERGY TO CONDOR OUT FILE
A,B = [-619.952559990907,-1.47859772642391]
with open(f'{args.name}_res_1A.log','r') as logfile:
    for line in logfile:
         t = line.split()
         if t[0] == "RISM-PSE3":
             err_tol = float(line.split()[4])
         if t[0] == 'total':
             out_energy = float(t[10])+float(t[8])*A+B
logfile.close()
print("G_SOLV")
print("FEP",out_energy)
print("EXP",G_solve)
print("ERR",err_tol)    
###################### GET RMSD
# placed from 3DRISM
placed  = md.Universe('out.gro')#.select_atoms('resname HOH').positions
# true waters
tru = md.Universe(f'full_{args.name}.pdb')
# waters placed by superwater
supr = md.Universe(f"combined_{args.name}.pdb")
#
num_orig_wat = len(tru.select_atoms("resname HOH"))
num_placed = len(placed.select_atoms("resname HOH"))
num_supr = len(supr.select_atoms("resname HOH"))

## set box sizes to be the same as precautionary

tru.dimensions = placed.dimensions
## select the backbone CA positions

# load into pymol for alignment
cmd.load(f"full_{args.name}.pdb")
cmd.load(f"out.gro")
cmd.load(f"combined_{args.name}.pdb")
print("# LOADED INTO PYMOL")
# Now align the RISM, and the supoerwater placed 
# WARNING! WILL FAIL FOR NUCLEOTIDE. but I am only doing prot. so this is a failsafe
cmd.align(f"polymer and name CA and (out)",f"polymer and name CA and (full_{args.name})",quiet=0,object="aln",reset=1)
cmd.align(f"polymer and name CA and (combined_{args.name})",f"polymer and name CA and (full_{args.name})",quiet=0,object="aln",reset=1)
print("# ALIGNED, SAVING...")
cmd.save(f"aligned_{args.name}.pdb",f"out")
cmd.save(f"aligned_combined_{args.name}.pdb",f"combined_{args.name}")
print("# SAVED")
cmd.reinitialize()
# ----------------------------------------------------------- SHELLS

tru = md.Universe(f'full_{args.name}.pdb')
# 3DRISM aligned
placed = md.Universe(f"aligned_{args.name}.pdb")
# SuperWater aligned
suprwtr = md.Universe(f"aligned_combined_{args.name}.pdb")

# for each solvation shell, around each residue count number of waters
print("# RESNAME LOW HIGH MU-NUMBER-OF-WATERS STDEV-NUMBER-OF-WATERS VAR-NUMBER-OF-WATERS")
SHELLS = [(0.0,1.0),(1.0,2.0),(2.0,3.0),(3.0,4.0)]
#
for u_name,u_i in zip(["TRU","3DP","SWP"],[tru,placed,suprwtr]):
    print("#",17*"-",u_name,18*"-")
    for resname_i in np.unique(u_i.select_atoms("protein").resnames):
        for (low_i,hi_i) in SHELLS:
            to_go = []
            for name_i in u_i.select_atoms(f"resname {resname_i}").resids:
                N_i = len(u_i.select_atoms(f"resname HOH and (sphlayer {low_i} {hi_i} resid {name_i})"))
                to_go.append(N_i)
            to_go = np.array(to_go)
            print(f"{u_name} {resname_i} {low_i} {hi_i} {to_go.mean():0.3f} {to_go.std():0.3f} {to_go.var():0.3f}")
# =========================================================== WAT SHELLS

print("# SHELL AROUND PLACED WATERS")
print("# NAME SHELLNAME LOW HI COUNT")
SHELLS = [(0.0,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.5,3.0),(3.5,4.0),(4.5,5.0),(5.0,5.5),(5.5,6.0)]
#
#
tru_wat_pos = tru.select_atoms("resname HOH").positions # select all waters
# SELECT WATER IN SPHERICAL SHELLS AROUND TRU WAT
for u_name,u_i in zip(["3DP","SWP"],[placed,suprwtr]):
    # make a new dictionary
    td = dict(zip([f'shell_{i}' for i in range(len(SHELLS))],[[i,int(),int()] for i in SHELLS]))
    print("#",17*"-",u_name,18*"-")
    if u_name == '3DP':
        sigma = SIGMA #1.52 # smearing function
    else:
        sigma = 0.0
    for i in u_i.select_atoms("resname HOH"):
        for tru_wat_pos_i in tru_wat_pos:
            per_wat_id = int() # for each water. if there are too many waters
            for key,((lo,hi),c,c1) in td.items():
                # take distance
                x = np.linalg.norm(i.position - tru_wat_pos_i) - sigma
                tmp_dist = np.piecewise(x,[x<=0,x>0],[0.1,x])
                if lo < tmp_dist <= hi:
                    td[key][1] += 1
                    per_wat_id += 1 # keep track of surrounding wat

                td[key][2] += per_wat_id
    # --------------- write to stdout 
    for key,((lo,hi),c,c1) in td.items():
        print(f"{u_name} {key} {lo} {hi} {c:0.3f} {c1:0.3f}")
# ------------------------------------------------------------------ PREC and RECALL
print("#\n# SHELL AROUND CRYSTAL WATERS")
print("# NAME SHELLNAME LOW HI COUNT")
SHELLS = [(0.0,0.5),(0.5,1.0),(1.0,1.5),(1.5,2.0),(2.5,3.0),(3.5,4.0),(4.5,5.0),(5.0,5.5),(5.5,6.0)]
#
#
tru = md.Universe(f'full_{args.name}.pdb')
tru_wat_pos = tru.select_atoms("resname HOH").positions # select all waters
# SELECT WATER IN SPHERICAL SHELLS AROUND TRU WAT
for u_name,u_i in zip(["3DP","SWP"],[placed,suprwtr]):
    # make a new dictionary
    print("#",17*"-",u_name,18*"-")
    td = dict(zip([f'shell_{i}' for i in range(len(SHELLS))],[[i,int(),int()] for i in SHELLS]))
    if u_name == '3DP':
        sigma = SIGMA # subtract uncertainty from 3drism grid
    else:
        sigma = 0.0
    for tru_wat_pos_i in tru_wat_pos:
        # now select the placed waters
        for i in u_i.select_atoms("resname HOH"):
            per_wat_id = int() # for each water.
            for key,((lo,hi),c,c1) in td.items():
                # take distance
                x = np.linalg.norm(i.position - tru_wat_pos_i) - sigma
                tmp_dist = np.piecewise(x,[x<=0,x>0],[0.1,x])
                if lo < tmp_dist <= hi:
                    td[key][1] += 1
                    per_wat_id += 1 # keep track of surrounding wat

                td[key][2] += per_wat_id
    # --------------- write to stdout 
    for key,((lo,hi),c,c1) in td.items():
        print(f"CWAT {u_name} {key} {lo} {hi} {c:0.3f} {c1:0.3f}")
# ---------------------------------------------------------- RMSD 
# select the protein atoms from the true PDB
tru_prot = tru.select_atoms("protein") # make atomgroup
print("# ABOUT TO COMPUTE DISTANCES")
#------------------- PASS NAME TO WATERS ONLY 
placed  = placed.select_atoms("resname HOH")
suprwtr = suprwtr.select_atoms("resname HOH")

# merge two selections together and select the waters
placed = md.Merge(tru_prot.atoms,placed.atoms).select_atoms("resname HOH").positions
suprwtr = md.Merge(tru_prot.atoms,suprwtr.atoms).select_atoms("resname HOH").positions

# select two copies
supr_tru = tru.select_atoms("resname HOH").positions
rism_tru = tru.select_atoms("resname HOH").positions

def match_points_and_get_distances(data1, data2):
    """
    from superwater .git
    """
    dist_matrix = cdist(data1, data2)
    row_ind, col_ind = linear_sum_assignment(dist_matrix)
    min_distances = np.array([dist_matrix[i, j] for i, j in zip(row_ind, col_ind)])
    return min_distances

print("# ================== START 3DRISM")
tt = match_points_and_get_distances(placed,rism_tru) 

print("ORIG",num_orig_wat)
print("PLACED_R",num_placed)

print(f"LSA RISM {np.mean(tt):0.3f}")
print(f"MUSQ_LSA_RISM {np.sqrt(np.mean(tt**2)):0.3f}")
tt = cdist(placed,rism_tru)
print(f"R_NEAREST_1 {tt.min(axis=0).mean():0.3f}")
print(f"R_NEAREST_2 {tt.min(axis=1).mean():0.3f}")
print(f"R_TRU_TRU {np.mean(tt):0.3f}")
print(f"R_TRU_TRU {np.mean(tt):0.3f}")
# ======================= PREC and RECALL no double count
print(f"R_FP_05 {np.float64(tt.min(axis=1) <= 0.5+SIGMA).mean() :0.3f}")
print(f"R_FP_1  {np.float64(tt.min(axis=1) <= 1+SIGMA).mean() :0.3f}")

print(f"R_TP_05 {np.float64(tt.min(axis=0) <= 0.5+SIGMA).mean() :0.3f}")
print(f"R_TP_1  {np.float64(tt.min(axis=0) <= 1+SIGMA).mean() :0.3f}")
# 
print("==================")
tt = cdist(placed,placed)
print(f"R_PLACED_PLACED {np.mean(tt):0.3f}")
print("================== END 3DRISM")

print("================== START SUPERWATER")
print("ORIG",num_orig_wat)
print("PLACED_S",num_supr)
tt = match_points_and_get_distances(suprwtr,supr_tru)
print(f"LSA SUPRWTR {np.mean(tt):0.3f}")
print(f"MUSQ_LSA_SUPRWTR {np.sqrt(np.mean(tt**2)):0.3f}")
print("==================")
tt = cdist(suprwtr,supr_tru)
print(f"S_TRU_TRU {np.mean(tt):0.3f}")
print(f"S_NEAREST_1 {tt.min(axis=0).mean():0.3f}")
print(f"S_NEAREST_2 {tt.min(axis=1).mean():0.3f}")
# ======================= PREC and RECALL no double count
print(f"S_FP_05 {np.float64(tt.min(axis=1) <= 0.5).mean() :0.3f}")
print(f"S_FP_1  {np.float64(tt.min(axis=1) <= 1).mean() :0.3f}")

print(f"S_TP_05 {np.float64(tt.min(axis=0) <= 0.5).mean() :0.3f}")
print(f"S_TP_1  {np.float64(tt.min(axis=0) <= 1).mean()  :0.3f}")
print("==================")
tt = cdist(suprwtr,suprwtr)
print(f"S_PLACED_PLACED {np.mean(tt):0.3f}")
print("================== END SUPERWATER")
###################################
print('Run completed, deleting solute, topology and dump files')
print("hope you have a good day!")

