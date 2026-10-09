# AD 6.17.2026 
import argparse
from requests import get
import os
import MDAnalysis as md
####
parser = argparse.ArgumentParser()
parser.add_argument("file",type=str,help="single column file with 4-letter PDBIDs")
args = parser.parse_args()
####
print("# USING FILE",args.file)

def download_rcsb(pdb_id: str,download_dir: str = "./") -> None:
    r = get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb")
    open(os.path.join(download_dir, f"{pdb_id.upper()}.pdb"), 'w').write(r.text)
    print(f"Downloaded {pdb_id}")
    return

ind = int()
with open(args.file,"r") as rr:
    for line in rr:
        try:
            download_rcsb(line.strip())
            # remove waters and overwrite downloaded pdb
            u = md.Universe(f"{line.strip()}.pdb")
            prot = u.select_atoms("protein and not resname HOH")
            prot.atoms.write(f"fixed_{line.strip().upper()}.pdb")
            ind += 1
            # get rid of old pdb
            os.remove(f"{line.strip()}.pdb")
        except Exception as exc:
            print("#FAILED FETCHING:",line,"AT INDEX:",ind)
            print("#REASON:",exc)
            continue
########################## COMPLETE
print("# DOWNLOAD COMPLETE: KYRIRE ELESON")

