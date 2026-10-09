# taken from stack exchange 
def download_rcsb(pdb_id: str,download_dir: str = "./") -> None:
    from requests import get
    import os
    r = get(f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb")
    open(os.path.join(download_dir, f"{pdb_id.upper()}.pdb"), 'w').write(r.text)
    print(f"Downloaded {pdb_id}")
    return

    
