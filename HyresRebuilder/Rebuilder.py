"""
this package is used to backmap HyRes model to atomistic model
Athour: Shanlong Li
"""

import sys
import argparse
import numpy as np
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.topology.guessers import guess_types
from .Rotamer import opt_side_chain


parser = argparse.ArgumentParser(
  prog='HyresRebuilder',
  description='Rebuild atomistic model from HyRes model.'
)
inp = parser.add_argument('input')
out = parser.add_argument('output')

def rebuild(inp, out):
    hyres = mda.Universe(inp)
    guessed_eles = guess_types(hyres.atoms.names)
    hyres.add_TopologyAttr('elements', guessed_eles)

    with mda.Writer(out, multiframe=False, reindex=False) as f:
        idx = 1
        for res in hyres.residues:
            res_name = res.resname
            if res_name in ['ARG','HIS','LYS','ASP','GLU','SER','THR','ASN','GLN','CYS','GLY','PRO','ALA','VAL','ILE','LEU','MET','PHE','TYR','TRP']:
                pdb = Path(__file__).parent / "../map/"+res_name+"_ideal.pdb"
                mobile = mda.Universe(pdb)
            else:
                print('Error: Unkown resname '+res_name)
                exit()
                
            segid = hyres.select_atoms("resid "+str(res.resid)).segids[0]
            chainID = hyres.select_atoms("resid "+str(res.resid)).chainIDs[0]
            for atom in mobile.atoms:
                atom.residue.resid = res.resid
                atom.segment.segid = segid
                atom.chainID = chainID
 
            #align.alignto(mobile.select_atoms("name N CA C O CB"), hyres.select_atoms("resid "+str(res.resid)+" and name N CA C O"), select='name N CA C O', match_atoms=False)
            align.alignto(mobile, hyres.select_atoms("resid "+str(res.resid)), select='name N CA C O', match_atoms=False)
 
            if res_name not in ['GLY', 'PRO', 'ALA']:
                refs = hyres.select_atoms("resid "+str(res.resid)+" and name CA CB CC CD CE CF")
                opt_side_chain(res_name, refs, mobile)
        
            if res_name != 'PRO':
                hyres_H = hyres.select_atoms("resid "+str(res.resid)+" and name H")
                for atom in hyres_H.atoms:
                    atom.id = idx
                    idx += 1
                f.write(hyres_H)
 
            mobile_sel = mobile.select_atoms("all")
            for atom in mobile_sel.atoms:
                atom.id = idx
                idx += 1
            #f.write(hyres_sel)
            f.write(mobile_sel)
