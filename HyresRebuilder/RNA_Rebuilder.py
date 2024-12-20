#!/usr/bin/python
from __future__ import division
import sys
import argparse
import numpy as np
from pathlib import Path
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.topology.guessers import guess_types

# This script is for backmapping CG RNA model to atomistic model
# Author: "Shanlong Li"
# Date: Sep-09-2024

def RNArebuild():
    ## input information
    parser = argparse.ArgumentParser(
      prog='RNARebuilder',
      description='Rebuild atomistic model from CG RNA model.'
    )
    parser.add_argument('-i', '--inp', help='input pdb file')
    parser.add_argument('-o', '--out', help='output pbd file')
    parser.add_argument('-n', '--num', default=1, type=int, help='the number of RNA chains')
    args = parser.parse_args()
    
    inp = args.inp
    out = args.out
    num = args.num
    
    ## load fragments
    pdb_tem = "map/ribose5.pdb"
    pdb = Path(__file__).parent / pdb_tem
    r5 = mda.Universe(pdb)
    
    pdb_tem = "map/ribose3.pdb"
    pdb = Path(__file__).parent / pdb_tem
    r3 = mda.Universe(pdb)
    
    pdb_tem = "map/bb.pdb"
    pdb = Path(__file__).parent / pdb_tem
    r0 = mda.Universe(pdb)
    
    pdb_tem = "map/po4.pdb"
    pdb = Path(__file__).parent / pdb_tem
    po4 = mda.Universe(pdb)
    
    pdb_tem = "map/baseA.pdb"
    pdb = Path(__file__).parent / pdb_tem
    baseA = mda.Universe(pdb)
    
    pdb_tem = "map/baseG.pdb"
    pdb = Path(__file__).parent / pdb_tem
    baseG = mda.Universe(pdb)
    
    pdb_tem = "map/baseC.pdb"
    pdb = Path(__file__).parent / pdb_tem
    baseC = mda.Universe(pdb)
    
    pdb_tem = "map/baseU.pdb"
    pdb = Path(__file__).parent / pdb_tem
    baseU = mda.Universe(pdb)
    
    ## main
    u = mda.Universe(inp)
    ## change the bead name to all atom name
    for atom in u.atoms:
        if atom.name == 'C1':
            atom.name = "C4'"
            atom.mass = 12.011
        elif atom.name == 'C2':
            atom.name = "C1'"
            atom.mass = 12.011
        elif atom.name == 'P':
            atom.mass = 30.974
        
        if atom.residue.resname in ['ADE', 'GUA']:
            if atom.name == 'NA':
                atom.name = 'N9'
                atom.mass = 14.007
            elif atom.name == 'NB':
                atom.name = 'N7'
                atom.mass = 14.007
            elif atom.name == 'NC':
                atom.name = 'N3'
                atom.mass = 14.007
            elif atom.name == 'ND':
                atom.name = 'N1'
                atom.mass = 14.007
        elif atom.residue.resname in ['CYT', 'URA']:
            if atom.name == 'NA':
                atom.name = 'N1'
                atom.mass = 14.007
            elif atom.name == 'NB':
                atom.name = 'C5'
                atom.mass = 12.011
            elif atom.name == 'NC':
                atom.name = 'N3'
                atom.mass = 14.007
    
    residues = u.residues
    n_res = int(len(residues)/num)
    segments = u.segments
    segids = segments.segids
    with mda.Writer(out, multiframe=False, reindex=False) as f:
        cnt = 0
        idx = 1
        for res in residues:
            segid = segids[cnt]
            cg = u.select_atoms("resid {} and segid {}".format(str(res.resid), segid))
            res_name = res.resname
            if res.resid == 1:
                mobile = r5
                if res_name in ['ADE', 'GUA']:
                    align.alignto(mobile, cg, select="name C4' C1' N9")
                else:
                    align.alignto(mobile, cg, select="name C4' C1' N1 N9")
                for atom in mobile.atoms:
                    atom.residue.resid = res.resid
                    atom.residue.segment.segid = segid
                    atom.residue.resname = res_name[0]
                    if atom.name != "N9":
                        atom.id = idx
                        idx += 1
                f.write(mobile.select_atoms("not name N9"))
    
            elif res.resid == n_res:
                mobile = r3
                if res_name in ['ADE', 'GUA']:
                    align.alignto(mobile, cg, select="name C4' C1' N9 P")
                else:
                    align.alignto(mobile, cg, select="name C4' C1' N1 N9 P")
                for atom in mobile.atoms:
                    atom.residue.resid = res.resid
                    atom.residue.segment.segid = segid
                    atom.residue.resname = res_name[0]
                    if atom.name != "N9":
                        atom.id = idx
                        idx += 1
                f.write(mobile.select_atoms("not name N9"))
                cnt += 1
    
            else:
                mobile = r0
                if res_name in ['ADE', 'GUA']:
                    align.alignto(mobile, cg, select="name C4' C1' N9 P")
                else:
                    align.alignto(mobile, cg, select="name C4' C1' N1 N9 P")
                for atom in mobile.atoms:
                    atom.residue.resid = res.resid
                    atom.residue.segment.segid = segid
                    atom.residue.resname = res_name[0]
                    if atom.name != "N9":
                        atom.id = idx
                        idx += 1
                f.write(mobile.select_atoms("not name N9"))
            
            # base groups
            if res.resname == 'ADE':
                mobile = baseA
            elif res.resname == 'GUA':
                mobile = baseG
            elif res.resname == 'CYT':
                mobile = baseC
            elif res.resname == 'URA':
                mobile = baseU
            else:
                print('Error: Unkown resname '+res_name)
                exit()
            
            if res.resname in ['ADE', 'GUA']:
                align.alignto(mobile, cg, select="name C1' N9 N7 N1 N3", strict=True)
            else:
                align.alignto(mobile, cg, select="name C1' N1 C5 N3", strict=True)
            for atom in mobile.atoms:
                atom.residue.resid = res.resid
                atom.residue.segment.segid = segid
                if atom.name != "C1'":
                        atom.id = idx
                        idx += 1
            f.write(mobile.atoms[1:])
    
    print('Done!')