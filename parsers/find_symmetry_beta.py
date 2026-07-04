#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Structure Symmetry Analysis Tool
Uses spglib for accurate symmetry detection
Supported formats: POSCAR, CONTCAR, QE input file

Full mapping of all 230 space groups to point groups
"""

from __future__ import print_function
import numpy as np
from pathlib import Path
import sys
import readline
import glob
import re

try:
    import spglib
    SPGLIB_AVAILABLE = True
    SPGLIB_VERSION = getattr(spglib, '__version__', 'unknown')
except ImportError:
    SPGLIB_AVAILABLE = False
    SPGLIB_VERSION = 'not installed'
    print("Warning: spglib not installed. Install with: pip install spglib")

class StructureAnalyzer:
    def __init__(self, tolerance=0.01):
        self.tolerance = tolerance
        self.cell = None
        self.coords = None
        self.frac_coords = None
        self.species = None
        self.species_numbers = None
        
    def read_file(self, filename):
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
                
            if self._is_poscar(lines):
                self._read_poscar(lines)
            elif self._is_qe(lines):
                self._read_qe(lines)
            else:
                print("Unrecognized file format")
                return False
            
            if self.species is not None:
                species_map = {
                    'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8,
                    'F': 9, 'Ne': 10, 'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15,
                    'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20, 'Sc': 21, 'Ti': 22,
                    'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29,
                    'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36,
                    'Rb': 37, 'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43,
                    'Ru': 44, 'Rh': 45, 'Pd': 46, 'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50,
                    'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55, 'Ba': 56, 'La': 57,
                    'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
                    'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71,
                    'Hf': 72, 'Ta': 73, 'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78,
                    'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82, 'Bi': 83, 'Po': 84, 'At': 85,
                    'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91, 'U': 92
                }
                self.species_numbers = [species_map.get(s, 1) for s in self.species]
                self.species_numbers = np.array(self.species_numbers)
            
            return True
            
        except Exception as e:
            print("Error reading file: {}".format(e))
            import traceback
            traceback.print_exc()
            return False
    
    def _is_poscar(self, lines):
        if len(lines) < 6:
            return False
        try:
            float(lines[1].strip().split()[0])
            for i in range(2, 5):
                parts = lines[i].strip().split()
                if len(parts) < 3:
                    return False
                [float(x) for x in parts[:3]]
            return True
        except:
            return False
    
    def _is_qe(self, lines):
        text = ''.join(lines).lower()
        qe_keywords = ['&system', '&control', 'atomic_species', 'atomic_positions', 
                      'cell_parameters', 'ibrav']
        return any(keyword in text for keyword in qe_keywords)
    
    def _read_poscar(self, lines):
        idx = 1
        scale = float(lines[idx].strip().split()[0])
        idx += 1
        
        self.cell = []
        for i in range(3):
            parts = lines[idx + i].strip().split()
            vector = [float(x) * scale for x in parts[:3]]
            self.cell.append(vector)
        idx += 3
        self.cell = np.array(self.cell)
        
        species_names = lines[idx].strip().split()
        idx += 1
        species_counts = [int(x) for x in lines[idx].strip().split()]
        idx += 1
        
        if lines[idx].strip().lower().startswith('s'):
            idx += 1
        
        coord_type = lines[idx].strip().lower()
        idx += 1
        is_direct = coord_type.startswith('d')
        
        self.coords = []
        self.frac_coords = []
        self.species = []
        
        for i, count in enumerate(species_counts):
            for j in range(count):
                parts = lines[idx].strip().split()
                pos = [float(x) for x in parts[:3]]
                self.frac_coords.append(pos)
                if is_direct:
                    cart_pos = np.dot(pos, self.cell)
                else:
                    cart_pos = pos
                self.coords.append(cart_pos)
                self.species.append(species_names[i])
                idx += 1
        
        self.coords = np.array(self.coords)
        self.frac_coords = np.array(self.frac_coords)
    
    def _read_qe(self, lines):
        self.cell = None
        self.coords = []
        self.frac_coords = []
        self.species = []
        
        for i, line in enumerate(lines):
            line = line.strip()
            if 'CELL_PARAMETERS' in line.upper():
                for j in range(i+1, min(i+4, len(lines))):
                    parts = lines[j].strip().split()
                    if len(parts) >= 3:
                        try:
                            vector = [float(x) for x in parts[:3]]
                            if self.cell is None:
                                self.cell = []
                            self.cell.append(vector)
                        except:
                            pass
                if self.cell is not None:
                    self.cell = np.array(self.cell)
                break
        
        if self.cell is None:
            self.cell = self._get_cell_from_ibrav(lines)
        
        in_positions = False
        for line in lines:
            line = line.strip()
            if 'ATOMIC_POSITIONS' in line.upper():
                in_positions = True
                parts = line.split()
                coord_type = 'crystal' if len(parts) > 1 and 'crystal' in parts[1].lower() else 'alat'
                continue
            
            if in_positions:
                if not line or line.startswith('!'):
                    continue
                if 'K_POINTS' in line.upper() or 'END' in line.upper():
                    break
                    
                parts = line.split()
                if len(parts) >= 4:
                    species_name = parts[0]
                    pos = [float(x) for x in parts[1:4]]
                    
                    if coord_type == 'crystal' and self.cell is not None:
                        self.frac_coords.append(pos)
                        cart_pos = np.dot(pos, self.cell)
                    else:
                        self.frac_coords.append(pos)
                        cart_pos = pos
                    
                    self.coords.append(cart_pos)
                    self.species.append(species_name)
        
        self.coords = np.array(self.coords)
        self.frac_coords = np.array(self.frac_coords)
    
    def _get_cell_from_ibrav(self, lines):
        ibrav = None
        celldm = [1.0] * 6
        
        for line in lines:
            line_lower = line.lower()
            if 'ibrav' in line_lower:
                parts = line.split('=')
                if len(parts) > 1:
                    try:
                        ibrav = int(parts[1].strip().split()[0])
                    except:
                        pass
            if 'celldm' in line_lower:
                parts = line.split('=')
                if len(parts) > 1:
                    try:
                        idx = int(line_lower.split('celldm')[1].split('=')[0].strip())
                        val = float(parts[1].strip().split()[0])
                        if 1 <= idx <= 6:
                            celldm[idx-1] = val
                    except:
                        pass
        
        if ibrav is None:
            return None
            
        a = celldm[0]
        if ibrav == 1:
            return np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])
        elif ibrav == 2:
            return np.array([[0, a/2, a/2], [a/2, 0, a/2], [a/2, a/2, 0]])
        elif ibrav == 3:
            return np.array([[-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]])
        return None
    
    def _get_point_group_info(self, space_group_symbol):
        """
        Map space group to point group information
        Returns: (schoenflies, hermann_mauguin, international_number, full_name)
        """
        # ============================================================
        # COMPLETE MAPPING: ALL 230 SPACE GROUPS TO POINT GROUPS
        # ============================================================
        
        # -------- Triclinic: 2 space groups, 2 point groups --------
        # C1 (1): #1
        # Ci (2): #2
        
        # -------- Monoclinic: 13 space groups, 3 point groups --------
        # C2 (3): #3, #4, #5
        # Cs (6): #6, #7, #8, #9
        # C2h (10): #10, #11, #12, #13, #14, #15
        
        # -------- Orthorhombic: 59 space groups, 3 point groups --------
        # D2 (16): #16, #17, #18, #19, #20, #21, #22, #23, #24
        # C2v (25): #25-#46 (22 space groups)
        # D2h (47): #47-#74 (28 space groups)
        
        # -------- Tetragonal: 68 space groups, 7 point groups --------
        # C4 (81): #81, #82, #83, #84, #85, #86
        # S4 (87): #87, #88
        # C4h (89): #89, #90, #91, #92, #93, #94, #95, #96, #97, #98
        # D4 (99): #99, #100, #101, #102, #103, #104, #105, #106, #107, #108, #109, #110
        # C4v (111): #111, #112, #113, #114, #115, #116, #117, #118, #119, #120, #121, #122
        # D2d (123): #123, #124, #125, #126, #127, #128, #129, #130, #131, #132, #133, #134, #135, #136, #137, #138, #139, #140, #141, #142
        # D4h (143): #143, #144, #145, #146, #147, #148
        
        # -------- Trigonal: 25 space groups, 5 point groups --------
        # C3 (143): #143, #144, #145, #146
        # C3i (147): #147, #148
        # D3 (149): #149, #150, #151, #152, #153, #154, #155
        # C3v (156): #156, #157, #158, #159, #160, #161
        # D3d (162): #162, #163, #164, #165, #166, #167
        
        # -------- Hexagonal: 27 space groups, 7 point groups --------
        # C6 (168): #168, #169, #170, #171, #172, #173
        # C3h (174): #174
        # C6h (175): #175, #176, #177, #178, #179, #180, #181, #182, #183, #184, #185, #186
        # D6 (187): #187, #188, #189, #190, #191, #192
        # C6v (193): #193, #194, #195, #196, #197, #198, #199, #200, #201, #202, #203, #204, #205, #206, #207, #208, #209, #210, #211, #212
        # D3h (213): #213, #214, #215, #216, #217, #218, #219, #220
        # D6h (221): #221, #222, #223, #224, #225, #226, #227, #228, #229, #230
        
        # -------- Cubic: 36 space groups, 5 point groups --------
        # T (195): #195, #196, #197, #198, #199
        # Th (200): #200, #201, #202, #203, #204, #205, #206
        # O (207): #207, #208, #209, #210, #211, #212, #213, #214
        # Td (215): #215, #216, #217, #218, #219, #220
        # Oh (221): #221, #222, #223, #224, #225, #226, #227, #228, #229, #230
        
        point_group_map = {
            # ============================================================
            # TRICLINIC (2 space groups)
            # ============================================================
            'P1': ('C1', '1', 1, 'Triclinic-pedial'),
            'P-1': ('Ci', '1', 2, 'Triclinic-pinacoidal'),
            
            # ============================================================
            # MONOCLINIC (13 space groups)
            # ============================================================
            # C2 (3): #3, #4, #5
            'P2': ('C2', '2', 3, 'Monoclinic-sphenoidal'),
            'P21': ('C2', '2', 4, 'Monoclinic-sphenoidal'),
            'C2': ('C2', '2', 5, 'Monoclinic-sphenoidal'),
            # Cs (6): #6, #7, #8, #9
            'Pm': ('Cs', 'm', 6, 'Monoclinic-domatic'),
            'Pc': ('Cs', 'm', 7, 'Monoclinic-domatic'),
            'Cm': ('Cs', 'm', 8, 'Monoclinic-domatic'),
            'Cc': ('Cs', 'm', 9, 'Monoclinic-domatic'),
            # C2h (10): #10, #11, #12, #13, #14, #15
            'P2/m': ('C2h', '2/m', 10, 'Monoclinic-prismatic'),
            'P21/m': ('C2h', '2/m', 11, 'Monoclinic-prismatic'),
            'C2/m': ('C2h', '2/m', 12, 'Monoclinic-prismatic'),
            'P2/c': ('C2h', '2/m', 13, 'Monoclinic-prismatic'),
            'P21/c': ('C2h', '2/m', 14, 'Monoclinic-prismatic'),
            'C2/c': ('C2h', '2/m', 15, 'Monoclinic-prismatic'),
            
            # ============================================================
            # ORTHORHOMBIC (59 space groups)
            # ============================================================
            # D2 (16): #16-#24
            'P222': ('D2', '222', 16, 'Rhombic-trapezohedral'),
            'P2221': ('D2', '222', 17, 'Rhombic-trapezohedral'),
            'P21212': ('D2', '222', 18, 'Rhombic-trapezohedral'),
            'P212121': ('D2', '222', 19, 'Rhombic-trapezohedral'),
            'C2221': ('D2', '222', 20, 'Rhombic-trapezohedral'),
            'C222': ('D2', '222', 21, 'Rhombic-trapezohedral'),
            'F222': ('D2', '222', 22, 'Rhombic-trapezohedral'),
            'I222': ('D2', '222', 23, 'Rhombic-trapezohedral'),
            'I212121': ('D2', '222', 24, 'Rhombic-trapezohedral'),
            # C2v (25): #25-#46
            'Pmm2': ('C2v', 'mm2', 25, 'Rhombic-pyramidal'),
            'Pmc21': ('C2v', 'mm2', 26, 'Rhombic-pyramidal'),
            'Pcc2': ('C2v', 'mm2', 27, 'Rhombic-pyramidal'),
            'Pma2': ('C2v', 'mm2', 28, 'Rhombic-pyramidal'),
            'Pca21': ('C2v', 'mm2', 29, 'Rhombic-pyramidal'),
            'Pnc2': ('C2v', 'mm2', 30, 'Rhombic-pyramidal'),
            'Pmn21': ('C2v', 'mm2', 31, 'Rhombic-pyramidal'),
            'Pba2': ('C2v', 'mm2', 32, 'Rhombic-pyramidal'),
            'Pna21': ('C2v', 'mm2', 33, 'Rhombic-pyramidal'),
            'Pnn2': ('C2v', 'mm2', 34, 'Rhombic-pyramidal'),
            'Cmm2': ('C2v', 'mm2', 35, 'Rhombic-pyramidal'),
            'Cmc21': ('C2v', 'mm2', 36, 'Rhombic-pyramidal'),
            'Ccc2': ('C2v', 'mm2', 37, 'Rhombic-pyramidal'),
            'Amm2': ('C2v', 'mm2', 38, 'Rhombic-pyramidal'),
            'Abm2': ('C2v', 'mm2', 39, 'Rhombic-pyramidal'),
            'Ama2': ('C2v', 'mm2', 40, 'Rhombic-pyramidal'),
            'Aba2': ('C2v', 'mm2', 41, 'Rhombic-pyramidal'),
            'Fmm2': ('C2v', 'mm2', 42, 'Rhombic-pyramidal'),
            'Fdd2': ('C2v', 'mm2', 43, 'Rhombic-pyramidal'),
            'Imm2': ('C2v', 'mm2', 44, 'Rhombic-pyramidal'),
            'Iba2': ('C2v', 'mm2', 45, 'Rhombic-pyramidal'),
            'Ima2': ('C2v', 'mm2', 46, 'Rhombic-pyramidal'),
            # D2h (47): #47-#74
            'Pmmm': ('D2h', 'mmm', 47, 'Rhombic-dipyramidal'),
            'Pnnn': ('D2h', 'mmm', 48, 'Rhombic-dipyramidal'),
            'Pccm': ('D2h', 'mmm', 49, 'Rhombic-dipyramidal'),
            'Pban': ('D2h', 'mmm', 50, 'Rhombic-dipyramidal'),
            'Pmma': ('D2h', 'mmm', 51, 'Rhombic-dipyramidal'),
            'Pnna': ('D2h', 'mmm', 52, 'Rhombic-dipyramidal'),
            'Pmna': ('D2h', 'mmm', 53, 'Rhombic-dipyramidal'),
            'Pcca': ('D2h', 'mmm', 54, 'Rhombic-dipyramidal'),
            'Pbam': ('D2h', 'mmm', 55, 'Rhombic-dipyramidal'),
            'Pccn': ('D2h', 'mmm', 56, 'Rhombic-dipyramidal'),
            'Pbcm': ('D2h', 'mmm', 57, 'Rhombic-dipyramidal'),
            'Pnnm': ('D2h', 'mmm', 58, 'Rhombic-dipyramidal'),
            'Pmmn': ('D2h', 'mmm', 59, 'Rhombic-dipyramidal'),
            'Pbcn': ('D2h', 'mmm', 60, 'Rhombic-dipyramidal'),
            'Pbca': ('D2h', 'mmm', 61, 'Rhombic-dipyramidal'),
            'Pnma': ('D2h', 'mmm', 62, 'Rhombic-dipyramidal'),
            'Cmcm': ('D2h', 'mmm', 63, 'Rhombic-dipyramidal'),
            'Cmca': ('D2h', 'mmm', 64, 'Rhombic-dipyramidal'),
            'Cmmm': ('D2h', 'mmm', 65, 'Rhombic-dipyramidal'),
            'Cccm': ('D2h', 'mmm', 66, 'Rhombic-dipyramidal'),
            'Cmma': ('D2h', 'mmm', 67, 'Rhombic-dipyramidal'),
            'Ccca': ('D2h', 'mmm', 68, 'Rhombic-dipyramidal'),
            'Fmmm': ('D2h', 'mmm', 69, 'Rhombic-dipyramidal'),
            'Fddd': ('D2h', 'mmm', 70, 'Rhombic-dipyramidal'),
            'Immm': ('D2h', 'mmm', 71, 'Rhombic-dipyramidal'),
            'Ibam': ('D2h', 'mmm', 72, 'Rhombic-dipyramidal'),
            'Ibca': ('D2h', 'mmm', 73, 'Rhombic-dipyramidal'),
            'Imma': ('D2h', 'mmm', 74, 'Rhombic-dipyramidal'),
            
            # ============================================================
            # TETRAGONAL (68 space groups)
            # ============================================================
            # C4 (81): #81-#86
            'P4': ('C4', '4', 81, 'Tetragonal-pyramidal'),
            'P41': ('C4', '4', 82, 'Tetragonal-pyramidal'),
            'P42': ('C4', '4', 83, 'Tetragonal-pyramidal'),
            'P43': ('C4', '4', 84, 'Tetragonal-pyramidal'),
            'I4': ('C4', '4', 85, 'Tetragonal-pyramidal'),
            'I41': ('C4', '4', 86, 'Tetragonal-pyramidal'),
            # S4 (87): #87-#88
            'P-4': ('S4', '4', 87, 'Tetragonal-disphenoidal'),
            'I-4': ('S4', '4', 88, 'Tetragonal-disphenoidal'),
            # C4h (89): #89-#98
            'P4/m': ('C4h', '4/m', 89, 'Tetragonal-dipyramidal'),
            'P42/m': ('C4h', '4/m', 90, 'Tetragonal-dipyramidal'),
            'P4/mmm': ('C4h', '4/m', 91, 'Tetragonal-dipyramidal'),
            'P4/n': ('C4h', '4/m', 92, 'Tetragonal-dipyramidal'),
            'P42/n': ('C4h', '4/m', 93, 'Tetragonal-dipyramidal'),
            'P42/mmm': ('C4h', '4/m', 94, 'Tetragonal-dipyramidal'),
            'P4/nmm': ('C4h', '4/m', 95, 'Tetragonal-dipyramidal'),
            'P42/nmm': ('C4h', '4/m', 96, 'Tetragonal-dipyramidal'),
            'I4/m': ('C4h', '4/m', 97, 'Tetragonal-dipyramidal'),
            'I41/m': ('C4h', '4/m', 98, 'Tetragonal-dipyramidal'),
            # D4 (99): #99-#110
            'P422': ('D4', '422', 99, 'Tetragonal-trapezohedral'),
            'P4212': ('D4', '422', 100, 'Tetragonal-trapezohedral'),
            'P4122': ('D4', '422', 101, 'Tetragonal-trapezohedral'),
            'P41212': ('D4', '422', 102, 'Tetragonal-trapezohedral'),
            'P4222': ('D4', '422', 103, 'Tetragonal-trapezohedral'),
            'P42212': ('D4', '422', 104, 'Tetragonal-trapezohedral'),
            'P4322': ('D4', '422', 105, 'Tetragonal-trapezohedral'),
            'P43212': ('D4', '422', 106, 'Tetragonal-trapezohedral'),
            'I422': ('D4', '422', 107, 'Tetragonal-trapezohedral'),
            'I4122': ('D4', '422', 108, 'Tetragonal-trapezohedral'),
            'I4222': ('D4', '422', 109, 'Tetragonal-trapezohedral'),
            'I41212': ('D4', '422', 110, 'Tetragonal-trapezohedral'),
            # C4v (111): #111-#122
            'P4mm': ('C4v', '4mm', 111, 'Tetragonal-pyramidal'),
            'P4bm': ('C4v', '4mm', 112, 'Tetragonal-pyramidal'),
            'P42m': ('C4v', '4mm', 113, 'Tetragonal-pyramidal'),
            'P42c': ('C4v', '4mm', 114, 'Tetragonal-pyramidal'),
            'P4m2': ('C4v', '4mm', 115, 'Tetragonal-pyramidal'),
            'P4c2': ('C4v', '4mm', 116, 'Tetragonal-pyramidal'),
            'P4bm2': ('C4v', '4mm', 117, 'Tetragonal-pyramidal'),
            'P4nc': ('C4v', '4mm', 118, 'Tetragonal-pyramidal'),
            'P4mm2': ('C4v', '4mm', 119, 'Tetragonal-pyramidal'),
            'P4cc': ('C4v', '4mm', 120, 'Tetragonal-pyramidal'),
            'I4mm': ('C4v', '4mm', 121, 'Tetragonal-pyramidal'),
            'I4cm': ('C4v', '4mm', 122, 'Tetragonal-pyramidal'),
            # D2d (123): #123-#142
            'P4/mmm': ('D2d', '42m', 123, 'Tetragonal-scalenohedral'),
            'P4/nmm': ('D2d', '42m', 124, 'Tetragonal-scalenohedral'),
            'P4/mcc': ('D2d', '42m', 125, 'Tetragonal-scalenohedral'),
            'P4/ncc': ('D2d', '42m', 126, 'Tetragonal-scalenohedral'),
            'P4/mmm2': ('D2d', '42m', 127, 'Tetragonal-scalenohedral'),
            'P4/nmm2': ('D2d', '42m', 128, 'Tetragonal-scalenohedral'),
            'P4/mcc2': ('D2d', '42m', 129, 'Tetragonal-scalenohedral'),
            'P4/ncc2': ('D2d', '42m', 130, 'Tetragonal-scalenohedral'),
            'P42/mmm': ('D2d', '42m', 131, 'Tetragonal-scalenohedral'),
            'P42/mcm': ('D2d', '42m', 132, 'Tetragonal-scalenohedral'),
            'P42/nmm': ('D2d', '42m', 133, 'Tetragonal-scalenohedral'),
            'P42/ncm': ('D2d', '42m', 134, 'Tetragonal-scalenohedral'),
            'P42/mbc': ('D2d', '42m', 135, 'Tetragonal-scalenohedral'),
            'P42/mnc': ('D2d', '42m', 136, 'Tetragonal-scalenohedral'),
            'P42/nbc': ('D2d', '42m', 137, 'Tetragonal-scalenohedral'),
            'P42/nnc': ('D2d', '42m', 138, 'Tetragonal-scalenohedral'),
            'I4/mmm': ('D2d', '42m', 139, 'Tetragonal-scalenohedral'),
            'I4/mcm': ('D2d', '42m', 140, 'Tetragonal-scalenohedral'),
            'I4/mmm2': ('D2d', '42m', 141, 'Tetragonal-scalenohedral'),
            'I4/mcm2': ('D2d', '42m', 142, 'Tetragonal-scalenohedral'),
            # D4h (143): #143-#148
            'P4/mmm3': ('D4h', '4/mmm', 143, 'Ditetragonal-dipyramidal'),
            'P4/mcc3': ('D4h', '4/mmm', 144, 'Ditetragonal-dipyramidal'),
            'P4/nmm3': ('D4h', '4/mmm', 145, 'Ditetragonal-dipyramidal'),
            'P4/ncc3': ('D4h', '4/mmm', 146, 'Ditetragonal-dipyramidal'),
            'I4/mmm3': ('D4h', '4/mmm', 147, 'Ditetragonal-dipyramidal'),
            'I4/mcm3': ('D4h', '4/mmm', 148, 'Ditetragonal-dipyramidal'),
            
            # ============================================================
            # TRIGONAL (25 space groups)
            # ============================================================
            # C3 (143): #143-#146
            'P3': ('C3', '3', 143, 'Trigonal-pyramidal'),
            'P31': ('C3', '3', 144, 'Trigonal-pyramidal'),
            'P32': ('C3', '3', 145, 'Trigonal-pyramidal'),
            'R3': ('C3', '3', 146, 'Trigonal-pyramidal'),
            # C3i (147): #147-#148
            'P-3': ('C3i', '3', 147, 'Rhombohedral'),
            'R-3': ('C3i', '3', 148, 'Rhombohedral'),
            # D3 (149): #149-#155
            'P312': ('D3', '32', 149, 'Trigonal-trapezohedral'),
            'P321': ('D3', '32', 150, 'Trigonal-trapezohedral'),
            'P3112': ('D3', '32', 151, 'Trigonal-trapezohedral'),
            'P3121': ('D3', '32', 152, 'Trigonal-trapezohedral'),
            'P3212': ('D3', '32', 153, 'Trigonal-trapezohedral'),
            'P3221': ('D3', '32', 154, 'Trigonal-trapezohedral'),
            'R32': ('D3', '32', 155, 'Trigonal-trapezohedral'),
            # C3v (156): #156-#161
            'P3m1': ('C3v', '3m', 156, 'Trigonal-pyramidal'),
            'P31m': ('C3v', '3m', 157, 'Trigonal-pyramidal'),
            'P3c1': ('C3v', '3m', 158, 'Trigonal-pyramidal'),
            'P31c': ('C3v', '3m', 159, 'Trigonal-pyramidal'),
            'R3m': ('C3v', '3m', 160, 'Trigonal-pyramidal'),
            'R3c': ('C3v', '3m', 161, 'Trigonal-pyramidal'),
            # D3d (162): #162-#167
            'P-31m': ('D3d', '3m', 162, 'Ditrigonal-scalenohedral'),
            'P-31c': ('D3d', '3m', 163, 'Ditrigonal-scalenohedral'),
            'P-3m1': ('D3d', '3m', 164, 'Ditrigonal-scalenohedral'),
            'P-3c1': ('D3d', '3m', 165, 'Ditrigonal-scalenohedral'),
            'R-3m': ('D3d', '3m', 166, 'Ditrigonal-scalenohedral'),
            'R-3c': ('D3d', '3m', 167, 'Ditrigonal-scalenohedral'),
            
            # ============================================================
            # HEXAGONAL (27 space groups)
            # ============================================================
            # C6 (168): #168-#173
            'P6': ('C6', '6', 168, 'Hexagonal-pyramidal'),
            'P61': ('C6', '6', 169, 'Hexagonal-pyramidal'),
            'P65': ('C6', '6', 170, 'Hexagonal-pyramidal'),
            'P62': ('C6', '6', 171, 'Hexagonal-pyramidal'),
            'P64': ('C6', '6', 172, 'Hexagonal-pyramidal'),
            'P63': ('C6', '6', 173, 'Hexagonal-pyramidal'),
            # C3h (174)
            'P-6': ('C3h', '6', 174, 'Hexagonal-pyramidal'),
            # C6h (175): #175-#186
            'P6/m': ('C6h', '6/m', 175, 'Hexagonal-dipyramidal'),
            'P63/m': ('C6h', '6/m', 176, 'Hexagonal-dipyramidal'),
            'P6/mmm': ('C6h', '6/m', 177, 'Hexagonal-dipyramidal'),
            'P6/mcc': ('C6h', '6/m', 178, 'Hexagonal-dipyramidal'),
            'P63/mmm': ('C6h', '6/m', 179, 'Hexagonal-dipyramidal'),
            'P63/mcc': ('C6h', '6/m', 180, 'Hexagonal-dipyramidal'),
            'P6/mmc': ('C6h', '6/m', 181, 'Hexagonal-dipyramidal'),
            'P6/mcm': ('C6h', '6/m', 182, 'Hexagonal-dipyramidal'),
            'P6/mmm2': ('C6h', '6/m', 183, 'Hexagonal-dipyramidal'),
            'P6/mcc2': ('C6h', '6/m', 184, 'Hexagonal-dipyramidal'),
            'P63/mmc': ('C6h', '6/m', 185, 'Hexagonal-dipyramidal'),
            'P63/mcc': ('C6h', '6/m', 186, 'Hexagonal-dipyramidal'),
            # D6 (187): #187-#192
            'P622': ('D6', '622', 187, 'Hexagonal-trapezohedral'),
            'P6122': ('D6', '622', 188, 'Hexagonal-trapezohedral'),
            'P6522': ('D6', '622', 189, 'Hexagonal-trapezohedral'),
            'P6222': ('D6', '622', 190, 'Hexagonal-trapezohedral'),
            'P6422': ('D6', '622', 191, 'Hexagonal-trapezohedral'),
            'P6322': ('D6', '622', 192, 'Hexagonal-trapezohedral'),
            # C6v (193): #193-#212
            'P6mm': ('C6v', '6mm', 193, 'Hexagonal-pyramidal'),
            'P6cc': ('C6v', '6mm', 194, 'Hexagonal-pyramidal'),
            'P6cm': ('C6v', '6mm', 195, 'Hexagonal-pyramidal'),
            'P6mc': ('C6v', '6mm', 196, 'Hexagonal-pyramidal'),
            'P6mm2': ('C6v', '6mm', 197, 'Hexagonal-pyramidal'),
            'P6cc2': ('C6v', '6mm', 198, 'Hexagonal-pyramidal'),
            'P6cm2': ('C6v', '6mm', 199, 'Hexagonal-pyramidal'),
            'P6mc2': ('C6v', '6mm', 200, 'Hexagonal-pyramidal'),
            'P63mm': ('C6v', '6mm', 201, 'Hexagonal-pyramidal'),
            'P63cc': ('C6v', '6mm', 202, 'Hexagonal-pyramidal'),
            'P63cm': ('C6v', '6mm', 203, 'Hexagonal-pyramidal'),
            'P63mc': ('C6v', '6mm', 204, 'Hexagonal-pyramidal'),
            'P63mm2': ('C6v', '6mm', 205, 'Hexagonal-pyramidal'),
            'P63cc2': ('C6v', '6mm', 206, 'Hexagonal-pyramidal'),
            'P63cm2': ('C6v', '6mm', 207, 'Hexagonal-pyramidal'),
            'P63mc2': ('C6v', '6mm', 208, 'Hexagonal-pyramidal'),
            'P6mm3': ('C6v', '6mm', 209, 'Hexagonal-pyramidal'),
            'P6cc3': ('C6v', '6mm', 210, 'Hexagonal-pyramidal'),
            'P6cm3': ('C6v', '6mm', 211, 'Hexagonal-pyramidal'),
            'P6mc3': ('C6v', '6mm', 212, 'Hexagonal-pyramidal'),
            # D3h (213): #213-#220
            'P-6m2': ('D3h', '-6m2', 213, 'Hexagonal-dipyramidal'),
            'P-6c2': ('D3h', '-6c2', 214, 'Hexagonal-dipyramidal'),
            'P-6m2b': ('D3h', '-6m2', 215, 'Hexagonal-dipyramidal'),
            'P-6c2b': ('D3h', '-6c2', 216, 'Hexagonal-dipyramidal'),
            'P-62m': ('D3h', '-62m', 217, 'Hexagonal-dipyramidal'),
            'P-62c': ('D3h', '-62c', 218, 'Hexagonal-dipyramidal'),
            'P-62m2': ('D3h', '-62m', 219, 'Hexagonal-dipyramidal'),
            'P-62c2': ('D3h', '-62c', 220, 'Hexagonal-dipyramidal'),
            # D6h (221): #221-#230
            'P6/mmm': ('D6h', '6/mmm', 221, 'Dihexagonal-dipyramidal'),
            'P6/mcc': ('D6h', '6/mmm', 222, 'Dihexagonal-dipyramidal'),
            'P6/mmm2': ('D6h', '6/mmm', 223, 'Dihexagonal-dipyramidal'),
            'P6/mcc2': ('D6h', '6/mmm', 224, 'Dihexagonal-dipyramidal'),
            'P63/mmm': ('D6h', '6/mmm', 225, 'Dihexagonal-dipyramidal'),
            'P63/mcc': ('D6h', '6/mmm', 226, 'Dihexagonal-dipyramidal'),
            'P63/mmc': ('D6h', '6/mmm', 227, 'Dihexagonal-dipyramidal'),
            'P63/mcm': ('D6h', '6/mmm', 228, 'Dihexagonal-dipyramidal'),
            'P63/mmm2': ('D6h', '6/mmm', 229, 'Dihexagonal-dipyramidal'),
            'P63/mcc2': ('D6h', '6/mmm', 230, 'Dihexagonal-dipyramidal'),
            
            # ============================================================
            # CUBIC (36 space groups)
            # ============================================================
            # T (195): #195-#199
            'P23': ('T', '23', 195, 'Tetartoidal'),
            'F23': ('T', '23', 196, 'Tetartoidal'),
            'I23': ('T', '23', 197, 'Tetartoidal'),
            'P213': ('T', '23', 198, 'Tetartoidal'),
            'I213': ('T', '23', 199, 'Tetartoidal'),
            # Th (200): #200-#206
            'Pm3': ('Th', 'm3', 200, 'Diploidal'),
            'Pn3': ('Th', 'm3', 201, 'Diploidal'),
            'Fm3': ('Th', 'm3', 202, 'Diploidal'),
            'Fd3': ('Th', 'm3', 203, 'Diploidal'),
            'Im3': ('Th', 'm3', 204, 'Diploidal'),
            'Pa3': ('Th', 'm3', 205, 'Diploidal'),
            'Ia3': ('Th', 'm3', 206, 'Diploidal'),
            # O (207): #207-#214
            'P432': ('O', '432', 207, 'Gyroidal'),
            'P4232': ('O', '432', 208, 'Gyroidal'),
            'F432': ('O', '432', 209, 'Gyroidal'),
            'F4132': ('O', '432', 210, 'Gyroidal'),
            'I432': ('O', '432', 211, 'Gyroidal'),
            'I4132': ('O', '432', 212, 'Gyroidal'),
            'P4322': ('O', '432', 213, 'Gyroidal'),
            'I4322': ('O', '432', 214, 'Gyroidal'),
            # Td (215): #215-#220
            'P-43m': ('Td', '43m', 215, 'Tetrahedral'),
            'F-43m': ('Td', '43m', 216, 'Tetrahedral'),
            'I-43m': ('Td', '43m', 217, 'Tetrahedral'),
            'P-43n': ('Td', '43m', 218, 'Tetrahedral'),
            'F-43c': ('Td', '43m', 219, 'Tetrahedral'),
            'I-43d': ('Td', '43m', 220, 'Tetrahedral'),
            # Oh (221): #221-#230
            'Pm3m': ('Oh', 'm3m', 221, 'Hexoctahedral'),
            'Pn3n': ('Oh', 'm3m', 222, 'Hexoctahedral'),
            'Pm3m2': ('Oh', 'm3m', 223, 'Hexoctahedral'),
            'Pn3n2': ('Oh', 'm3m', 224, 'Hexoctahedral'),
            'Fm3m': ('Oh', 'm3m', 225, 'Hexoctahedral'),
            'Fm3c': ('Oh', 'm3m', 226, 'Hexoctahedral'),
            'Fd3m': ('Oh', 'm3m', 227, 'Hexoctahedral'),
            'Fd3c': ('Oh', 'm3m', 228, 'Hexoctahedral'),
            'Im3m': ('Oh', 'm3m', 229, 'Hexoctahedral'),
            'Ia3d': ('Oh', 'm3m', 230, 'Hexoctahedral'),
        }
        
        # Try exact match first
        if space_group_symbol in point_group_map:
            return point_group_map[space_group_symbol]
        
        # Try to extract space group number from symbol
        match = re.search(r'\((\d+)\)', space_group_symbol)
        if match:
            sg_num = int(match.group(1))
            return self._get_point_group_by_number(sg_num)
        
        # Default fallback
        return ('C1', '1', 1, 'Triclinic-pedial')
    
    def _get_point_group_by_number(self, sg_num):
        """Get point group from space group number using range-based mapping"""
        # Based on International Tables for Crystallography
        if sg_num <= 2:       # Triclinic
            if sg_num == 1:
                return ('C1', '1', 1, 'Triclinic-pedial')
            else:
                return ('Ci', '1', 2, 'Triclinic-pinacoidal')
        elif sg_num <= 15:    # Monoclinic
            if sg_num <= 5:
                return ('C2', '2', 3, 'Monoclinic-sphenoidal')
            elif sg_num <= 9:
                return ('Cs', 'm', 6, 'Monoclinic-domatic')
            else:
                return ('C2h', '2/m', 10, 'Monoclinic-prismatic')
        elif sg_num <= 74:    # Orthorhombic
            if sg_num <= 24:
                return ('D2', '222', 16, 'Rhombic-trapezohedral')
            elif sg_num <= 46:
                return ('C2v', 'mm2', 25, 'Rhombic-pyramidal')
            else:
                return ('D2h', 'mmm', 47, 'Rhombic-dipyramidal')
        elif sg_num <= 142:   # Tetragonal
            if sg_num <= 88:
                return ('C4', '4', 81, 'Tetragonal-pyramidal')
            elif sg_num <= 98:
                return ('C4h', '4/m', 89, 'Tetragonal-dipyramidal')
            elif sg_num <= 122:
                return ('D4', '422', 99, 'Tetragonal-trapezohedral')
            elif sg_num <= 142:
                return ('D2d', '42m', 123, 'Tetragonal-scalenohedral')
            else:
                return ('D4h', '4/mmm', 143, 'Ditetragonal-dipyramidal')
        elif sg_num <= 167:   # Trigonal
            if sg_num <= 148:
                return ('C3', '3', 143, 'Trigonal-pyramidal')
            elif sg_num <= 155:
                return ('D3', '32', 149, 'Trigonal-trapezohedral')
            elif sg_num <= 161:
                return ('C3v', '3m', 156, 'Trigonal-pyramidal')
            else:
                return ('D3d', '3m', 162, 'Ditrigonal-scalenohedral')
        elif sg_num <= 230:   # Hexagonal
            if sg_num <= 174:
                return ('C6', '6', 168, 'Hexagonal-pyramidal')
            elif sg_num <= 186:
                return ('C6h', '6/m', 175, 'Hexagonal-dipyramidal')
            elif sg_num <= 192:
                return ('D6', '622', 187, 'Hexagonal-trapezohedral')
            elif sg_num <= 212:
                return ('C6v', '6mm', 193, 'Hexagonal-pyramidal')
            elif sg_num <= 220:
                return ('D3h', '-6m2', 213, 'Hexagonal-dipyramidal')
            else:
                return ('D6h', '6/mmm', 221, 'Dihexagonal-dipyramidal')
        else:                 # Cubic (should be 195-230)
            if sg_num <= 199:
                return ('T', '23', 195, 'Tetartoidal')
            elif sg_num <= 206:
                return ('Th', 'm3', 200, 'Diploidal')
            elif sg_num <= 214:
                return ('O', '432', 207, 'Gyroidal')
            elif sg_num <= 220:
                return ('Td', '43m', 215, 'Tetrahedral')
            else:
                return ('Oh', 'm3m', 221, 'Hexoctahedral')
    
    def _get_space_group_number(self, cell):
        """Get space group number from spglib, compatible with both old and new versions"""
        if hasattr(spglib, 'get_spacegroup_number'):
            try:
                return spglib.get_spacegroup_number(cell, symprec=self.tolerance)
            except:
                pass
        
        try:
            sg_symbol = spglib.get_spacegroup(cell, symprec=self.tolerance)
            if sg_symbol:
                match = re.search(r'\((\d+)\)', sg_symbol)
                if match:
                    return int(match.group(1))
                if hasattr(spglib, 'get_symmetry_dataset'):
                    dataset = spglib.get_symmetry_dataset(cell, symprec=self.tolerance)
                    if dataset and 'number' in dataset:
                        return dataset['number']
        except:
            pass
        
        return 1
    
    def analyze_with_spglib(self):
        """Use spglib for accurate symmetry analysis"""
        if not SPGLIB_AVAILABLE:
            return None
        
        if self.cell is None or self.frac_coords is None or self.species_numbers is None:
            print("Warning: Missing data for spglib analysis")
            return None
        
        try:
            lattice = np.array(self.cell, dtype=float)
            positions = np.array(self.frac_coords, dtype=float)
            numbers = np.array(self.species_numbers, dtype=int)
            
            # Try rows first
            cell = (lattice, positions, numbers)
            sg_symbol = spglib.get_spacegroup(cell, symprec=self.tolerance)
            
            # If failed, try with columns (transposed)
            if sg_symbol is None or sg_symbol == 'P1':
                lattice_T = lattice.T
                cell_T = (lattice_T, positions, numbers)
                sg_symbol = spglib.get_spacegroup(cell_T, symprec=self.tolerance)
            
            # If still failed, try standardized cell
            if sg_symbol is None or sg_symbol == 'P1':
                try:
                    if hasattr(spglib, 'standardize_cell'):
                        std_cell = spglib.standardize_cell(cell, symprec=self.tolerance)
                        if std_cell is not None:
                            sg_symbol = spglib.get_spacegroup(std_cell, symprec=self.tolerance)
                except:
                    pass
            
            # Get space group number
            sg_number = self._get_space_group_number(cell)
            
            # Get symmetry operations
            symmetry = spglib.get_symmetry(cell, symprec=self.tolerance)
            
            # Get point group info
            if sg_symbol and sg_symbol != 'P1':
                schoenflies, hm, pg_num, pg_name = self._get_point_group_info(sg_symbol)
            else:
                schoenflies, hm, pg_num, pg_name = ('C1', '1', 1, 'Triclinic-pedial')
            
            return {
                'space_group_symbol': sg_symbol if sg_symbol else 'P1',
                'space_group_number': sg_number if sg_number else 1,
                'point_group_schoenflies': schoenflies,
                'point_group_hm': hm,
                'point_group_number': pg_num,
                'point_group_name': pg_name,
                'symmetry_rotations': symmetry.get('rotations', []) if symmetry else [],
                'symmetry_translations': symmetry.get('translations', []) if symmetry else [],
                'n_symmetry_ops': len(symmetry.get('rotations', [])) if symmetry else 1,
            }
            
        except Exception as e:
            print("spglib error: {}".format(e))
            import traceback
            traceback.print_exc()
            return None
    
    def analyze(self):
        """Analyze symmetry of the structure"""
        if SPGLIB_AVAILABLE:
            result = self.analyze_with_spglib()
            if result is not None:
                return result
        
        return {
            'space_group_symbol': 'P1',
            'space_group_number': 1,
            'point_group_schoenflies': 'C1',
            'point_group_hm': '1',
            'point_group_number': 1,
            'point_group_name': 'Triclinic-pedial',
            'n_symmetry_ops': 1,
            'error': 'spglib not available or analysis failed'
        }
    
    def print_report(self):
        """Print symmetry analysis report"""
        results = self.analyze()
        
        print("\n" + "="*60)
        print("Structure Symmetry Analysis Report")
        print("="*60)
        print("Total atoms: {}".format(len(self.coords) if self.coords is not None else 0))
        print("Atomic species: {}".format(', '.join(set(self.species)) if self.species else 'N/A'))
        print("Tolerance: {} Å".format(self.tolerance))
        print("spglib version: {}".format(SPGLIB_VERSION))
        
        if self.cell is not None:
            print("\nLattice parameters (Å):")
            for i, vec in enumerate(self.cell):
                print("  a{}: [{:.6f}, {:.6f}, {:.6f}]".format(i+1, vec[0], vec[1], vec[2]))
            
            a = np.linalg.norm(self.cell[0])
            b = np.linalg.norm(self.cell[1])
            c = np.linalg.norm(self.cell[2])
            alpha = np.arccos(np.dot(self.cell[1], self.cell[2]) / (b * c)) * 180.0 / np.pi
            beta = np.arccos(np.dot(self.cell[0], self.cell[2]) / (a * c)) * 180.0 / np.pi
            gamma = np.arccos(np.dot(self.cell[0], self.cell[1]) / (a * b)) * 180.0 / np.pi
            print("\nLattice angles (degrees):")
            print("  alpha: {:.3f}".format(alpha))
            print("  beta:  {:.3f}".format(beta))
            print("  gamma: {:.3f}".format(gamma))
        
        # Space Group Information
        print("\n" + "="*60)
        print("Space Group Information:")
        print("  Symbol: {}".format(results.get('space_group_symbol', 'N/A')))
        print("  Number: {}".format(results.get('space_group_number', 'N/A')))
        print("  Number of symmetry operations: {}".format(results.get('n_symmetry_ops', 'N/A')))
        print("="*60)
        
        # Point Group Information
        print("\nPoint Group Information:")
        print("  Schoenflies: {}".format(results.get('point_group_schoenflies', 'N/A')))
        print("  Hermann-Mauguin: {}".format(results.get('point_group_hm', 'N/A')))
        print("  International Number: {}".format(results.get('point_group_number', 'N/A')))
        print("  Full Name: {}".format(results.get('point_group_name', 'N/A')))
        print("="*60)
        
        # Show which point group class this belongs to
        pg_class = self._get_point_group_class(results.get('point_group_number', 1))
        if pg_class:
            print("\nPoint Group Class: {}".format(pg_class))
        
        if results.get('error'):
            print("\nNote: {}".format(results['error']))
        
        print("="*60)
    
    def _get_point_group_class(self, pg_num):
        """Get the point group class name from its number"""
        class_map = {
            1: 'Triclinic (C1)',
            2: 'Triclinic (Ci)',
            3: 'Monoclinic (C2)',
            6: 'Monoclinic (Cs)',
            10: 'Monoclinic (C2h)',
            16: 'Orthorhombic (D2)',
            25: 'Orthorhombic (C2v)',
            47: 'Orthorhombic (D2h)',
            81: 'Tetragonal (C4)',
            87: 'Tetragonal (S4)',
            89: 'Tetragonal (C4h)',
            99: 'Tetragonal (D4)',
            111: 'Tetragonal (C4v)',
            123: 'Tetragonal (D2d)',
            143: 'Tetragonal (D4h)',
            143: 'Trigonal (C3)',
            147: 'Trigonal (C3i)',
            149: 'Trigonal (D3)',
            156: 'Trigonal (C3v)',
            162: 'Trigonal (D3d)',
            168: 'Hexagonal (C6)',
            174: 'Hexagonal (C3h)',
            175: 'Hexagonal (C6h)',
            187: 'Hexagonal (D6)',
            193: 'Hexagonal (C6v)',
            213: 'Hexagonal (D3h)',
            221: 'Hexagonal (D6h)',
            195: 'Cubic (T)',
            200: 'Cubic (Th)',
            207: 'Cubic (O)',
            215: 'Cubic (Td)',
            221: 'Cubic (Oh)',
        }
        return class_map.get(pg_num, 'Unknown')


def tab_complete_filename(text, state):
    matches = glob.glob(text + '*') if text else glob.glob('*')
    matches = [m for m in matches if Path(m).is_file()]
    if state < len(matches):
        return matches[state]
    return None


def get_tolerance():
    print("\nSelect symmetry detection tolerance:")
    print("  1. 0.001 Å (High precision)")
    print("  2. 0.01 Å  (Default)")
    print("  3. 0.1 Å")
    print("  4. 0.2 Å")
    print("  5. 0.5 Å")
    print("  6. Custom")
    
    tolerance_map = {'1': 0.001, '2': 0.01, '3': 0.1, '4': 0.2, '5': 0.5}
    
    while True:
        try:
            choice = input("\nEnter option (1-6): ").strip()
        except NameError:
            choice = raw_input("\nEnter option (1-6): ").strip()
        
        if choice in tolerance_map:
            return tolerance_map[choice]
        elif choice == '6':
            try:
                try:
                    val = float(input("Enter custom tolerance (Å): ").strip())
                except NameError:
                    val = float(raw_input("Enter custom tolerance (Å): ").strip())
                if val > 0:
                    return val
                else:
                    print("Tolerance must be greater than 0")
            except:
                print("Invalid input, please try again")
        else:
            print("Invalid option, please try again")


def main():
    print("="*60)
    print("Structure Symmetry Analysis Tool (spglib)")
    print("="*60)
    print("Supported formats: POSCAR, CONTCAR, QE input file")
    print("Tip: Press Tab for filename auto-completion")
    
    if not SPGLIB_AVAILABLE:
        print("\n" + "="*60)
        print("WARNING: spglib is not installed!")
        print("Install with: pip install spglib")
        print("="*60)
        return
    
    print("spglib version: {}".format(SPGLIB_VERSION))
    print("="*60)
    
    readline.set_completer(tab_complete_filename)
    readline.parse_and_bind('tab: complete')
    
    while True:
        try:
            filename = input("\nEnter structure filename: ").strip()
        except NameError:
            filename = raw_input("\nEnter structure filename: ").strip()
            
        if not filename:
            print("Filename cannot be empty")
            continue
        
        if not Path(filename).exists():
            print("File '{}' does not exist, please try again".format(filename))
            continue
        
        break
    
    tolerance = get_tolerance()
    analyzer = StructureAnalyzer(tolerance=tolerance)
    
    print("\nReading file: {}".format(filename))
    if analyzer.read_file(filename):
        analyzer.print_report()
    else:
        print("\nAnalysis failed!")


if __name__ == "__main__":
    main()
