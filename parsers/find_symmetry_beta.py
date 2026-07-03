#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Structure Symmetry Analysis Tool
Supported formats: POSCAR, CONTCAR, QE input file
"""

from __future__ import print_function
import numpy as np
from pathlib import Path
import sys
import readline
import glob
from typing import Tuple, List, Optional, Dict
from scipy.spatial import KDTree

class StructureAnalyzer:
    def __init__(self, tolerance=0.01):
        self.tolerance = tolerance
        self.cell = None
        self.coords = None
        self.frac_coords = None
        self.species = None
        self.dim = None
        self.structure_type = None
        self.kd_tree = None
        self.center = None
        self.inv_cell = None
        
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
                
            self._determine_dimension()
            self._build_kd_tree()
            
            # Pre-compute inverse cell for fractional coordinate transformations
            if self.cell is not None:
                self.inv_cell = np.linalg.inv(self.cell)
            
            return True
            
        except Exception as e:
            print("Error reading file: {}".format(e))
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
    
    def _determine_dimension(self):
        if self.cell is None or len(self.cell) != 3:
            self.dim = 3
            self.structure_type = "3D (Bulk)"
            return
            
        lengths = np.linalg.norm(self.cell, axis=1)
        periodic = lengths > self.tolerance
        non_periodic_dims = np.sum(~periodic)
        
        if non_periodic_dims == 3:
            self.dim = 0
            self.structure_type = "0D (Cluster/Molecule)"
        elif non_periodic_dims == 2:
            self.dim = 1
            self.structure_type = "1D (Chain/Nanowire)"
        elif non_periodic_dims == 1:
            self.dim = 2
            self.structure_type = "2D (Film/Layered)"
        else:
            self.dim = 3
            self.structure_type = "3D (Bulk)"
    
    def _build_kd_tree(self):
        if self.coords is not None and len(self.coords) > 0:
            self.kd_tree = KDTree(self.coords)
            self.center = np.mean(self.coords, axis=0)
    
    def _find_nearest_atom(self, position):
        if self.kd_tree is None:
            return -1, float('inf')
        dist, idx = self.kd_tree.query(position)
        if dist < self.tolerance:
            return idx, dist
        return -1, float('inf')
    
    def _wrap_frac(self, frac_pos):
        """Wrap fractional coordinates to [0, 1)"""
        return frac_pos - np.floor(frac_pos)
    
    def _check_rotation_frac(self, axis, angle):
        """Check rotation symmetry using fractional coordinates"""
        if self.cell is None or self.inv_cell is None:
            return False
            
        # Build rotation matrix in Cartesian
        c = np.cos(angle)
        s = np.sin(angle)
        axis = np.array(axis) / np.linalg.norm(axis)
        
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[0], axis[1], 0]])
        R_cart = np.eye(3) + s * K + (1 - c) * np.dot(K, K)
        
        # Convert rotation to fractional coordinates: R_frac = inv_cell @ R_cart @ cell
        R_frac = np.dot(self.inv_cell, np.dot(R_cart, self.cell))
        
        # For each atom, check if rotation maps it to another atom
        for frac_pos in self.frac_coords:
            # Rotate in fractional coordinates
            rotated_frac = np.dot(R_frac, frac_pos)
            # Wrap to [0, 1)
            rotated_frac = self._wrap_frac(rotated_frac)
            
            # Check all periodic images
            found = False
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        shift = np.array([i, j, k])
                        test_frac = rotated_frac + shift
                        test_cart = np.dot(test_frac, self.cell)
                        idx, dist = self._find_nearest_atom(test_cart)
                        if idx >= 0 and dist < self.tolerance:
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if not found:
                return False
        return True
    
    def _check_mirror_frac(self, normal):
        """Check mirror symmetry using fractional coordinates"""
        if self.cell is None or self.inv_cell is None:
            return False
            
        normal = np.array(normal) / np.linalg.norm(normal)
        # Reflection in Cartesian: R = I - 2*n*n^T
        R_cart = np.eye(3) - 2 * np.outer(normal, normal)
        
        # Convert to fractional
        R_frac = np.dot(self.inv_cell, np.dot(R_cart, self.cell))
        
        for frac_pos in self.frac_coords:
            reflected_frac = np.dot(R_frac, frac_pos)
            reflected_frac = self._wrap_frac(reflected_frac)
            
            found = False
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        shift = np.array([i, j, k])
                        test_frac = reflected_frac + shift
                        test_cart = np.dot(test_frac, self.cell)
                        idx, dist = self._find_nearest_atom(test_cart)
                        if idx >= 0 and dist < self.tolerance:
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if not found:
                return False
        return True
    
    def _check_inversion_frac(self):
        """Check inversion symmetry using fractional coordinates"""
        if self.cell is None or self.inv_cell is None or self.center is None:
            return False
            
        # Inversion center in fractional coordinates
        center_frac = np.dot(self.inv_cell, self.center)
        
        for frac_pos in self.frac_coords:
            inv_frac = 2 * center_frac - frac_pos
            inv_frac = self._wrap_frac(inv_frac)
            
            found = False
            for i in [-1, 0, 1]:
                for j in [-1, 0, 1]:
                    for k in [-1, 0, 1]:
                        shift = np.array([i, j, k])
                        test_frac = inv_frac + shift
                        test_cart = np.dot(test_frac, self.cell)
                        idx, dist = self._find_nearest_atom(test_cart)
                        if idx >= 0 and dist < self.tolerance:
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if not found:
                return False
        return True
    
    def analyze(self):
        """Analyze symmetry of the structure"""
        if self.coords is None or len(self.coords) == 0:
            return {}
        
        results = {
            'inversion': False,
            'rotations': [],
            'mirror_planes': [],
            'point_group': None
        }
        
        # Use fractional coordinate analysis if cell is available
        if self.cell is not None and self.inv_cell is not None:
            # Check inversion
            results['inversion'] = self._check_inversion_frac()
            
            # Check rotations along z-axis (c-axis)
            # For hexagonal systems, check 6-fold first
            for n in [6, 4, 3, 2]:
                angle = 2 * np.pi / n
                if self._check_rotation_frac([0, 0, 1], angle):
                    results['rotations'].append((n, 'c'))
                    # If 6-fold found, don't check lower orders along same axis
                    if n == 6:
                        break
            
            # Check rotations along a and b axes
            for axis, name in [([1, 0, 0], 'a'), ([0, 1, 0], 'b')]:
                for n in [2]:
                    angle = 2 * np.pi / n
                    if self._check_rotation_frac(axis, angle):
                        results['rotations'].append((n, name))
            
            # Check mirror planes
            # Perpendicular to z-axis (c)
            if self._check_mirror_frac([0, 0, 1]):
                results['mirror_planes'].append('c')
            
            # Perpendicular to a and b axes
            for normal, name in [([1, 0, 0], 'a'), ([0, 1, 0], 'b')]:
                if self._check_mirror_frac(normal):
                    results['mirror_planes'].append(name)
            
        else:
            # Fallback: Cartesian analysis
            results['inversion'] = self._check_inversion()
            
            for n in [6, 4, 3, 2]:
                angle = 2 * np.pi / n
                if self._check_rotation([0, 0, 1], angle):
                    results['rotations'].append((n, 'c'))
                    if n == 6:
                        break
            
            if self._check_mirror([0, 0, 1]):
                results['mirror_planes'].append('c')
        
        # Determine point group
        results['point_group'] = self._determine_point_group(results)
        
        return results
    
    def _check_rotation(self, axis, angle):
        """Check rotation symmetry in Cartesian coordinates (fallback)"""
        c = np.cos(angle)
        s = np.sin(angle)
        axis = np.array(axis) / np.linalg.norm(axis)
        
        K = np.array([[0, -axis[2], axis[1]],
                      [axis[2], 0, -axis[0]],
                      [-axis[0], axis[1], 0]])
        R = np.eye(3) + s * K + (1 - c) * np.dot(K, K)
        
        for pos in self.coords:
            rotated = np.dot(R, pos - self.center) + self.center
            idx, dist = self._find_nearest_atom(rotated)
            if idx < 0 or dist > self.tolerance:
                return False
        return True
    
    def _check_mirror(self, normal):
        """Check mirror symmetry in Cartesian coordinates (fallback)"""
        normal = np.array(normal) / np.linalg.norm(normal)
        R = np.eye(3) - 2 * np.outer(normal, normal)
        
        for pos in self.coords:
            reflected = np.dot(R, pos - self.center) + self.center
            idx, dist = self._find_nearest_atom(reflected)
            if idx < 0 or dist > self.tolerance:
                return False
        return True
    
    def _check_inversion(self):
        """Check inversion symmetry in Cartesian coordinates (fallback)"""
        if self.center is None:
            return False
        for pos in self.coords:
            inv_pos = 2 * self.center - pos
            idx, dist = self._find_nearest_atom(inv_pos)
            if idx < 0 or dist > self.tolerance:
                return False
        return True
    
    def _determine_point_group(self, results):
        """Determine point group from detected symmetries"""
        has_inv = results['inversion']
        rotations = results['rotations']
        mirrors = results['mirror_planes']
        
        rot_orders = [r[0] for r in rotations]
        has_6fold = 6 in rot_orders
        has_4fold = 4 in rot_orders
        has_3fold = 3 in rot_orders
        has_2fold = 2 in rot_orders
        
        has_mirror_c = 'c' in mirrors
        has_mirror_a = 'a' in mirrors
        has_mirror_b = 'b' in mirrors
        n_mirrors = len(mirrors)
        
        # Debug output
        print("\n[Debug] Detected symmetry elements:")
        print("  Rotations: {}".format(rotations))
        print("  Mirror planes: {}".format(mirrors))
        print("  Inversion: {}".format(has_inv))
        
        # Hexagonal system - C6v (6mm)
        if has_6fold:
            if has_inv:
                if has_mirror_c and has_mirror_a and has_mirror_b:
                    return ('D6h', '6/mmm', 191, 'Dihexagonal-dipyramidal')
                elif has_mirror_c:
                    return ('C6h', '6/m', 175, 'Hexagonal-dipyramidal')
                else:
                    return ('C6h', '6/m', 175, 'Hexagonal-dipyramidal')
            else:
                if has_mirror_c:
                    return ('C6v', '6mm', 183, 'Hexagonal-pyramidal')
                else:
                    return ('C6', '6', 169, 'Hexagonal-pyramidal')
        
        # Trigonal system
        elif has_3fold:
            if has_inv:
                if has_mirror_c and has_mirror_a:
                    return ('D3d', '3m', 164, 'Ditrigonal-scalenohedral')
                else:
                    return ('C3i', '3', 147, 'Rhombohedral')
            else:
                if has_mirror_c:
                    return ('C3v', '3m', 156, 'Trigonal-pyramidal')
                else:
                    return ('C3', '3', 143, 'Trigonal-pyramidal')
        
        # Tetragonal system
        elif has_4fold:
            if has_inv:
                if has_mirror_c and has_mirror_a:
                    return ('D4h', '4/mmm', 123, 'Ditetragonal-dipyramidal')
                elif has_mirror_c:
                    return ('C4h', '4/m', 87, 'Tetragonal-dipyramidal')
                else:
                    return ('C4h', '4/m', 87, 'Tetragonal-dipyramidal')
            else:
                if has_mirror_c:
                    return ('C4v', '4mm', 99, 'Tetragonal-pyramidal')
                else:
                    return ('C4', '4', 81, 'Tetragonal-pyramidal')
        
        # Orthorhombic system
        elif has_2fold and n_mirrors >= 2:
            if has_inv:
                return ('D2h', 'mmm', 47, 'Rhombic-dipyramidal')
            else:
                if has_mirror_c:
                    return ('C2v', 'mm2', 25, 'Rhombic-pyramidal')
                else:
                    return ('D2', '222', 16, 'Rhombic-trapezohedral')
        
        # Monoclinic system
        elif has_2fold:
            if has_inv:
                return ('C2h', '2/m', 10, 'Monoclinic-prismatic')
            else:
                if n_mirrors >= 1:
                    return ('C2v', 'mm2', 25, 'Rhombic-pyramidal')
                else:
                    return ('C2', '2', 3, 'Monoclinic-sphenoidal')
        
        # Triclinic system
        else:
            if has_inv:
                return ('Ci', '1', 2, 'Triclinic-pinacoidal')
            else:
                if n_mirrors >= 1:
                    return ('Cs', 'm', 6, 'Monoclinic-domatic')
                else:
                    return ('C1', '1', 1, 'Triclinic-pedial')
    
    def print_report(self):
        """Print symmetry analysis report"""
        results = self.analyze()
        
        print("\n" + "="*60)
        print("Structure Symmetry Analysis Report")
        print("="*60)
        print("Structure type: {}".format(self.structure_type))
        print("Dimension number: {}".format(self.dim))
        print("Total atoms: {}".format(len(self.coords)))
        print("Atomic species: {}".format(', '.join(set(self.species))))
        print("Tolerance: {} Å".format(self.tolerance))
        
        if self.cell is not None:
            print("\nLattice parameters (Å):")
            for i, vec in enumerate(self.cell):
                print("  a{}: [{:.6f}, {:.6f}, {:.6f}]".format(i+1, vec[0], vec[1], vec[2]))
        
        # Point group
        pg = results.get('point_group')
        if pg:
            schoenflies, hm, num, name = pg
            print("\n" + "="*60)
            print("Point Group Information:")
            print("  Schoenflies: {}".format(schoenflies))
            print("  Hermann-Mauguin: {}".format(hm))
            print("  International Number: {}".format(num))
            print("  Full Name: {}".format(name))
            print("="*60)
        
        print("\nSymmetry Analysis:")
        print("  Inversion symmetry: {}".format('Yes' if results.get('inversion') else 'No'))
        
        # Rotations
        if results.get('rotations'):
            print("\n  Rotational symmetries:")
            for n, axis in results['rotations']:
                print("    - {}fold rotation about {}-axis".format(n, axis))
        else:
            print("  Rotational symmetries: None detected")
        
        # Mirror planes
        if results.get('mirror_planes'):
            print("\n  Mirror planes:")
            for plane in results['mirror_planes']:
                print("    - Perpendicular to {}-axis".format(plane))
        else:
            print("  Mirror planes: None detected")
        
        print("="*60)


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
    
    tolerance_map = {
        '1': 0.001,
        '2': 0.01,
        '3': 0.1,
        '4': 0.2,
        '5': 0.5,
    }
    
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
    print("Structure Symmetry Analysis Tool")
    print("="*60)
    print("Supported formats: POSCAR, CONTCAR, QE input file")
    print("Tip: Press Tab for filename auto-completion")
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
