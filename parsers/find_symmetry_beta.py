#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Structure Symmetry Analysis Tool
Supported formats: POSCAR, CONTCAR, QE input file
"""

import numpy as np
from pathlib import Path
import sys
import readline
import glob
from typing import Tuple, List, Optional
import itertools

class StructureAnalyzer:
    def __init__(self, tolerance=0.01):
        """
        Initialize structure analyzer
        
        Parameters:
        -----------
        tolerance : float
            Tolerance for symmetry detection (Angstrom)
        """
        self.tolerance = tolerance
        self.cell = None
        self.coords = None
        self.species = None
        self.dim = None
        self.structure_type = None
        
    def read_file(self, filename):
        # type: (str) -> bool
        """
        Read structure file, auto-detect format
        
        Returns:
        --------
        bool : True if read successfully
        """
        try:
            with open(filename, 'r') as f:
                lines = f.readlines()
                
            # Detect file format
            if self._is_poscar(lines):
                self._read_poscar(lines)
            elif self._is_qe(lines):
                self._read_qe(lines)
            else:
                print("Unrecognized file format")
                return False
                
            # Determine dimension
            self._determine_dimension()
            return True
            
        except Exception as e:
            print("Error reading file: {}".format(e))
            return False
    
    def _is_poscar(self, lines):
        # type: (List[str]) -> bool
        """Check if the file is in POSCAR/CONTCAR format"""
        if len(lines) < 6:
            return False
        # POSCAR format characteristics: line 1 comment, line 2 scale factor, lines 3-5 lattice vectors
        try:
            # Try to parse line 2 (scale factor)
            float(lines[1].strip().split()[0])
            # Check if lines 3-5 are three lattice vectors
            for i in range(2, 5):
                parts = lines[i].strip().split()
                if len(parts) < 3:
                    return False
                [float(x) for x in parts[:3]]
            return True
        except:
            return False
    
    def _is_qe(self, lines):
        # type: (List[str]) -> bool
        """Check if the file is in QE input format"""
        text = ''.join(lines).lower()
        qe_keywords = ['&system', '&control', 'atomic_species', 'atomic_positions', 
                      'cell_parameters', 'ibrav']
        return any(keyword in text for keyword in qe_keywords)
    
    def _read_poscar(self, lines):
        # type: (List[str]) -> None
        """Read POSCAR/CONTCAR file"""
        # Skip first comment line
        idx = 1
        # Scale factor
        scale = float(lines[idx].strip().split()[0])
        idx += 1
        
        # Read lattice vectors
        self.cell = []
        for i in range(3):
            parts = lines[idx + i].strip().split()
            vector = [float(x) * scale for x in parts[:3]]
            self.cell.append(vector)
        idx += 3
        self.cell = np.array(self.cell)
        
        # Read species names and counts
        species_names = lines[idx].strip().split()
        idx += 1
        species_counts = [int(x) for x in lines[idx].strip().split()]
        idx += 1
        
        # Check for selective dynamics
        if lines[idx].strip().lower().startswith('s'):
            idx += 1
        
        # Determine coordinate type (Direct/Cartesian)
        coord_type = lines[idx].strip().lower()
        idx += 1
        is_direct = coord_type.startswith('d')
        
        # Read coordinates
        total_atoms = sum(species_counts)
        self.coords = []
        self.species = []
        
        for i, count in enumerate(species_counts):
            for j in range(count):
                parts = lines[idx].strip().split()
                pos = [float(x) for x in parts[:3]]
                if is_direct:
                    # Convert fractional to Cartesian coordinates
                    pos = np.dot(pos, self.cell)
                self.coords.append(pos)
                self.species.append(species_names[i])
                idx += 1
        
        self.coords = np.array(self.coords)
    
    def _read_qe(self, lines):
        # type: (List[str]) -> None
        """Read QE input file"""
        self.cell = None
        self.coords = []
        self.species = []
        
        # Get species mapping
        species_map = {}
        in_atomic_species = False
        
        # First read atomic species
        for line in lines:
            line = line.strip()
            if 'ATOMIC_SPECIES' in line.upper():
                in_atomic_species = True
                continue
            if in_atomic_species:
                if not line or line.startswith('!'):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    species_map[parts[0]] = parts[0]
                if 'ATOMIC_POSITIONS' in line.upper() or 'K_POINTS' in line.upper():
                    break
        
        # Read lattice parameters
        for i, line in enumerate(lines):
            line = line.strip()
            if 'CELL_PARAMETERS' in line.upper():
                # Find lattice parameter block
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
        
        # If lattice parameters not found, try to get from ibrav
        if self.cell is None:
            self.cell = self._get_cell_from_ibrav(lines)
        
        # Read atomic positions
        in_positions = False
        for line in lines:
            line = line.strip()
            if 'ATOMIC_POSITIONS' in line.upper():
                in_positions = True
                parts = line.split()
                # Check coordinate type
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
                        pos = np.dot(pos, self.cell)
                    elif coord_type == 'alat' and self.cell is not None:
                        pos = np.dot(pos, self.cell)
                    
                    self.coords.append(pos)
                    self.species.append(species_name)
        
        self.coords = np.array(self.coords)
    
    def _get_cell_from_ibrav(self, lines):
        # type: (List[str]) -> Optional[np.ndarray]
        """Get lattice parameters from ibrav parameter"""
        ibrav = None
        celldm = [1.0] * 6  # Default values
        
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
            
        # Generate lattice matrix based on ibrav (common cases)
        a = celldm[0]
        if ibrav == 1:  # Cubic
            return np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])
        elif ibrav == 2:  # FCC
            return np.array([[0, a/2, a/2], [a/2, 0, a/2], [a/2, a/2, 0]])
        elif ibrav == 3:  # BCC
            return np.array([[-a/2, a/2, a/2], [a/2, -a/2, a/2], [a/2, a/2, -a/2]])
        # More ibrav values can be added here
        
        return None
    
    def _determine_dimension(self):
        # type: () -> None
        """Determine the dimension of the structure"""
        if self.cell is None or len(self.cell) != 3:
            self.dim = 3
            self.structure_type = "3D (Bulk)"
            return
            
        # Check lattice vector lengths and angles
        lengths = np.linalg.norm(self.cell, axis=1)
        
        # Check which dimensions are periodic (length > tolerance)
        periodic = lengths > self.tolerance
        
        # Check non-periodic dimensions
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
    
    def check_translation_symmetry(self, points):
        # type: (np.ndarray) -> bool
        """Check translation symmetry"""
        if len(points) == 0:
            return False
            
        # Check for identical points (considering periodic boundary conditions)
        for i in range(len(points)):
            for j in range(i+1, len(points)):
                diff = points[i] - points[j]
                # Consider periodicity
                for k in range(self.dim):
                    diff[k] = diff[k] - np.round(diff[k] / self.cell[k][k]) * self.cell[k][k]
                
                if np.linalg.norm(diff) < self.tolerance:
                    return True
        return False
    
    def find_symmetry_operations(self):
        # type: () -> dict
        """
        Find symmetry operations (simplified version)
        Mainly detects inversion center, mirror planes, rotations, etc.
        """
        if self.coords is None or len(self.coords) == 0:
            return {}
            
        operations = {}
        center = np.mean(self.coords, axis=0)
        
        # 1. Check inversion symmetry
        is_centrosymmetric = True
        for pos in self.coords:
            inv_pos = 2 * center - pos
            # Check if corresponding inversion position exists
            found = False
            for other in self.coords:
                if np.linalg.norm(inv_pos - other) < self.tolerance:
                    found = True
                    break
            if not found:
                is_centrosymmetric = False
                break
        operations['Inversion'] = is_centrosymmetric
        
        # 2. Check mirror symmetry along principal axes (simplified)
        mirror_planes = []
        axes = ['x', 'y', 'z']
        for i, axis in enumerate(axes[:self.dim]):
            is_mirror = True
            for pos in self.coords:
                mirrored = pos.copy()
                mirrored[i] = -pos[i] + 2 * center[i]
                found = False
                for other in self.coords:
                    if np.linalg.norm(mirrored - other) < self.tolerance:
                        found = True
                        break
                if not found:
                    is_mirror = False
                    break
            if is_mirror:
                mirror_planes.append("Mirror plane (perpendicular to {}-axis)".format(axis))
        operations['Mirror planes'] = mirror_planes
        
        # 3. Check rotational symmetry (2-fold, 3-fold, 4-fold, 6-fold)
        rotations = self._find_rotational_symmetry(center)
        if rotations:
            operations['Rotational symmetry'] = rotations
        
        return operations
    
    def _find_rotational_symmetry(self, center):
        # type: (np.ndarray) -> List[str]
        """Find rotational symmetry axes"""
        rotations = []
        
        # Check for 2-fold, 3-fold, 4-fold, 6-fold rotation about z-axis
        for n in [2, 3, 4, 6]:
            angle = 2 * np.pi / n
            rot_matrix = np.array([
                [np.cos(angle), -np.sin(angle), 0],
                [np.sin(angle), np.cos(angle), 0],
                [0, 0, 1]
            ])
            
            if self._check_rotation_symmetry(rot_matrix, center):
                rotations.append("{}fold rotation (about z-axis)".format(n))
        
        # Check for rotation about x-axis
        for n in [2]:
            angle = 2 * np.pi / n
            rot_matrix = np.array([
                [1, 0, 0],
                [0, np.cos(angle), -np.sin(angle)],
                [0, np.sin(angle), np.cos(angle)]
            ])
            
            if self._check_rotation_symmetry(rot_matrix, center):
                rotations.append("{}fold rotation (about x-axis)".format(n))
        
        # Check for rotation about y-axis
        for n in [2]:
            angle = 2 * np.pi / n
            rot_matrix = np.array([
                [np.cos(angle), 0, np.sin(angle)],
                [0, 1, 0],
                [-np.sin(angle), 0, np.cos(angle)]
            ])
            
            if self._check_rotation_symmetry(rot_matrix, center):
                rotations.append("{}fold rotation (about y-axis)".format(n))
        
        return rotations
    
    def _check_rotation_symmetry(self, rot_matrix, center):
        # type: (np.ndarray, np.ndarray) -> bool
        """Check if structure has a given rotational symmetry"""
        for pos in self.coords:
            # Rotate position around center
            rotated = np.dot(rot_matrix, pos - center) + center
            # Check if rotated position exists
            found = False
            for other in self.coords:
                if np.linalg.norm(rotated - other) < self.tolerance:
                    found = True
                    break
            if not found:
                return False
        return True
    
    def analyze(self):
        # type: () -> dict
        """Comprehensive structure symmetry analysis"""
        result = {
            'Dimension': self.structure_type,
            'Dimension number': self.dim,
            'Total atoms': len(self.coords),
            'Atomic species': list(set(self.species)),
            'Lattice parameters': self.cell.tolist() if self.cell is not None else None,
            'Tolerance': self.tolerance,
        }
        
        # Find symmetry operations
        sym_ops = self.find_symmetry_operations()
        result.update(sym_ops)
        
        return result
    
    def print_report(self):
        # type: () -> None
        """Print symmetry analysis report"""
        result = self.analyze()
        
        print("\n" + "="*60)
        print("Structure Symmetry Analysis Report")
        print("="*60)
        print("Structure type: {}".format(result['Dimension']))
        print("Dimension number: {}".format(result['Dimension number']))
        print("Total atoms: {}".format(result['Total atoms']))
        print("Atomic species: {}".format(', '.join(result['Atomic species'])))
        print("Tolerance: {} Å".format(result['Tolerance']))
        
        if result['Lattice parameters']:
            print("\nLattice parameters (Å):")
            for i, vec in enumerate(result['Lattice parameters']):
                print("  a{}: [{:.6f}, {:.6f}, {:.6f}]".format(i+1, vec[0], vec[1], vec[2]))
        
        print("\nSymmetry Analysis:")
        print("  Inversion symmetry: {}".format('Yes' if result.get('Inversion', False) else 'No'))
        
        if result.get('Mirror planes'):
            print("  Mirror planes:")
            for plane in result['Mirror planes']:
                print("    - {}".format(plane))
        else:
            print("  Mirror planes: None detected")
        
        if result.get('Rotational symmetry'):
            print("  Rotational symmetry:")
            for rot in result['Rotational symmetry']:
                print("    - {}".format(rot))
        else:
            print("  Rotational symmetry: None detected")
        
        # Point group summary (simplified)
        self._suggest_point_group(result)
        print("="*60)
    
    def _suggest_point_group(self, result):
        # type: (dict) -> None
        """Suggest possible point group (simplified)"""
        has_inversion = result.get('Inversion', False)
        mirror_planes = result.get('Mirror planes', [])
        rotations = result.get('Rotational symmetry', [])
        
        if has_inversion and len(mirror_planes) == 3 and len(rotations) == 0:
            print("\nPossible point group: Ci (S2)")
        elif has_inversion and len(mirror_planes) == 3 and len(rotations) == 3:
            print("\nPossible point group: Oh (Octahedral)")
        elif has_inversion and len(mirror_planes) == 1 and len(rotations) == 0:
            print("\nPossible point group: C2h")
        elif has_inversion and len(mirror_planes) == 0 and len(rotations) == 0:
            print("\nPossible point group: Ci")
        elif len(mirror_planes) == 1 and len(rotations) == 0:
            print("\nPossible point group: Cs")
        elif len(rotations) == 1 and "2-fold" in rotations[0]:
            if len(mirror_planes) > 0:
                print("\nPossible point group: C2v")
            else:
                print("\nPossible point group: C2")
        elif len(rotations) == 1 and "3-fold" in rotations[0]:
            if len(mirror_planes) > 0:
                print("\nPossible point group: C3v")
            else:
                print("\nPossible point group: C3")
        elif len(rotations) == 1 and "4-fold" in rotations[0]:
            if len(mirror_planes) > 0:
                print("\nPossible point group: C4v")
            else:
                print("\nPossible point group: C4")
        else:
            print("\nPoint group: Not determined (may be C1 or lower symmetry)")


def tab_complete_filename(text, state):
    """Tab completion for filenames"""
    matches = glob.glob(text + '*') if text else glob.glob('*')
    # Return only files (exclude directories)
    matches = [m for m in matches if Path(m).is_file()]
    if state < len(matches):
        return matches[state]
    return None


def get_tolerance():
    # type: () -> float
    """Interactive tolerance selection"""
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
        choice = raw_input("\nEnter option (1-6): ").strip()
        if choice in tolerance_map:
            return tolerance_map[choice]
        elif choice == '6':
            try:
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
    # type: () -> None
    """Main program"""
    print("="*60)
    print("Structure Symmetry Analysis Tool")
    print("="*60)
    print("Supported formats: POSCAR, CONTCAR, QE input file")
    print("Tip: Press Tab for filename auto-completion")
    print("="*60)
    
    # Setup tab completion
    readline.set_completer(tab_complete_filename)
    readline.parse_and_bind('tab: complete')
    
    # Get filename
    while True:
        filename = raw_input("\nEnter structure filename: ").strip()
        if not filename:
            print("Filename cannot be empty")
            continue
        
        if not Path(filename).exists():
            print("File '{}' does not exist, please try again".format(filename))
            continue
        
        break
    
    # Get tolerance
    tolerance = get_tolerance()
    
    # Analyze structure
    analyzer = StructureAnalyzer(tolerance=tolerance)
    
    print("\nReading file: {}".format(filename))
    if analyzer.read_file(filename):
        analyzer.print_report()
    else:
        print("\nAnalysis failed!")


if __name__ == "__main__":
    main()
