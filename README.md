# PRIMoRDiA Auxiliary tools

Laboratório de Química Quântica e Computcional Molecular modelling tools is a application that reunites a set of functionalities to prepare and extract information of and from molecular simulations. Also is the designed tool to prepare all quantum mechanical simulations of biological systems for the PRIMoRDiA softwre.

---

## Table of Contents
1. [Building the Software](#building-the-software)
2. [Usage Guide](#usage-guide)
3. [Available Functions](#available-functions)
4. [Examples](#examples)

---

## Building the Software

### Prerequisites
- CMake 3.10 or higher
- C++17 compiler
- Standard development tools (make, gcc/clang)

### Build Instructions

```bash
cd /path/to/PRIMoRDiA_aux
mkdir build
cd build
cmake ..
make
```

The compiled executable will be located at `build/primordia_aux`.

---

## Usage Guide

The PRIMoRDiA Auxiliary software is command-line based. General syntax:

```bash
./primordia_aux [FLAG] [ARGUMENTS]
```

For help information:
```bash
./primordia_aux -h
# or
./primordia_aux --help
```

---

## Available Functions

### 1. **-rmsd_rg**: RMSD and Radius of Gyration Analysis
Calculate Root Mean Square Deviation (RMSD) and Radius of Gyration (RG) for a molecular dynamics trajectory using MDTraj.

**Usage:**
```bash
./primordia_aux -rmsd_rg <trajectory_file>
```

**Arguments:**
- `<trajectory_file>`: Path to the trajectory file (supported formats: XTC, DCD, etc.)

**Output:** 
- RMSD and RG values calculated and written to output files

---

### 2. **-atom_BID**: Distance Analysis Between Atom Pairs
Analyze distances between specific atoms in a PDB trajectory file.

**Usage:**
```bash
./primordia_aux -atom_BID <trajectory_file> <atom_id_1> [atom_id_2] [atom_id_3] ...
```

**Arguments:**
- `<trajectory_file>`: Path to the PDB trajectory file
- `<atom_id_N>`: Atom indices to calculate distances (one or more)

**Output:**
- Distance values between specified atoms over the trajectory

---

### 3. **-ext_frame_dcd**: Extract Sampled Frames from DCD Trajectory
Extract frames at regular intervals from a DCD trajectory file and save as PDB.

**Usage:**
```bash
./primordia_aux -ext_frame_dcd <trajectory_file> <topology_file> <interval>
```

**Arguments:**
- `<trajectory_file>`: Path to the DCD trajectory file
- `<topology_file>`: Path to the topology/structure file
- `<interval>`: Sampling interval (extract every Nth frame)

**Output:**
- `sampled.pdb`: PDB file containing sampled frames

---

### 4. **-ext_frame**: Extract Single or Multiple Frames from Trajectory
Extract specific frame(s) from a trajectory file in any supported format.

**Usage:**

**Extract single frame:**
```bash
./primordia_aux -ext_frame <trajectory_file> <frame_number>
```

**Extract multiple frames:**
```bash
./primordia_aux -ext_frame <trajectory_file> <start_frame> <end_frame>
```

**Arguments:**
- `<trajectory_file>`: Path to the trajectory file
- `<frame_number>`: Frame index to extract (single frame mode)
- `<start_frame>`: Starting frame index (multiple frames mode)
- `<end_frame>`: Ending frame index (multiple frames mode, use 0 for single frame)

**Output:**
- Extracted frame(s) saved as PDB file(s)

---

### 5. **-spherical_prune**: Spherical Pruning
Apply spherical pruning to a structure file based on a specified radius.

**Usage:**
```bash
./primordia_aux -spherical_prune <structure_file> <radius>
```

**Arguments:**
- `<structure_file>`: Path to the PDB structure file
- `<radius>`: Pruning radius (in Angstroms)

**Output:**
- Pruned structure file

---

### 6. **-analysis_mol_from_center**: Molecular Analysis from Center
Analyze molecules within a specified radius from a central residue across a trajectory.

**Usage:**
```bash
./primordia_aux -analysis_mol_from_center <trajectory_file> <topology_file> <residue_index> <molecule_name> <radius> [prune_distance]
```

**Arguments:**
- `<trajectory_file>`: Path to the trajectory file (XTC, DCD, etc.)
- `<topology_file>`: Path to the topology/GRO file
- `<residue_index>`: Index of the central residue
- `<molecule_name>`: Name of the molecule to analyze (e.g., "CO2", "H2O")
- `<radius>`: Search radius in Angstroms
- `[prune_distance]` (optional): Distance threshold for pruning

**Output:**
- Analysis results of molecular interactions around the central residue

---

### 7. **-QCP_inp**: Prepare Quantum Chemistry Package Input Files
Generate input files for quantum chemistry calculations from geometry files.

**Usage:**
```bash
./primordia_aux -QCP_inp <package> [OPTIONS]
```

**Arguments:**
- `<package>`: Target quantum chemistry package
  - `gamess`: GAMESS package
  - `orca`: ORCA package (requires `-NP <num_processes>`)
  - `mopac`: MOPAC package

**Options:**
- `-gExt <extension>`: Geometry file extension (default: .xyz)
- `-QMm <method>`: QM method (default: am1)
- `-basis <basis_set>`: Basis set (default: am1)
- `-FD`: Use finite difference mode
- `-chg <charge>`: Molecular charge (default: 0)
- `-bmulti <multiplicity>`: Spin multiplicity (default: 1)
- `-runOP <operation>`: Type of calculation (default: Energy)
  - `Energy`: Single point energy calculation
  - `Frequency`: Vibrational frequency analysis
  - `Optimization`: Geometry optimization
- `-chgDiff <value>`: Charge difference (default: 1)
- `-markCharge`: Mark specific charge on residue
- `-markResidue <residue_name>`: Residue name to mark charge
- `-topology <file>`: Topology file for marked charges
- `-NP <processes>`: Number of processes for ORCA (required for ORCA)

**Output:**
- Input files for the specified quantum chemistry package

**Example:**
```bash
./primordia_aux -QCP_inp orca -QMm DFT -basis 6-31G -chg 0 -bmulti 1 -NP 4
```

---

### 8. **-test**: Run Software Tests
Execute basic functionality tests for the software.

**Usage:**
```bash
./primordia_aux -test
```

**Output:**
- Test results and diagnostic information

---

### 9. **-rH2O_traj**: Remove Waters from Trajectory (Beta)
Remove water molecules from a trajectory. Currently supports processing PDB files in the current folder.

**Usage:**
```bash
./primordia_aux -rH2O_traj pdbs_folder <parameter1> <parameter2>
```

**Arguments:**
- `pdbs_folder`: Keyword indicating to process PDB files in current directory
- `<parameter1>`: First parameter for water removal
- `<parameter2>`: Second parameter for water removal

**Output:**
- Trajectory files with water molecules removed

---

## Under Development

The following features are still in development:
- `-MDprep`: Prepare files for molecular dynamics simulations
- `-check_QCP`: Check quantum chemistry package output files

---

## Examples

### Example 1: Calculate RMSD and RG for a trajectory
```bash
./primordia_aux -rmsd_rg trajectory.xtc
```

### Example 2: Analyze distances between specific atoms
```bash
./primordia_aux -atom_BID trajectory.pdb 100 150 200
```

### Example 3: Extract frames 0-100 from trajectory
```bash
./primordia_aux -ext_frame trajectory.dcd 0 100
```

### Example 4: Prepare ORCA input files with DFT/6-31G
```bash
./primordia_aux -QCP_inp orca -QMm DFT -basis 6-31G -chg 0 -bmulti 1 -NP 8
```

### Example 5: Analyze CO2 molecules around residue 466
```bash
./primordia_aux -analysis_mol_from_center trajectory.xtc topology.gro 466 CO2 4.2 30.0
```

---

## Author Information

**Created by:** Igor Barden Grillo  
**Institution:** Federal University of Paraíba (UFPB) / Pontifícia Universidade Católica do Rio Grande do Sul (PUCRS)  
**Contact:** 
- Personal: barden.igor@gmail.com
- Academic: igor.grillo@acad.pucrs.br  
**Group Site:** quantum-chem.pro.br  
**GitHub:** IgorChem

---

## License

This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.  
If a copy of the MPL was not distributed with this file, you can obtain one at https://mozilla.org/MPL/2.0/.
