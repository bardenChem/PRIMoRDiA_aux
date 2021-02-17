//system.h

/*********************************************************************/
/* This source code file is part of LQQCMMtools software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
/*********************************************************************/

#ifndef SYSTEM
#define SYSTEM

#include <vector>
#include <string>

/*********************************************************************/
class atom{
	public:
		double xc;
		double yc;
		double zc;
		float pCharge;
		double aMass;
		std::string element;
		unsigned int nmb;
		atom();
		~atom();
		atom(double xc, double yc, double zc, std::string type);
		atom(const atom& rhs);
		atom& operator=(const atom& rhs);
		atom(atom&& rhs) noexcept;
		atom& operator=(atom&& rhs) noexcept;
		void set_element( std::string type);
		void set_coord(double xc, double yc, double zc);
		double get_distance(const atom& a2);
};
/********************************************************************/
class system{
	public:
		unsigned int nAtoms;
		unsigned int nElectrons;
		std::string name;
		std::string type;
		std::vector<atom> atoms;
		int fCharge;
		double molar_mass;
		double ver_inf[3];
		double ver_sup[3];
		system();
		~system();
		system(const system& rhs);
		system& operator=(const system& rhs);
		system(system&& rhs) noexcept;
		system& operator=(system&& rhs) noexcept;
		
};
/*********************************************************************/
class pdbAtom{
	public:
		std::string atom_name;
		unsigned int indx;
		std::string res_name;
		std::string chain_name;
		double occupancy;
		double b_factor;
		bool sideC;
		double xc;
		double yc;
		double zc;
		pdbAtom();
		~pdbAtom();
		pdbAtom(const pdbAtom& rhs);
		pdbAtom& operator=(const pdbAtom& rhs);
		pdbAtom(pdbAtom&& rhs) noexcept;
		pdbAtom& operator=(pdbAtom&& rhs) noexcept;
		bool pdbAtom operator!=(const pdbAtom& rhs);
};
/*********************************************************************/
class residue{
	public:
		std::string type;
		bool ligand;
		bool terminal;
		bool WAT;
		bool ION;
		int fCharge;
		unsigned int nHydrogens;
		unsigned int nAtoms;
		std::vector<pdbAtom> r_atoms;		
		residue();
		~residue();
		residue(const residue& rhs);
		residue& operator=(const residue& rhs);
		residue( residue&& rhs) noexcept;
		residue& operator=( residue&& rhs) noexcept;
		void set_charge();
};
/*********************************************************************/
class pdbModel{
	public:
		std::string index;
		std::string remark;
		std::string title;
		std::vector<residue> monomers;
		unsigned int model;
		unsigned int nChains;
		unsigned int nResidues;
		unsigned int nAtoms;
		pdbModel();
		pdbModel(std::vector<residue> residues);
		~pdbModel();
		pdbModel(const pdbModel& rhs);
		pdbModel& operator=(const pdbModel& rhs);
		pdbModel(pdbModel&& rhs) noexcept;
		pdbModel& operator=(pdbModel&& rhs) noexcept;
		void write_model(std::string out_name);
		void prune_atoms();
		void remove_waters();
		void remove_waters(double radius);
		void remove_ions();
		void remove_residue(unsigned int i);
		void update_residues();
};
/*********************************************************************/
class PDB{
	public:
		unsigned int nModels;
		std::string PDB_ID;
		std::vector<pdbModel> models;
		bool Traj;
		bool NMR;
		PDB();
		PDB(std::vector<pdbModel> mods);
		PDB(const char* pdb_file);
		PDB(std::vector< std::string> > pdb_file_list);
		~PDB();
		PDB(const PDB& rhs);
		PDB& operator=(const PDB& rhs);
		PDB(PDB&& rhs) noexcept;
		PDB& operator=(PDB&& rhs) noexcept;
		void split_models_in_files();
		void read_models_from_file(const char* pdb_file);
		void read_model_from_file(const char* pdb_file,unsigned int modN);
		void add_model(pdbModel model);
		void remove_model(unsigned int model);
		void write_pdb;
		system get_system_from_model(unsigned int model);
		std::vector<system> get_systems();
		
};
/*********************************************************************/
#endif