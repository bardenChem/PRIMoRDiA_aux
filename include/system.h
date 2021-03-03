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

enum AAres{
	OTH=-1,
	ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE,
	PRO, SER, THR, TRP, TYR, VAL	
};
/*********************************************************************/
enum res_type { 
	UNK,
	WAT,
	HOH,
	ION,
	AA,
	DNA 
};
	
/*********************************************************************/
class atom{
	public:
		double xc;
		double yc;
		double zc;
		float pCharge;
		float aMass;
		unsigned int aNmb;
		std::string element;
		atom();
		~atom();
		atom(double x, double y, double z, std::string type);
		atom(const atom& rhs);
		atom& operator=(const atom& rhs);
		atom(atom&& rhs) noexcept;
		atom& operator=(atom&& rhs) noexcept;
		void set_element( std::string Type);
		void set_coord(double x, double y, double z);
		double get_distance(const atom& a2);
		void set_pCharge(double chg);
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
		system(std::vector<atom> ats, std::string nme);
		~system();
		system(const system& rhs);
		system& operator=(const system& rhs);
		system(system&& rhs) noexcept;
		system& operator=(system&& rhs) noexcept;
		void add_atom(atom a);
		void add_atom(double x,double y, double z, std::string el);
		void remove_atom(unsigned int i);
		
};
/*********************************************************************/
class pdbAtom{
	public:
		std::string atom_name;
		unsigned int indx;
		std::string res_name;
		int res_indx;
		std::string chain_name;
		float occupancy;
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
		pdbAtom(std::string& pdb_line);
		bool operator==(const pdbAtom& rhs);
		bool is_hydrogen();
		
};
/*********************************************************************/
class residue{
	public:
		std::string res3n;
		std::string res1n;
		res_type type;
		AAres AAname;
		bool ligand;
		bool terminal;
		bool first;
		int fCharge;
		int pdb_index;
		unsigned int nHydrogens;
		unsigned int nAtoms;
		std::vector<pdbAtom> r_atoms;		
		residue();
		residue(std::vector<pdbAtom> resAtoms, int resTyp, int resMon);
		~residue();
		residue(const residue& rhs);
		residue& operator=(const residue& rhs);
		residue( residue&& rhs) noexcept;
		residue& operator=( residue&& rhs) noexcept;
		void set_charge();
		bool is_ion();
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
		pdbModel()
		pdbModel(std::vector<residue> residues);
		pdbModel(const char* pdb_file, int mdl);
		~pdbModel();
		pdbModel(const pdbModel& rhs);
		pdbModel& operator=(const pdbModel& rhs);
		pdbModel(pdbModel&& rhs) noexcept;
		pdbModel& operator=(pdbModel&& rhs) noexcept;
		void write_model(std::string out_name);
		void prune_atoms();
		void remove_waters();
		void remove_waters(double radius, unsigned int res);
		void remove_ions();
		void remove_atom(unsigned int res, unsigned int at);
		void remove_residue(unsigned int i);
		void split_complex(std::string mol);
		void built_complex(const char* pdb_mol);
};
/*********************************************************************/
#endif