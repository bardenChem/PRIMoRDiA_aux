//geometry.h

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


class system;
class pdbAtom;
class pdbModel;

enum units{
	Ang,
	Bohr
};
/*********************************************************************/
class XYZ{
	public:
		unsigned int nAtoms;
		std::vector<double> xc;
		std::vector<double> yc;
		std::vector<double> zc;
		std::vector<std::string> elements;
		std::string commentary;
		XYZ();
		~XYZ();
		XYZ(const XYZ& rhs);
		XYZ& operator=(const XYZ& rhs);
		XYZ(XYZ&& rhs) noexcept;
		XYZ& operator=(XYZ&& rhs) noexcept;
		void write_xyz(std::string out_name);
		system get_molecule();
		void init_from_system(const system& molecule);
};
/*********************************************************************/
class PDB{
	public:
		unsigned int nModels;
		std::string PDB_ID;
		std::string basename;
		std::vector<pdbModel> models;
		bool Traj;
		bool NMR;
		PDB();
		~PDB();
		PDB(const PDB& rhs);
		PDB& operator=(const PDB& rhs);
		PDB(PDB&& rhs) noexcept;
		PDB& operator=(PDB&& rhs) noexcept;
		void split_models_in_files(std::string ref_name);
		void read_models_from_file(const char* pdb_file);
		void read_model_from_file(const char* pdb_file,unsigned int modN);
		void add_model(pdbModel model);
		void cat_pdbs(std::string file_list);
		void remove_model(unsigned int model);
		void write_pdb(std::string out_name);
		void init_from_system(const system& molecule);
		system get_system_from_model(unsigned int model);
		std::vector<system> get_systems();
		
};
/*********************************************************************/
class MOL2{
	public:
		MOL2();
		~MOL2();
		MOL2(const MOL2& rhs);
		MOL2& operator=(const MOL2& rhs);
		MOL2(MOL2&& rhs) noexcept;
		MOL2& operator=(MOL2&& rhs) noexcept;
};
/*********************************************************************/
class geometry{
	public:
		std::string type;
		XYZ xyz;
		PDB pdb;
		MOL2 mol2;
		system molecule;
		units cUnit;
		geometry();
		~geometry();
		geometry(const geometry& rhs);
		geometry& operator=(const geometry& rhs);
		geometry(geometry&& rhs) noexcept;
		geometry& operator=(geometry&& rhs) noexcept;
		void init_from_file(const char* file_name, std:string Type);
		void read_QCPinput(const char* file_name, std:string program);
		void read_QCPoutput(const char* file_name, std:string program,bool last);
		void convert_to_ang();
		void convert_to_bohr();
		void write_to_file(std::string out_name,std::string format);
		void center_coord();
};