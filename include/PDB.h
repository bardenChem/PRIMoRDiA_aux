//PDB.h

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

#ifndef PDB_H
#define PDB_H

#include <vector>
#include <string>
#include <fstream>

class pdbModel;
class molecule;

class PDB{
	public:
		unsigned int nModels;
		std::string PDB_ID;
		std::string basename;
		std::vector<pdbModel> models;
		bool MULTI;
		PDB();
		~PDB();
		PDB(const char* pdb_file);
		PDB(const PDB& rhs);
		PDB& operator=(const PDB& rhs);
		PDB(PDB&& rhs) noexcept;
		PDB& operator=(PDB&& rhs) noexcept;
		void split_models_in_files();
		void add_model(pdbModel model);
		void cat_pdbs(std::vector<std::string> file_list);
		void remove_model(unsigned int model);
		void write_pdb(std::string out_name);
		void init_from_system(const molecule& mol);
		molecule get_system_from_model(unsigned int model);
		std::vector<molecule> get_systems();
		friend std::ostream& operator<<(std::ostream& out, const PDB& obj);
		void print();

};
/******************************************************************************/
void UnitTest_PDB();

#endif