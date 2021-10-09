//pdbModel.h

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

#ifndef PDB_MODEL_H
#define PDB_MODEL_H

#include <string>
#include <vector>

class pdbAtom;
class residue;
/***********************************************************/
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
		void remove_chain( unsigned nchain );
		std::vector<unsigned> spherical_selection(unsigned center_atom, double radius, bool within );
		void split_complex( std::string mol );
		void built_complex( const char* pdb_mol );
		double get_distance( int a1, int a2 );
};
/***********************************************************/
#endif 