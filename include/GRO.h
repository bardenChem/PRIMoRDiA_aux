//GRO.h

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

#ifndef GRO_H
#define GRO_H

#include <vector>
#include <string>
#include <fstream>

class molecule;
class pdbAtom;
class pdbModel;

class GRO{
	public:
		unsigned int nAtoms;
		std::vector<double> sides;
		std::vector<pdbAtom> atoms;
		std::string basename;
		GRO();
		~GRO();
		GRO(const char* gro_file);
		GRO(const GRO& rhs);
		GRO& operator=(const GRO& rhs);
		GRO(GRO&& rhs) noexcept;
		GRO& operator=(GRO&& rhs) noexcept;
		void write_gro(std::string out_name);
		void init_from_system(const molecule& mol);
		molecule get_system_from_gro(unsigned int model);
		std::vector<molecule> get_systems();
		pdbModel get_pdb_from_gro();
		friend std::ostream& operator<<(std::ostream& out, const GRO& obj);
		void print();
};
/******************************************************************************/
void UnitTest_GRO();

#endif