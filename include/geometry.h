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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <string>

#include "../include/PDB.h"
#include "../include/XYZ.h"
#include "../include/molecule.h"


//===============================================
enum units	{ Ang, Bohr };
//----------------------------------
enum GEO_file{ INVALID, xyz_, mol2_	, GRO_, pdb_ };
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
		GEO_file Type;
		XYZ xyz;
		PDB pdb;
		MOL2 mol2;
		molecule Molecule;
		units cUnit;
		geometry();
		geometry(const char* file_name);
		~geometry();
		geometry(const geometry& rhs);
		geometry& operator=(const geometry& rhs);
		geometry(geometry&& rhs) noexcept;
		geometry& operator=(geometry&& rhs) noexcept;
		void read_QCPinput(const char* file_name, std::string program);
		void read_QCPoutput(const char* file_name, std::string program,bool last);
		void convert_to_ang();
		void convert_to_bohr();
		void write_to_file(std::string out_name,std::string format);
		void center_coord();
};
///////////////////////////////////////////////////////////
#endif