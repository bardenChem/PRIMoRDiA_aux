//molecule.h

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

#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <string>

class atom;

/********************************************************************/
class molecule{
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
		molecule();
		molecule(std::vector<atom> ats, std::string nme);
		~molecule();
		molecule(const molecule& rhs);
		molecule& operator=(const molecule& rhs);
		molecule(molecule&& rhs) noexcept;
		molecule& operator=(molecule&& rhs) noexcept;
		void add_atom(atom a);
		void add_atom(double x,double y, double z, std::string el);
		void remove_atom(unsigned int i);
};
/********************************************************************/
#endif