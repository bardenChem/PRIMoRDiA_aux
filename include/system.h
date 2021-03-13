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

class atom;

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
/********************************************************************/
#endif