//atom.h

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

#ifndef ATOM_H
#define ATOM_H

#include <string>

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
		friend std::ostream& operator<<(std::ostream& out, const atom& obj);
		void print();
};		



/**************************************************************/
void UnitTest_atom();


#endif