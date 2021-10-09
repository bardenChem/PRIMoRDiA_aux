//unit_test.h

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

#ifndef UNIT_TEST
#define UNIT_TEST

#include <fstream>

class UnitTests{
	public:
		UnitTests();
		UnitTests(const UnitTests& rhs) = delete;
		UnitTests& operator=(const UnitTests& rhs) = delete;
		~UnitTests();
		void run_unit_tests();
};

#endif