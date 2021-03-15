//Line.h
//Declarations for functions and source file with constant values. 

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

#ifndef LINE_H
#define LINE_H

#include <string>
#include <vector> 

class Line{
	public:
		std::vector<std::string> words;
		unsigned int line_len;
		Line();
		~Line();
		Line(std::string line);
		Line(const char* line);
		Line(const Line& rhs);
		Line& operator=(const Line& rhs);
		Line(Line&& rhs) noexcept;
		Line& operator=(Line&& rhs) noexcept;
		double get_double(int pos);
		int get_int(int pos);
};

#endif