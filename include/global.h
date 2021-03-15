//global.h
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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <vector> 
#include <experimental/filesystem>


//========================================================
bool IF_file(const char* name); 
bool IF_file(std::experimental::filesystem::path& name); 
bool check_file_ext(std::string ext,const char* file_name);
std::string remove_extension(const char* file_name);
std::string change_extension(const char* file_name,std::string new_ext);
void rename_file(const char* file_name,std::string new_file_name);
std::string get_file_name(const char* path_name);
std::string str_array(std::string& line, int in, int fin);
double mean_dvec(std::vector<double>& vec);
double max_dvec(std::vector<double>& vec);
double min_dvec(std::vector<double>& vec);
double sum_dvec(std::vector<double>& vec);
double sd_dvec(std::vector<double>& vec);
float get_atom_mass(std::string sym);
int get_atomic_number(std::string sym);
std::vector<double> scale_dvec(std::vector<double>& vec);
std::string get_res1n( int i );
std::string get_res3n( int i );
int get_AAnHy(int i);
//------------------------------------------------------------
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
