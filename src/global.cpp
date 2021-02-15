//global.h
//Declarations for functions and source file with constant values. 

/*********************************************************************/
/* This source code file is part of OOCCuPy++ software project created 
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

//Including header from the c++
#include <iostream>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <experimental/filesystem>

#include "../include/global.h"
//=============================================================================
using std::string;
using std::vector;
namespace fs = std::experimental::filesystem;

/*********************************************************************************/
bool IF_file(const char* name){
	fs::path file_name(name);
	return fs::exists(file_name);
}
/*********************************************************************************/
bool IF_file(fs::path& name){
	return fs::exists(name);
}
/*********************************************************************************/
bool check_file_ext(string ext,const char* file_name){
	fs::path f_name(file_name);
	if ( f_name.extension() == ext ) return true;
	else return false;
}
/*********************************************************************************/
string get_file_name(const char* path_name){
	fs::path f_name(path_name);
	string resultS( f_name.filename() );
	return resultS;
}
/*********************************************************************************/
string remove_extension(const char* file_name){
	string file_name_string  = file_name;
	int point = 0;
	char dot = '.';
	for(unsigned int i = 0;i<file_name_string.size();i++){
		char character = file_name_string[i];
		if ( character == dot){	point = i; }
	}
	int pos = file_name_string.size()-point;
	return file_name_string.substr(0,file_name_string.size()-pos);
}
/*********************************************************************************/
string change_extension(const char* file_name,string new_ext){
	string name_wth_ext = remove_extension(file_name);
	name_wth_ext += new_ext;
	return name_wth_ext;
}
/*********************************************************************************/
void rename_file(const char* file_name,string new_file_name){
	fs::path f_name(file_name);
	fs::path nf_name( fs::current_path() / new_file_name );
	fs::rename(f_name,nf_name);
}
/*********************************************************************************/
string str_array(string& line, int in, int fin){
	string result = "";
	result.resize(fin-in+2);
	int cnt = 0;
	for(int i=in;i<fin;i++){ 
		result[cnt++] = line[i];
	}
	return result;
}
/*********************************************************************************/
double mean_dvec(vector<double>& vec){
	double result = 0.0;
	for(unsigned int i=0;i<vec.size();i++){
		result += vec[i];
	}	
	return result/vec.size();
}
/********************************************************************************/
double max_dvec(vector<double>& vec){
	std::sort( vec.begin(), vec.end() );
	if ( vec.size() > 0) return vec[vec.size()-1];
	else return 0.0;
}
/********************************************************************************/
double min_dvec(vector<double>& vec){
	std::sort( vec.begin(), vec.end() );
	if ( vec.size() > 0) return vec[0];
	else return 0.0;
}
/********************************************************************************/
double sum_dvec(std::vector<double>& vec){
	double result = 0.0;
	for(unsigned int i=0;i<vec.size();i++){
		result += vec[i];
	}	
	return result;
}
/******************************************************************************/
double sd_dvec(vector<double>& vec){
	double mean = mean_dvec(vec);
	double sd = 0;
	for(unsigned int i=0;i<vec.size();i++){
		sd += (vec[i]-mean)*(vec[i]-mean);
	}
	sd/=vec.size();
	sd = sqrt(sd);
}
/******************************************************************************/
vector<double> scale_dvec(vector<double>& vec){
	vector<double> result( vec.size() );
	double mean = mean_dvec(vec);
	for(unsigned int i=0;i<vec.size();i++){
		result[i] = vec[i] - mean;
	}
	double sd = sd_dvec(result);
	for(unsigned int i=0;i<result.size();i++){
		result[i]/=sd;
	}
	return result;
}

//==================================================================================