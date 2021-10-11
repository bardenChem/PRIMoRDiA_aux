//Line.cpp
//Definitions for functions and source file with constant values. 

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
#include <string>
#include <vector> 
#include <sstream>
#include <iostream>

#include "../include/Line.h"

using std::string;
using std::vector;

////////////////////////////////////////////////////
Line::Line()	:
	line_len(0)	{
}
/*****************************************/
Line::~Line(){}
/*****************************************/
Line::Line(string line):
	line_len(0)			  {
		
	words.reserve(30);
	std::stringstream stream(line);
	string word;
	word.reserve(20);
	
	while (stream >> word){
		words.emplace_back( word );
		line_len++;
	}
}
/*****************************************/
Line::Line(const char* line):
	line_len(0)				{
		
	words.reserve(30);
	std::stringstream stream(line);
	string word;
	word.reserve(20);
	
	while (stream >> word){
		words.emplace_back( word );
		line_len++;
	}
}
/*****************************************/
Line::Line(const Line& rhs)	:
	words(rhs.words)		,
	line_len(rhs.line_len)	{
}
/*****************************************/
Line& Line::operator=(const Line& rhs){
	if ( this != &rhs ){
		words 		= rhs.words;
		line_len 	= rhs.line_len;
	}
	return *this;
}
/*****************************************/
Line::Line(Line&& rhs) noexcept	:
	words( move(rhs.words) )	,
	line_len(rhs.line_len)		{	
}
/*****************************************/
Line& Line::operator=(Line&& rhs) noexcept{
	if ( this != &rhs ){
		words = move(rhs.words);
		line_len = rhs.line_len;
	}
	return *this;
}
/*****************************************/
double Line::get_double(int pos){
	double value = 0.0000;
	if ( pos > words.size() ){
		std::cout << "Impossible to convert position beyond the size of line!" << std::endl;
		return value;
	}
	std::stringstream stream( words[pos] );
	try{
		stream >> value;
	}catch( const std::invalid_argument& ){
		std::cout << "Impossible to convert to double!" << std::endl;
	}
	return value;
}
/*****************************************/
int Line::get_int(int pos){
	int value = 0;
	if ( pos > words.size() ){
		std::cout << "Impossible to convert position beyond the size of line!" << std::endl;
		return value;
	}
	std::stringstream stream( words[pos] );
	try{
		stream >> value;
	}catch( const std::invalid_argument& ){
		std::cout << "Impossible to convert to double!" << std::endl;
	}
	return value;
}
///////////////////////////////////////////