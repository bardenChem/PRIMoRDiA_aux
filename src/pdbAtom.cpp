//pdbAtom.cpp

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

#include "../include/global.h"
#include "../include/pdbAtom.h"

#include <string>
#include <cmath>
#include <iostream>
#include<algorithm>
using std::endl;
using std::string;

/******************************************************************/
pdbAtom::pdbAtom()		:
	atom_name("nonamed"),
	atom_type("H")		,
	indx(1)				,
	res_name("UNK")		,
	res_indx(0)			,
	chain_name(" ")		,
	occupancy(0.0)		,
	b_factor(0.00)		,
	sideC(false)		,
	xc(0.00)			,
	yc(0.00)			,
	zc(0.00)			{
}
/*********************************************************/
pdbAtom::pdbAtom(std::string _res_name, unsigned _res_num, std::string element, double _xc, double _yc, double _zc):
	atom_name(element)	,
	atom_type(element)  ,
	indx(1)				,
	res_name(_res_name)	,
	res_indx(_res_num)	,
	chain_name(" ")		,
	occupancy(0.0)		,
	b_factor(0.00)		,
	sideC(false)		,
	xc(_xc)			    ,
	yc(_yc)			    ,
	zc(_zc)			    {

	 if ( this->is_hydrogen() ){
		 atom_type = "H";
	 }else{
		 atom_type = element[0];
	 }
}
/*********************************************************/
pdbAtom::~pdbAtom(){}

/*********************************************************/
pdbAtom::pdbAtom(const pdbAtom& rhs):
	atom_name(rhs.atom_name)		,
	atom_type(rhs.atom_type)		,
	indx(rhs.indx)					,
	res_name(rhs.res_name)			,
	res_indx(rhs.res_indx)			,
	chain_name(rhs.chain_name)		,
	occupancy(rhs.occupancy)		,
	b_factor(rhs.b_factor)			,
	sideC(rhs.sideC)				,
	xc(rhs.xc)						,
	yc(rhs.yc)						,
	zc(rhs.zc)						{
}
/*********************************************************/
pdbAtom& pdbAtom::operator=(const pdbAtom& rhs){
	if ( this != &rhs ){
		atom_name	= rhs.atom_name;
		atom_type	= rhs.atom_type;
		indx		= rhs.indx;
		res_name	= rhs.res_name;
		res_indx	= rhs.res_indx;
		chain_name	= rhs.chain_name;
		occupancy	= rhs.occupancy;
		b_factor	= rhs.b_factor;
		sideC		= rhs.sideC;
		xc			= rhs.xc;
		yc			= rhs.yc;
		zc			= rhs.zc;
	}
	return *this;
}
/*********************************************************/
pdbAtom::pdbAtom(pdbAtom&& rhs) noexcept:
	atom_name( move(rhs.atom_name) )	,
	atom_type( move(rhs.atom_type) )	,
	indx(rhs.indx)						,
	res_name( move(rhs.res_name) )		,
	res_indx( rhs.res_indx )			,
	chain_name( move(rhs.chain_name) )  ,
	occupancy(rhs.occupancy)			,
	b_factor(rhs.b_factor)  			,
	sideC(rhs.sideC) 					,
	xc(rhs.xc)							,
	yc(rhs.yc)							,
	zc(rhs.zc)							{
	
}
/*********************************************************/
pdbAtom& pdbAtom::operator=(pdbAtom&& rhs) noexcept{
	if ( this != &rhs ){
		atom_name	= move(rhs.atom_name);
		atom_type	= move(rhs.atom_type);
		indx		= rhs.indx;
		res_name	= move(rhs.res_name);
		res_indx	= rhs.res_indx; 
		chain_name	= move(rhs.chain_name);
		occupancy	= rhs.occupancy;
		b_factor	= rhs.b_factor;
		sideC		= rhs.sideC;
		xc			= rhs.xc;
		yc			= rhs.yc;
		zc			= rhs.zc;
	}
	return *this;
}
/*********************************************************/
pdbAtom::pdbAtom(std::string& pdb_line)	:
	atom_name(pdb_line,12,4)			,
	res_name(pdb_line,17,3)				,
	chain_name(pdb_line,21,1)			{
	
	string atom_symbol = atom_name;
	std::remove(atom_symbol.begin(), atom_symbol.end(), ' ');
	atom_symbol = atom_symbol.substr(0,2);
	if ( atom_symbol == "Mg") {
		atom_type = "Mg";
	}else if( atom_symbol == "Cl"){
		atom_type = "Cl";
	}else if( atom_symbol == "Zn"){
		atom_type = "Zn";
	}else if( atom_symbol == "Ca"){
		atom_type = "Ca";
	}else if( atom_symbol == "Na"){
		atom_type = "Na";
	}else{
		atom_type = atom_symbol.substr(0,1);
	}
	
	string tmp_rsin(pdb_line,6,5);
	indx	 = stoi(tmp_rsin);
	string tmp_rsi(pdb_line,22,4);
	res_indx = stoi(tmp_rsi);
	string temp_xc(pdb_line,31,7);
	string temp_yc(pdb_line,39,7);
	string temp_zc(pdb_line,47,7);
	string occ(pdb_line,56,3);
	string bfactor(pdb_line,61,5);
	xc = stod(temp_xc);
	yc = stod(temp_yc);
	zc = stod(temp_zc);
	occupancy = stod(occ);
	b_factor = stod(bfactor);
}
/*********************************************************/
bool pdbAtom::is_hydrogen(){
	
	string tmp1(atom_name,0,2);
	
	if 		( tmp1 == " H" ) return true;
	else if ( tmp1 == "1H" ) return true;
	else if ( tmp1 == "2H" ) return true;
	else if ( tmp1 == "3H" ) return true;
	else return false;	
}
/*********************************************************/
double pdbAtom::get_distance(const pdbAtom& a2){
	double dist = 0.0;
	dist = (xc - a2.xc)*(xc - a2.xc);
	dist+= (yc - a2.yc)*(yc - a2.yc);
	dist+= (zc - a2.zc)*(zc - a2.zc);
	return sqrt(dist);
}
/*********************************************************/
std::ostream& operator<<(std::ostream& out, const pdbAtom& obj){
	out << "Outputting information about object of a class 'pdbAtom'"
		<< "\nAtom name: "				<< obj.atom_name
		<< "\n\t index: "				<< obj.indx
		<< "\n\t residue "				<< obj.res_name
		<< "\n\t residue index "		<< obj.res_indx
		<< "\n\t chain name "			<< obj.chain_name
		<< "\n\t b-factor "				<< obj.b_factor
		<< "\n\t occupancy "			<< obj.occupancy
		<< "\n\t x-coordinate "			<< obj.xc
		<< "\n\t y-coordinate "			<< obj.yc
		<< "\n\t z-coordinate "			<< obj.zc;
		
		return out;
}
/*********************************************************/
void pdbAtom::print(){
	std::cout << *this << std::endl;
}
/*********************************************************/
void UnitTest_pdbAtom(){
	string pdb_line_A = "ATOM      1  N   ALA  1         53.325  28.151  42.519  0.00  0.00          N  ";
	string pdb_line_B = "HETATM 3732  C02 LIG  248       41.143  37.015  22.775  0.00  0.00          C  ";
	
	ut_log.input_line("========================================");
	ut_log.input_line("Starting unit test for 'pdbAtom' class");
	ut_log.input_line("Default constructor:");
	pdbAtom _atom_a;
	ut_log.data << _atom_a << endl;
	
	ut_log.input_line("String line constructor:");
	pdbAtom _atom_b(pdb_line_A);
	ut_log.data << _atom_b << endl;
	
	ut_log.input_line("copy constructor:");
	pdbAtom _atom_c(_atom_b);
	ut_log.data << _atom_c << endl;
	
	ut_log.input_line("Assign operator overloading:");
	pdbAtom _atom_d = _atom_c;
	ut_log.data << _atom_d << endl;
	
	ut_log.input_line("move constructor:");
	pdbAtom _atom_e( std::move(_atom_b) );
	ut_log.input_line("\tmoved atom:");
	ut_log.data << _atom_b << endl;
	ut_log.input_line("\tnew atom:");
	ut_log.data << _atom_e << endl;
	
	ut_log.input_line("move Assign operator overloading:");
	pdbAtom _atom_f = std::move(_atom_e);
	ut_log.data << _atom_f << endl;
	
	pdbAtom _atom_g(pdb_line_B);
	double distance = _atom_f.get_distance(_atom_g);
	ut_log.input_line("calculated distance between two atoms: ");
	ut_log.data << distance << endl;
	
	ut_log.input_line("Finished the unit test of the 'pdbAtom' class!\n");

}

////////////////////////////////////////////////////////////
