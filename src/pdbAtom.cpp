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

#include "../include/pdbAtom.h"

#include <string>
#include <cmath>
using std::string;

/******************************************************************/
pdbAtom::pdbAtom()		:
	atom_name("nonamed"),
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
pdbAtom::~pdbAtom(){}
/*********************************************************/
pdbAtom::pdbAtom(const pdbAtom& rhs):
	atom_name(rhs.atom_name)		,
	indx(rhs.indx)					,
	res_name(rhs.res_name)			,
	res_indx(0)						,
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
	chain_name(pdb_line,21,2)			{
	
	string tmp_rsi(pdb_line,23,3);
	res_indx = stoi(tmp_rsi);
	string temp_xc(pdb_line,31,6);
	string temp_yc(pdb_line,39,6);
	string temp_zc(pdb_line,47,6);
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
double pdbAtom::get_distance(const pdbAtom& a2){
	double dist = 0.0;
	dist = (xc - a2.xc)*(xc - a2.xc);
	dist+= (yc - a2.yc)*(yc - a2.yc);
	dist+= (zc - a2.zc)*(zc - a2.zc);
	return sqrt(dist);
}////////////////////////////////////////////////////////////
