//XYZ.cpp

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

#include <vector>
#include <string>

#include "../include/atom.h"
#include "../include/molecule.h"
#include "../include/XYZ.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;

/*********************************************************************/
XYZ::XYZ()									:
	nAtoms(0)								,
	name("noname")							,
	commentary("Created by LQQCMM tools!")	{		
}
/*********************************************************************/
XYZ::XYZ(const char* xyz_file)			:
	name( remove_extension(xyz_file) )	{
	char tmp_line[30];
	int line = 0;
	double tmp_crdx	= 0.0;
	double tmp_crdy	= 0.0;
	double tmp_crdz	= 0.0;
	string tmp_str;

	if ( IF_file( xyz_file ) ){
		std::ifstream buf(_file);
		while( !buf.eof() ){
			buf.getline(tmp_line,30);
			std::sstream _stream(tmp_line);
			if ( line == 0 ){
				_stream >> nAtoms;
			}else if ( line == 1 ){
				_stream >> commentary;
			}else{
				_stream >> tmp_str;
				_stream >> tmp_crdx;
				_stream >> tmp_crdy;
				_stream >> tmp_crdz;
				atoms.emplace_back(tmp_crdx,tmp_crdy,tmp_crdz,tmp_str);
			}			
		}
	}
}
/*********************************************************************/
XYZ::~XYZ(){}
/*********************************************************************/
XYZ::XYZ(const molecule& mol)					:
	name(mol.name)							,
	nAtoms(mol.nAtoms)						,
	atoms(mol.atoms)						,
	commentary("Created by LQQCMM tools!")	{	
}
/*********************************************************************/
XYZ::XYZ(const XYZ& rhs)		:
	nAtoms(rhs.nAtoms)			,
	atoms(rhs.atoms)			,
	commentary(rhs.commentary)	{	
}
/*********************************************************************/
XYZ& XYZ::operator=(const XYZ& rhs){
	if ( this != &rhs ){
		nAtoms		= rhs.nAtoms;
		atoms 		= rhs.atoms;
		commentary	=rhs.commentary;		
	}
	return *this;	
}
/*********************************************************************/
XYZ::XYZ(XYZ&& rhs) noexcept			:
	nAtoms(rhs.nAtoms)					,
	atoms( move(rhs.atoms) )			,
	commentary( move(rhs.commentary) )	{
}
/*********************************************************************/
XYZ& XYZ::operator=(XYZ&& rhs) noexcept{
	if ( this != &rhs ){
		nAtoms		= rhs.nAtoms;
		atoms 		= rhs.atoms;
		commentary	= move(rhs.commentary);		
	}
	return *this;		
}
/*********************************************************************/
void XYZ::write_xyz(std::string out_name){
	std::ofstream out_file;
	out_file.open( out_name.c_str() );
	out_file << std::fixed; 
	out_file << std::precision(4);
	
	out_file << nAtoms << endl;
	out_file << commentary << endl;
	for(int i=0;i<nAtoms;i++){
		out_file << atoms[i].element	<< " "
				 << atoms[i].xc			<< " "
				 << atoms[i].yc 		<< " "
				 << atoms[i].zc 		<< "\n"
	}
}
/*********************************************************************/
molecule XYZ::get_molecule(){
	molecule mol(atoms);
	return mol;
}
//////////////////////////////////////////////////////////////////////