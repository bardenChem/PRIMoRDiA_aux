//pdbModel.cpp

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
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>


#include "../include/global.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"
#include "../include/pdbModel.h"

using std::string;
using std::vector;
using std::move;
using std::endl;

//=========================================================

pdbModel::pdbModel():
	index(0)		,
	remark("")		,
	title("")		,
	model(0)		,
	nChains(0)		,
	nResidues(0)	,
	nAtoms(0)		{	
}
/*********************************************************/
pdbModel::pdbModel(std::vector<residue> residues):
	index(0)					,
	remark("")					,
	title("")					,
	model(0)					,
	nChains(0)					,
	nResidues(residues.size())	,
	nAtoms(0)					,
	monomers(residues)			{
}
/*********************************************************/
pdbModel::pdbModel(const char* pdb_file, int mdl){
	if ( !check_file_ext(".pdb",pdb_file) )	{
		std::cout << "Warning! The file has wrong extension name!" << std::endl;
	}
	
	vector<pdbAtom> tmp_atoms;
	
	char pdb_line[100];
	int line = 0;
	string old_res = "0";
	string curr_res = "_";	
	string pdb_lineS;
	
	if ( IF_file( pdb_file ) ){
		std::ifstream buf(pdb_file);
		while( !buf.eof() ){
			if ( line >= mdl ) {
				buf.getline(pdb_line,100);
				string word(pdb_line,0,6);
				if ( word == "ATOM  " || word == "HETATM" ) {
					pdb_lineS = pdb_line;
					pdbAtom _atom(pdb_lineS);					
					if ( old_res == "0" ) { old_res = _atom.res_name; }
					curr_res = _atom.res_name;
					tmp_atoms.emplace_back( move(_atom) );
					if ( curr_res != old_res ){
						residue _residue(tmp_atoms);
						monomers.emplace_back(_residue);
						old_res = curr_res;
					}
				}
				else if ( 	word == "TER   " ||
							word == "ENDMDL" || 
							word == "END   " ){
					break;
				}
			}
			line++;
		}		
	}
}
/*********************************************************/
pdbModel::~pdbModel(){}
/*********************************************************/
pdbModel::pdbModel(const pdbModel& rhs):
	index(rhs.index)					,
	remark(rhs.remark)					,
	title(rhs.title)					,
	model(rhs.model)					,
	nChains(rhs.nChains)				,
	nResidues(rhs.nResidues)			,
	nAtoms(rhs.nAtoms)					,
	monomers(rhs.monomers)				{
}
/*********************************************************/
pdbModel& pdbModel::operator=(const pdbModel& rhs){
	if ( this != &rhs ){
		index		= rhs.index;					
		remark		= rhs.remark;					
		title		= rhs.title;					
		model		= rhs.model;					
		nChains		= rhs.nChains;				
		nResidues	= rhs.nResidues;			
		nAtoms		= rhs.nAtoms;					
		monomers	= rhs.monomers;			
	}
	return *this;
}
/*********************************************************/
pdbModel::pdbModel(pdbModel&& rhs) noexcept:
	index(rhs.index)					,
	remark( move(rhs.remark) )			,
	title( move(rhs.title) )			,
	model(rhs.model)					,
	nChains(rhs.nChains)				,
	nResidues(rhs.nResidues)			,
	nAtoms(rhs.nAtoms)					,
	monomers( move(rhs.monomers) )		{
}
/*********************************************************/
pdbModel& pdbModel::operator=(pdbModel&& rhs) noexcept{
	if ( this != &rhs ){
		index		= rhs.index;					
		remark		= move(rhs.remark);					
		title		= move(rhs.title);					
		model		= rhs.model;					
		nChains		= rhs.nChains;				
		nResidues	= rhs.nResidues;			
		nAtoms		= rhs.nAtoms;					
		monomers	= move(rhs.monomers);			
	}
	return *this;
	
}
/*********************************************************/
void pdbModel::write_model(std::string out_name){
	std::ofstream pdb_file;
	pdb_file.open( out_name.c_str() );
	pdb_file << std::fixed;
	pdb_file.precision(3);
	pdb_file << "PDB file written by LQQCMMtools created by barden.igor@gmail.com" << endl;
	pdb_file << title << endl;
	pdb_file << remark << endl;
	
	unsigned int i,j,cont;
	for(i=0;i<nResidues;i++){
		for(j=0;j<monomers[i].nAtoms;j++){
			pdb_file<< std::setw(6) << std::left  << "ATOM" 
					<< " "
					<< std::setw(4) << std::right  << (cont+1) 
					<< " "
					<< std::setw(4) << monomers[i].r_atoms[j].atom_name
					<< " "
					<< std::left << std::setw(4) << monomers[i].r_atoms[0].res_name
					<< " "
					<< std::right << std::setw(4) << (i+1)
					<< std::setw(5) << " "
					<< std::setw(7) << monomers[i].r_atoms[j].xc
					<< " "
					<< std::setw(7) << monomers[i].r_atoms[j].yc
					<< " "
					<< std::setw(7) << monomers[i].r_atoms[j].zc
					<< " "
					<< std::setw(5)  << "1.00"
					<< " "
					<< std::setw(5) << monomers[i].r_atoms[j].b_factor
					<< "\n";
					cont++;
		}
	}
	pdb_file << "ENDMDL" << endl;
	pdb_file.close();
}
/*********************************************************/
void pdbModel::remove_atom(unsigned int res, unsigned int at){	
	
	monomers[res].r_atoms.erase( monomers[res].r_atoms.begin()+at);
	monomers[res].nAtoms--;
	nAtoms--;
	
}
/*********************************************************/
void pdbModel::remove_residue(unsigned int i){
	nAtoms -= monomers[i].nAtoms;
	monomers.erase( monomers.begin() +i);
	nResidues--;
}
/*********************************************************/
void pdbModel::prune_atoms(){
	unsigned int i,j;	
	for(i=0;i<nResidues;i++){
		for(j=0;j<monomers[i].nAtoms;j++){
			if ( monomers[i].r_atoms[j].res_name.substr(0,1) == "B" ) 
				this->remove_atom( i, j ); 	
		}
	}
}
/*********************************************************/
void pdbModel::remove_waters(){
	unsigned int i;
	for(i=0;i<nResidues;i++){
		if ( monomers[i].type == WAT ) this->remove_residue(i);
	}
}
/*********************************************************/
void pdbModel::remove_waters(double radius, unsigned int res){
	unsigned int i;
	double refXC, refYC, refZC, distTemp = 0.000;
	refXC = monomers[res].r_atoms[0].xc;
	refYC = monomers[res].r_atoms[0].yc;
	refZC = monomers[res].r_atoms[0].zc;
	
	for(i=0;i<nResidues;i++){
		if ( monomers[i].type == WAT ) {
			distTemp =  (monomers[i].r_atoms[0].xc - refXC)*(monomers[i].r_atoms[0].xc - refXC);
			distTemp += (monomers[i].r_atoms[0].yc - refYC)*(monomers[i].r_atoms[0].yc - refYC);
			distTemp += (monomers[i].r_atoms[0].zc - refZC)*(monomers[i].r_atoms[0].zc - refZC);
			distTemp = sqrt(distTemp);
			if ( distTemp > radius )
				this->remove_residue(i);
		}
	}
}
/*********************************************************/
void pdbModel::remove_ions(){
	unsigned int i = 0;
	for(i;i<nResidues;i++){
		if ( monomers[i].type == ION ) this->remove_residue(i);
	}
}
/*********************************************************/
void pdbModel::split_complex(std::string mol){
	unsigned int i;
	for(i=nResidues-1 ;i>0;i--){
		if ( mol == monomers[i].r_atoms[0].res_name ) {
			pdbModel ligand;
			ligand.monomers.emplace_back(monomers[i]);
			ligand.title  = "Ligand " + mol;
			ligand.remark = "Ligand splitted from complex files.";
			ligand.nAtoms = monomers[i].nAtoms;
			ligand.write_model( mol +".pdb" );
		}
	}
}
/*********************************************************/
void pdbModel::built_complex(const char* pdb_mol){
	pdbModel temp(pdb_mol,0);
	monomers.emplace_back( temp.monomers[0] );
	nResidues++;
	nAtoms += temp.monomers[0].nAtoms;
}
//////////////////////////////////////////////////////////