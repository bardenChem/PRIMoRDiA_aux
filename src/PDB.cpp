//PDB.cpp

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
#include <fstream>
#include <iomanip>

#include "../include/global.h"
#include "../include/atom.h"
#include "../include/molecule.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"
#include "../include/pdbModel.h"
#include "../include/PDB.h"

using std::vector;
using std::string;
using std::move;
using std::endl; 

////////////////////////////////////////////////////////////////////////
PDB::PDB()				:
	nModels(0)			,
	PDB_ID("0XXX")		,
	basename("noname")	,
	MULTI(false)		{	
}
/*********************************************************************/
PDB::PDB(const char* pdb_file){
	vector<int> mod_lines;	
	
	basename = remove_extension(pdb_file);
	char pdb_line[10];
	string word = "#";
	word.reserve(10);
	int line = 0;
	
	if ( IF_file( pdb_file ) ){
		std::ifstream buf(pdb_file);
		while( !buf.eof() ){
			buf.getline(pdb_line,10);
			word = pdb_line;
			word = word.substr(0,6);
			if ( word == "MODEL " ) { 
				mod_lines.emplace_back(line);
			}
			line++;
		}
	}	
	
	if ( mod_lines.size() > 0 ) {
		MULTI = true;
		for(unsigned int j=0;j<mod_lines.size();j++){
			models.emplace_back(pdb_file,mod_lines[j]);
		}
	}else{
		models.emplace_back(pdb_file,0);
		nModels++;
	}
		
}
/*********************************************************************/
PDB::~PDB(){}
/*********************************************************************/
PDB::PDB(const PDB& rhs)	:
	nModels(rhs.nModels)	,
	PDB_ID(rhs.PDB_ID)		,
	basename(rhs.basename)	,
	MULTI(rhs.MULTI)		,
	models(rhs.models)		{
}
/*********************************************************************/
PDB& PDB::operator=(const PDB& rhs){
	if ( this != &rhs ){
		nModels	= rhs.nModels;
		PDB_ID	= rhs.PDB_ID;
		basename= rhs.basename;
		MULTI	= rhs.MULTI;
		models	= rhs.models;
	}
	return *this;
}
/*********************************************************************/
PDB::PDB(PDB&& rhs) noexcept		:
	nModels(rhs.nModels)			,
	PDB_ID( move(rhs.PDB_ID) )		,
	basename( move(rhs.basename) )	,
	MULTI(rhs.MULTI)				,
	models( move(rhs.models) )		{	
}
/*********************************************************************/
PDB& PDB::operator=(PDB&& rhs) noexcept{
	if ( this != &rhs ){
		nModels	= rhs.nModels;
		PDB_ID	= move(rhs.PDB_ID);
		basename= move(rhs.basename);
		MULTI	= rhs.MULTI;
		models	= move(rhs.models);
	}
	return *this;	
}
/*********************************************************************/
void PDB::split_models_in_files(){
	string name_ = "";
	for(unsigned int i=0;i<nModels;i++){
		name_ = basename;
		name_ +="_";
		name_ +=std::to_string(i);
		name_ +="_.pdb";
		
		models[i].write_model(name_);
	}
}
/*********************************************************************/
void PDB::add_model(pdbModel model){
	models.emplace_back(model);
	nModels++;
}
/*********************************************************************/
void PDB::cat_pdbs(vector<string> file_list){
	for(unsigned int i=0;i<nModels;i++){
		models.emplace_back( file_list[i].c_str(),0 );
		nModels++;
	}
}
/*********************************************************************/
void PDB::remove_model(unsigned int model){
	models.erase( models.begin()+model );
	nModels--;
}
/*********************************************************************/
void PDB::write_pdb(std::string out_name){
	std::ofstream pdb_file;
	pdb_file.open( out_name.c_str() );
	pdb_file << std::fixed;
	pdb_file.precision(3);
	pdb_file << "PDB file written by LQQCMMtools created by barden.igor@gmail.com" << endl;
	pdb_file << models[0].title << endl;
	pdb_file << models[0].remark << endl;
	
	unsigned int i,j,k,cont = 0;
	for(k=0;k<nModels;k++){
		for(i=0;i<models[k].nResidues;i++){
			for(j=0;j<models[k].monomers[i].nAtoms;j++){
				pdb_file<< std::setw(6) << std::left  << "ATOM" 
						<< " "
						<< std::setw(4) << std::right  << (cont+1) 
						<< " "
						<< std::setw(4) << models[k].monomers[i].r_atoms[j].atom_name
						<< " "
						<< std::left << std::setw(4) << models[k].monomers[i].r_atoms[0].res_name 
						<< " "
						<< std::right << std::setw(4) << (i+1)
						<< std::setw(5) << " "
						<< std::setw(7) << models[k].monomers[i].r_atoms[j].xc
						<< " "
						<< std::setw(7) << models[k].monomers[i].r_atoms[j].yc
						<< " "
						<< std::setw(7) << models[k].monomers[i].r_atoms[j].zc
						<< " "
						<< std::setw(5)  << "1.00"
						<< " "
						<< std::setw(5) << models[k].monomers[i].r_atoms[j].b_factor
						<< "\n";
						cont++;
			}
		}	
		pdb_file << "ENDMDL" << endl;
	}
	pdb_file.close();
}
/*********************************************************************/
void PDB::init_from_system(const molecule& mol){
	vector<pdbAtom> _atoms;
	
	for(int i=0;i<mol.nAtoms;i++){
		pdbAtom _atom;
		_atom.atom_name = mol.atoms[i].element;
		_atom.xc		= mol.atoms[i].xc;
		_atom.yc		= mol.atoms[i].yc;
		_atom.zc		= mol.atoms[i].zc;
		_atom.b_factor	= 0.000;
		_atom.chain_name= "X";
		_atom.indx		= models[0].nAtoms + 1;
		_atom.occupancy = 1.00;
		_atom.res_indx	= models[0].nResidues + 1;
		_atom.res_name	= mol.name.substr(0,3);		
		_atoms.emplace_back(_atom);
	}
	vector<residue> _res;
	_res.emplace_back(_atoms);
	models.emplace_back(_res);
}
/*********************************************************************/
molecule PDB::get_system_from_model(unsigned int model){
	vector<atom> _atoms;
	for(int i=0;i<models[model].nResidues;i++){
		for(int j=0;j<models[model].monomers[i].nAtoms;j++){		
			_atoms.emplace_back(models[model].monomers[i].r_atoms[j].xc			,
								models[model].monomers[i].r_atoms[j].yc			,
								models[model].monomers[i].r_atoms[j].zc			,
								models[model].monomers[i].r_atoms[j].atom_name	);
		}
	}
	molecule mol(_atoms,basename);
	return mol;
}
/*********************************************************************/
vector<molecule> PDB::get_systems(){
	vector<molecule> systems;
	for(unsigned int i=0;i<nModels;i++){
		systems.emplace_back( get_system_from_model(i) );
	}
	return systems;
}
//////////////////////////////////////////////////////
