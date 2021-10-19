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
#include <iostream>

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
PDB::PDB(const char* pdb_file)	:
	nModels(0)					{
	vector<int> mod_lines;
	
	basename = remove_extension(pdb_file);
	char pdb_line[100];
	string word = "#";
	word.reserve(100);
	int line = 0;
	
	if ( IF_file( pdb_file ) ){
		std::ifstream buf(pdb_file);
		while( !buf.eof() ){
			buf.getline(pdb_line,100);
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
			nModels++;
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
	for(unsigned int i=0;i<file_list.size(); i++ ){
		models.emplace_back( file_list[i].c_str(), 0 );
		nModels++;
	}
	basename = remove_extension( file_list[0].c_str() );
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
	
	unsigned int cont = 0;
	for( unsigned k=0; k<models.size(); k++ ){
		pdb_file << "MODEL\n";
		for(unsigned i=0; i<models[k].monomers.size(); i++ ){
			for( unsigned j=0; j<models[k].monomers[i].r_atoms.size(); j++ ){
				pdb_file<< std::setw(6) << std::left  << "ATOM" 
						<< " "
						<< std::setw(4) << std::right  << models[k].monomers[i].r_atoms[j].indx 
						<< " "
						<< std::setw(4) << models[k].monomers[i].r_atoms[j].atom_name
						<< " "
						<< std::left << std::setw(4) << models[k].monomers[i].r_atoms[j].res_name 
						<< " "
						<< std::right << std::setw(4) << models[k].monomers[i].r_atoms[j].res_indx
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
			}
		}
		pdb_file << "ENDMDL" << endl;
	}
	pdb_file.close();
}
/*********************************************************************/
void PDB::init_from_system(const molecule& mol){
	vector<pdbAtom> _atoms;
	
	for(unsigned i=0; i<mol.nAtoms; i++ ){
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
void PDB::iterate_models(std::string func_call, std::vector<std::string>& parameters ){
	if ( func_call == "remove_waters"){
		if ( parameters.size() == 0 ){
			for( unsigned i=0; i<models.size(); i++){
				models[i].remove_waters();
			}
		}else if ( parameters.size() == 2 ){
			double rad = std::stod(parameters[0]);
			double res = std::stoi(parameters[1]);
			for( unsigned i=0; i<models.size(); i++){
				models[i].remove_waters(rad,res);
			}
		}
	}else if ( func_call == "prune_atoms"){
		if ( parameters.size() == 0 ){
			for( unsigned i=0; i<models.size(); i++) models[i].prune_atoms();
		}else if ( parameters.size() == 4 ){
			unsigned catom = std::stoi(parameters[0]);
			unsigned rad = std::stod(parameters[1]);
			bool Within = false;
			bool ByRes = true;
			if ( parameters[2] == "within" ){ Within = true;}
			if ( parameters[3] == "byAtom" ){ ByRes = false; }
			if ( ByRes ){
				for( unsigned i=0; i<models.size(); i++) {
					vector<unsigned> selection = models[i].spherical_selection(catom,rad,Within,ByRes);
					models[i] = models[i].prune_atoms_by_residue( selection );
				}
			}else{
				for( unsigned i=0; i<models.size(); i++) {
					vector<unsigned> selection = models[i].spherical_selection(catom,rad,Within,ByRes);
					models[i] = models[i].prune_atoms( selection );
				}
			}
		}
	}else if ( func_call == "prune_atoms_ByRes"){ 
		if ( parameters.size() == 3 ){
			unsigned res = std::stoi(parameters[0]);
			unsigned rad = std::stod(parameters[1]);
			bool Within = false;
			if ( parameters[2] == "within" ){ Within = true;}
			for( unsigned i=0; i<models.size(); i++) {
				vector<unsigned> selection = models[i].spherical_selection(res,rad,Within,true);
				models[i] = models[i].prune_atoms_by_residue(selection);
			}
		}else if ( parameters.size() == 4 ){
			unsigned catom = std::stoi(parameters[0]);
			unsigned rad = std::stod(parameters[1]);
			bool Within = false;
			bool ByRes = true;
			if ( parameters[2] == "within" ){ Within = true; }
			if ( parameters[3] == "byAtom" ){	ByRes = false; }
			for( unsigned i=0; i<models.size(); i++) {
				vector<unsigned> selection = models[i].spherical_selection(catom,rad,Within,ByRes);
				models[i] = models[i].prune_atoms( selection );
			}
		}
	}
}
/*********************************************************************/
molecule PDB::get_system_from_model(unsigned int model){
	vector<atom> _atoms;
	for(int i=0;i<models[model].nResidues;i++){
		for(int j=0;j<models[model].monomers[i].nAtoms;j++){		
			_atoms.emplace_back(models[model].monomers[i].r_atoms[j].xc			,
								models[model].monomers[i].r_atoms[j].yc			,
								models[model].monomers[i].r_atoms[j].zc			,
								models[model].monomers[i].r_atoms[j].atom_type	);
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
/*********************************************************************/
std::ostream& operator<<(std::ostream& out, const PDB& obj){
	out << obj.basename	<< "\n"
		<< obj.nModels	<< "\n"
		<< obj.PDB_ID	<< "\n";
	if ( obj.models.size() > 0 ){
		out << "printing info of first model in PDB:" << std::endl;
		out << obj.models[0] << std::endl;
	}
	return out;
}
/*********************************************************************/
void PDB::print(){
	std::cout << *this << std::endl;
}
/*********************************************************************/
void UnitTest_PDB(){
	
	const char* pdb_file  = "/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/1l2y.pdb";
	const char* pdb_1dnk  = "/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/1dnk.pdb";
	vector<string> files_list;
	
	files_list.push_back("/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/frame0.pdb");
	files_list.push_back("/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/frame1.pdb");
	files_list.push_back("/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/frame2.pdb");
	files_list.push_back("/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/structure/frame3.pdb");
	
	ut_log.input_line("=======================================");
	ut_log.input_line("Starting unit test for 'PDB class!");
	ut_log.input_line("Default constructor: ");	
	PDB _PDB_A;
	ut_log.data << _PDB_A << endl;
	
	ut_log.input_line("File constructor!");
	PDB _PDB_B(pdb_file);
	ut_log.data << _PDB_B << endl;
	
	_PDB_B.split_models_in_files();
	
	ut_log.input_line("Testing get system methods.");
	vector<molecule> mols = _PDB_B.get_systems();
	
	ut_log.input_line("Testing concatenating pdbs.");
	
	PDB _PDB_C;
	_PDB_C.cat_pdbs(files_list);
	PDB _PDB_D = _PDB_C;
	PDB _PDB_E = _PDB_C;
	PDB _PDB_F = _PDB_C;
	_PDB_C.write_pdb("cat_pdbs.pdb");
	system("pymol cat_pdbs.pdb");
	
	vector<string> pars;
	_PDB_C.iterate_models("remove_waters",pars);
	_PDB_C.write_pdb("frames_without.pdb");
	system("pymol frames_without.pdb");
	
	pars.push_back("10");
	pars.push_back("94");
	
	_PDB_D.iterate_models("remove_waters",pars);
	_PDB_D.write_pdb("frames_without_water_10.pdb");
	system("pymol frames_without_water_10.pdb");
	
	pars[0]	= "3732";
	pars[1] = "10.0";
	pars.push_back("within");
	pars.push_back("ByRes");
	
	_PDB_E.iterate_models("prune_atoms",pars);
	_PDB_E.write_pdb("frames_pruned.pdb");
	system("pymol frames_pruned.pdb");
	
	
	ut_log.input_line("Finished the unit test of the 'pdbModel' class!\n");

}

//////////////////////////////////////////////////////
