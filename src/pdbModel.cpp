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
	index("0")		,
	remark("REMARK"),
	title("TITLE")	,
	model(0)		,
	nChains(0)		,
	nResidues(0)	,
	nAtoms(0)		{
}
/*********************************************************/
pdbModel::pdbModel(std::vector<residue> residues):
	index("0")					,
	remark("REMARK")			,
	title("TITLE")				,
	model(0)					,
	nChains(0)					,
	nResidues(residues.size())	,
	nAtoms(0)					,
	monomers(residues)			{
}
/*********************************************************/
pdbModel::pdbModel(const char* pdb_file, int mdl)	:
	index("0")										,
	remark("REMARK")								,
	title("TITLE")									,
	model(0)										,
	nChains(0)										,
	nResidues(0)									,
	nAtoms(0)										{
	
	if ( !check_file_ext(".pdb",pdb_file) )	{
		std::cout << "Warning! The file has wrong extension name!" << std::endl;
	}
	
	vector<pdbAtom> tmp_atoms;
	
	char pdb_line[100];
	int line = 0;
	unsigned old_res = -1;
	unsigned curr_res = 0;
	string old_chain = "a";
	string curr_chain = "_";
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
					nAtoms++;
					if ( old_res == -1 ) { 
						old_res = _atom.res_indx;
					}
					curr_res = _atom.res_indx;
					if ( curr_res != old_res ){
						residue _residue(tmp_atoms);
						monomers.emplace_back(_residue);
						old_res = curr_res;
						nResidues++;
						tmp_atoms.clear();
						tmp_atoms.emplace_back( move(_atom) );
					}else{
						tmp_atoms.emplace_back( move(_atom) );
					}
					/**********************************/
					if ( old_chain == "a" ){
						old_chain == _atom.chain_name;
					}
					curr_chain = _atom.chain_name;
					if ( curr_chain != old_chain ){
						old_chain = curr_chain;
						nChains++;
					}
					/**********************************/
				}
				else if (	word == "TER   " || word == "ENDMDL" || word == "END   " ){
					if ( nResidues > 0 ){
						residue _residue(tmp_atoms);
						monomers.emplace_back(_residue);
						old_res = curr_res;
						nResidues++;
					}
					break;
				}
			}
			line++;
		}
		residue _residue(tmp_atoms);
		monomers.emplace_back(_residue);
		nResidues++;
	}else{
		std::cout << "PDB file not open!" << std::endl;
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
					<< std::setw(4) << std::right  << monomers[i].r_atoms[j].indx
					<< " "
					<< std::setw(4) << monomers[i].r_atoms[j].atom_name
					<< " "
					<< std::left << std::setw(4) << monomers[i].r_atoms[0].res_name
					<< " "
					<< std::right << std::setw(4) << monomers[i].r_atoms[0].res_indx
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
pdbAtom& pdbModel::pick_atom(unsigned i){
	for(unsigned i=0;i<monomers.size();i++){
		for(unsigned j=0;j<monomers[i].r_atoms.size();j++){
			if ( monomers[i].r_atoms[j].indx == i){
				return monomers[i].r_atoms[j];
			}
		}
	}
	pdbAtom empty;
	return empty;
}
/****************************************************************************/
vector<unsigned> pdbModel::spherical_selection(unsigned center_atom			,
												double radius				,
												bool within					, 
												bool byres 					){
	vector<unsigned> selection;
	vector<unsigned> unselection;
	
	pdbAtom c_atom = this->pick_atom(center_atom);
	
	double dist_ref = 0.0;
	double dist_calc= 0.0;
	unsigned count = 0;
	
	if ( byres ){
		for(unsigned i=0;i<monomers.size();i++){
			for(unsigned j=0;j<monomers[i].r_atoms.size();j++){
				dist_calc = c_atom.get_distance( monomers[i].r_atoms[j]);
				if ( dist_calc > radius ){
					unselection.push_back(i);
				}else{
					selection.push_back(i);
				}
			}
		}
	}else{
		dist_ref = c_atom.get_distance( monomers[0].r_atoms[0]);
		for(unsigned i=0;i<monomers.size();i++){
			for(unsigned j=0;j<monomers[i].r_atoms.size();j++){
				dist_calc = c_atom.get_distance( monomers[i].r_atoms[j]);
				if ( dist_calc > radius ){
					unselection.push_back(count++);
				}else{
					selection.push_back(count++);
				}
			}
		}
	}
	
	
	if ( within ){
		return selection;
	}else{
		return unselection;
	}
}
/*********************************************************/
pdbModel pdbModel::prune_atoms( std::vector<unsigned> selection ){
	pdbModel pruned = *this;
	
	unsigned count = nAtoms;
	for(unsigned i=monomers.size();i>0;i--){
		for(unsigned j=monomers[i].r_atoms.size();j>0;j--){
			for( unsigned k=selection.size();k>0;k--){
				if ( count == selection[count] ){
					pruned.remove_atom(i,j);
				}
			}
			count--;
		}
	}	
	return pruned;
}
/*********************************************************/
pdbModel pdbModel::prune_atoms_by_residue( unsigned res, double radius ){
	pdbModel pruned = *this;
	return pruned;
}
/*********************************************************/
pdbModel pdbModel::prune_atoms_beyond_residue( unsigned res, double radius ){
	pdbModel pruned = *this;
	return pruned;
}
/*********************************************************/
void pdbModel::remove_atom(unsigned int res, unsigned int at){	
	
	monomers[res].r_atoms.erase( monomers[res].r_atoms.begin()+at);
	monomers[res].nAtoms--;
	nAtoms--;
	if ( monomers[res].r_atoms.size() == 0 ){
		this->remove_residue(res);
	}
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
	for(unsigned i=monomers.size()-1;i>0;i--){
		if ( monomers[i].type == WAT ) 
			this->remove_residue(i);
	}
}
/*********************************************************/
void pdbModel::remove_waters(double radius, unsigned int res){
	double refXC, refYC, refZC, distTemp = 0.000;
	refXC = monomers[res].r_atoms[0].xc;
	refYC = monomers[res].r_atoms[0].yc;
	refZC = monomers[res].r_atoms[0].zc;
	
	for(unsigned i=monomers.size()-1;i>0;i--){
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
	for(unsigned i=monomers.size()-1;i>0;i--){
		if ( monomers[i].type == ION ) 
			this->remove_residue(i);
	}
}
/*********************************************************/
void pdbModel::split_complex(std::string mol){
	for(unsigned i=nResidues-1;i>0;i--){
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
/*********************************************************/
double pdbModel::atom_distance(unsigned a1, unsigned a2){
	int count_a1	= 0;
	int count_a2	= 0;
	int resnmb_a1	= 0;
	int resnmb_a2	= 0;
	pdbAtom atom1;
	pdbAtom atom2;
	
	for( int i=0;i<monomers.size();i++){
		if ( count_a1 < a1 ){
			count_a1 += monomers[i].nAtoms;
		}else{
			resnmb_a1 = i;
		}
		if ( count_a2 < a2 ){
			count_a2 += monomers[i].nAtoms;
		}else{
			resnmb_a2 = i;
		}
	}
	
	for(int j=0;j<monomers[resnmb_a1].nAtoms;j++){
		if ( monomers[resnmb_a1].r_atoms[j].indx == a1 ){
			atom1 = monomers[resnmb_a1].r_atoms[j];
		}
	}
	
	for(int j=0;j<monomers[resnmb_a2].nAtoms;j++){
		if ( monomers[resnmb_a2].r_atoms[j].indx == a2 ){
			atom2 = monomers[resnmb_a2].r_atoms[j];
		}
	}
	return atom1.get_distance(atom2);
}
/*********************************************************/
double pdbModel::atom_angle(unsigned a1, unsigned a2, unsigned a3){
	
}
/*********************************************************/
std::ostream& operator<<(std::ostream& out, const pdbModel& obj){
	out << "Outputting information for model in a PDB file\n"
		<< obj.title 
		<< "\n"
		<< obj.remark
		<< "\n\tNumber of atoms in the model "	<< obj.nAtoms
		<< "\n\tNumber of residues "			<< obj.nResidues
		<< "\n\tNumber of chains "				<< obj.nChains
		<< "\n\tmodel index	"					<< obj.index;
	return out;
}
/*********************************************************/
void pdbModel::print(){
	std::cout << *this << endl;
}
/*********************************************************/
void UnitTest_pdbModel(){
	ut_log.input_line("=======================================");
	ut_log.input_line("Starting unit test for 'pdbModel class!");
	ut_log.input_line("Default constructor: ");
	pdbModel _model_A;
	ut_log.data << _model_A << endl;
	
	const char* pdb_file  = "/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/tim.pdb";
	pdbModel _model_B(pdb_file,0);
	ut_log.input_line("File opening constructor: ");
	ut_log.data << _model_B << endl;
	
	ut_log.input_line("Writting method testing. \nWill open Pymol! ");
	_model_B.write_model("test.pdb");
	
	system("pymol test.pdb");
	
	ut_log.input_line("Testing copy constructor!");
	pdbModel _model_C(_model_B);
	ut_log.data << _model_C << endl;
	
	ut_log.input_line("Testing assing operator overloading!");
	pdbModel _model_D = _model_C; 
	ut_log.data << _model_C << endl;
	
	ut_log.input_line("Testing move constructor!");
	pdbModel _model_E( std::move(_model_D) );
	ut_log.data << _model_E << endl;
	
	ut_log.input_line("Testing move constructor!");
	pdbModel _model_F = std::move(_model_E);
	ut_log.data << _model_F << endl;
	
	ut_log.input_line("Testing the water removal!");
	_model_F.remove_waters();
	_model_F.write_model("test_wwater.pdb");
	system("pymol test_wwater.pdb");
	
	ut_log.input_line("Testing the water removal from a given radius of a residue!");
	_model_C.remove_waters(7.0,93);
	_model_C.write_model("test_wwater_r.pdb");
	system("pymol test_wwater_r.pdb");
	
	ut_log.input_line("Finished the unit test of the 'pdbModel' class!\n");

}	

//////////////////////////////////////////////////////////