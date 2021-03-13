//geometry.cpp

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
#include <sstream>

#include "../include/global.h"
#include "../include/system.h"
#include "../include/geometry.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;

//////////////////////////////////////////////////////////////////////
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
	double tmp_coord;
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
				elements.emplace_back(tmp_str);
				_stream >> tmp_coord;
				xc.emplace_back(tmp_coord);
				_stream >> tmp_coord;
				yc.emplace_back(tmp_coord);
				_stream >> tmp_coord;
				zc.emplace_back(tmp_coord);
			}			
		}
	}
}
/*********************************************************************/
XYZ::~XYZ(){}
/*********************************************************************/
XYZ::XYZ(const system& mol)					:
	name(mol.name)							,
	nAtoms(mol.nAtoms)						,
	commentary("Created by LQQCMM tools!")	{
	
	for(int i=0;i<nAtoms;i++){
		elements.emplace_bacl( mol.atoms[i].element );
		xc.emplace_back( mol.atoms[i].xc );
		yc.emplace_back( mol.atoms[i].yc );
		zc.emplace_back( mol.atoms[i].zc );

	}
}
/*********************************************************************/
XYZ::XYZ(const XYZ& rhs)		:
	nAtoms(rhs.nAtoms)			,
	xc(rhs.xc)					,
	yc(rhs.yc)					,
	zc(rhs.zc)					,
	elements(rhs.elements)		,
	commentary(rhs.commentary)	{	
}
/*********************************************************************/
XYZ& XYZ::operator=(const XYZ& rhs){
	if ( this != &rhs ){
		nAtoms		= rhs.nAtoms;
		xc 			= rhs.xc;
		yc 			= rhs.yc;
		zc 			= rhs.zc;
		elements	= rhs.elements;
		commentary	=rhs.commentary;		
	}
	return *this;	
}
/*********************************************************************/
XYZ::XYZ(XYZ&& rhs) noexcept			:
	nAtoms(rhs.nAtoms)					,
	xc( move(rhs.xc) )					,
	yc( move(rhs.yc) )					,
	zc( move(rhs.zc) )					,
	elements( move(rhs.elements) )		,
	commentary( move(rhs.commentary) )	{
}
/*********************************************************************/
XYZ& XYZ::operator=(XYZ&& rhs) noexcept{
	if ( this != &rhs ){
		nAtoms		= rhs.nAtoms;
		xc 			= move(rhs.xc);
		yc 			= move(rhs.yc);
		zc 			= move(rhs.zc);
		elements	= move(rhs.elements);
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
		out_file << elements[i] << " "
				 << xc[i] 		<< " "
				 << yc[i] 		<< " "
				 << zc[i] 		<< "\n"
	}
}
/*********************************************************************/
system XYZ::get_molecule(){
	std::vector<atom> tmp_atoms;
	for(unsigned int i=0;i<nAtoms;i++){
		tmp_atoms.emplace_back(xc[i],yc[i],zc[i],elements[i]);
	}
	system mol(tmp_atoms);
	return mol;
}
////////////////////////////////////////////////////////////////////////
PDB::PDB()			:
	nModels(0)		,
	PDB_ID("0XXX"	,
	basename(" ")	,
	MULTI(false)	{	
}
/*********************************************************************/
PDB::PDB(const char* pdb_file){
	vector<int> mod_lines;	
	
	basename = remove_extension(pdb_file);
	char pdb_line[10];
	string word = "#";
	word.reserve(10);
	int line = 0;
	
	if ( IF_file( pdb_name ) ){
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
		for(int j=0;j<mod_lines.size();j++){
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
	nModels(rhs.models)		,
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
	nModels(rhs.models)				,
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
	for( int i=0;i<nModels;i++){
		name_ = basename + "_" + string(i) +"_.pdb";
		models[i].write_model(name_);
	}
}
/*********************************************************************/
void PDB::add_model(pdbModel model){
	models.eḿplace_back(model);
	nModels++;
}
/*********************************************************************/
void PDB::cat_pdbs(vector<string> file_list){
	for(int i=0;i<){
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
	for(k,k<nModels;k++){
		for(i;i<models[k].nResidues;i++){
			for(j;j<models[k].monomers[i].nAtoms;j++){
				pdb_file<< std::setw(6) << std::left  << "ATOM" 
						<< " "
						<< std::setw(4) << std::right  << (cont+1) 
						<< " "
						<< std::setw(4) << models[k].monomers[i].r_atoms[j].atom_name; 
						<< " "
						<< std::left << std::setw(4) << models[k].monomers[i].res3n; 
						<< " "
						<< std::right << std::setw(4) << (i+1)
						<< std::setw(5) << " "
						<< std::setw(7) << models[k].monomers[i].r_atoms[j].xc; 
						<< " "
						<< std::setw(7) << models[k].monomers[i].r_atoms[j].yc;
						<< " "
						<< std::setw(7) << models[k].monomers[i].r_atoms[j].zc;
						<< " "
						<< std::setw(5)  << "1.00"
						<< " "
						<< std::setw(5) << models[k].monomers[i].r_atoms[j].b_factor;
						<< "\n";
						cont++;
			}
		}	
		pdb_file << "ENDMDL" << endl;
	}
	pdb_file.close();
}
/*********************************************************************/
void PDB::init_from_system(const system& molecule){
	vector<pdbAtom> _atoms;
	
	for(int i=0;i<molecule.nAtoms;i++){
		pdbAtom _atom;
		_atom.atom_name = molecule[i].element;
		_atom.xc		= molecule[i].xc;
		_atom.yc		= molecule[i].yc;
		_atom.zc		= molecule[i].zc;
		_atom.b_factor	= 0.000;
		_atom.chain_name= "X";
		_atom.indx		= models[0].nAtoms + 1;
		_atom.occupancy = 1.00;
		_atom.res_indx	= models[0].nResidues + 1;
		_atom.res_name	= molecule.name.substr(0,3);		
		_atoms.emplace_back(_atom);
	}
	vector<residue> _res;
	_res.emplace_back(_atoms);
	models.emplace_back(_res);
}
/*********************************************************************/
system PDB::get_system_from_model(unsigned int model){
	vector<atom> _atoms;
	for(int i=0;i<models[model].nResidues;i++){
		for(int j=0;j<models[model].monomers[i];j++){		
			_atoms.emplace_back(models[model].monomers[i].r_atoms[j].element,
								models[model].monomers[i].r_atoms[j].xc		,
								models[model].monomers[i].r_atoms[j].yc		,
								models[model].monomers[i].r_atoms[j].zc		);
		}
	}
	system mol(_atoms,basename);
	return mol;
}
/*********************************************************************/
vector<system> PDB::get_systems(){
	vector<system> systems;
	for(int i=0;i<nModels;i++){
		systems.emplace_back( get_system_from_model(i) );
	}
	return systems;
}
///////////////////////////////////////////////////////////////////////
MOL2::MOL2(){}
/*********************************************************************/
MOL2::~MOL2(){}
/*********************************************************************/
MOL2::MOL2(const MOL2& rhs){}
/*********************************************************************/
MOL2& MOL2::operator=(const MOL2& rhs){}
/*********************************************************************/
MOL2::MOL2(MOL2&& rhs) noexcept{}
/*********************************************************************/
MOL2& MOL2::operator=(MOL2&& rhs) noexcept{}
///////////////////////////////////////////////////////////////////////
geometry::geometry():
	type(INVALID)	,
	units(Ang)		{	
}
/*********************************************************************/
geometry::geometry(	const char* file_name)	:
	type(INVALID)							,
	cUnit(Ang)								{	
	
	if ( check_file_ext(".xyz",file_name ) ){
		type = xyz_;
	}else if ( check_file_ext(".pdb",file_name ) ){
		type = pdb_;
	}
	switch( type ){
		case xyz_:
			xyz = XYZ(file_name);
			molecule = xyz.get_molecule();
			break;
		case pdb_:
			pdb = PDB(file_name);
			molecule = pdb.get_system_from_model(0);
			break;
	}
		
}
/*********************************************************************/
geometry::~geometry(){}
/*********************************************************************/
geometry::geometry(const geometry& rhs) :
	type(rhs.type)						,
	xyz(rhs.xyz)						,
	pdb(rhs.pdb)						,
	mol2(rhs.mol2)						,
	molecule(rhs.molecule)				,
	cUnit(rhs.cUnit)					{		
}
/*********************************************************************/
geometry& geometry::operator=(const geometry& rhs){
	if ( this != &rhs ){
		type	= rhs.type;
		xyz		= rhs.xyz;
		pdb		= rhs.pdb;
		mol2	= rhs.mol2;
		molecule= rhs.molecule;
		cUnit	= rhs.cUnit;
	}
	return *this;
}
/*********************************************************************/
geometry::geometry(geometry&& rhs) noexcept:
	type(rhs.type)						,
	xyz( move(rhs.xyz) )				,
	pdb( move(rhs.pdb) )				,
	mol2( move(rhs.mol2) )				,
	molecule( move(rhs.molecule) )		,
	cUnit( move(rhs.cUnit) )			{
	
}
/*********************************************************************/
geometry& geometry::operator=(geometry&& rhs) noexcept{
	if ( this != &rhs ){
		type	= rhs.type;
		xyz		= move(rhs.xyz);
		pdb		= move(rhs.pdb);
		mol2	= move(rhs.mol2);
		molecule= move(rhs.molecule);
		cUnit	= move(rhs.cUnit);
	}
	return *this;	
}
/*********************************************************************/
void geometry::read_QCPinput(const char* file_name, std:string program){
	
}
/*********************************************************************/
void geometry::read_QCPoutput(const char* file_name, std:string program,bool last){
	
}
/*********************************************************************/
void geometry::convert_to_ang(){
	const double bohrtoang = 0.52917726;
	if ( cUnit == Bohr ){
		for(int i=0;i<molecule.nAtoms;i++){
			molecule.atoms[i].x *= bohrtoang;
			molecule.atoms[i].y *= bohrtoang;
			molecule.atoms[i].z *= bohrtoang;
		}
	}

}
/*********************************************************************/
void geometry::convert_to_bohr(){
	const double angtobohr = 1.0/0.52917726;
	if ( cUnit == Ang ){
		for(int i=0;i<molecule.nAtoms;i++){
			molecule.atoms[i].x *= angtobohr;
			molecule.atoms[i].y *= angtobohr;
			molecule.atoms[i].z *= angtobohr;
		}
	}
	
}
/*********************************************************************/
void geometry::write_to_file(std::string out_name,std::string format){
	if ( format == "xyz" ){
		if ( type == xyz_){
			xyz.write_xyz(out_name);
		}else{
			xyz = XYZ(molecule);
			xyz.write_xyz(out_name);
		}		
	}else if ( format == "pdb"){
		if ( type == pdb_){
			pdb.write_pdb(out_name);
		}else{
			pdb.init_from_system(molecule);
			pdb.write_pdb(out_name);
		}
	}
}
/*********************************************************************/
void geometry::center_coord(){
	double center_x = molecule.atoms[0].xc;
	double center_y = molecule.atoms[0].yc;
	double center_z = molecule.atoms[0].zc;
	for(int i=0;i<molecule.nAtoms;i++){
		molecule.atoms[i].x -= center_x;
		molecule.atoms[i].y -= center_y;
		molecule.atoms[i].z -= center_z;
	}
}
///////////////////////////////////////////////////////////////////////