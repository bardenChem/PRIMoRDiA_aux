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

#include <"../include/global.h"
#include <"../include/system.h"
#include <"../include/geometry.h"

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
	if ( IF_file( pdb_name ) ){
		std::ifstream buf(pdb_file);
		while( !buf.eof() ){
			
		}
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
void PDB::split_models_in_files(std::string ref_name){
	
}
/*********************************************************************/
void PDB::read_models_from_file(const char* pdb_file){
	
}
/*********************************************************************/
void PDB::read_model_from_file(const char* pdb_file,unsigned int modN){
	
}
/*********************************************************************/
void PDB::add_model(pdbModel model){
	
}
/*********************************************************************/
void PDB::cat_pdbs(std::string file_list){
	
}
/*********************************************************************/
void PDB::remove_model(unsigned int model){
	
}
/*********************************************************************/
void PDB::write_pdb(std::string out_name){
	
}
/*********************************************************************/
void PDB::init_from_system(const system& molecule){
	
}
/*********************************************************************/
system PDB::get_system_from_model(unsigned int model){
	
}
/*********************************************************************/
std::vector<system> PDB::get_systems(){
	
}
///////////////////////////////////////////////////////////////////////
MOL2::MOL2(){
	
}
/*********************************************************************/
MOL2::~MOL2(){
	
}
/*********************************************************************/
MOL2::MOL2(const MOL2& rhs){
	
}
/*********************************************************************/
MOL2& MOL2::operator=(const MOL2& rhs){
	
}
/*********************************************************************/
MOL2::MOL2(MOL2&& rhs) noexcept{
	
}
/*********************************************************************/
MOL2& MOL2::operator=(MOL2&& rhs) noexcept{
	
}
///////////////////////////////////////////////////////////////////////
geometry::geometry(){
	
}
/*********************************************************************/
geometry::~geometry(){
	
}
/*********************************************************************/
geometry::geometry(const geometry& rhs){
	
}
/*********************************************************************/
geometry& geometry::operator=(const geometry& rhs){
	
}
/*********************************************************************/
geometry::geometry(geometry&& rhs) noexcept{
	
}
/*********************************************************************/
geometry& geometry::operator=(geometry&& rhs) noexcept{
	
}
/*********************************************************************/
void geometry::init_from_file(const char* file_name, std:string Type){
	
}
/*********************************************************************/
void geometry::read_QCPinput(const char* file_name, std:string program){
	
}
/*********************************************************************/
void geometry::read_QCPoutput(const char* file_name, std:string program,bool last){
	
}
/*********************************************************************/
void geometry::convert_to_ang(){
	
}
/*********************************************************************/
void geometry::convert_to_bohr(){
	
}
/*********************************************************************/
void geometry::write_to_file(std::string out_name,std::string format){
	
}
/*********************************************************************/
void geometry::center_coord(){
	
}
///////////////////////////////////////////////////////////////////////