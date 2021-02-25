//system.cpp

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

#include <"../include/global.h"
#include <"../include/system.h"

using std::vector;
using std::string;
using std::move;

/**********************************************************/
atom::atom()	:
	xc(0.00)	,
	yc(0.00)	,
	zc(0.00)	,
	pCharge(0.0),
	aMass(1.0)	,
	element("H"){	
}
/**********************************************************/
atom::~atom(){}
/**********************************************************/
atom::atom(double x			,
		   double y			, 
		   double z			, 
		   std::string type):
		   xc(x)			,
		   yc(y)			,
		   zc(z)			,
		   element(type)	{
	
	aMass = get_atom_mass(element);
}
/**********************************************************/
atom::atom(const atom& rhs)	:
	xc(rhs.xc)				,
	yc(rhs.yc)				,
	zc(rhs.zc)				,
	element(rhs.element)	,
	pCharge(rhs.pCharge)	,
	aMass(rhs.aMass)		{
	
}
/**********************************************************/
atom& atom::operator=(const atom& rhs){
	if ( this != &rhs ){
		xc		= rhs.xc;
		yc		= rhs.yc;
		zc		= rhs.zc;
		element	= rhs.element;
		pCharge	= rhs.pCharge;
		aMass	= rhs.aMass;	
	}
	return *this;
}
/**********************************************************/
atom::atom(atom&& rhs) noexcept:
	xc(rhs.xc)					,
	yc(rhs.yc)					,
	zc(rhs.zc)					,
	element( move(rhs.element) ),
	pCharge(rhs.pCharge)		,
	aMass(rhs.aMass)			{	
}
/**********************************************************/
atom& atom::operator=(atom&& rhs) noexcept{
	if ( this != &rhs ){
		xc		= rhs.xc;
		yc		= rhs.yc;
		zc		= rhs.zc;
		element	= move(rhs.element);
		pCharge	= rhs.pCharge;
		aMass	= rhs.aMass;
	}
	return *this;
}
/**********************************************************/
void atom::set_element( std::string Type ){
	element = Type;
}
/**********************************************************/
void atom::set_coord(double x, double y, double z){
	xc = x;
	yc = y;
	zc = z;	
}
/**********************************************************/
double atom::get_distance(const atom& a2){
	double dist = 0.0;
	dist = (xc - a2.xc)*(xc - a2.xc);
	dist+= (yc - a2.yc)*(yc - a2.yc);
	dist+= (zc - a2.zc)*(zc - a2.zc);
	return sqrt(dist);
}
////////////////////////////////////////////////////////////
system::system()	:
	nAtoms(0)		,
	nElectrons(0)	,
	nHydrogens(0)	,
	name("nonamed")	,
	type("none")	,
	fCharge(0.0)	,
	molar_mass(0.0)	{

	for(int i=0;i<3;i++){
		ver_inf[i] = 0.000;
		ver_sup[i] = 0.000;
	}
}
/**********************************************************/
system::~system(){}
/**********************************************************/
system::system(const system& rhs):
	nAtoms(rhs.nAtoms)			,
	nElectrons(rhs.nElectrons)	,
	nHydrogens(rhs.nHydrogens)	,
	name(rhs.name)				,
	type(rhs.type)				,
	atoms(rhs.atoms)			,
	molar_mass(rhs.molar_mass)  {
		
	for(int i=0;i<3;i++){
		ver_inf[i] = rhs.ver_inf[i];
		ver_sup[i] = rhs.ver_sup[i];
	}
	
}
/**********************************************************/
system& system::operator=(const system& rhs){
	if( this != &rhs ){
		nAtoms 		= rhs.nAtoms;
		nElectrons 	= rhs.nElectrons;	
		nHydrogens 	= rhs.nHydrogens;	
		name 		= rhs.name;		
		type 		= rhs.type;			
		atoms 		= rhs.atoms;			
		molar_mass 	= rhs.molar_mass;  		
		
		for(int i=0;i<3;i++){
			ver_inf[i] = rhs.ver_inf[i];
			ver_sup[i] = rhs.ver_sup[i];
		}
	}
	return *this;
}
/*********************************************************/
system::system(system&& rhs) noexcept:
	nAtoms(rhs.nAtoms)				,
	nElectrons(rhs.nElectrons)		,
	nHydrogens(rhs.nHydrogens)		,
	name( move(rhs.name) )			,
	type( move(rhs.type) )			,
	atoms( move(rhs.atoms) ) 		,
	molar_mass(rhs.molar_mass)		{
	
	for(int i=0;i<3;i++){
		ver_inf[i] = rhs.ver_inf[i];
		ver_sup[i] = rhs.ver_sup[i];
	}
}
/*********************************************************/
system& system::operator=(system&& rhs) noexcept{
	if( this != &rhs ){
		nAtoms 		= rhs.nAtoms;
		nElectrons 	= rhs.nElectrons;	
		nHydrogens 	= rhs.nHydrogens;	
		name 		= move(rhs.name);		
		type 		= move(rhs.type);			
		atoms 		= move(rhs.atoms);			
		molar_mass 	= rhs.molar_mass;  		
		
		for(int i=0;i<3;i++){
			ver_inf[i] = rhs.ver_inf[i];
			ver_sup[i] = rhs.ver_sup[i];
		}
	}
	return *this;
}
/*********************************************************/
void system::add_atom(atom a){
	atoms.emplace_back( move(a) );
	nAtoms++;
}
/*********************************************************/
void system::add_atom(double x		,
					  double y		, 
					  double z		, 
					  std::string el){
						  
	atoms.emplace(x,y,z,el);
	nAtoms++;
}
/*********************************************************/
void system::remove_atom(unsigned int i){
	atoms.erase( atoms.begin()+i );
	nAtoms--;
}
////////////////////////////////////////////////////////////
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
	b_factor(rhs.b_factor)				,
	sideC( move(rhs.sideC) )			,
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
		sideC		= move(rhs.sideC);
		xc			= rhs.xc;
		yc			= rhs.yc;					
		zc			= rhs.zc;
	}
	return *this;
}
/*********************************************************/
bool pdbAtom::is_hydrogen(){
	if ( atom_name[0] == "H" ) return true;
	else if ( atom_name.substr( atom_name.begin(),2 ) == "1H" ) return true;
	else if ( atom_name.substr( atom_name.begin(),2 ) == "2H" ) return true;
	else if ( atom_name.substr( atom_name.begin(),2 ) == "3H" ) return true;
	else false;
}
////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////
residue::residue()	:
	res1n("n")		,
	res3n("UNK")	,
	type(UNK)		,
	AAname(OTH)     ,
	ligand(false)	,
	terminal(false)	,
	first(false)    ,
	nHydrogens(0)	,
	fCharge(0)		,
	pdb_index(0)	,
	nAtoms(0)		{
}
/*********************************************************/
residue::residue( vector<pdbAtom> resAtoms	,
				int resType					,
				int resMon 					):
				
	r_atoms(resAtoms)	,
	type(resType)		,
	AAname(OTH)			,
	ligand(false)		,
	terminal(false)		,
	first(false)        ,
	pdb_index			,
	fCharge(0)			{
		
	res1n 	= get_res1n(resMon);
	res3n 	= get_res3n(resMon);
	AAname 	= resMon;
	
	
	for(int=0;i<resAtoms.size();i++){
		if ( resAtoms[i].is_hydrogen() ) nHydrogens++;
	}
	 
}
/*********************************************************/
residue::~residue(){}
/*********************************************************/
residue::residue(const residue& rhs):
	res1n(rhs.res1n)				,
	res3n(rhs.res3n)				,
	type(rhs.type)					,
	AAname(rhs.AAname)				,
	ligand(rhs.ligand)				,
	terminal(rhs.terminal)			,
	nHydrogens(rhs.nHydrogens)		,
	fCharge(rhs.fCharge)			,
	pdb_index(rhs.pdb_index)        ,
	nAtoms(rhs.nAtoms)				,
	r_atoms(rhs.r_atoms)			{	
}
/*********************************************************/
residue& residue::operator=(const residue& rhs){
	if ( this != &rhs ){
		res1n 		= rhs.res1n;				
		res3n 		= rhs.res3n;			
		type 		= rhs.type;
		AAname		= rhs.AAname;
		ligand		= rhs.ligand;				
		terminal	= rhs.terminal;	
		first		= rhs.first;
		pdb_index	= rhs.pdb_index;
		nHydrogens	= rhs.nHydrogens;		
		fCharge 	= rhs.fCharge;			
		nAtoms 		= rhs.nAtoms;				
		r_atoms 	= rhs.r_atoms;
	}
	return *this;
}
/*********************************************************/
residue::residue( residue&& rhs) noexcept:
	res1n( move(rhs.res1n) )			,
	res3n( move(rhs.res3n) )			,
	type( rhs.type )					,
	AAname(rhs.AAname)					,
	ligand(rhs.ligand)					,
	terminal(rhs.terminal)				,
	first(rhs.first)					,
	pdb_index(rhs.pdb_index )			,
	nHydrogens(rhs.nHydrogens)			,
	fCharge(rhs.fCharge)				,
	nAtoms( move(rhs.nAtoms) )			,
	r_atoms( move(rhs.r_atoms) )		{
	
}
/*********************************************************/
residue& residue::operator=( residue&& rhs) noexcept{
	if ( this != &rhs ){
		res1n 		= move(rhs.res1n);				
		res3n 		= move(rhs.res3n);			
		type 		= rhs.type;
		AAname		= rhs.AAname;
		ligand		= rhs.ligand;				
		terminal	= rhs.terminal;	
		first		= rhs.first;
		pdb_index	= rhs.pdb_index;
		nHydrogens	= rhs.nHydrogens;		
		fCharge 	= rhs.fCharge;			
		nAtoms 		= rhs.nAtoms;				
		r_atoms 	= move(rhs.r_atoms);
	}
	return *this;	
}
/*********************************************************/
bool residue::is_ion(){
	if ( r_atoms[0].res_name == "Cl-" ||
		 r_atoms[0].res_name == "Na+" ||
		 r_atoms[0].res_name == "K+"  ||
		 r_atoms[0].res_name == "Mg+" ||
		 r_atoms[0].res_name == "Zn+" ||
		 r_atoms[0].res_name == "SO4" ) {
	
		type = ION;
		return true;
	}
}
/*********************************************************/
void residue::set_charge(){
	
	int base_HN = get_AAnHy(AAname);
	
	if ( type == AA ){
		switch ( AAname ){
			case ASP:
			case GLU:
				if ( first || terminal ){
					fCharge = nHydrogens - base_HN - 2;				
				}else{
					fCharge = nHydrogens - base_HN  -1;
				} 
			break;
			/-------
			case ARG:
			case LYS:
				if ( first || terminal ) {
					fCharge = nHydrogens - base_HN;					
				}
				else {
					fCharge = nHydrogens - base_HN +1;
				}
			break;
			/-------
			case HIS:
				if ( first || terminal ) {
					fCharge = nHydrogens - base_HN;
				}
				else {
					fCharge = nHydrogens - base_HN +1;
				}
			break;
			/-------
			case CYS:
				fCharge = nHydrogens - base_HN - 1;
			break;
			/-------
			case default:
				if ( first || terminal ){
					fCharge = nHydrogens - base_HN -1;
				}
				else fCharge = 0;
			break;
			
		}		
	}
	else if ( type == ION ){
		if ( r_atoms[0].res_name == "Cl-" )
			fCharge = -1;
		else if  ( 	r_atoms[0].res_name == "Na+" ||
					r_atoms[0].res_name == "K+"  )
			fCharge = +1;
		else if (   r_atoms[0].res_name == "Mg+" ||
					r_atoms[0].res_name == "Zn+" )
			fCharge = +2;
		else if (  r_atoms[0].res_name == "SO4" ) 
			fCharge = -2;
		else fCharge = 0;
		
	}	
}
////////////////////////////////////////////////////////////
pdbModel::pdbModel(){
	
}
/*********************************************************/
pdbModel::pdbModel(std::vector<residue> residues){
	
}
/*********************************************************/
pdbModel::~pdbModel(){}
/*********************************************************/
pdbModel::pdbModel(const pdbModel& rhs){
	
}
/*********************************************************/
pdbModel& pdbModel::operator=(const pdbModel& rhs){
	
}
/*********************************************************/
pdbModel::pdbModel(pdbModel&& rhs) noexcept:
{
	
}
/*********************************************************/
pdbModel& pdbModel::operator=(pdbModel&& rhs) noexcept:
{
	
}
/*********************************************************/
void pdbModel::write_model(std::string out_name){
	
}
/*********************************************************/
void pdbModel::prune_atoms(){
	
}
/*********************************************************/
void pdbModel::remove_waters(){
	
}
/*********************************************************/
void pdbModel::remove_waters(double radius){
	
}
/*********************************************************/
void pdbModel::remove_ions(){
	
}
/*********************************************************/
void pdbModel::remove_residue(unsigned int i){
	
}
/*********************************************************/
void pdbModel::update_residues(){
	
}
/*********************************************************/
void pdbModel::split_complex(std::string mol){
	
}
/*********************************************************/
void pdbModel::built_complex(const char* pdb_mol){
	
}
//////////////////////////////////////////////////////////