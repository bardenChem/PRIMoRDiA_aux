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
using std::cout;
using std::endl;
using std::stoi;
using std::stod;

/**********************************************************/
atom::atom()	:
	xc(0.00)	,
	yc(0.00)	,
	zc(0.00)	,
	pCharge(0.0),
	aMass(1.0)	,
	aNmb(1)		,
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
	aNmb  = get_atomic_number(element);
}
/**********************************************************/
atom::atom(const atom& rhs)	:
	xc(rhs.xc)				,
	yc(rhs.yc)				,
	zc(rhs.zc)				,
	element(rhs.element)	,
	pCharge(rhs.pCharge)	,
	aMass(rhs.aMass)		,
	aNmb(rhs.aNmb0			{
	
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
		aNmb	= rhs.aNmb;
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
	aMass(rhs.aMass)			,
	aNmb(rhs.aNmb)				{	
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
		aNmb	= rhs.aNmb;
	}
	return *this;
}
/**********************************************************/
void atom::set_element( std::string Type ){
	element = Type;
	aMass 	= get_atom_mass(element);
	aNmb  	= get_atomic_number(element);
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
/**********************************************************/
double atom::set_pCharge(double chg){
	pCharge = chg;
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

	for(unsigned int i=0;i<3;i++){
		ver_inf[i] = 0.000;
		ver_sup[i] = 0.000;
	}
}
/**********************************************************/
system::system(std::vector<atom> ats,
				std::string nme		):
				
	nAtoms( ats.size() ),
	nElectrons(0)		,
	nHydrogens(0)		,
	name(nme)			,
	type("none")		,
	fCharge(0.0)		,
	molar_mass(0.0)		,
	atoms(ats)			{
		
	for(unsigned int i=0;i<nAtoms;i++){
		nElectrons	+= atoms[i].aNmb;
		molar_mass	+= atoms[i].molar_mass;
		fCharge		+= atoms[i].pCharge;
		if ( atoms[i].aNmb == 1 ) nHydrogens++;
		if( i == 0 ){
			ver_inf[0] = atoms[i].xc;
			ver_inf[0] = atoms[i].yc;
			ver_inf[0] = atoms[i].zc;
			ver_sup[0] = atoms[i].xc;
			ver_sup[0] = atoms[i].yc;
			ver_sup[0] = atoms[i].zc;
		}
		else{
			if ( ver_inf[0] > atoms[i].xc ) ver_inf[0] = atoms[i].xc;
			if ( ver_inf[1] > atoms[i].yc ) ver_inf[1] = atoms[i].yc;
			if ( ver_inf[2] > atoms[i].zc ) ver_inf[2] = atoms[i].zc;
			if ( ver_sup[0] < atoms[i].xc ) ver_sup[0] = atoms[i].xc;
			if ( ver_sup[1] < atoms[i].yc ) ver_sup[1] = atoms[i].yc;
			if ( ver_sup[2] < atoms[i].zc ) ver_sup[2] = atoms[i].zc;
		}		
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
	if ( atom_name[0] == "H" ) return true;
	else if ( atom_name.substr( atom_name.begin(),2 ) == "1H" ) return true;
	else if ( atom_name.substr( atom_name.begin(),2 ) == "2H" ) return true;
	else if ( atom_name.substr( atom_name.begin(),2 ) == "3H" ) return true;
	else false;
}
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
pdbModel::pdbModel(const char* pdb_file, int& mdl){
	if ( !check_file_ext(".pdb",pdb_name) )	{
		cout << "Warning! The file has wrong extension name!" << endl;
	}
	
	vector<pdbAtom> tmp_atoms;
	
	char pdb_line[100];
	int line = 0;
	string old_res = "0";
	string curr_res = "_";
	string pdb_line = "";
	
	
	if ( IF_file( pdb_name ) ){
		std::ifstream buf(pdb_file);
		while( !buf.eof() ){
			if line >= mdl{
				buf.getline(pdb_line,100);
				string word(pdb_line,0,6);
				if ( word == "ATOM  " || word == "HETATM" ) {					
					pdbAtom _atom(pdb_line);					
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
	
	unsigned int i,j,cont = 0;
	for(i;i<nResidues;i++){
		for(j;j<monomers[i].nAtoms;j++){
			pdb_file<< std::setw(6) << std::left  << "ATOM" 
					<< " "
					<< std::setw(4) << std::right  << (cont+1) 
					<< " "
					<< std::setw(4) << monomers[i].r_atoms[j].atom_name; 
					<< " "
					<< std::left << std::setw(4) << monomers[i].res3n; 
					<< " "
					<< std::right << std::setw(4) << (i+1)
					<< std::setw(5) << " "
					<< std::setw(7) << monomers[i].r_atoms[j].xc; 
					<< " "
					<< std::setw(7) << monomers[i].r_atoms[j].yc;
					<< " "
					<< std::setw(7) << monomers[i].r_atoms[j].zc;
					<< " "
					<< std::setw(5)  << "1.00"
					<< " "
					<< std::setw(5) << monomers[i].r_atoms[j].b_factor;
					<< "\n";
					cont++;
		}
	}
	pdb_file << "ENDMDL" << endl;
	pdb_file.close();
}
/*********************************************************/
void pdbModel::remove_atom(unsigned int res, unsigned int at){	
	
	monomers[res].r_atoms.erase( monomers[res].r_atoms.begin(), at);
	monomers[res].nAtoms--;
	nAtoms--;
	
}
/*********************************************************/
void pdbModel::remove_residue(unsigned int i){
	nAtoms -= monomers[i].nAtoms;
	monomers.erase( monomers.begin(), i);
	nResidues--;
}
/*********************************************************/
void pdbModel::prune_atoms(){
	unsigned int i,j = 0;	
	for(i;i<nResidues;i++){
		for(j;j<monomers[i].nAtoms;j++){
			if ( monomers[i].r_atoms[j].res_name[0] == "B" ) 
				this->remove_atom( i, j ); 	
		}
	}
}
/*********************************************************/
void pdbModel::remove_waters(){
	unsigned int i = 0;
	for(i;i<nResidues;i++){
		if ( monomers[i].type == WAT ) this->remove_residue(i);
	}
}
/*********************************************************/
void pdbModel::remove_waters(double radius, unsigned int res){
	unsigned int i = 0;
	double refXC, refYC, refZC, distTemp = 0.000;
	refXC = monomers[res].r_atoms[0].xc;
	refYC = monomers[res].r_atoms[0].yc;
	refZC = monomers[res].r_atoms[0].zc;
	
	for(i;i<nResidues;i++){
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
		if ( monomers[i].type == ION ) this->remove_residue(i)
}
/*********************************************************/
void pdbModel::split_complex(std::string mol){
	unsigned int i = nResidues-1;
	for(i;i>0;i--){
		if ( mol == monomers[i].res3n ) {
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
	monomers.emplace_black( temp.monomers[0] );
	nResidues++;
	nAtoms += temp.monomers[0].nAtoms;
}
//////////////////////////////////////////////////////////