//gamessInput.cpp

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
#include <experimental/filesystem>

#include <"../include/global.h"
#include <"../include/system.h"
#include <"../include/geometry.h"
#include <"../include/gamessInput.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;
namespace fs = std::experimental::filesystem;

string gEnd = " $end\n";
//==================================================
GMS_basis::GMS_basis()	:
	type(MINI)			,
	ngauss(0)			,
	npfunc(0)			,
	ndfunc(0)			,
	ndiffuseS(false)	,
	ndiffuseP(false)	,
	gaussBasis("mini")	,
	path_to_dftb("/")	{
	
}
/**************************************************/
GMS_basis::GMS_basis(GMS_BasisSet bs):
	type(bs)						,
	ngauss(0)						,
	npfunc(0)						,
	ndfunc(0)						,
	ndiffuseS(false)					,
	ndiffuseP(false)					,
	gaussBasis("mini")				,
	path_to_dftb("/")				{
		
	switch( type ){
		case STO3G:
		gaussBasis	= "STO";
		ngauss		= 3;
		break;
		case x321G:
		gaussBasis	= "N21";
		ngauss		= 3;
		break;
		case x631G;
		gaussBasis	= "N31";
		ngauss		= 6;
		break;
		case x631Gd:
		gaussBasis	= "N31";
		ngauss		= 6;
		ndfunc		= 1;
		break;
		case x6311Gd:
		gaussBasis	= "N311";
		ngauss		= 6;
		ndfunc		= 1;
		break;
		case x6311Gdp:
		gaussBasis	= "N311";
		ngauss		= 6;
		ndfunc		= 1;
		npfunc		= 1;
		break;
		case x6311GdpD:
		gaussBasis	= "N311";
		ngauss		= 6;
		ndfunc		= 1;
		npfunc		= 1;
		ndiffusep	= true;
		break;
		case x6311Gdp2D:
		gaussBasis	= "N311";
		ngauss		= 6;
		ndfunc		= 1;
		npfunc		= 1;
		ndiffusep	= true;
		ndiffuseS	= true;
		break;
	}	
}
/**************************************************/
GMS_basis::~GMS_basis(){}
/**************************************************/
GMS_basis::GMS_basis(const GMS_basis& rhs)	:
	type(rhs.type)							,
	ngauss(rhs.ngauss)						,
	npfunc(rhs.npfunc)						,
	ndfunc(rhs.ndfunc)						,
	ndiffuseS(rhs.ndiffuseS)				,
	ndiffuseP(rhs.ndiffuseP)				,
	gaussBasis(rhs.gaussBasis)				,
	path_to_dftb(rhs.path_to_dftb)			{
}
/**************************************************/
GMS_basis::GMS_basis& operator=(const GMS_basis& rhs){
	if( this != &rhs ){
		type 		=rhs.type;
		ngauss		=rhs.ngauss;
		npfunc		=rhs.npfunc;
		ndfunc		=rhs.ndfunc;
		ndiffuseS	=rhs.ndiffuseS;
		ndiffuseP	=rhs.ndiffuseP;
		gaussBasis	=rhs.gaussBasis;
		path_to_dftb=rhs.path_to_dftb;
	}
	return *this;	
}
//==================================================
gms_group::gms_group():
	group(OTHER)	,
	grpName("")		{
	
		
	switch( group ){
		case CTRL:
			grpName	= "$contrl ";
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("runtyp=");
			inp_text.emplace_back("energy ");	// variable 2
			inp_text.emplace_back("icharg=");
			inp_text.emplace_back("0 ");		// variable 4
			inp_text.emplace_back("multi=");
			inp_text.emplace_back("1 ");		// variable 6
			inp_text.emplace_back(gEnd);
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("maxit=");
			inp_text.emplace_back("150 ");		// variable 10 
			inp_text.emplace_back("nprint=");
			inp_text.emplace_back("7 ");		// variable 12 
			inp_text.emplace_back("scftyp=");
			inp_text.emplace_back("rhf ");		// variable 14
			inp_text.emplace_back(gEnd);		
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("swoff=");
			inp_text.emplace_back("1e-03 ");	// variable 18
			inp_text.emplace_back(gEnd);
		break;
		case SYS:
			grpName = "$sys ";
			inp_text +=grpName;
		break;
		case SCF:
			grpName = "$scf ";
			inp_text +=grpName;
		break;
		case BASIS:
			grpName = "$basis ";
			inp_text +=grpName;
		break;
		case DFTB:
			grpName = "$dftb";
			inp_text +=grpName;
		break;
		case OPT:
			grpName = "$statpt ";
			inp_text +=grpName;
		break;
		case PCM:
			grpName = "$pcm ";
			inp_text +=grpName;
		break;
		case GUESS:
			grpName = "";
			inp_text +=grpName;
		break;
		case DATA:
			grpName = "";
			inp_text +=grpName;
		break;
		case VEC:
			grpName = "";
			inp_text +=grpName;
		break;
	}
	
}
/**************************************************/
gms_group::gms_group( GMS_Group grp_type )	:
	group(grp_type)							,
	grpName("")								{	
}
/**************************************************/
gms_group::~gms_group(){}
/**************************************************/
gms_group::gms_group(const gms_group& rhs)	:
	group(rhs.group)						,
	inp_text(rhs.inp_text)					,
	grpName(rhs.grpName)					{
}
/**************************************************/
gms_group& gms_group::operator=(const gms_group& rhs){
	if( this != &rhs ){
		group	= rhs.group;						
		inp_text= rhs.inp_text;					
		grpName	= rhs.grpName;		
	}
	return *this;
}
/**************************************************/
gms_group::gms_group(gms_group&& rhs) noexcept	:
	group(rhs.group)							,
	inp_text( move(rhs.inp_text) )				,
	grpName( move(rhs.grpName) )				{	
}
/**************************************************/
gms_group& operator=(gms_group&& rhs) noexcept{
	if( this != &rhs ){
		group	= move(rhs.group);						
		inp_text= move(rhs.inp_text);					
		grpName	= move(rhs.grpName);		
	}
	return *this;
}
/**************************************************/
friend std::ostream& gms_group::operator<<(std::ostream& out, gms_group& grp){
	return out << inp_text << endl;
}
//====================================================
gms_input::gms_input()	:
	multi(1)			,
	charge(0)			,
	gbasis(MINI)		,
	RunType(Energy)		,
	QM_method("HF")		,
	dfttyp("b3lyp")		,
	nprint(7)			,
	maxit(150)			,
	mwords(200)			,
	npunch(2)			,
	nsteps(10)			,
	copt(SOSCF)			,
	copt2(None)			,
	QMlevel(DFT)		,
	scftyp(RHF)			,
	solvent("none")		,
	guess("huckel")		,
	opttol(0.0001)		,
	disp("UFF")			{	
}
/**************************************************/
gms_input::~gms_input(){}
/**************************************************/
void gms_input::init(	int chg				, 
						unsigned int mpcty	,
						std::string rtype	,
						std::string method	,		
						std::string basis)	{
	
	bool def = this->load_default_options("/home/igorchem/LQQCMMtools");
	//creating the basic groups
	if ( method == "DFT " || "B3LYP "){
		QMlevel = GMS_TheoryLevel::DFT;
	}
	
						

	groups.emplace_back( GMS_Group::CTRL );
	if ( def ){
		groups[0].inp_text[2]	= rtype;
		groups[0].inp_text[4]	= std::int_to_string(charge);
		groups[0].inp_text[6]	= std::int_to_string(multi);
		groups[0].inp_text[10]	= std::int_to_string(maxit);
		groups[0].inp_text[12]	= std::int_to_string(nprint);
		groups[0].inp_text[14]	= std::int_to_string(maxit);
		groups[0].inp_text[18]	= std::int_to_string(swoff);		
	}
	
	switch ( QMlevel ){
		case DFT:
			groups[0].inp_text.emplace_back( groups[0].grp_name );
			groups[0].inp_text.emplace_back( " dfttyp=");
			groups[0].inp_text.emplace_back( dfttyp );
		break;
		case MP2:
			groups[0].inp_text.emplace_back( groups[0].grp_name );
			groups[0].inp_text.emplace_back( " mplevl=" );
			groups[0].inp_text.emplace_back( std::int_to_string(2) );
	}
	
}	
/**************************************************/		
bool gms_input::load_default_options(string path_dir){
	
	bool ndef = false;
	
	if ( charge != 0 || multi != 1 || RunType != Energy) { ndef = true; } 
	
	//reading from user file
	string file_p = path_dir + "/" + "gms_data";
	int lin = 0;
	char gms_l[50];
	bool change_p = false;
	std::ifstream gms_file(file_p.c_str());
	if ( IF_file(file_p.c_str() ) ){
		while( !gms_file.eof() ){
			gms_file.getline(gms_l,50);
			Line line(gms_file);
			if ( lin == 1 && line.words[2] == "TRUE" ){
				change_p true;
			}
			if ( change_p ){
				if ( lin == 2 )
					nprint = line.get_double(2);
				else if ( lin == 3 )
					maxit	= line.get_int(2);
				else if ( lin == 4 )
					mwords	= line.get_int(2);
				else if ( lin == 5 )
					npunch	= line.get_int(2);
				else if ( lin == 6 )
					nsteps	= line.get_int(2);
				else if ( lin == 7 )
					opttol	= line.get_int(2);
				else if ( lin == 8 )
					disp	= line.words[2];
				else if ( lin == 9 )
					swoff	= line.get_double[2];
				else if (lin == 10 )
					dfttyp	= line.words[2];
				else if (lin == 11 ){
					if ( line.words[2] == "UHF") scftyp = UHF;
					else if ( line.words[2] == "ROHF") scftyp = ROHF;
				}
			}
			lin++;
		}
	}
	if ( change_p || ndef )	
		return true;
	else 
		return false;
	
}
/**************************************************/
gms_group gms_input::load_molecule_info( system& molecule ){
	
}
/**************************************************/
void gms_input::read_input(const char* file_name){
	
}
/**************************************************/
void gms_input::restart_input(	const char* inp_name,
								const char* vec_data){
						
}
/**************************************************/
void gms_input::write_input(std::string out_name){
	
}
/**************************************************/
void gms_input::clear_directory(){
	
}
//==================================================