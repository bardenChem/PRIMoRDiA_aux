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

//==================================================
GMS_basis::GMS_basis()	:
	type(MINI)			,
	ngauss(0)			,
	npfunc(0)			,
	ndfunc(0)			,
	ndiffuseS(0)		,
	ndiffuseP(0)		,
	gaussBasis("mini")	,
	path_to_dftb("/")	{
	
}
/**************************************************/
GMS_basis::GMS_basis(GMS_BasisSet bs):
	type(bs)						,
	ngauss(0)						,
	npfunc(0)						,
	ndfunc(0)						,
	ndiffuseS(0)					,
	ndiffuseP(0)					,
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
		ndiffusep	= 1;
		break;
		case x6311Gdp2D:
		gaussBasis	= "N311";
		ngauss		= 6;
		ndfunc		= 1;
		npfunc		= 1;
		ndiffusep	= 1;
		ndiffuseS	= 1;
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
	inp_text(" ")	,
	grpName("")		{
	
}
/**************************************************/
gms_group::gms_group( GMS_Group grp_type )	:
	group(grp_type)							,
	inp_text(" ")							,
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
	nprint(7)			,
	maxit(150)			,
	mwords(200)			,
	npunch(2)			,
	nsteps(10)			,
	copt(SOSCF)			,
	copt2(None)			,
	QMlevel(DFT)		,
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
						std::string method	,
						std::string basis)	{
							

}	
/**************************************************/		
void gms_input::load_default_options(string path_dir){
	string file_p = path_dir + "/" + "gms_data";
	int lin = 0;
	std::ifstream gms_file(file_p.c_str());
	if ( IF_file(file_p.c_str() ) ){
		while( !gms_file.eof() ){
			
		}
	}
	
}
/**************************************************/
void gms_input::load_molecule_info( system& molecule ){
	
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