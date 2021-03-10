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
	ndiffuseS(false)				,
	ndiffuseP(false)				,
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
			grpName	= " $contrl ";
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
			grpName = " $sys ";
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("mwords=");
			inp_text.emplace_back("200 ");		// variable 2 
			inp_text.emplace_back(gEnd);			
		break;
		case SCF:
			grpName = " $scf ";
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("dirscf=.t.");
			inp_text.emplace_back("npunch="); 
			inp_text.emplace_back("2 ");		// variable 3 
			inp_text.emplace_back(gEnd);
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("soscf=");
			inp_text.emplace_back(".t. ");		// variable 7
			inp_text.emplace_back(gEnd);
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("diis=");
			inp_text.emplace_back(".t. ");		// variable 11
			inp_text.emplace_back("ethrsh=");
			inp_text.emplace_back(".t. ");		// variable 13
			inp_text.emplace_back(gEnd);
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("shift=");
			inp_text.emplace_back(".f. ");		// variable 19
			inp_text.emplace_back("rstrct=");
			inp_text.emplace_back(".f. ");		// variable 21
			inp_text.emplace_back("damp=");
			inp_text.emplace_back(".f. ");		// variable 23
			inp_text.emplace_back(gEnd);
		break;
		case DFTB:
			grpName = " $dftb ";
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("scc=.true.");
			inp_text.emplace_back("ndftb=");
			inp_text.emplace_back("2 ");		// variable 3
			inp_text.emplace_back("dampxh=.t.");
			inp_text.emplace_back("dampex=");
			inp_text.emplace_back("4.2 ");		// variable 6 
			inp_text.emplace_back("itypmx=");
			inp_text.emplace_back("-1");		// variable 8 
			inp_text.emplace_back(" $dftbsk \n");
			inp_text.emplace_back( "   C H "+ path_to_dftb +"/C-H.skf\n");
			inp_text.emplace_back( "   C O "+ path_to_dftb +"/C-O.skf\n");
			inp_text.emplace_back( "   C N "+ path_to_dftb +"/C-N.skf\n");
			inp_text.emplace_back( "   H C "+ path_to_dftb +"/H-C.skf\n");
			inp_text.emplace_back( "   H H "+ path_to_dftb +"/H-H.skf\n");
			inp_text.emplace_back( "   H O "+ path_to_dftb +"/H-O.skf\n");
			inp_text.emplace_back( "   H N "+ path_to_dftb +"/H-N.skf\n");
			inp_text.emplace_back( "   O C "+ path_to_dftb +"/O-C.skf\n");
			inp_text.emplace_back( "   O H "+ path_to_dftb +"/O-H.skf\n");
			inp_text.emplace_back( "   O N "+ path_to_dftb +"/O-N.skf\n");
			inp_text.emplace_back( "   O O "+ path_to_dftb +"/O-O.skf\n");
			inp_text.emplace_back( "   N C "+ path_to_dftb +"/N-C.skf\n");
			inp_text.emplace_back( "   N H "+ path_to_dftb +"/N-H.skf\n");
			inp_text.emplace_back( "   N O "+ path_to_dftb +"/N-O.skf\n");
			inp_text.emplace_back( "   N N "+ path_to_dftb +"/N-N.skf\n");
			inp_text.emplace_back( "   S S "+ path_to_dftb +"/S-S.skf\n");
			inp_text.emplace_back( "   S H "+ path_to_dftb +"/S-H.skf\n");
			inp_text.emplace_back( "   S C "+ path_to_dftb +"/S-C.skf\n");
			inp_text.emplace_back( "   S O "+ path_to_dftb +"/S-O.skf\n");
			inp_text.emplace_back( "   S N "+ path_to_dftb +"/S-N.skf\n");
			inp_text.emplace_back( "   N S "+ path_to_dftb +"/N-S.skf\n");
			inp_text.emplace_back( "   O S "+ path_to_dftb +"/O-S.skf\n");
			inp_text.emplace_back( "   H S "+ path_to_dftb +"/H-S.skf\n");
			inp_text.emplace_back( "   C S "+ path_to_dftb +"/C-S.skf\n");
			inp_text.emplace_back(gEnd)
		break;
		case OPT:
			grpName = " $statpt ";
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("opttol=");
			inp_text.emplace_back("0.0001 ");	// variable 2
			inp_text.emplace_back("nstep=");
			inp_text.emplace_back("150 ");		// variable 4 
			inp_text.emplace_back("method=");
			inp_text.emplace_back("QA ");		// variable 6 
			inp_text.emplace_back(gEnd);					
		break;
		case PCM:
			grpName = " $pcm ";
			inp_text.emplace_back(grpName);
			inp_text.emplace_bacl(" solvnt=");
			inp_text.emplace_bacl("h2o "); 		// variable 2
			inp_text.emplace_bacl(gEnd);
		break;
		case GUESS:
			grpName = " $guess ";
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("guess=");
			inp_text.emplace_back("huckel ");
			inp_text.emplace_back(gEnd);			
		break;
		case DATA:
			grpName = " $data\n";
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("Molecule Specification\n");
			inp_text.emplace_back("C1\n");
		break;
		case BASIS:
			grpName = " $basis ";
			inp_text.emplace_back(grpName); 
			inp_text.emplace_back("gbasis=");
			inp_text.emplace_back("sto"); 			// variable 2 
			inp_text.emplace_back(" polar=");
			inp_text.emplace_back("");			// varible 4 		
			inp_text.emplace_back(gEnd);
			inp_text.emplace_back(grpName);
			inp_text.emplace_back("ngauss=");
			inp_text.emplace_back("0 ");			// varible 8
			inp_text.emplace_back("ndfunc=");
			inp_text.emplace_back("0 ");			// varible 10 
			inp_text.emplace_back("ndpunc=");
			inp_text.emplace_back("0 ");			// varible 12
			inp_text.emplace_back("diffp=");
			inp_text.emplace_back(".f.");		// varible 14
			inp_text.emplace_back("diffs");		
			inp_text.emplace_back(".f.");		// varible 16
			inp_text.emplace_back(gEnd);			
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
	swoff(0.001)		,
	swdiis(0.005)		,
	ethrsh(2.0)			,
	alg("QA")			,
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
	}else if ( method == "AM1" ){
		QMlevel = GMS_TheoryLevel::Semi;
	}else if ( method == "DFTB2" ){
		QMlevel = GMS_TheoryLevel::DFTB2;
		gbasis	= GMS_BasisSet::smDFTB;
	}else if ( method == "DFTB3" ){
		QMlevel = GMS_TheoryLevel::DFTB3;
		gbasis	= GMS_BasisSet::smDFTB;
	}
						
	// creating the basic groups
	groups.emplace_back( GMS_Group::CTRL );
	groups.emplace_back( GMS_Group::SYS  );
	GMS_basis basis_set ( this->init_basis(basis) );
	groups.emplace_back( GMS_Group::BASIS );
	
	if ( def ){
		groups[0].inp_text[2]	= rtype;
		groups[0].inp_text[4]	= std::int_to_string(charge);
		groups[0].inp_text[6]	= std::int_to_string(multi);
		groups[0].inp_text[10]	= std::int_to_string(maxit);
		groups[0].inp_text[12]	= std::int_to_string(nprint);
		groups[0].inp_text[14]	= std::int_to_string(maxit);
		groups[0].inp_text[18]	= std::int_to_string(swoff);		
	}
	
	
	groups[2].inp_text.emplace_back("gbasis=");
	groups[2].inp_text.emplace_back(basis_set.gaussBasis); 		
	groups[2].inp_text.emplace_back(gEnd);
	
	
	groups[2].inp_text.emplace_back(grpName);
	groups[2].inp_text.emplace_back("ngauss=");
	groups[2].inp_text.emplace_back("0 ");			
	groups[2].inp_text.emplace_back("ndfunc=");
	groups[2].inp_text.emplace_back("0 ");			
	groups[2].inp_text.emplace_back("ndpunc=");
	groups[2].inp_text.emplace_back("0 ");			
	groups[2].inp_text.emplace_back("diffp=");
	groups[2].inp_text.emplace_back(".f.");		
	groups[2].inp_text.emplace_back("diffs");		
	groups[2].inp_text.emplace_back(".f.");		
	groups[2].inp_text.emplace_back(gEnd);		
	
	
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
		case DFTB2:
		break;
		case DFTB3:
		break;
	}
	
	switch ( QMlevel ){
		case HF:
		case MP2:
		case DFT:
			groups[2].inp_text.emplace_back(grpName);
			groups[2].inp_text.emplace_back("ngauss=");
			groups[2].inp_text.emplace_back( int_to_string(basis_set.ngauss) );
			groups[2].inp_text.emplace_back("ndfunc=");
			groups[2].inp_text.emplace_back("0 ");			
			groups[2].inp_text.emplace_back("ndpunc=");
			groups[2].inp_text.emplace_back("0 ");			
			groups[2].inp_text.emplace_back("diffp=");
			groups[2].inp_text.emplace_back(".f.");		
			groups[2].inp_text.emplace_back("diffs");		
			groups[2].inp_text.emplace_back(".f.");		
			groups[2].inp_text.emplace_back(gEnd);
		break;
		
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
				else if ( lin == 10 )
					dfttyp	= line.words[2];
				else if ( lin == 11 ){
					if ( line.words[2] == "UHF") scftyp = UHF;
					else if ( line.words[2] == "ROHF") scftyp = ROHF;
				}
				else if ( lin == 12 )
					ethrsh = line.get_double[2];
				else if ( lin == 13 )
					swdiis = line.get_double[2];
				else if ( lin == 14 )
					guess  = line.words[2];
				else if ( lin == 15 )
					damph  = line.get_double[2];
				else if ( lin == 16 )
					alg    = line.words[2];
				else if ( lin == 17 )
					solvent= line.words[2];				
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
void gms_input::load_molecule_info( system& molecule ){
	gms_group data_group(GMS_Group::DATA);
	for(int i=0;i<molecule.nAtoms;i++){
		data_group.inp_text.emplace_back( molecule.atoms[i].element );
		data_group.inp_text.emplace_back( " " );
		data_group.inp_text.emplace_back( molecule.atoms[i].aNmb );
		data_group.inp_text.emplace_back( " " );
		data_group.inp_text.emplace_back( molecule.atoms[i].xc );
		data_group.inp_text.emplace_back( " " );
		data_group.inp_text.emplace_back( molecule.atoms[i].yc );
		data_group.inp_text.emplace_back( " " );
		data_group.inp_text.emplace_back( molecule.atoms[i].zc );
		data_group.inp_text.emplace_back( " " );	
	}
	data_group.inp_text.emplace_back(gEnd);
	groups.emplace_back( move (data_group ) );
	
}
/**************************************************/
GMS_basis gms_input::init_basis(std::string& bsis_nm){
	if 		( bsis_nm == "MINI"   		) gbasis = MINI;
	else if ( bsis_nm == "STO3G"  		) gbasis = STO3G;
	else if ( bsis_nm == "321G"   		) gbasis = x321G;
	else if ( bsis_nm == "631G"   		) gbasis = x631G;
	else if ( bsis_nm == "631G*"  		) gbasis = x631Gd;
	else if ( bsis_nm == "6311G*" 		) gbasis = x6311Gd;
	else if ( bsis_nm == "6311G**"		) gbasis = x6311Gdp;
	else if ( bsis_nm == "6311G(2d)+"	) gbasis = x6311GdpD;
	else if ( bsis_nm == "6311G(2d)++"	) gbasis = x6311Gdp2D;
	
	GMS_basis basisset(gbasis);
	return basisset;
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