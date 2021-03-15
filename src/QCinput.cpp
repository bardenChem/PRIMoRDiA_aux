//QCinput.cpp

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
#include <algorithm>
#include <experimental/filesystem>

#include "../include/global.h"
#include "../include/molecule.h" 
#include "../include/geometry.h"
#include "../include/gamessInput.h"
#include "../include/orcaInput.h"
#include "../include/mopacInput.h"
#include "../include/QCinput.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;
namespace fs = std::experimental::filesystem;

/************************************************************/
QCPinput::QCPinput()	:
	program(0)			,
	base_multi(1)		,
	base_charge(0)		,
	QMmethod("HF")		,
	base_basis("3-21G")	,
	runtype("Energy")	,
	geo_ext("xyz")		{
		
}
/************************************************************/
QCPinput::QCPinput(	string _geo_ext	  ,
					string _runtype	  ,
					string _basis	  ,
					string _method   ):
	program(0)							,
	base_basis(_basis)					,
	base_charge(0)						,
	base_multi(1)						,
	geo_ext(_geo_ext)					,
	runtype(_runtype)					,
	QMmethod(_method)					{
}
/************************************************************/
QCPinput::~QCPinput(){}
/************************************************************/
void QCPinput::make_input_from_folder(	int QCP				,
										unsigned int bs_mlt	,
										int bs_chg)			{
											
	program 	= QCP;
	base_charge = bs_chg;
	base_multi 	= bs_mlt;
										
	fs::path c_path = fs::current_path();
	std::vector<string> fnames; 
		
	for ( const auto & entry : fs::directory_iterator(c_path) ){
		string tmp_name = entry.path();
		if ( check_file_ext( geo_ext,tmp_name.c_str() ) ){
			fnames.push_back( tmp_name );
		}
	}
	//cout << c_path << endl;
	
	std::sort( fnames.begin(),fnames.end() );
	
	for(unsigned int i=0;i<fnames.size();i++){
		geometry geo_file( fnames[i].c_str() );
		string oname = remove_extension( fnames[i].c_str() );
		switch( program ){
			case package::GAMESS:
				gms_input g_input();
				g_input.init(base_charge,base_multi,runtype,QMmethod, base_basis);
				g_input.load_molecule_info(geo_file.molecule);
				g_input.write_input(oname);
			break;
			case package::MOPAC:
				mopac_input mpc_input();
				mpc_input.init(base_charge,base_multi,runtype,QMmethod,base_basis);
				mpc_input.write_file(geo_file.molecule,oname);
			break;
			case package::ORCA:
				orcaInput orc_input();
				orc_input.init(geo_file.molecule,base_charge,base_multi,QMmethod,runtype,1,base_basis);
				orc_input.write_input_file(oname);
			break;
		}
	}
}
/************************************************************/								
void QCPinput::make_input_from_folder_FD(int QCP			,
										 unsigned int bs_mlt,
										 int bs_chg			,
										 int chg_diff		){
	
	program 	= QCP;
	base_charge = bs_chg;
	base_multi 	= bs_mlt;
	
	int cat_charge = base_charge + chg_diff;
	int an_charge = base_charge - chg_diff;
	int ion_multi = base_multi;
	
	if ( chg_diff % 2 != 0 ) {
		ion_multi = 2;
	}
										
	fs::path c_path = fs::current_path();
	std::vector<string> fnames; 
		
	for ( const auto & entry : fs::directory_iterator(c_path) ){
		string tmp_name = entry.path();
		if ( check_file_ext( geo_ext,tmp_name.c_str() ) ){
			fnames.push_back( tmp_name );
		}
	}
	//cout << c_path << endl;
	
	std::sort( fnames.begin(),fnames.end() );
	
	for(unsigned int i=0;i<fnames.size();i++){
		geometry geo_file( fnames[i].c_str() );
		string oname = remove_extension( fnames[i].c_str() );
		switch( program ){
			case package::GAMESS:
				gms_input g_input_neu();
				gms_input g_input_cat();
				gms_input g_input_an();
				g_input_neu.init(base_charge,base_multi,runtype,QMmethod, base_basis);
				g_input_neu.load_molecule_info(geo_file.molecule);
				g_input_neu.write_input(oname);
				g_input_cat.init(cat_charge,ion_multi,runtype,QMmethod, base_basis);
				g_input_cat.load_molecule_info(geo_file.molecule);
				g_input_cat.write_input(oname);
				g_input_an.init(an_charge,ion_multi,runtype,QMmethod, base_basis);
				g_input_an.load_molecule_info(geo_file.molecule);
				g_input_an.write_input(oname);			
			break;
			case package::MOPAC:
				mopac_input mpc_input_neu();
				mopac_input mpc_input_cat();
				mopac_input mpc_input_an();
				mpc_input_neu.init(base_charge,base_multi,runtype,QMmethod,base_basis);
				mpc_input_neu.write_file(geo_file.molecule,oname);
				mpc_input_cat.init(cat_charge,ion_multi,runtype,QMmethod,base_basis);
				mpc_input_cat.write_file(geo_file.molecule,oname);
				mpc_input_an.init(an_charge,ion_multi,runtype,QMmethod,base_basis);
				mpc_input_an.write_file(geo_file.molecule,oname);				
			break;
			case package::ORCA:
				orcaInput orc_input();
				orcaInput orc_input_cat();
				orcaInput orc_input_an();
				orc_input.init(geo_file.molecule,base_charge,base_multi,QMmethod,runtype,1,base_basis);
				orc_input.write_input_file(oname);
				orc_input.init(geo_file.molecule,cat_charge,ion_multi,QMmethod,runtype,1,base_basis);
				orc_input.write_input_file(oname);
				orc_input.init(geo_file.molecule,an_charge,ion_multi,QMmethod,runtype,1,base_basis);
				orc_input.write_input_file(oname);
			break;
		}
	}										 
}
//////////////////////////////////////////////////////////////