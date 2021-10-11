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
#include <iostream>

#include "../include/global.h"
#include "../include/molecule.h" 
#include "../include/geometry.h"
#include "../include/gamessInput.h"
#include "../include/orcaInput.h"
#include "../include/mopacInput.h"
#include "../include/QCPinput.h"

using std::string;
using std::vector;
using std::cout;
using std::endl;
using std::move;
namespace fs = std::experimental::filesystem;

/************************************************************/
QCPinput::QCPinput()		:
	program(package::MOPAC)	,
	base_multi(1)			,
	base_charge(0)			,
	QMmethod("HF")			,
	base_basis("3-21G")		,
	runtype("Energy")		,
	geo_ext("xyz")			{
		
}
/************************************************************/
QCPinput::QCPinput(	string _geo_ext	  ,
					string _runtype	  ,
					string _basis	  ,
					string _method   ):
	program(package::MOPAC)				,
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
void QCPinput::make_input_from_folder(	package QCP			,
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
			{
				gms_input g_input;
				g_input.init(base_charge,base_multi,runtype,QMmethod, base_basis);
				g_input.load_molecule_info(geo_file.Molecule);
				g_input.write_input(oname);
				break;			
			}
			case package::MOPAC:
			{
				mopac_input mpc_input;
				mpc_input.init(base_charge,base_multi,runtype,QMmethod,base_basis);
				mpc_input.write_file(geo_file.Molecule,oname);
				break;
			}
			case package::ORCA:
			{
				orcaInput orc_input(base_charge,base_multi,m_NumOfProcess,runtype);
				orc_input.write_inp(geo_file.Molecule,QMmethod,base_basis,oname);
			}
			break;
		}
	}
}
/************************************************************/								
void QCPinput::make_input_from_folder_FD(package QCP			,
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
		string oname_cat = oname + "_cat";
		string oname_an = oname + "_an";		
		switch( program ){
			case package::GAMESS:
			{
				gms_input g_input_neu;
				gms_input g_input_cat;
				gms_input g_input_an;
				g_input_neu.init(base_charge,base_multi,runtype,QMmethod, base_basis);
				g_input_neu.load_molecule_info(geo_file.Molecule);
				g_input_neu.write_input(oname);
				g_input_cat.init(cat_charge,ion_multi,runtype,QMmethod, base_basis);
				g_input_cat.load_molecule_info(geo_file.Molecule);
				g_input_cat.write_input(oname_cat);
				g_input_an.init(an_charge,ion_multi,runtype,QMmethod, base_basis);
				g_input_an.load_molecule_info(geo_file.Molecule);
				g_input_an.write_input(oname_an);			
				break;			
			}
			case package::MOPAC:
			{
				mopac_input mpc_input_neu;
				mopac_input mpc_input_cat;
				mopac_input mpc_input_an;
				mpc_input_neu.init(base_charge,base_multi,runtype,QMmethod,base_basis);
				mpc_input_neu.write_file(geo_file.Molecule,oname);
				mpc_input_cat.init(cat_charge,ion_multi,runtype,QMmethod,base_basis);
				mpc_input_cat.write_file(geo_file.Molecule,oname);
				mpc_input_an.init(an_charge,ion_multi,runtype,QMmethod,base_basis);
				mpc_input_an.write_file(geo_file.Molecule,oname);				
				break;				
			}
			case package::ORCA:
			{
				//cout << an_charge << " " << cat_charge <<  " " << ion_multi << endl;
				orcaInput orc_input(base_charge,base_multi,m_NumOfProcess,runtype);
				orc_input.write_inp(geo_file.Molecule,QMmethod,base_basis,oname);
				orcaInput orc_input_cat(cat_charge,ion_multi,m_NumOfProcess,runtype);
				orc_input_cat.write_inp(geo_file.Molecule,QMmethod,base_basis,oname_cat);
				orcaInput orc_input_an(an_charge,ion_multi,m_NumOfProcess,runtype);
				orc_input_an.write_inp(geo_file.Molecule,QMmethod,base_basis,oname_an);	
				break;
			}
		}
	}
	this->make_sh(fnames);
}
/************************************************************/
void QCPinput::make_sh(vector<string> _fnames ){
	
	string run_txt_1 = " ";
	string run_txt_2 = " ";
	
	vector<string> inps( _fnames.size() );
	vector<string> outs( _fnames.size() );
	
	if ( program == package::GAMESS ){	
		for(int i=0;i<inps.size();i++){
			inps[i] = remove_extension( _fnames[i].c_str() );
			outs[i] = inps[i]+".log";
		}
	}else if ( program == package::ORCA ){
		for(int i=0;i<inps.size();i++){
			inps[i] = change_extension( _fnames[i].c_str(), ".inp");
			outs[i] = change_extension( inps[i].c_str(),".out" );
		}
	}else if (program == package::MOPAC ){
		for(int i=0;i<inps.size();i++){
			inps[i] = change_extension( _fnames[i].c_str(), ".mop");
			outs[i] = change_extension( inps[i].c_str(),".log" );
		}
	}
	
	if ( program == package::ORCA ){
		run_txt_1 = "orca ";
		run_txt_2 = " > ";
	}else if ( program == package::GAMESS ){
		run_txt_1 = "gms ";
		run_txt_2 = " 00 ";
		run_txt_2 += std::to_string(m_NumOfProcess);
		run_txt_2 += " > ";
	}else if ( program == package::MOPAC ){
		run_txt_1 = "/opt/mopac/MOPAC2016.exe ";
		run_txt_2 = " > ";
	}
	
	string lgc_ee = "&&\n";
	
	
	std::ofstream sh_file;
	sh_file.open("run_QCP_inp.sh");
	sh_file << "#!/bin/sh\n";
	for(int i=0;i<_fnames.size();i++){
		if ( i == _fnames.size()-1 ){
			lgc_ee = " ";
		}
		sh_file << run_txt_1
				<< inps[i]
				<< run_txt_2
				<< outs[i]
				<< lgc_ee;
	}
	cout << "writting file" << endl;
	sh_file.close();
}
/*************************************************************/
std::ostream& operator<<(std::ostream& out, const QCPinput& obj){
	
}
/*************************************************************/
void QCPinput::print(){
	
}
/*************************************************************/
void UnitTest_QCPinput(){
	
}
//////////////////////////////////////////////////////////////