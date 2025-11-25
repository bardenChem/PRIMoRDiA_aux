//interface.cpp

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

#include <iostream>
#include <string>
#include <vector>

#include "../include/global.h"
#include "../include/interface.h"
#include "../include/traj_analysis.h"
#include "../include/read_traj.h"
#include "../include/QCPinput.h"
#include "../include/PDB.h"
#include "../include/GRO.h"

#include "../include/unit_test.h"

#include <experimental/filesystem>

using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::to_string;
using std::stoi;

namespace fs = std::experimental::filesystem;
//=====================================================================
interface::interface(){}
/**********************************************************************/
interface::~interface(){}
/**********************************************************************/
interface::interface(int argc, char* argv[]):
	m_argc(argc)							{
	
	for(int i =0;i<m_argc;i++){
		m_argv.emplace_back(argv[i]);
	}	
	cout << "Starting PRIMoRDIA Auxialliary software! "	<< endl;
}
/**********************************************************************/
void interface::run(){
	
	/***************************************************************/
	//Trajectory analysis with mdtraj
	if ( m_argv[1] == "-rmsd_rg"){
		cout << "Selected: RMSD and RG calculations using MDTraj" << endl;
		traj_an Trajectory( m_argv[2] );
		Trajectory.mdtraj_geo();
	}
	/***************************************************************/
	//Distance analysis from PDB traj
	else if ( m_argv[1] == "-atom_BID" ){
		vector<int> ats; 
		for(int i=3;i<m_argc;i++){
			ats.push_back( std::stoi(m_argv[i]) );
		}
		traj_an Trajectory(m_argv[2],ats);
		Trajectory.calc_distances( m_argv[2].c_str() );
	}
	/****************************************************************/
	//Extract frame from a trajectory
	else if ( m_argv[1] == "-ext_frame_dcd"){
		int interval = stoi(m_argv[4]);
		ReadTraj dcd( m_argv[2].c_str(), m_argv[3].c_str() );
		dcd.parse();
		PDB _Samp = dcd.sample( interval );
		_Samp.write_pdb( "sampled.pdb" );
	}
	/****************************************************************/
	else if ( m_argv[1] == "-ext_frame" ){
		int arg3 = stoi(m_argv[3]);
		int arg4 = stoi(m_argv[4]);
		traj_an Trajectory;
		if ( arg4 > 0 ){
			Trajectory.extract_frames(m_argv[2].c_str(),arg3,arg4);
		}else{
			Trajectory.extract_frame(m_argv[2].c_str(),arg3);
		}
	}
	else if ( m_argv[1] == "-spherical_prune"){
		double arg4 = stod(m_argv[3]); // tamanho do raio
		PDB( m_argv[2].c_str() );
	}
	else if ( m_argv[1] == "-analysis_mol_from_center"){
		ReadTraj TRJ( m_argv[2].c_str(), m_argv[3].c_str() );  //file and topology 
		TRJ.parse();
		unsigned res_indx    = stoi(m_argv[4]);
		std::string mol_name = m_argv[5];
		double radius        = stod(m_argv[6]);
		double prune         = 0.0;
		if ( m_argc > 6 ) { prune = stod(m_argv[7]); }
		TRJ.analysis_ac_from_molecules(res_indx,mol_name,radius,prune);
	} 
	/*****************************************************************/
	//creat input from QCP calculations
	else if ( m_argv[1] == "-QCP_inp" ){
		this->input_QM();
	}
	/*****************************************************************/
	else if ( m_argv[1] == "-MDprep" ){}
	/*****************************************************************/
	else if ( m_argv[1] == "-check_QCP" ){}
	/*****************************************************************/
	else if ( m_argv[1] == "-test" ) {
		//this->UnitTest();
		this->test();
	}
	else if ( m_argv[1] == "-rH2O_traj" ){
		if ( m_argv[2] == "pdbs_folder" ){
			fs::path c_path = fs::current_path();
			std::vector<string> fnames; 
		
			for ( const auto & entry : fs::directory_iterator(c_path) ){
				string tmp_name = entry.path();
				if ( check_file_ext( ".pdb",tmp_name.c_str() ) ){
					fnames.push_back( tmp_name );
				}
			}
			PDB _pdb;
			_pdb.cat_pdbs(fnames);
			std::vector<std::string> pars = { m_argv[3], m_argv[4] };
			_pdb.iterate_models("remove_waters",pars);
			_pdb.split_models_in_files();
		}
	}
	/****************************************************************/
	else{
		cout << "None valid flag was selected!\n Exiting PRIMoRDiA AUX\n" << endl;
		exit(-1);
	}
}
//======================================================================
//SETTING QM INPUT 
/**********************************************************************/
void interface::input_QM(){
	
	cout << "Quantum Chemistry Input Option preparation Selected!" << endl;
	
	package prog;
	if		( m_argv[2] == "gamess" ){
		prog = package::GAMESS;	
	}else if( m_argv[2] == "orca" ){
		prog = package::ORCA;
		this->set_nprocs();
	}else if( m_argv[2] == "mopac" ){
		prog = package::MOPAC;
	}
	
	string geo_ext	= ".xyz";
	string QM_		= "am1";
	string basis_	= "am1";
	string mode		= "Normal"; // or FD
	int bcharge		= 0;
	int bmulti		= 1;
	string rtype	= "Energy";
	int charge_diff	= 1;
	bool markCharge = false; 
	int  Mcharge    = 0;
	string _residueM= "UNK";
	string _topology= "";
	
	for( int i =0; i<m_argc; i++ ){ 
		if ( m_argv[i] == "-gExt" ){
			geo_ext = m_argv[i+1];
		}
		else if ( m_argv[i] == "-QMm" )		   { QM_         = m_argv[i+1];			}
		else if ( m_argv[i] == "-basis" )	   { basis_      = m_argv[i+1];			}
		else if ( m_argv[i] == "-FD" )		   { mode        = "FD";				}
		else if ( m_argv[i] == "-chg" )		   { bcharge     = stoi( m_argv[i+1] );	}
		else if ( m_argv[i] == "-bmulti" )	   { bmulti      = stoi( m_argv[i+1] );	}
		else if ( m_argv[i] == "-runOP" )	   { rtype       = m_argv[i+1];			}
		else if ( m_argv[i] == "-chgDiff" )	   { charge_diff = stoi( m_argv[i+1] );	}
		else if ( m_argv[i] == "-markCharge" ) { markCharge  = true; Mcharge =  stoi(m_argv[i+1]) ; }
		else if ( m_argv[i] == "-markResidue" ){ _residueM   = m_argv[i+1] ; }
		else if ( m_argv[i] == "-topology" )   { _topology   = m_argv[i+1] ; }
	}
	
	//this->print_options();
	QCPinput Input(geo_ext,rtype,basis_,QM_);
	if ( markCharge ) {
		Input.make_input_mopac_marked(prog,_topology.c_str(),Mcharge,_residueM,bmulti,bcharge);
	}
	else if ( mode == "FD" ){
		Input.make_input_from_folder_FD(prog,bmulti,bcharge,charge_diff);
	}else{
		Input.make_input_from_folder(prog,bmulti,bcharge);
	}
	
}
//=================================================================================

//==================================================================================
// AUXILIARY INTERFACE FUNCTIONS
/**********************************************************************/
void interface::print_options(){
	for(unsigned int i=0; i<m_argc; i++ ){
		cout << "#" << to_string(i+1) << " " << m_argv[i] << endl;
	}
}
/**********************************************************************/
void interface::help(){
	cout<< "PRIMoRDiA Auxilliary Help info!\n\n"
		<<"-h,--help: Print this option\n " 
		<<"-QCP_inp: Prepare input for quantum chemistry package\n"
		<<"-MDprep: Prepare files for molecular dynamics simulations\n"
		<<"-check_QCP_: Check the quantum chemistry package output files in the folder\n"
		<<"-extract_frame: Extract one frame from PDB trajectory file \n"
		<<"-extract_frames: Extract several frames from PDB trajectory file\n"
		<<"-traj_geo: Geometric analysis of a molecular dynamics trajectory\n"
		<<"-dist_an: Distance Analysis of pair of atoms of a molecular dynamics trajectory\n"
		<<"-------------------------------------------------------------------------------\n";
}
/**********************************************************************/
void interface::test(){
	/*
	std::string test_path = "/home/igorchem/primordia-code/PRIMoRDiA_aux/test_data/";
	std::string Top = test_path + "tim.pdb";
	std::string DCD = test_path + "scan1d.dcd";
	
	ReadTraj TestTraj( DCD.c_str(), Top.c_str() );
	TestTraj.parse();
	 
	vector<int> ats = {4672, 5325, 5326, 5357 };
	
	std::string fn = "md.pdb";
	traj_an Trajectory(fn,ats);
	Trajectory.calc_distances( fn.c_str() );
	*/
	
	const char* xtc_file = "/home/igorchem/CCDIR/PETROBRAS-F3/AC_DM/AC/traj.xtc";
	const char* gro_file = "/home/igorchem/CCDIR/PETROBRAS-F3/AC_DM/AC/traj.gro";
	
	ReadTraj _xtc_test(xtc_file,gro_file);
	_xtc_test.parse();
	_xtc_test.analysis_ac_from_molecules(466,"CO2",4.2,30.0);
	
}
/**********************************************************************/
void interface::UnitTest(){
	UnitTests tests;
	tests.run_unit_tests();
}
/**********************************************************************/
void interface::set_nprocs(){
	for(unsigned int i=0;i<m_argc;i++){
		if ( m_argv[i] == "-NP" ){
			m_NumOfProcess = stoi( m_argv[i+1]  );
		}
	}
}
////////////////////////////////////////////////////////////////////////