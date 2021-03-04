//gamessInput.h

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

#ifndef GAMESS_INPUT
#define GAMESS_INPUT

#include <vector>
#include <string>
#include <iostream>


/*********************************************************************/
enum GMS_Group{	
	OTHER=0, CTRL, SYS, SCF, BASIS,	DFTB, OPT, PCM,
	GUESS,	DATA, VEC, ELSDEN, GRID, ELSPOT, 	
	NumOfGMSgroups
};
//----------------------------
enum GMS_TheoryLevel { 	
	HF, DFT, DFTB2, DFTB3,
	FMO , P_HF, Semi,						
	NumOfTheory
};
//----------------------------
enum SCF_TYPE { 
	Default, UHF, ROHF,				
	NumOfSCFoptions
};

//----------------------------
enum GMS_Run_Type { 
	InvalidRunType	,
	Energy			,
	GradientRun		,
	HessianRun		,
	OptimizeRun		,
	SadPointRun		,
	IRCRun			,
	SurfaceRun		,
	FMOEnergy		,
	FMOoptimize		,
	FMOdynamics		,					
	NumOfRunTypeOptions					
};
//----------------------------
enum GMS_Conv_OPT { 
	DIIS  , SOSCF, MIXED ,
	
	NumOfConvOpt
};
//----------------------------		
enum GMS_Conv_OPT2{ 
	None  ,
	SHIFT ,
	DAMP  ,
	RSTRCT,
	
	NumOfConvOpt2
};
//----------------------------
enum GMS_BasisSet {
	
};
/*********************************************************************/
class GMS_basis{
	public:
		GMS_BasisSet type;
		std::string name;
		unsigned int ngauss;
		unsigned int npfunc;
		unsigned int ndfunc;
		unsigned int ndiffuseS;
		unsigned int ndiffuseP;
		std::string gaussBasis;
		std::string path_to_dftb;
		GMS_basis();
		GMS_basis(int typ)
		~GMS_basis();
		GMS_basis(const GMS_basis& rhs);
		GMS_basis& operator=(const GMS_basis& rhs);
		void set_type( int typ );
};
/*********************************************************************/
class FMO_fragment{
	
};
/*********************************************************************/
class FMO_options{
	
};
/*********************************************************************/
class gms_group{
	public:
		GMS_Group group;
		std::vector<std::string> inp_text;
		std::string grpName;
		gms_group();
		gms_group( int grp_type );
		~gms_group();
		gms_group(const gms_group& rhs);
		gms_group& operator=(const gms_group& rhs);
		gms_group(gms_group&& rhs) noexcept;
		gms_group& operator=(gms_group&& rhs) noexcept;
		friend std::ostream& operator<<(std::ostream& out, gms_group& grp);
};
/*********************************************************************/
class gms_input{
	public:
		//basic data
		unsigned int multi;
		int charge;		
		std::string QM_method;
		GMS_basis gbasis;
		std::string solvent ;
		GMS_Run_Type RunType;
		//control data
		unsigned int nprint;
		unsigned int maxit;
		unsigned int mwords;
		unsigned int npunch;
		unsigned int nsteps;
		std::string solvent;
		GMS_Conv_OPT copt;
		GMS_Conv_OPT2 copt2;
		GMS_TheoryLevel QMlevel;
		std::vector<double> tolerances;
		//data to write
		std::vector<gms_group> groups;
		std::ostream fl_data;		
		
		//constructors and methods
		gms_input();
		~gms_input();
		gms_input(const gms_input& rhs) = delete;
		gms_input& operator=(gms_input& rhs) = delete;
		void init(int chg, unsigned int mpcty, std::string method, std::string basis);
		void load_default_options();
		void load_molecule_info( system& molecule );
		void read_input(const char* file_name);
		void restart_input(const char* inp_name, const char* vec_data);
		void write_input(std::string out_name);
		void clear_directory();

		
};

#endif