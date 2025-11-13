//read_traj.h

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


#ifndef READ_TRAJ
#define READ_TRAJ

#include <vector>
#include <string>
#include <fstream>

#include "../include/geometry.h"
#include "../include/pdbModel.h"
/*****************************************************/
class PDB;
class pdbModel;
/*****************************************************/
enum traj_type{ pdb, dcd, xtc, gro, xyz, mol2, NONE };


class TrjCRD{
	public:
		std::vector< std::vector<double> >  xc;
		std::vector< std::vector<double> >  yc;
		std::vector< std::vector<double> >  zc;
		pdbModel topology;
		TrjCRD();
		TrjCRD( pdbModel _topology);
		TrjCRD( const TrjCRD& rhs );
		TrjCRD& operator=( const TrjCRD& rhs );
		TrjCRD(  TrjCRD&& rhs ) noexcept;
		TrjCRD& operator=( TrjCRD&& rhs ) noexcept;
		~TrjCRD();
		pdbModel create_pdb( unsigned _frame); 
		void init_frames(unsigned _nFrames);
};

/*****************************************************/
class ReadTraj{
	public:
		std::string traj_file;
		unsigned int natoms;
		unsigned int nframes;
		TrjCRD Positions;
		traj_type Type;
		ReadTraj();
		ReadTraj( const char* file_name );
		ReadTraj( const char* file_name, const char* topol_file );
		ReadTraj( const ReadTraj& rhs );
		ReadTraj& operator=( const ReadTraj& rhs );
		ReadTraj(  ReadTraj&& rhs ) noexcept;
		ReadTraj& operator=( ReadTraj&& rhs ) noexcept;
		void parse();
		PDB sample(unsigned interval);
		PDB sample_chunk(unsigned _init, unsigned _final );
		void analysis_ac_from_molecules(unsigned _res_indx, std::string _res_name);
		friend std::ostream& operator<<(std::ostream& out, const ReadTraj& obj);
		void print();
};
/*****************************************************/
void UnitTest_positions3D();
void UnitTest_ReadTraj();
/*****************************************************/
#endif