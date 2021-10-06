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

#include <chemfiles.hpp>

class geometry;

enum traj_type{ pdb, dcd, xtc, gro, xyz, mol2 };

class positions3D{
	public:
		unsigned int natoms;
		unsigned int nframes;
		std::vector< std::vector<double> > xc;
		std::vector< std::vector<double> > yc;
		std::vector< std::vector<double> > zc;
		positions3D();
		~positions3D();
		positions3D( const positions3D& rhs );
		positions3D& operator=( const positions3D& rhs);
		positions3D( positions3D&& rhs ) noexcept;
		positions3D& operator=( positions3D&& rhs ) noexcept;
		double get_distance(unsigned atom1, unsigned atom2, unsigned frame);
		std::vector<double> get_distances(unsigned atom1, unsigned atom2);
		double avg_distance(unsigned atom1, unsigned atom2);
}

/*****************************************************/
class ReadTraj{
	public:
		std::string traj_file;
		unsigned int natoms;
		unsigned int nframes;
		geometry topology;
		positions3D coordinates;
		traj_type Type;
		ReadTraj();
		ReadTraj( const char* file_name );
		ReadTraj( const char* file_name, const char* topol_file );
		ReadTraj( const ReadTraj& rhs );
		ReadTraj& operator=( ReadTraj&& rhs ) noexcept;
		ReadTraj( const ReadTraj& rhs );
		ReadTraj& operator=( ReadTraj&& rhs ) noexcept;
		void parse();
};

#endif