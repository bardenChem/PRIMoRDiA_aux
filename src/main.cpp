//main.cpp

/*********************************************************************/
/* This source code file is part of LQQCMModels software project created 
 * by Igor Barden Grillo at Federal University of Paraíba. 
 * barden.igor@gmail.com ( Personal e-mail ) 
 * igor.grillo@acad.pucrs.br ( Academic e-mail )
 * quantum-chem.pro.br ( group site )
 * IgorChem ( Git Hub account )
 */ 

/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
/********************************************************/


/********************************************************
 * 
 * QM = Quantum Mechanics
 * MM = Molecular Mechanics
 * MD = Molecular Dynamics 
 *
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 +	***Functionalities of the PRIMoRDiA AUX 0.1***		+ 
 +														+
 +	1. Automatization of MDtraj to obtain RMSD and RG	+
 +	2. Atomic Distances Distribution Analysis from Traj	+
 +	3. Extract Frames from trajectories					+
 +	4. Check Quantum Chemistry simulation outputs files	+
 +	5. Make input for Quantum Chemistry Packages		+
 +	6. Make QM/MM inputs for energy refinement in Mopac +
 +	7. Automate the preparation of simple MD Simulations+
 +	8. PDB file manipulations							+
 +														+
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++

*********************************************************/

#include <iostream>
#include "../include/interface.h"

/***********************************/
int main(int argc, char** argv){
	interface program(argc,argv);
	program.run();
	return 0;
}
//======================================= 