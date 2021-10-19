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

#include "../include/global.h"
#include "../include/unit_test.h"
#include "../include/atom.h"
#include "../include/molecule.h"
#include "../include/pdbAtom.h"
#include "../include/residue.h"
#include "../include/pdbModel.h"
#include "../include/PDB.h"
#include "../include/read_traj.h"

/*********************************************************************/
UnitTests::UnitTests(){}
/*********************************************************************/
UnitTests::~UnitTests(){}
/*********************************************************************/
void UnitTests::run_unit_tests(){
	ut_log.open();
	//UnitTest_atom();
	//UnitTest_molecule();
	//UnitTest_pdbAtom();
	//UnitTest_residue();
	//UnitTest_pdbModel();
	//UnitTest_PDB();
	UnitTest_ReadTraj();
}
///////////////////////////////////////////////////////////////////////