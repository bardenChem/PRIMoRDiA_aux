//MDprep.h

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

#ifndef MD_PREP_H
#define MD_PREP_H

#include <string>
#include <vector>

class pdbModel;

//-----------------
enum mdPKG{
	undef, AMBER,
	GROMACS, NumOfPKGs
}
//-----------------
enum mdTask{
	undefined, parametrization,	topology,
	minimization, equilibration, production,
	NumofTasks
}
/****************************************/
class MDprep{
	public:
		mdPKG program;
		mdTask job;
		pdbModel bio_sys_ref;
		pdbModel bio_sys_md;
		unsigned int net_chg;
		unsigned int ref_temp;
		unsigned int prod_time;
		std::string top_file;
		std::string crd_file;
		std::string sys_basename;
		std::vector<std::string> ligands_n;
		std::vector<int> liganf_chg;
		MDprep();
		MDprep(std::string bs_name,mdPKG pkg, mdTask Job);
		MDprep(const MDprep& rhs) = delete;
		MDprep& operator=(const MDprep& rhs) = delete ;
		~MDprep();
		void prepare_ligand(int nligand, bool ambTools);
		void prepare_complex(int nligand);
		void built_topology();
		void prepare_minimization();
		void prepare_equilibration();
		void prepare_production();
		void organize_dir_files();
}
/////////////////////////////////////////

#endif