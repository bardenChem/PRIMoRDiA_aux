//residue.h

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


#ifndef RESIDUE_H
#define RESIDUE_H

#include <string>
#include <vector>

///////////////////////////////////////////////////////////////////////
enum AAres{
	OTH=-1,
	ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE,
	PRO, SER, THR, TRP, TYR, VAL	
};
/*********************************************************************/
enum DNAres{
	OTHd=-1,
	DC, DG, DT, DA
};
/*********************************************************************/
enum IONres{
	OTHi=-1,
	K,Ca,Cl,Mg,SO4
};
/*********************************************************************/
enum res_type {
	UNK, WAT, ION, AA, DNA, RNA, LIG
};
/*********************************************************************/

//------------------------------------------------
class pdbAtom; // foward declaration

class residue{
	public:
		res_type type;
		AAres AAname;
		IONres IONname;
		DNAres DNAname;
		bool ligand;
		bool terminal;
		bool first;
		int fCharge;
		int pdb_index;
		unsigned int nHydrogens;
		unsigned int nAtoms;
		std::vector<pdbAtom> r_atoms;		
		residue();
		residue(std::vector<pdbAtom> resAtoms);
		~residue();
		residue(const residue& rhs);
		residue& operator=(const residue& rhs);
		residue( residue&& rhs) noexcept;
		residue& operator=( residue&& rhs) noexcept;
		res_type get_type();
		void set_charge();
};

#endif