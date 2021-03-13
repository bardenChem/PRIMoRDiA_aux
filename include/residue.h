

#ifndef RESIDUE_H
#define RESIDUE_H

#include <string>
#include <vector>

enum AAres{
	OTH=-1,
	ALA, ARG, ASN, ASP, CYS, GLN, GLU, GLY, HIS, ILE, LEU, LYS, MET, PHE,
	PRO, SER, THR, TRP, TYR, VAL	
};
/*********************************************************************/
enum res_type {
	UNK, WAT, HOH, ION, AA, DNA 
};

//------------------------------------------------
class pdbAtom; // foward declaration

class residue{
	public:
		std::string res3n;
		std::string res1n;
		res_type type;
		AAres AAname;
		bool ligand;
		bool terminal;
		bool first;
		int fCharge;
		int pdb_index;
		unsigned int nHydrogens;
		unsigned int nAtoms;
		std::vector<pdbAtom> r_atoms;		
		residue();
		residue(std::vector<pdbAtom> resAtoms, int resTyp, int resMon);
		~residue();
		residue(const residue& rhs);
		residue& operator=(const residue& rhs);
		residue( residue&& rhs) noexcept;
		residue& operator=( residue&& rhs) noexcept;
		void set_charge();
		bool is_ion();
};

#endif