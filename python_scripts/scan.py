#Python Class for automatization of relaxed scan


from __future__ import print_function

import os 
from pDynamo_generic import *
import pymp

#=======================================================================
class SCAN:
		
	# . Class to automate relaxed scan for pDynamo package
		
	def __init__(self,molecule,name):
		
		# . Default Constructor
		
		self.base_name 			= name
		self.molecule  			= molecule 
		self.ndim	   			= 0
		self.rcs       			= []
		self.atoms	   			= [] # indices of the atoms 
		self.nprocs    			= 1
		self.text	   			= " "
		self.energies  			= []
		self.DMINIMUM  			= [ 0, 0 ]
		self.DINCREMENT			= [ 0, 0 ]
		self.forceC				= 4000
		self.mass_constraint	= True
		self.multiple_d			= False
		self.nsteps	   			= [ 1, 1 ]
		self.maxIt				= 30
		self.rmsGT				= 0.1
		self.minimumAlg			= "cgrad" 
		
		self.LOG  = open(self.base_name+"_scan.log","w")
		self.text = "Information of pDynamo SCAN \n"
		
	#---------------------------------------------------
	def set_rcs(self,na,atoms,dincre):
		
		# . Set reaction coordinates 
		
		ndim = self.ndim				
		self.ndim += 1	
			
		self.atoms.append(atoms)
		self.DINCREMENT[ndim] = dincre
		sigma_a1_a3 = 1 
		sigma_a3_a1 = -1
		
		if self.mass_constraint:
			atomic_n1 = self.molecule.atoms.items[ self.atoms[ndim][0] ].atomicNumber
			atomic_n3 = self.molecule.atoms.items[ self.atoms[ndim][2] ].atomicNumber
			mass_a1 = atomic_mass[ atomic_n1 ]
			mass_a3 = atomic_mass[ atomic_n3 ]
			sigma_a1_a3 = mass_a1 /(mass_a1+mass_a3)
			sigma_a3_a1 = mass_a3 /(mass_a1+mass_a3)
			sigma_a3_a1 = sigma_a3_a1*-1
		
		if na == 3:			
			self.multiple_d = True
			dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )
			dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[ndim][1], self.atoms[ndim][2] )
			if self.mass_constraint:
				self.DMINIMUM[ndim] = (sigma_a1_a3 * dist_a1_a2) - ( sigma_a3_a1 * dist_a2_a3*-1)
			else:
				self.DMINIMUM[ndim] = dist_a1_a2 - dist_a2_a3								
		elif na == 2:
			self.DMINIMUM[ndim] = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )			
		else:
			print( "Neither '3' or '2' are the size of atoms container") 
	
	#----------------------------------------------------
	def run_1D_Scan(self,nsx,step_s,):
		pass
	
	#----------------------------------------------------
	def run_2D_Scan(self,nsx,nsy,step_sx,step_sy):		
		pass
		
	#----------------------------------------------------
	def save_traj(self):
		pass
		
	#----------------------------------------------------	
	def change_parameters(self,rmsTopt=0.1	,
							   fc=4000		,
							   maxit=30		,
							   massC=False  ,
							   minAlg="lfbgs"):
		self.rmsGT 	= 0.1
		self.forceC	=fc
		self.maxIt	=maxit
		self.mass_constraint = massC
		self.minimumAlg = minAlg
		
	#----------------------------------------------------	
	def chang_nprocs(self,nc):
		self.nprocs = nc
		
	#----------------------------------------------------
	def write_log(self):		
		pass
		
	#----------------------------------------------------	
	def Print(self):
		
		atoms_r1 = [] 
		atoms_r2 = [] 
		
		lbl = self.molecule.atoms.items[ self.atoms[0][0] ].label
		atoms_r1.append(lbl)
		lbl = self.molecule.atoms.items[ self.atoms[0][1] ].label
		atoms_r1.append(lbl)
		
		if len(self.atoms[0]) == 3:
			lbl = self.molecule.atoms.items[ self.atoms[0][2] ].label
			atoms_r1.append(lbl)
			
		if self.ndim > 1:
			lbl = self.molecule.atoms.items[ self.atoms[1][0] ].label
			atoms_r2.append(lbl)
			lbl = self.molecule.atoms.items[ self.atoms[1][1] ].label
			atoms_r2.append(lbl)
			if len(self.atoms[1]) == 3:
				lbl = self.molecule.atoms.items[ self.atoms[1][2] ].label
				atoms_r2.append(lbl)
		
		
		#. Print basic parameters for checking 
		text = "Print Scan parameters\n"
		text+= "Num of RCs:  {}\n".format(self.ndim) 
		text+= "RC{}: dmin = {}; increment = {} \n".format(1,self.DMINIMUM[0],self.DINCREMENT[0])
		if len(self.atoms[0]) == 3:
			text+= "Atoms: {}--{}--{}\n".format(atoms_r1[0],atoms_r1[1],atoms_r1[2])
		else:
			text+= "Atoms: {}--{}\n".format(atoms_r1[0],atoms_r1[1])		
		print(text)


#=======================================================================

if __name__ == "__main__":	
	filename = os.path.join("pDynamo_tests","TIM_qcmm_opt.pkl")
	tim_qcmm= Unpickle( filename )
	C02 	= AtomSelection.FromAtomPattern(tim_qcmm,"*:LIG.248:C02")
	proton	= AtomSelection.FromAtomPattern(tim_qcmm,"*:LIG.248:H02")
	oxygen	= AtomSelection.FromAtomPattern(tim_qcmm,"*:GLU.164:OE2")
	
	atom1	= C02.selection.pop()
	atom2	= proton.selection.pop()
	atom3	= oxygen.selection.pop()
	
	rc1_atoms	= [atom1,atom2,atom3]	
	scan1d		= SCAN(tim_qcmm,"TIM")	
	scan1d.set_rcs(3,rc1_atoms,0.1)
	scan1d.Print()
