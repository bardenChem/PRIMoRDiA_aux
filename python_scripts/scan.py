#Python Class for automatization of relaxed scan


from __future__ import print_function

import os 
import numpy as np
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
		self.forceC				= 4000.0
		self.mass_constraint	= True
		self.multiple_d			= [False,False]
		self.nsteps	   			= [ 1, 1 ]
		self.maxIt				= 30
		self.rmsGT				= 0.1
		self.minimumAlg			= "cgrad" 
		self.sigma_a1_a3		= [0,0]
		self.sigma_a3_a1		= [0,0]
		self.real_distance_1	= []
		self.real_distance_2	= []
		
		self.LOG  = open(self.base_name+"_scan.log","w")
		self.text = "Information of pDynamo SCAN \n"
	
		self.folder = os.path.join( self.base_name+"pDynamo_pj", "scan_trajectory" )
		self.folder_trj = os.path.join( self.base_name+"pDynamo_pj", "trj" )
		if not os.path.exists(self.folder):
			os.makedirs( self.folder )
		if not os.path.exists(self.folder_trj):
			os.makedirs( self.folder_trj )
		
	#.===========================================================
	def set_rcs(self,na,atoms,dincre):
		
		# . Set reaction coordinates 
		
		ndim = self.ndim				
		self.ndim += 1	
			
		self.atoms.append(atoms)
		self.DINCREMENT[ndim] = dincre
		self.sigma_a1_a3[ndim] = 1 
		self.sigma_a3_a1[ndim] = -1		
		
		if self.mass_constraint:
			atomic_n1 = self.molecule.atoms.items[ self.atoms[ndim][0] ].atomicNumber
			atomic_n3 = self.molecule.atoms.items[ self.atoms[ndim][2] ].atomicNumber
			mass_a1 = atomic_mass[ atomic_n1 ]
			mass_a3 = atomic_mass[ atomic_n3 ]
			self.sigma_a1_a3[ndim] = mass_a1 /(mass_a1+mass_a3)
			self.sigma_a3_a1[ndim] = mass_a3 /(mass_a1+mass_a3)
			self.sigma_a3_a1[ndim] = self.sigma_a3_a1[ndim]*-1
		
		if na == 3:			
			self.multiple_d[ndim] = True
			dist_a1_a2 = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )
			dist_a2_a3 = self.molecule.coordinates3.Distance( self.atoms[ndim][1], self.atoms[ndim][2] )
			if self.mass_constraint:
				self.DMINIMUM[ndim] = (self.sigma_a1_a3[ndim] * dist_a1_a2) - ( self.sigma_a3_a1[ndim] * dist_a2_a3*-1)
			else:
				self.DMINIMUM[ndim] = dist_a1_a2 - dist_a2_a3
		elif na == 2:
			self.DMINIMUM[ndim] = self.molecule.coordinates3.Distance( self.atoms[ndim][0], self.atoms[ndim][1] )
		else:
			print( "Neither '3' or '2' are the size of atoms container") 
	
	#.==============================================================
	
	def run_1D_Scan(self,step_s):
		self.nsteps[0] = step_s
		
		distance = 0.0
		constraints = SoftConstraintContainer( )
		self.molecule.DefineSoftConstraints( constraints )
			
		self.text += "Step RC Energy\n"	
		
		# .---------------------------------------------------------
		if self.multiple_d[0]:
			ATOM1 = self.atoms[0][0]
			ATOM2 = self.atoms[0][1]
			ATOM3 = self.atoms[0][2]
			SIGMA13 = self.sigma_a1_a3[0]
			SIGMA31 = self.sigma_a3_a1[0]
			
			# .-----------------------------------------------------	
			for i in range (self.nsteps[0]):				
				distance = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]								
				scModel  = SoftConstraintEnergyModelHarmonic ( distance, self.forceC )
				constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
				constraints["ReactionCoord"] = constraint
				
				out_name = os.path.join( self.folder,"frame{}.pdb".format(i) )
				Geometry_Optimization(self.molecule, out_name, self.rmsGT, traj=False )
				self.energies.append( self.molecule.Energy() )
				self.rcs.append(  self.molecule.coordinates3.Distance( ATOM1, ATOM2 ) - self.molecule.coordinates3.Distance( ATOM2, ATOM3 ) ) 
				self.text += "{} {} {}\n".format(i,self.rcs[i],self.energies[i])
				Pickle( os.path.join( self.folder_trj , "frame{}.pkl".format(i) ), self.molecule.coordinates3 ) 
				
		# .----------------------------------------------------------------------------
		else:
			ATOM1 = self.atoms[0][0]
			ATOM2 = self.atoms[0][1]
			for i in range ( self.nsteps[0] ):
				distance = self.DINCREMENT * float(i) + self.DMINIMUM
				scModel  = SoftConstraintEnergyModelHarmonic ( distance, self.forceC )
				constraint = SoftConstraintDistance ( int(ATOM1), int(ATOM2), scModel )
				constraints["ReactionCoord"] = constraint	
				
				out_name = os.path.join( self.folder,"frame{}.pdb".format(i) )
				Geometry_Optimization(self.molecule, out_name, self.rmsGT, traj=False )
				self.energies.append( self.molecule.Energy() )
				self.real_distance_1.append( self.molecule.coordinates3.Distance( ATOM1, ATOM2 ) )
				self.text += "{} {} {}\n".format(i,self.real_distance_1[i], self.energies[i] - self.energies[0])
				Pickle( os.path.join( self.folder_trj , "frame{}.pkl".format(i) ), self.molecule.coordinates3 ) 
				
		# .------------------------------------------------------------------------------

		DCDTrajectory_FromSystemGeometryTrajectory( os.path.join(self.folder, "scan1d.dcd"), self.folder_trj , self.molecule )
		self.molecule.DefineSoftConstraints ( None )
			
	#----------------------------------------------------
	
	def run_2D_Scan(self,step_sx,step_sy):
		self.nsteps[0] = step_sx
		self.nsteps[1] = step_sy
		M = self.nsteps[0]
		N = self.nsteps[1]
		
		total_steps = self.nsteps[0]*self.nsteps[1]
		
		distance_1		= 0.0
		distance_2		= 0.0
		self.real_distance_1	= [0] * total_steps
		self.real_distance_2	= [0] * total_steps
		self.energies			= [0] * total_steps
		
		ATOM1 = self.atoms[0][0]
		ATOM2 = self.atoms[0][1]
		ATOM3 = self.atoms[0][2]
		ATOM4 = self.atoms[1][0]
		ATOM5 = self.atoms[1][1]
		ATOM6 = self.atoms[1][2]
		SIGMA13 = self.sigma_a1_a3[0]
		SIGMA31 = self.sigma_a3_a1[0]
		SIGMA46 = self.sigma_a1_a3[1]
		SIGMA64 = self.sigma_a3_a1[1]
		
		constraints = SoftConstraintContainer( )
		self.molecule.DefineSoftConstraints( constraints )
		
		for i in range (self.nsteps[0] ):
			#.----
				
			distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
			scModel  = SoftConstraintEnergyModelHarmonic ( distance_1, self.forceC )
			constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
			constraints["ReactionCoord"] = constraint
			for j in range ( self.nsteps[1] ):		
			#.----
				distance_2	= self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
				scModel		= SoftConstraintEnergyModelHarmonic ( distance_2, self.forceC )
				constraint2	= SoftConstraintMultipleDistance ( [ [ATOM5, ATOM4, SIGMA46], [ATOM5, ATOM6, SIGMA64] ], scModel )
				constraints["ReactionCoord2"] = constraint2
						
				out_name = os.path.join( self.folder,"frame{}_{}.pdb".format(i,j) )
				Geometry_Optimization( self.molecule, out_name, self.rmsGT, traj=False )
				self.energies[ i*M +j ]			=	self.molecule.Energy() 
				self.real_distance_1[ i*M +j ]	=	self.molecule.coordinates3.Distance( ATOM1, ATOM2 ) - self.molecule.coordinates3.Distance( ATOM2, ATOM3 )   
				self.real_distance_2[ i*M +j ]	=	self.molecule.coordinates3.Distance( ATOM4, ATOM5 ) - self.molecule.coordinates3.Distance( ATOM5, ATOM6 ) 
				Pickle( os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j ) ), self.molecule.coordinates3 ) 
				
		# .----------------------------------------------------
				
		self.molecule.DefineSoftConstraints ( None )
		
		count = 0 
		for i in range(self.nsteps[0]):
			for j in range(self.nsteps[1]):
				self.text += "{} {} {} {} {}\n".format( i,j,real_distance_1[i*M+j],real_distance_2[i*M+j],energies[i*M+j]-energies[0] )
		

	#.----------------------------------------------
	
	def run_2D_Scan_parallel(self,step_sx,step_sy):
		self.nsteps[0] = step_sx
		self.nsteps[1] = step_sy
		
		M = self.nsteps[0]
		N = self.nsteps[1]
		
		total_steps = self.nsteps[0]*self.nsteps[1]
		
		ATOM1 = self.atoms[0][0]
		ATOM2 = self.atoms[0][1]
		ATOM3 = self.atoms[0][2]
		ATOM4 = self.atoms[1][0]
		ATOM5 = self.atoms[1][1]
		ATOM6 = self.atoms[1][2]
		
		SIGMA13 = self.sigma_a1_a3[0]
		SIGMA31 = self.sigma_a3_a1[0]
		SIGMA46 = self.sigma_a1_a3[1]
		SIGMA64 = self.sigma_a3_a1[1]
		
		distance_1			= 0.0
		distance_2			= 0.0
		self.real_distance_1= [0] * total_steps
		self.real_distance_2= [0] * total_steps
		self.energies		= [0] * total_steps
		
		constraints = SoftConstraintContainer( )
		self.molecule.DefineSoftConstraints( constraints )
		self.text += "Step1 Step2 RC1 RC2 Energy\n"	
				
		with pymp.Parallel(self.nprocs) as p:
			for i in p.range ( 0, M ):	
				#.----				
				distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
				scModel  = SoftConstraintEnergyModelHarmonic ( distance_1, self.forceC )
				constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
				constraints["ReactionCoord"] = constraint
				
				j=0
				distance_2	= self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
				scModel		= SoftConstraintEnergyModelHarmonic ( distance_2, self.forceC )
				constraint2	= SoftConstraintMultipleDistance ( [ [ATOM5, ATOM4, SIGMA46], [ATOM5, ATOM6, SIGMA64] ], scModel )
				constraints["ReactionCoord2"] = constraint2		
									
				if i > 0:
					init_coordinate_file = os.path.join( self.folder_trj,"frame{}_{}.pkl".format(0,0) )	
					self.molecule.coordinates3 = Unpickle(init_coordinate_file)
				
				coordinate_file = os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j ) )
				out_name = os.path.join( self.folder,"frame{}_{}.pdb".format(i,j) )	
				Geometry_Optimization( self.molecule, out_name, self.rmsGT, traj=False )
				Pickle( coordinate_file, self.molecule.coordinates3 ) 
		
		with pymp.Parallel(self.nprocs) as p:
			for i in p.range ( 0, M ):
				distance_1 = self.DINCREMENT[0] * float(i) + self.DMINIMUM[0]
				scModel  = SoftConstraintEnergyModelHarmonic ( distance_1, self.forceC )
				constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
				constraints["ReactionCoord"] = constraint
				
				for j in range( 1,N ):
					distance_2	= self.DINCREMENT[1] * float(j) + self.DMINIMUM[1]
					scModel		= SoftConstraintEnergyModelHarmonic ( distance_2, self.forceC )
					constraint2	= SoftConstraintMultipleDistance ( [ [ATOM5, ATOM4, SIGMA46], [ATOM5, ATOM6, SIGMA64] ], scModel )
					constraints["ReactionCoord2"] = constraint2		
						
					out_name = os.path.join( self.folder,"frame{}_{}.pdb".format(i,j) )
					init_coordinate_file = ""
					if  j==0:
						init_coordinate_file = os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, 0 ) )
					elif j>0:
						init_coordinate_file = os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j-1 ) )
			
					coordinate_file = os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j ) )
					self.molecule.coordinates3 = Unpickle(init_coordinate_file)
					Geometry_Optimization( self.molecule, out_name, self.rmsGT, traj=False )
					Pickle( coordinate_file, self.molecule.coordinates3 ) 
		
		#.-----------------
		self.molecule.DefineSoftConstraints ( None )
		
		for i in range(self.nsteps[0]):
			for j in range(self.nsteps[1]):
				molecule = self.molecule
				molecule.coordinates3 = Unpickle( os.path.join( self.folder_trj , "frame{}_{}.pkl".format( i, j ) ) )
				self.energies[i*M+j] = molecule.Energy()
				self.real_distance_1[ i*M +j ]	=	molecule.coordinates3.Distance( ATOM1, ATOM2 ) - molecule.coordinates3.Distance( ATOM2, ATOM3 )   
				self.real_distance_2[ i*M +j ]	=	molecule.coordinates3.Distance( ATOM4, ATOM5 ) - molecule.coordinates3.Distance( ATOM5, ATOM6 )
				self.text += "{} {} {} {} {}\n".format( i,j,self.real_distance_1[i*M+j],self.real_distance_2[i*M+j],self.energies[i*M+j]-self.energies[0] )
	
	#-------------------------------------------------------------------
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
		
	#-------------------------------------------------------------------
	def chang_nprocs(self,nc):
		self.nprocs = nc
		
	#-------------------------------------------------------------------
	def write_log(self):		
		self.LOG.write(self.text)
		self.LOG.close()
		
	#-------------------------------------------------------------------
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
	'''
	filename = os.path.join("pDynamo_tests","TIM_qcmm_opt.pkl")
	tim_qcmm= Unpickle( filename )
	
	C02 	= AtomSelection.FromAtomPattern(tim_qcmm,"*:LIG.248:C02")
	proton	= AtomSelection.FromAtomPattern(tim_qcmm,"*:LIG.248:H02")
	oxygen	= AtomSelection.FromAtomPattern(tim_qcmm,"*:GLU.164:OE2")
	
	NE2		= AtomSelection.FromAtomPattern(tim_qcmm,"*:HIE.94:NE2")
	HE2		= AtomSelection.FromAtomPattern(tim_qcmm,"*:HIE.94:HE2")
	O06		= AtomSelection.FromAtomPattern(tim_qcmm,"*:LIG.248:O06")
	
	atom1	= C02.selection.pop()
	atom2	= proton.selection.pop()
	atom3	= oxygen.selection.pop()
	
	atom4	= NE2.selection.pop()
	atom5	= HE2.selection.pop()
	atom6	= O06.selection.pop()
	
	rc1_atoms	= [atom1,atom2,atom3]
	rc2_atoms	= [atom4,atom5,atom6]

	scan1d		= SCAN(tim_qcmm,"TIM1D")
	scan1d.set_rcs(3,rc1_atoms,0.09)
	scan1d.Print()
	scan1d.run_1D_Scan(12)
	scan1d.write_log()
	
	scan2d = SCAN(tim_qcmm,"TIM2D_seq")
	scan2d.set_rcs(3,rc1_atoms,0.1)
	scan2d.set_rcs(3,rc2_atoms,0.1)
	scan2d.run_2D_Scan(15,15)
	scan2d.write_log()
	
	scan2d_p = SCAN(tim_qcmm,"TIM2D_par")
	scan2d_p.nprocs = 8
	scan2d_p.set_rcs(3,rc1_atoms,0.07)
	scan2d_p.set_rcs(3,rc2_atoms,0.07)
	scan2d_p.run_2D_Scan_parallel(20,20)
	scan2d_p.write_log()
	'''
	
