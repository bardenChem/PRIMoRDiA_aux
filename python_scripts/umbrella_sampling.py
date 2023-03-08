#Python class to automate Umbrella Sampling in pDynamo

from __future__ import print_function

import os 
import numpy as np
from pDynamo_generic import *
import pymp

#=======================================================================
class MD:
	'''
	'''
	def __init__(self, molecule, output_path, constrained ):
		self.molecule		= molecule
		self.temperature	= 300.15
		self.time 		 	= 1 	# 1 ns
		self.time_step	 	= 0.001	# 0.1 ps 
		self.nsteps_equi	= 5000 # 5000 * 0.001 ps = 5ps
		self.nsteps_prod	= 20000	# 20000 * 0.001 ps = 20ps
		self.energies	 	= [0] * self.nsteps_prod
		self.base_name		= output_path
		self.out_traj		= False
		self.sample_traj	= 0
		self.trajectory		= None
		self.constrained	= constrained 
		self.samplen		= 0
		
		
		self.folder = os.path.join( self.base_name+"_trj" )
		self.folder_sc = os.path.join( self.base_name+"_softc_trj" )
		if not os.path.exists(self.folder):
			os.makedirs( self.folder )
		if not os.path.exists(self.folder_sc):
			if constrained:
				os.makedirs( self.folder_sc )
		
	
	#.------------------------------------------------------------------
	def Run(self):
		
		randomNumberGenerator  = RandomNumberGenerator ( )
		normalDeviateGenerator = NormalDeviateGenerator.WithRandomNumberGenerator ( randomNumberGenerator )
		
		randomNumberGenerator.SetSeed ( 292731  )
		log_equi_path = os.path.join(self.folder,"equi_log.txt")
		log_equi  = TextLogFileWriter( log_equi_path )
		
		self.molecule.energyModel.qcModel.converger.SetOptions( maximumSCFCycles = 10000 )
		self.molecule.energyModel.qcModel.converger.SetOptions( densityTolerance = 2.0e-6 )
		self.molecule.energyModel.qcModel.converger.SetOptions( energyTolerance = 2.0e-3 )
		
		#. Leap Frog equilibration
		LeapFrogDynamics_SystemGeometry(self.molecule									,
										logFrequency			= self.nsteps_equi/500	,
										log 					= log_equi				,
										normalDeviateGenerator	= normalDeviateGenerator,
										steps					= self.nsteps_equi		,
										temperature				= self.temperature		,
										temperatureControl		= True					,
										temperatureCoupling		= 0.1					,
										timeStep				= 0.001					)
									
        # . Data-collection.
        
		log_prod_path = os.path.join(self.folder,"prod_log.txt")
		log_prod  = TextLogFileWriter ( log_prod_path )
        
		if self.constrained:
			if self.samplen >0:
				print("sampling:{}".format( self.nsteps_prod/self.samplen ) )
				trajectoryS	= SystemSoftConstraintTrajectory ( self.folder_sc , self.molecule, mode = "w" )
				trajectory	= SystemGeometryTrajectory ( self.folder , self.molecule, mode = "w" )
				LeapFrogDynamics_SystemGeometry ( self.molecule						,
                                      logFrequency			= self.nsteps_prod/500	,
                                      log					= log_prod				,
                                      steps					= self.nsteps_prod		,
                                      temperature			= self.temperature		,
                                      temperatureControl	= True					,
                                      temperatureCoupling	= 0.1					,
                                      timeStep				= 0.001					,
                                      trajectories			= [ ( trajectoryS, 1 ), ( trajectory, self.samplen) ] )
               
			else:
				print("Not Sampling")
				trajectoryS	= SystemSoftConstraintTrajectory ( self.folder_sc , self.molecule, mode = "w" )
				LeapFrogDynamics_SystemGeometry ( self.molecule						,
                                      logFrequency			=	self.nsteps_prod/500,
                                      log					= log_prod				,
                                      steps					= self.nsteps_prod		,
                                      temperature			= self.temperature		,
                                      temperatureControl	= True					,
                                      temperatureCoupling	=	0.1					,
                                      timeStep				= 0.001					,
                                      trajectories			= [ ( trajectoryS, 1 ) ] )
		else:
			trajectory	= SystemGeometryTrajectory ( self.folder , self.molecule, mode = "w" )
			LeapFrogDynamics_SystemGeometry ( self.molecule						,
                                      logFrequency			= self.nsteps_prod/500	,
                                      log					= log_prod				,
                                      steps					= self.nsteps_prod		,
                                      temperature			= self.temperature		,
                                      temperatureControl	= True					,
                                      temperatureCoupling	= 0.1					,
                                      timeStep				= 0.001					,
                                      trajectories			= [ (trajectory, self.samplen) ] )
			

                                      
	#.------------------------------------------------------------------
	def change_parameters(self,param):
		self.temperature	= param[0]
		self.time 		 	= param[1]
		self.time_step	 	= param[2]
		self.nsteps_equi	= param[3]
		self.nsteps_prod	= param[4]
		

#=======================================================================
class umbrella_Sampling:
	
	def __init__(self,molecule,name,traj_folder):
		self.base_name			= name
		self.traj_folder		= traj_folder 
		self.molecule 			= molecule 
		self.ndim				= 0
		self.atoms				= [] # indices of the atoms 
		self.nprocs				= 1
		self.text				= " "
		self.energies			= []
		self.forceC				= 2500.0
		self.nsteps				= [ 1, 1 ]
		self.maxIt				= 30
		self.rmsGT				= 0.1
		self.minimumAlg			= "cgrad" 
		self.sigma_a1_a3		= [0,0]
		self.sigma_a3_a1		= [0,0]
		self.all_files			= []
		self.completed			= []

		self.folder = os.path.join( self.base_name )
		if not os.path.exists(self.folder):
			os.makedirs( self.folder )
	
	#.------------------------------------------------------------------
	def set_modes(self,atoms):
		# . Set reaction coordinates 
		
		ndim = self.ndim
		self.ndim += 1	
			
		self.atoms.append(atoms)
		self.sigma_a1_a3[ndim] = 1 
		self.sigma_a3_a1[ndim] = -1
		
		atomic_n1 = self.molecule.atoms.items[ self.atoms[ndim][0] ].atomicNumber
		atomic_n3 = self.molecule.atoms.items[ self.atoms[ndim][2] ].atomicNumber
		mass_a1 = atomic_mass[ atomic_n1 ]
		mass_a3 = atomic_mass[ atomic_n3 ]
		self.sigma_a1_a3[ndim] = mass_a1 /(mass_a1+mass_a3)
		self.sigma_a3_a1[ndim] = mass_a3 /(mass_a1+mass_a3)
		self.sigma_a3_a1[ndim] = self.sigma_a3_a1[ndim]*-1
		
	
	#.------------------------------------------------------------------
	def Sample1D(self):
		pat = os.path.join(self.traj_folder,"trj","")
		files_list = glob.glob(pat+"*.pkl")
		
		ATOM1 = self.atoms[0][0]
		ATOM2 = self.atoms[0][1]
		ATOM3 = self.atoms[0][2]
		SIGMA13 = self.sigma_a1_a3[0]
		SIGMA31 = self.sigma_a3_a1[0]
		self.all_files = files_list
		
		for i in range( len(files_list) ):
			molecule = self.molecule
			molecule.coordinates3 = Unpickle( files_list[i] )
				
			onstraints = SoftConstraintContainer( )
			molecule.DefineSoftConstraints( constraints )
			
			distance = molecule.coordinates3.Distance( ATOM1, ATOM2 ) - molecule.coordinates3.Distance( ATOM2, ATOM3 ) 
			scModel  = SoftConstraintEnergyModelHarmonic ( distance, self.forceC )
			constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
			constraints["ReactionCoord"] = constraint
			
			temp = files_list[i][:-4]
			temp = os.path.basename(temp)
			md_path	= os.path.join( self.folder, temp )
			md_run	= MD(molecule,md_path,True)
			md_run.sample_traj = 20
			md_run.Run()
			self.completed.append( files_list[i] )
	
	#.------------------------------------------------------------------
	def Sample1D_parallel(self):
		pat = os.path.join(self.traj_folder,"trj","")
		files_list = glob.glob(pat+"*.pkl")
		self.all_files = files_list
		
		ATOM1 = self.atoms[0][0]
		ATOM2 = self.atoms[0][1]
		ATOM3 = self.atoms[0][2]
		SIGMA13 = self.sigma_a1_a3[0]
		SIGMA31 = self.sigma_a3_a1[0]
		
		with pymp.Parallel(self.nprocs) as p:
			for i in p.range(0, len(files_list) ):
			
				molecule = self.molecule
				molecule.coordinates3 = Unpickle( files_list[i] )
				
				constraints = SoftConstraintContainer( )
				molecule.DefineSoftConstraints( constraints )
				
				distance = molecule.coordinates3.Distance( ATOM1, ATOM2 ) - molecule.coordinates3.Distance( ATOM2, ATOM3 ) 
				scModel  = SoftConstraintEnergyModelHarmonic ( distance, self.forceC )
				constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
				constraints["ReactionCoord"] = constraint
				
				temp = files_list[i][:-4]
				temp = os.path.basename(temp)
				md_path	= os.path.join( self.folder, temp )
				md_run	= MD(molecule,md_path,True)
				md_run.sample_traj = 20
				md_run.Run()
				self.completed.append( files_list[i] )
				
	#.------------------------------------------------------------------
	
	def Sample2D(self,x,y):
		pat = os.path.join(self.traj_folder,"trj","")
		files_list = glob.glob(pat+"*.pkl")
		files_list.sort()
		self.all_files = files_list
		
		ATOM1 = self.atoms[0][0]        #0
		ATOM2 = self.atoms[0][1]		#1
		ATOM3 = self.atoms[0][2]		#2
		ATOM4 = self.atoms[1][0]		#3
		ATOM5 = self.atoms[1][1]		#4
		ATOM6 = self.atoms[1][2]		#5
		SIGMA13 = self.sigma_a1_a3[0]	#6
		SIGMA31 = self.sigma_a3_a1[0]	#7
		SIGMA46 = self.sigma_a1_a3[1]	#8
		SIGMA64 = self.sigma_a3_a1[1]	#9
		
		constraints = SoftConstraintContainer( )
		self.molecule.DefineSoftConstraints( constraints )
		
		for i in range( x ):
			for j in range ( y ):
				
				coordinate_file = pat +"frame{}_{}.pkl".format( i, j )
				self.molecule.coordinates3 = Unpickle( coordinate_file )
				
				distance_1 = self.molecule.coordinates3.Distance( ATOM1, ATOM2 ) - self.molecule.coordinates3.Distance( ATOM2, ATOM3 )
				scModel  = SoftConstraintEnergyModelHarmonic ( distance_1, self.forceC )
				constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
				constraints["ReactionCoord"] = constraint
				
				distance_2 = self.molecule.coordinates3.Distance( ATOM4, ATOM5 ) - self.molecule.coordinates3.Distance( ATOM5, ATOM6 )
				scModel2  = SoftConstraintEnergyModelHarmonic ( distance_2, self.forceC )
				constraint2=SoftConstraintMultipleDistance ( [ [ATOM5, ATOM4, SIGMA46], [ATOM5, ATOM6, SIGMA64] ], scModel2 )
				constraints["ReactionCoord2"] = constraint2
				
				temp = coordinate_file[:-4]
				temp = os.path.basename(temp)
				md_path	= os.path.join( self.folder, temp )
				md_run	= MD(self.molecule,md_path,True)
				md_run.Run()
	
	#.------------------------------------------------------------------
	def Sample2D_parallel(self, sample):
		
		pat = os.path.join(self.traj_folder,"trj","")
		files_list = glob.glob(pat+"*.pkl")
		files_list.sort()
		self.all_files = files_list
		
		ATOM1 = self.atoms[0][0]		#0
		ATOM2 = self.atoms[0][1]		#1
		ATOM3 = self.atoms[0][2]		#2
		ATOM4 = self.atoms[1][0]		#3
		ATOM5 = self.atoms[1][1]		#4
		ATOM6 = self.atoms[1][2]		#5
		SIGMA13 = self.sigma_a1_a3[0]	#6
		SIGMA31 = self.sigma_a3_a1[0]	#7
		SIGMA46 = self.sigma_a1_a3[1]	#8
		SIGMA64 = self.sigma_a3_a1[1]	#9
		
		constraints = SoftConstraintContainer( )
		self.molecule.DefineSoftConstraints( constraints )
		
		with pymp.Parallel(self.nprocs) as p:
			for i in p.range( len(files_list) ):
								
				coordinate_file = files_list[i]

				temp = coordinate_file[:-4]
				temp = os.path.basename(temp)
				md_path	= os.path.join( self.folder, temp )
				if  not os.path.exists( md_path +"_softc_trj" ):
					self.molecule.coordinates3 = Unpickle( coordinate_file )
									
					distance_1 = self.molecule.coordinates3.Distance( ATOM1, ATOM2 ) - self.molecule.coordinates3.Distance( ATOM2, ATOM3 )
					scModel  = SoftConstraintEnergyModelHarmonic ( distance_1, self.forceC )
					constraint=SoftConstraintMultipleDistance ( [ [ATOM2, ATOM1, SIGMA13], [ATOM2, ATOM3, SIGMA31] ], scModel )
					constraints["ReactionCoord"] = constraint
				
					distance_2 = self.molecule.coordinates3.Distance( ATOM4, ATOM5 ) - self.molecule.coordinates3.Distance( ATOM5, ATOM6 )
					scModel2  = SoftConstraintEnergyModelHarmonic ( distance_2, self.forceC )
					constraint2=SoftConstraintMultipleDistance ( [ [ATOM5, ATOM4, SIGMA46], [ATOM5, ATOM6, SIGMA64] ], scModel2 )
					constraints["ReactionCoord2"] = constraint2
					md_run	= MD(self.molecule,md_path,True)
					if sample:
						md_run.samplen = 1000
						print("set sampling")
					md_run.Run()
				else:
					continue
		
#=======================================================================
class PMF:

	def __init__(self,name,molecule,folder,ndim):
		self.base_name	= name
		self.molecule	= molecule
		self.ndim		= ndim
		self.src_folder	= folder
		self.fileNames	= []
		
		self.text		= ""
		self.LOG		= open(self.base_name+"_FE.log","w")
		pat = os.path.join( folder,"" )
		
		self.fileNames = glob.glob ( pat+"*_softc_trj" )
		self.fileNames.sort()
		
	#.------------------------------------------------------------------
	def calculate_1D(self,bins):
		
		state = WHAM_ConjugateGradientMinimize ( self.fileNames			,
                                         bins                 = [ bins ],
                                         logFrequency         =      1  ,
                                         maximumIterations    =   1000  ,
                                         rmsGradientTolerance = 1.0e-3  ,
                                         temperature          = 300.15   )

		histogram = state["Histogram"]
		pmf       = state["PMF"      ]
		FE		  = state["Free Energies"]
		histogram.ToTextFileWithData ( self.base_name+".dat" , [ pmf ], format = "{:20.3f} {:20.3f}\n" )
		
		self.text += "X FE\n"
		for i in range( len(FE) ):
			a = os.path.basename( self.fileNames[i] )
			try:
				b = int( a[5:6] )
			except:
				b = int( a[5:7] )

			self.text += "{}  {}\n".format( x,  FE[i] )
			
			self.LOG.write( self.text )
			self.LOG.close( )
		
	#.------------------------------------------------------------------
	def calculate_2D(self,bins_x,bins_y):
		
		state = WHAM_ConjugateGradientMinimize ( self.fileNames						,
                                         bins                 = [ bins_x, bins_y ]	,
                                         logFrequency         =      1				,
                                         maximumIterations    =   1000				,
                                         rmsGradientTolerance = 1.0e-3				,
                                         temperature          = 300.15				)

		histogram 	= state["Histogram"]
		pmf       	= state["PMF"      ]
		FE			= state["Free Energies"]
		
		self.text += "X Y FE\n"
		x = -1
		y = -1
		for i in range( len(FE) ):
			x = -1
			y = -1
			a = os.path.basename( self.fileNames[i] )
			b = a[5:-10]
			if len(b) == 5:
				x = int( b[0:2] )
				y = int( b[3:5] )
			elif len(b) == 3:
				x = int( b[0:1] )
				y = int( b[2:3] )
			elif len(b) == 4:
				try: 
					x = int( b[0:1] )
					y = int( b[2:4] )
				except:
					x = int( b[0:2] )
					y = int( b[3:4] )
			if x == -1:
				self.text += "{} {}\n".format( i, FE[i]-FE[0] )
			else:
				self.text += "{} {} {}\n".format( x,y, FE[i]-FE[0] )
		
		histogram.ToTextFileWithData ( self.base_name+".dat" , [ pmf ], format = "{:20.3f} {:20.3f} {:20.3f}\n" )
		
		self.LOG.write( self.text )
		self.LOG.close( )

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
	
	simulation  = umbrella_Sampling(tim_qcmm,"tim_US_2Dpar","TIM2D_seqpDynamo_pj")
	simulation.nprocs = 5
	simulation.set_modes(rc1_atoms)
	simulation.set_modes(rc2_atoms)
	simulation.Sample2D_parallel(15,15)
	
	pmf = PMF("tim_US_2Dpar",tim_qcmm,"tim_US_2Dpar",2)
	pmf.calculate_2D(60,60)
	'''
