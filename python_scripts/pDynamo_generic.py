#pdynamo simulation functions definitions, tests and functionalities

#=======================================================================
#-----------------------------------------------------------------------
#import python libraries
from __future__ import print_function
import glob, math, os, numpy
import multiprocessing as mp

import pymp

#-----------------------------------------------------------------------
#import pDynamo libraries
from pCore import *
from pBabel import *
from pMolecule import *
from pMoleculeScripts import *
   
atomic_mass = {1   :1.007940  ,
2   :4.002602  ,
3   :6.941000  ,
4   :9.012182  ,
5   :10.811000 ,
6   :12.010700 ,
7   :14.006700 ,
8   :15.999400 ,
9   :18.998404 ,
10  :20.179701 ,
11  :22.989771 ,
12  :24.305000 ,
13  :26.981539 ,
14  :28.085501 ,
15  :30.973761 ,
16  :32.064999 ,
17  :35.452999 ,
18  :39.948002 ,
19  :39.098301 ,
20  :40.077999 ,
21  :44.955910 ,
22  :47.867001 ,
23  :50.941502 ,
24  :51.996101 ,
25  :54.938049 ,
26  :55.845001 ,
27  :58.933201 ,
28  :58.693401 ,
29  :63.546001 ,
30  :65.408997 ,
31  :69.723000 ,
32  :72.639999 ,
33  :74.921600 ,
34  :78.959999 ,
35  :79.903999 ,
36  :83.797997 ,
37  :85.467796 ,
38  :87.620003 ,
39  :88.905853 ,
40  :91.223999 ,
41  :92.906380 ,
42  :95.940002 ,
43  :98.000000 ,
44  :101.070000,
45  :102.905502,
46  :106.419998,
47  :107.868202,
48  :112.411003,
49  :114.818001,
50  :118.709999,
51  :121.760002,
52  :127.599998,
53  :126.904472,
54  :131.292999,
55  :132.905457,
56  :137.326996,
57  :138.905502,
58  :140.115997,
59  :140.907654,
60  :144.240005,
61  :145.000000,
62  :150.360001,
63  :151.964005,
64  :157.250000,
65  :158.925339,
66  :162.500000,
67  :164.930313,
68  :167.259003,
69  :168.934204,
70  :173.039993,
71  :174.966995,
72  :178.490005,
73  :180.947906,
74  :183.839996,
75  :186.207001,
76  :190.229996,
77  :192.216995,
78  :195.078003,
79  :196.966553,
80  :200.589996,
81  :204.383301,
82  :207.199997,
83  :208.980377,
84  :209.000000,
85  :210.000000,
86  :222.000000,
87  :223.000000,
88  :226.000000,
89  :227.000000,
90  :232.038101,
91  :231.035873,
92  :238.028915,
93  :237.000000,
94  :244.000000,
95  :243.000000,
96  :247.000000,
97  :247.000000,
98  :251.000000,
99  :252.000000,
100 :257.000000,
101 :258.000000,
102 :259.000000,
103 :262.000000,
104 :261.000000,
105 :262.000000,
106 :263.000000,
107 :264.000000,
108 :265.000000,
109 :268.000000,
0   :0.000000  ,
0   :0.000000  }
#===============================================================================
#-----------------------------------------------------------------------
#pre-defining QC energy models

             
PM6 = QCModelMNDO ( "pm6"  )
PM3 = QCModelMNDO ( "pm3"  )
RM1 = QCModelMNDO ( "rm1"  )
AM1 = QCModelMNDO ( "am1"  )
AM1dPhot = QCModelMNDO ( "am1dphot" )

nbmodel = NBModelABFS ( )
             
#-----------------------------------------------------------------------
#Prepare from gromacs simulations
def load_gromacs_top(base_name):
	"""	
		Function to load coordinates and topologies.
		base_name = base name for the gromacs topologies and coordinates
		system	  = instance of System class
	"""
	
	parameters = GromacsParameters_ToParameters ( "{}.top".format(base_name) )
	system     = GromacsDefinitionsReader.ToSystem   ( "{}.top".format(base_name) ,parameters = parameters )
	system.coordinates3 = ImportCoordinates3 ( "{}.gro".format(base_name) )
	
	return (system)
	
#-----------------------------------------------------------------------	

def load_system_from_coord(file_name):
	
	filename, file_extension = os.path.splitext(file_name)
		
	if file_extension == ".xyz":
		mol = XYZFile_ToSystem ( file_name )
		return mol
	elif file_extension == ".pdb":
		mol = PDBFile_ToSystem ( file_name )
		return mol
	
#-----------------------------------------------------------------------

def calc_molecule_energies(molecule=None, models = None):
	
	energies = []
	for model in models:		
		molecule.DefineQCModel ( model )
		molecule.Summary ( )
		energy  = molecule.Energy ( )
		energies.append( energy )
	return energies

#-----------------------------------------------------------------------

def Geometry_Optimization(molecule=None,out_path=None,rmsT=0.1,alg="cgrad",traj=False):
	
	
	if traj:
		t_opt = SystemGeometryTrajectory ( out_path[:-4] + ".trj" , molecule, mode = "w" )	
		if alg == "cgrad":
			ConjugateGradientMinimize_SystemGeometry(molecule				,
												logFrequency=10			,
												maximumIterations=5000	,
												trajectories =  [( t_opt, 10 ) ],
												rmsGradientTolerance=rmsT )
		elif alg == "steep":
			SteepestDescentMinimize_SystemGeometry(molecule						,
                                               logFrequency  = 10			,
                                               maximumIterations  = 5000	,
                                               trajectories =  [( t_opt, 10 ) ],
                                               rmsGradientTolerance =rmsT 	)
		elif alg == "lbfgs":
			LBFGSMinimize_SystemGeometry(molecule					,
                                     logFrequency = 50			,
									 maximumIterations= 5000	,
									 trajectories =  [( t_opt, 10 ) ],
									 rmsGradientTolerance = rmsT)
	else:
		if alg == "cgrad":
			ConjugateGradientMinimize_SystemGeometry(molecule				,
												logFrequency=10			,
												maximumIterations=5000	,
												rmsGradientTolerance=rmsT )
		elif alg == "steep":
			SteepestDescentMinimize_SystemGeometry(molecule						,
                                               logFrequency  = 10			,
                                               maximumIterations  = 5000	,
                                               rmsGradientTolerance =rmsT 	)
		elif alg == "lbfgs":
			LBFGSMinimize_SystemGeometry(molecule					,
                                     logFrequency = 50			,
									 maximumIterations= 5000	,
									 rmsGradientTolerance = rmsT)		
	
	
	PDBFile_FromSystem(out_path,molecule)
	if traj:
		DCDTrajectory_FromSystemGeometryTrajectory(out_path[:-4] + ".dcd",out_path[:-4] + ".trj",molecule)


#-----------------------------------------------------------------------
def Get_total_charge(molecule):
	totalCharge = 0	
	Charge=molecule.energyModel.mmAtoms.AtomicCharges()
	for i in range(len(Charge)):
		totalCharge+=Charge[i]
	return (totalCharge)

#-----------------------------------------------------------------------
def reescale_charges(molecule, tc):
	scaled_system = Clone(molecule)
	
	Charges = scaled_system.energyModel.mmAtoms.AtomicCharges()
	pTC = Get_total_charge(scaled_system)
	
	print ("Old Charges Sum:",pTC)
	
	nAtoms = len(scaled_system.atoms.items)
	
	new_charges = []
	new_tc = 0
	frac = (tc - pTC)/nAtoms
	
	for i in range( len(Charges) ):
		new_charges.append(Charges[i] + frac)
	
	for i in range( len(new_charges) ):
		Charges[i] = new_charges[i]
		
	
	scaled_system.energyModel.mmAtoms.SetAtomicCharges(Charges)
	
	new_tc = Get_total_charge(scaled_system)
	print ("New Charges Sum:",new_tc)	
	
	return(scaled_system)	
	
#-----------------------------------------------------------------------

def amber12_to_amber11_topology_converter (filein, fileout):
	filein = open(filein, 'r')
	text   = []
	print_line = True

	for line in filein:
		line2 = line.split()
		try:
			if line2[0] == '%FLAG':
				if   line2[1] == 'ATOMIC_NUMBER':
					print ('excluding flag:', line)
					print_line = False

				elif   line2[1] == 'SCEE_SCALE_FACTOR':
					print ('excluding flag:', line)
					print_line = False

				elif   line2[1] == "SCNB_SCALE_FACTOR":
					print ('excluding flag:', line)
					print_line = False			

				elif   line2[1] == 'IPOL':
					print ('excluding flag:', line)
					print_line = False
		
				else:
					print_line = True	
					#print print_line
		except:
			a= None
		if print_line == True:
			text.append(line)

	fileout = open(fileout, 'w')
	fileout.writelines(text)
	fileout.close()
	
#=======================================================================
class SCAN:
	'''	
	
	'''	
	def __init__(self,molecule,name):
		'''
		'''
		self.base_name 			= name
		self.molecule  			= molecule 
		self.ndim	   			= 0
		self.rcs       			= []
		self.atoms	   			= [] # indices of the atoms 
		self.nprocs    			= 1
		self.log 	   			= open(self.base_name + "scan.log","w")
		self.text	   			= " "
		self.energies  			= []
		self.DMINIMUM  			= [ 0, 0 ]
		self.DINCREMENT			= [ 0, 0 ]
		self.mass_constraint	= True
		self.multiple_d			= False
		self.nsteps	   			= [ 1, 1 ]
		self.maxIt				= 30
		self.rmsGT				= 0.1
		self.minimumAlg			= "cgrad" 
		
		self.text = "Information of pDynamo SCAN \n"
		
	#---------------------------------------------------
	def set_rcs(self,na,atoms,dincre):
		
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
	def change_parameters(self):
		
		pass
		
	#----------------------------------------------------
	def write_log(self):
		
		pass
		
	#----------------------------------------------------	
	def Print(self):
		text = "Print Scan parameters\n"
		text+= "Num of RCs:  {}".format(self.ndim) 
		print(text)
		pass
#=======================================================================
class NEB:
	pass	
#=======================================================================	
class umbrella_Sampling:
	pass 



if __name__ == "__main__":	
	tim_qcmm = Unpickle("TIM_qcmm_opt.pkl")
	C02 =  AtomSelection.FromAtomPattern(tim_qcmm,"*:LIG.248:C02")
	proton = AtomSelection.FromAtomPattern(tim_qcmm,"*:LIG.248:H02")
	oxygen = AtomSelection.FromAtomPattern(tim_qcmm,"*:GLU.164:OE2")
	
	atom1 = C02.selection.pop()
	atom2 = proton.selection.pop()
	atom3 = oxygen.selection.pop()
	
	rc1_atoms = [atom1,atom2,atom3]	
	scan1d = SCAN(tim_qcmm,"TIM")	
	scan1d.set_rcs(3,rc1_atoms,0.1)
	#print(scan1d.DMINIMUM[0])



