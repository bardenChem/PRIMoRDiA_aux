#Python Class for automatization of Nudget Elastic Band 


from pDynamo_generic import *

#=======================================================================

class NEB:
	def __init__(self,molecule,name):
		self.molecule = molecule
		self.base_name= name
	
	def run(self,reac,prod,bins):
		reactants = Unpickle(reac)
		products  = Unpickle(prod)
		
		trajectoryPath = os.path.join ( self.base_name, "NEB_trj" )
		GrowingStringInitialPath ( self.molecule, bins, reactants, products, trajectoryPath )
		
		trajectory = SystemGeometryTrajectory ( trajectoryPath, self.molecule, mode = "a+" )
		ChainOfStatesOptimizePath_SystemGeometry (	self.molecule				, 
													trajectory					,
													logFrequency         = 1	,
													maximumIterations    = 100	,
													rmsGradientTolerance = 0.1	)
	
#=======================================================================
if __name__ == "__main__":	
	pass
