'''
Class implementation for .hist files 
Variable names follow the naming scheme used in the histogram files with the exception of a few
WARNING: This file is pretty dependent on how the .hist file is outputted so if the .hist file is of an older version, this code might not work in getting all
the correct header information especially if spacing is changed
'''

#===============================================================================
# IMPORTS
#===============================================================================


#===============================================================================
# Histogram Data Class
#===============================================================================

class histData():

    def __init__(self, generationDate=[], eulerAngles=[], lambdaRange=[], betaRange=[], lambdaBinSize=[], betaBinSize=[], nbody=[], evolveBackwardTime=[], evolveFowardTime=[], bestLikelihood=[], timestep=[], criterion=[], sunGCDist=[], theta=[], quadrupoleMoments=[], eps=[], barTimeError=[], barAngleError=[], spherical=[], primaryDisk=[], secondaryDisk=[], halo=[], n=[], massPerParticle=[], totalSimulated=[], lambdaBins=[], betaBins=[], columnHeaders=[], useBin=[], lamb=[], beta=[], normalizedCounts=[], countError=[], betaDispersion=[], betaDispersionError=[], losVelocityDispersion=[], velocityDispersionError=[], losVelocity=[], losVelocityError=[], betaAverage=[], betaAverageError=[], distance=[], distanceError=[], orbitFitting = False):
        
        #Header portion of Histogram 
        self.generationDate = generationDate
        self.eulerAngles = eulerAngles #degrees
        self.lambdaRange = lambdaRange 
        self.betaRange = betaRange
        self.lambdaBinSize = lambdaBinSize
        self.betaBinSize = betaBinSize

        self.nbody = nbody #total number of bodies in n-body simulation
        self.evolveBackwardTime = evolveBackwardTime #Gyr
        self.evolveFowardTime = evolveFowardTime #Gyr
        self.bestLikelihood = bestLikelihood
        self.timestep = timestep #Gyr
        self.sunGCDist = sunGCDist #kpc
        self.criterion = criterion
        self.theta = theta
        self.quadrupoleMoments = quadrupoleMoments
        self.eps = eps #kpc
        self.barTimeError = barTimeError #Gyr
        self.barAngleError = barAngleError #radians
        
        #Potential portion of Histogram
        self.spherical = spherical #buldge potential
        self.primaryDisk = primaryDisk
        self.secondaryDisk = secondaryDisk
        self.halo = halo

        #More histogram infomation
        self.n = n #number of bodies within histogram range
        self.massPerParticle = massPerParticle #structure mass units
        self.totalSimulated = totalSimulated #total number of baryons in simulation
        self.lambdaBins = lambdaBins
        self.betaBins = betaBins

        #Column Headers
        self.columnHeaders = columnHeaders 

        #Data Portion of Histogram
        self.useBin = useBin
        self.lamb = lamb #lamb instead of lambda since lambda is a keyword in python 
        self.beta = beta
        self.normalizedCounts = normalizedCounts
        self.countError = countError
        self.betaDispersion = betaDispersion
        self.betaDispersionError = betaDispersionError
        self.losVelocityDispersion = losVelocityDispersion
        self.velocityDispersionError = velocityDispersionError
        self.losVelocity = losVelocity
        self.losVelocityError = losVelocityError
        self.betaAverage = betaAverage
        self.betaAverageError = betaAverageError
        self.distance = distance
        self.distanceError = distanceError

        #Marker for if orbit fitting is included in histogram
        self.orbitFitting = orbitFitting


        self.array_dict = {'generationDate':self.generationDate, 'eulerAngles':self.eulerAngles, 'lambdaRange':self.lambdaRange, 'betaRange':self.betaRange, \
                            'lambdaBinSize':self.lambdaBinSize, 'betaBinSize':self.betaBinSize, 'nbody':self.nbody, 'evolveBackwardTime':self.evolveBackwardTime, \
                            'evolveFowardTime':self.evolveFowardTime, 'bestLikelihood':self.bestLikelihood, 'timestep':self.timestep, 'sundGCDist':self.sunGCDist, \
                            'criterion':self.criterion, 'theta':self.theta, 'quadrupoleMoments':self.quadrupoleMoments, 'eps':self.eps, 'barTimeError':self.barTimeError, \
                            'barAngleError':self.barAngleError, 'spherical':self.spherical, \
                            'primaryDisk':self.primaryDisk, 'secondaryDisk':self.secondaryDisk, 'halo':self.halo, 'n':self.n, 'massPerParticle':self.massPerParticle,\
                            'columnheaders':self.columnHeaders, 'totalSimulated':self.totalSimulated, 'lambdaBins':self.lambdaBins, 'betaBins':self.betaBins, \
                            'useBin':self.useBin, 'lamb':self.lamb, 'beta':self.beta, 'normalizedCounts':self.normalizedCounts, 'countError':self.countError,  \
                            'betaDispersion':self.betaDispersion, 'betaDispersionError':self.betaDispersionError, 'losVelocityDispersion':self.losVelocityDispersion, \
                            'velocityDispersionError':self.velocityDispersionError, 'losVelocity':self.losVelocity, 'losVelocityError':self.losVelocityError, \
                            'betaAverage':self.betaAverage, 'betaAverageError':self.betaAverageError, \
                            'distance':self.distance, 'distanceError':self.distanceError, 'orbitFitting':self.orbitFitting}


