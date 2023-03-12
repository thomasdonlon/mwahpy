'''
Import data from milkyway .hist files
Variable names follow the naming scheme used in the histogram files
WARNING: This file is pretty dependent on how the .hist file is outputted so if the .hist file is of an older version, this code might not work in getting all
the correct header information
'''

#===============================================================================
# IMPORTS
#===============================================================================


#===============================================================================
# Histogram Data Class
#===============================================================================

class histData():

    def __init__(self, generationDate=[], eulerAngles=[], lambdaRange=[], betaRange=[], lambdaBinSize=[], betaBinSize=[], nbody=[], evolveBackwardTime=[], evolveFowardTime=[], bestLikelihood=[], timestep=[], criterion=[], sunGCDist=[], theta=[], quadrupoleMoments=[], eps=[], barTimeError=[], barAngleError=[], spherical=[], primaryDisk=[], secondaryDisk=[], halo=[], n=[], massPerParticle=[], totalSimulated=[], lambdaBins=[], betaBins=[], columnHeaders=[], useBin=[], lamb=[], beta=[], normalizedCounts=[], countError=[], betaDispersion=[], betaDispersionError=[], losVelocityDispersion=[], velocityDispersionError=[], losVelocity=[], losVelocityError=[], betaAverage=[], betaAverageError=[], distance=[], distanceError=[]):
        
        #Header portion of Histogram 
        self.generationDate = generationDate
        self.eulerAngles = eulerAngles
        self.lambdaRange = lambdaRange
        self.betaRange = betaRange
        self.lambdaBinSize = lambdaBinSize
        self.betaBinSize = betaBinSize

        self.nbody = nbody
        self.evolveBackwardTime = evolveBackwardTime
        self.evolveFowardTime = evolveFowardTime
        self.bestLikelihood = bestLikelihood
        self.timestep = timestep
        self.sunGCDist = sunGCDist
        self.criterion = criterion
        self.theta = theta
        self.quadrupoleMoments = quadrupoleMoments
        self.eps = eps
        self.barTimeError = barTimeError
        self.barAngleError = barAngleError
        
        #Potential portion of Histogram
        self.spherical = spherical
        self.primaryDisk = primaryDisk
        self.secondaryDisk = secondaryDisk
        self.halo = halo

        #More histogram infomation
        self.n = n
        self.massPerParticle = massPerParticle
        self.totalSimulated = totalSimulated
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


        self.array_dict = {'date':self.generationDate, 'eulerAngles':self.eulerAngles, 'lambdaRange':self.lambdaRange, 'betaRange':self.betaRange, \
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
                            'distance':self.distance, 'distanceError':self.distanceError}

#===============================================================================
# FUNCTIONS
#===============================================================================

#Reads in histogram file path as a string and creates histData class
def read_histogram(hist_file_path):
    
    #initialising histData class instance
    hist = histData()

    with open(str(hist_file_path)) as f:
        line = f.readline()
        lineNumber = 1
        while (line):
            line = line.strip()
            line = line.split(' ')
            #print(line)
            if (len(line)> 1):
                #Date histogram was generated
                if (line[1]=='Generated'):
                    generationDate = ''
                    for i in range(len(line) - 2):
                        generationDate = generationDate + line[i+2]+ ' '
                    hist.generationDate = generationDate
                    #print(hist.generationDate)
                #Euler Angles
                if (line[1]=='(phi,'):
                    eulerPhi = float(line[5].strip('(').strip(','))
                    eulerTheta = float(line[6].strip(','))
                    eulerPsi = float(line[7].strip(')').strip(','))
                    hist.eulerAngles = [eulerPhi, eulerTheta, eulerPsi]
                    #print(hist.eulerAngles)
                #Lambda Rnage
                if (line[1]=='(lambdaStart,'):
                    lambdaStart = float(line[5].strip('(').strip(','))
                    lambdaCenter = float(line[6].strip(','))
                    lambdaEnd = float(line[7].strip(')').strip(','))
                    hist.lambdaRange = [lambdaStart, lambdaCenter, lambdaEnd]
                    #print(hist.lambdaRange)
                #Beta Range
                if (line[1]=='(betaStart,'):
                    betaStart = float(line[5].strip('(').strip(','))
                    betaCenter = float(line[6].strip(','))
                    betaEnd = float(line[7].strip(')').strip(','))
                    hist.betaRange = [betaStart, betaCenter, betaEnd]
                    #print(hist.betaRange)
                #Lambda Bin Size
                if (line[1]=='Lambda'):
                    hist.lambdaBinSize = float(line[5])
                    #print(hist.lambdaBinSize)
                #Beta Bin Size
                if (line[1]=='Beta'):
                    hist.betaBinSize = float(line[5])
                    #print(hist.betaBinSize)
                #Nbody (total numbers of body in Nbody simulation)
                if (line[1]=='Nbody'):
                    hist.nbody = int(line[3])
                    #print(hist.nbody)
                #Evolve backward time and evolve foward time
                if (line[1]=='Evolve'):
                    if (line[2]=='backward'):
                        hist.evolveBackwardTime = float(line[5])
                        #print(hist.evolveBackwardTime)
                    elif (line[2]=='forward'):
                        hist.evolveFowardTime = float(line[5])
                        #print(hist.evolveFowardTime)
                #Best Likelihood (absolute value of the likelihood score)
                if (line[1]=='Best'):
                    hist.bestLikelihood = float(line[4])
                    #print(hist.bestLikelihood)
                #Timestep
                if (line[1]=='Timestep'):
                    hist.timestep = float(line[3])
                    #print(hist.timestep)
                #Sun GC Distance 
                if (line[1]=='Sun'):
                    hist.sunGCDist = float(line[5])
                    #print(hist.sunGCDist)
                #Criterion
                if (line[1]=='Criterion'):
                    hist.criterion = line[3]
                    #print(hist.criterion)
                #Theta (not sure what this is for)
                if (line[1]=='Theta'): 
                    hist.theta = float(line[3])
                    #print(hist.theta)
                #Quadrupole Moments (not sure what this is for)
                if (line[1]=='Quadrupole'):
                    hist.quadrupoleMoments = line[4]
                    #print(hist.quadrupoleMoments)
                #Eps (not sure what this is for)
                if (line[1]=='Eps'):
                    hist.eps = float(line[3])
                    #print(hist.eps)
                #Bar time error and bar angle error in radians
                if (line[1]=='Bar'):
                    if (line[2]=='Time'):
                        hist.barTimeError = float(line[5])
                        #print(hist.barTimeError)
                    if (line[2]=='Angle'):
                        hist.barAngleError = float(line[5])
                        #print(hist.barAngleError)
                #Spherical Potential Type (just saves the type of potential and not the mass and such)
                if (line[1]=='Spherical:'):
                    hist.spherical = line[2]
                    #print(hist.spherical)
                #Primary Disk Potential Type  (just saves the type of potential and not the mass and such)
                if (line[1]=='Primary'):
                    hist.primaryDisk = line[3]
                    #print(hist.primaryDisk)
                #Secondary Disk Potential Type (just saves the type of potential and not the mass and such)
                if (line[1]=='Secondary'):
                    hist.secondaryDisk = line[3]
                    #print(hist.secondaryDisk)
                #Halo Potential Type (just saves the type of potential and not the mass and such)
                if (line[1]=='Halo:'):
                    hist.halo = line [2]
                    #print(hist.halo)
                #Column headers
                if (line[1]=='UseBin,'):
                    columnHeaders = ''
                    for i in range(len(line) - 1):
                        columnHeaders += line[i + 1] + ' '
                    hist.columnHeaders = columnHeaders
                    #print(hist.columnHeaders)
                #n (not sure where this number comes from)
                if (line[0]=='n'):
                    hist.n = int(line[2])
                    #print(hist.n)
                #Mass per particle
                if (line[0]=='massPerParticle'):
                    hist.massPerParticle = float(line[2])
                    #print(hist.massPerParticle)
                if (line[0]=='totalSimulated'):
                    hist.totalSimulated = int(line[2])
                    #print(hist.totalSimulated)
                #Number of lambda bins
                if (line[0]=='lambdaBins'):
                    hist.lambdaBins = int(line[2])
                    #print(hist.lambdaBins)
                #number of Beta Bins
                if (line[0]=='betaBins'):
                    hist.betaBins = int(line[2])
                    #print(hist.betaBins)
                #Getting Data (Not sure if useBin can have any other value other than 1)
                if (line[0]=='1'):
                    try:
                        hist.useBin.append(float(line[0]))
                        hist.lamb.append(float(line[1]))
                        hist.beta.append(float(line[2]))
                        hist.normalizedCounts.append(float(line[3]))
                        hist.countError.append(float(line[4]))
                        hist.betaDispersion.append(float(line[5]))
                        hist.betaDispersionError.append(float(line[6]))
                        hist.losVelocityDispersion.append(float(line[7]))
                        hist.velocityDispersionError.append(float(line[8]))
                        hist.losVelocity.append(float(line[9]))
                        hist.losVelocityError.append(float(line[1]))
                        hist.betaAverage.append(float(line[11]))
                        hist.betaAverageError.append(float(line[12]))
                        hist.distance.append(float(line[13]))
                        hist.distanceError.append(float(line[12]))
                    except ValueError:
                        print('Invalid histagram data entry: Non-numerical value in histogram data on line ' + str(lineNumber))
                
            line = f.readline()
            lineNumber += 1
    f.close()
    return hist


