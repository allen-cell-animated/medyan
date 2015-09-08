from numpy import *
import pylab as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import scipy.sparse.linalg as la
import math

#A filament snapshot
class FilamentSnapshot:
	def __init__(self):
		self.id=None
		self.length=None
		self.delta_left=None	
		self.delta_right=None	
		self.coords={}

#A cylinder snapshot
class CylinderSnapshot:
	def __init__(self):
		self.id = None
		self.coords={}

#A linker snapshot
class LinkerSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.coords={}

#A motor snapshot
class MotorSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.coords={}

#A brancher snapshot
class BrancherSnapshot:
	def __init__(self):
		self.id=None
		self.type=None;
		self.coords={}

#A complete snapshot
class Snapshot:
	def __init__(self):
		self.step=None
		self.time=None
		self.n_filaments=None
		self.n_linkers=None
		self.n_motors=None
		self.n_branchers=None
		self.filaments={}
		self.cylinders={}
		self.linkers={}
		self.motors={}
		self.branchers={}

#A chemical snapshot
class ChemSnapshot:
	def __init__(self):
		self.step=None
		self.time=None
		self.diffusingSpecies={}
		self.bulkSpecies={}
		self.filamentSpecies={}
		self.plusEndSpecies={}
		self.minusEndSpecies={}
		self.linkerSpecies={}
		self.motorSpecies={}
		self.brancherSpecies={}

########################################################
#
#			         PARSING FUNCTIONS
#
########################################################


#Read all data from a trajectory file
#returns a list of Snapshots with all data
def readTrajectory(filename):

	print "Reading " + filename + "..."

	#Open the traj file
	traj_file=open(filename)

	SnapshotList=[]
	first_Snapshot_line=True
	first_line=True
	reading_filament=False
	reading_linker=False
	reading_motor=False
	reading_brancher=False
	n_beads_filament=0
	n_beads_linker=0
	n_beads_motor=0
	n_beads_brancher=0;

	line_number = 0;

	for line in traj_file:
		line = line.strip()

		#get Snapshot data

		if(first_Snapshot_line):
			F=Snapshot()
			F.step, F.time, F.n_filaments, F.n_linkers, \
			F.n_motors, F.n_branchers = map(double,line.split())
			F.n_filaments = int(F.n_filaments)
			F.n_linkers = int(F.n_linkers)
			F.n_motors = int(F.n_motors)
			F.n_branchers = int(F.n_branchers)
			first_Snapshot_line=False
			line_number+=1
			continue

		if(len(line)==0):
			first_Snapshot_line=True
			first_filament_line=True
			SnapshotList.append(F)
			assert F.n_filaments == len(F.filaments)
			assert F.n_linkers   == len(F.linkers)
			assert F.n_motors    == len(F.motors)
			assert F.n_branchers == len(F.branchers)
			n_beads_filament=0
			n_beads_linker=0
			n_beads_motor=0
			n_beads_brancher=0;
			line_number+=1
			continue
				
		#Read the line
		if(first_line):

			line_split = line.split()

			if(line[0] == "F"):
				FS=FilamentSnapshot()
				FS.id, FS.length, FS.delta_left, FS.delta_right = map(int,line_split[1:])
				first_line=False
				F.filaments[FS.id]=FS
				n_beads_filament+=FS.length
				reading_filament=True
				line_number+=1
				continue

			if(line[0] == "L"):
				LS = LinkerSnapshot()
				LS.id, LS.type = map(int, line_split[1:])
				first_line = False
				F.linkers[LS.id] = LS
				n_beads_linker+=2
				reading_linker=True
				line_number+=1
				continue

			if(line[0] == "M"):
				MS = MotorSnapshot()
				MS.id, MS.type = map(int, line_split[1:])
				first_line = False
				F.motors[MS.id] = MS
				n_beads_motor+=2
				reading_motor=True
				line_number+=1
				continue

			if(line[0] == "B"):
				BS = BrancherSnapshot()
				BS.id, BS.type = map(int, line_split[1:])
				first_line = False
				F.branchers[BS.id] = BS
				n_beads_brancher+=1
				reading_brancher=True
				line_number+=1
				continue

		#Format coordinates
		if(reading_filament):
			FS.coords=array(line.split(),'f')
			N=len(FS.coords)
			FS.coords = FS.coords.reshape(N/3,3)

			#create cylinder snapshots
			for i in xrange(0, len(FS.coords) - 1):
				CS = CylinderSnapshot()
				#create unique id using cylinder position and fil id
				CS.id = (FS.id + i) * (FS.id + i + 1) / 2 + FS.id
				CS.coords = [FS.coords[i], FS.coords[i+1]]
				F.cylinders[CS.id] = CS

			first_line=True
			reading_filament=False

		if(reading_linker):
			LS.coords=array(line.split(),'f')
			N=len(LS.coords)
			LS.coords = LS.coords.reshape(N/3,3)
			first_line=True
			reading_linker=False

		if(reading_motor):
			MS.coords=array(line.split(),'f')
			N=len(MS.coords)
			MS.coords = MS.coords.reshape(N/3,3)
			first_line=True
			reading_motor=False

		if(reading_brancher):
			BS.coords=array(line.split(),'f')
			N=len(BS.coords)
			BS.coords = BS.coords.reshape(N/3,3)
			first_line=True
			reading_brancher=False

		line_number+=1

	#return the final Snapshot list
	print str(len(SnapshotList)) + " snapshots read."

	return SnapshotList	

#Read a number of trajectories, return array of Snapshot lists
def readTrajectories(fileNames):

	SnapshotLists = []

	for trajFile in fileNames:
		print "Reading trajectory..."
		SnapshotLists.append(readTrajectory(trajFile))

	return SnapshotLists

#Read all data from a chemical output file
#return a list of chemical snapshots
def readChemistryOutput(fileName):

	print "Reading " + fileName + "..."

	ChemSnapshotList = []

	#Open the traj file
	chem_file=open(fileName)

	SnapshotList=[]
	first_Snapshot_line=True
	first_line=True
	line_number=0

	for line in chem_file:
		line = line.strip()

		if(first_Snapshot_line):
			C=ChemSnapshot()
			C.step, C.time = map(double,line.split())
			first_Snapshot_line=False
			line_number+=1
			continue

		elif(len(line)==0):
			first_Snapshot_line=True
			ChemSnapshotList.append(C)
			line_number+=1
			continue

		else:
			line_split = line.split()
			
			species_split = line_split[0].split(":")

			species_name = species_split[0]
			species_type = species_split[1]

			if(species_type == "DIFFUSING"):
				C.diffusingSpecies[species_name] = line_split[1]

			if(species_type == "BULK"):
				C.bulkSpecies[species_name] = line_split[1]

			if(species_type == "FILAMENT"):
				C.filamentSpecies[species_name] = line_split[1]

			if(species_type == "PLUSEND"):
				C.plusEndSpecies[species_name] = line_split[1]

			if(species_type == "MINUSEND"):
				C.minusEndSpecies[species_name] = line_split[1]

			if(species_type == "LINKER"):
				C.linkerSpecies[species_name] = line_split[1]

			if(species_type == "MOTOR"):
				C.motorSpecies[species_name] = line_split[1]

			if(species_type == "BRANCHER"):
				C.brancherSpecies[species_name] = line_split[1]

			line_number+=1

	print str(len(ChemSnapshotList)) + " snapshots read."

	return ChemSnapshotList


#Read a number of chem outputs, return array of chem lists
def readChemistryOutputs(fileNames):

	ChemSnapshotLists = []

	for chemFile in fileNames:
		print "Reading chemistry output..."
		ChemLists.append(readChemistryOutput(chemFile))

	return ChemSnapshotLists

########################################################
#
#
#			     ANALYSIS FUNCTIONS
#
#	contains: 
#		-average filament end to end distance
#		-average filament length
#		-minimum sphere enclosing volume
#		-Nematic order parameter
#	    -contractile velocity
#		-density of entire system
#		
#
########################################################


#Calculate average filament end to end distance at each step
#returns a list of time and distance pairs
def avgEndToEndDistance(SnapshotList):

	avgDists = []

	for Snapshot in SnapshotList:

		totalDist = 0;

		for filamentID, FS in Snapshot.filaments.iteritems():

			#get first and last bead
			b1 = FS.coords[0]
			b2 = FS.coords[len(FS.coords) - 1]

			#calculate distance
			dist = sqrt((b2 - b1).dot(b2 - b1))
			totalDist += dist

		#average
		totalDist /= len(Snapshot.filaments)
		totalDist *= 0.001 #to um

		#add to list
		avgDists.append((Snapshot.time, totalDist))

	return avgDists

#Calculate average filament end to end distance from a number of trajectories
#plot the data over time. If file is provided, save to that file
def avgEndToEndDistances(SnapshotLists, saveFile=''):

	avgDistsList = []

	#get all data
	for SnapshotList in SnapshotLists:

		avgDists = avgEndToEndDistance(SnapshotList)
		avgDistsList.append(avgDists)

	#now combine data
	finalAvgDists = []

	for i in xrange(0, len(avgDists)):

		finalTime = 0
		finalAvgDist = 0

		avgDistsSet = []

		for avgDists in avgDistsList:

			finalTime += avgDists[i][0] 
			finalAvgDist += avgDists[i][1]
			avgDistsSet.append(avgDists[i][1])

		#average
		finalTime /= len(SnapshotLists)
		finalAvgDist /= len(SnapshotLists)
		error = numpy.std(avgDistsSet)

		finalAvgDists.append((finalTime, finalAvgDist, error))

	#plot
	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	times  = [x[0] for x in finalAvgDists]
	dists  = [x[1] for x in finalAvgDists]
	errors = [x[2] for x in finalAvgDists]

	plt.errorbar(times, dists, yerr=errors)

	plt.title("Average end to end distance of filaments")
	plt.xlabel("Time (s)")
	plt.ylabel("End to end distance (um)")

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


#Calculate average filament length at each step
#returns a list of time and distance pairs
def avgLength(SnapshotList):

	avgLengths = []

	for Snapshot in SnapshotList:

		totalLength = 0;

		for filamentID, FS in Snapshot.filaments.iteritems():

			filLength = 0;

			for i in xrange(0, len(FS.coords) - 1):

				#get first and last bead
				b1 = FS.coords[i]
				b2 = FS.coords[i+1]

				#calculate distance
				filLength += sqrt((b2 - b1).dot(b2 - b1))
			
			totalLength += filLength

		#average
		totalLength /= len(Snapshot.filaments)

		#add to list
		avgLengths.append((Snapshot.time,totalLength))

	return avgLengths

#Calculate average lengths from a number of trajectories
#plot the data over time. If file is provided, save to that file
def avgLengths(SnapshotLists, saveFile=''):

	avgLengthsList = []

	#get all data
	for SnapshotList in SnapshotLists:

		avgLengths = avgLength(SnapshotList)
		avgLengthsList.append(avgLengths)

	#now combine data
	finalAvgLengths = []

	for i in xrange(0, len(avgLengths)):

		finalTime = 0
		finalAvgLength = 0

		avgLengthsSet = []

		for avgLengths in avgLengthsList:

			finalTime += avgLengths[i][0]
			finalAvgLength += avgLengths[i][1]

			avgLengthsSet.append(avgLengths[i][1])

		#average
		finalTime /= len(SnapshotLists)
		finalAvgLength /= len(SnapshotLists)

		error = numpy.std(avgLengthsSet)

		finalAvgLengths.append((finalTime, finalAvgLength, error))

	#plot
	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	times  = [x[0] for x in finalAvgLengths]
	lens   = [x[1] for x in finalAvgLengths]
	errors = [x[2] for x in finalAvgLengths]

	plt.errorbar(times, lens, yerr=errors)

	plt.title("Average length of filaments")
	plt.xlabel("Time (s)")
	plt.ylabel("Length (um)")

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)

#CREDIT TO: BlenderArtists.org

#minimum covering sphere candidate (mcsc)
#given candidate_points find minimum covering sphere
#find coordinates of point(candidate center) that is equidistant from candidate points
#if a coordinate is negative, point(candidate center) is outside of simplex,
#and the candidate point with the smallest value is removed

def mcsc(points, candidate_indices):

    candidate_count = len(candidate_indices)
    candidate_center = np.array([0.0,0.0,0.0])

    if candidate_count == 1:

        candidate_center = points[candidate_indices[0]]

    if candidate_count == 2:

        q0 = points[candidate_indices[0]]
        q1 = points[candidate_indices[1]]
        #l0 = l1 = 0.5
        #candidate_center then is midpoint formula
        candidate_center = (q0 + q1)/2.0

    if candidate_count == 3:

        q0 = points[candidate_indices[0]]
        q1 = points[candidate_indices[1]]
        q2 = points[candidate_indices[2]]

        a00 = (q0 - q2).dot(q0 - q2)
        a01 = (q0 - q2).dot(q1 - q2)
        a10 = (q1 - q2).dot(q0 - q2)
        a11 = (q1 - q2).dot(q1 - q2)
  
        A = np.matrix([[a00,a01],[a10,a11]])
        b = np.array([a00/2.0,a11/2.0])

        if A.determinant() == 0:

            candidate_indices.pop()
            candidate_center, candidate_indices = mcsc(points, candidate_indices)

        else:
            l = b*A.inverted()
  
            l_sum = l[0] + l[1]
            l2 = 1.0 - l_sum
        
            l_list = [l[0], l[1], l2]
  
            drop_index = None
            minimum = 0
  
            for index, number in enumerate(l_list):
                if number < minimum:
                    drop_index = index
                    minimum = number
     
            if drop_index != None:
                candidate_indices.pop(drop_index)
                candidate_center, candidate_indices = mcsc(points, candidate_indices)
            else:
                candidate_center = l[0]*q0 + l[1]*q1 + l2*q2  
 
    if candidate_count == 4:

        q0 = points[candidate_indices[0]]
        q1 = points[candidate_indices[1]]
        q2 = points[candidate_indices[2]]
        q3 = points[candidate_indices[3]]

        a00 = (q0 - q3).dot(q0 - q3)
        a01 = (q0 - q3).dot(q1 - q3)
        a02 = (q0 - q3).dot(q2 - q3)
        a10 = (q1 - q3).dot(q0 - q3)
        a11 = (q1 - q3).dot(q1 - q3)
        a12 = (q1 - q3).dot(q2 - q3)
        a20 = (q2 - q3).dot(q0 - q3)
        a21 = (q2 - q3).dot(q1 - q3)
        a22 = (q2 - q3).dot(q2 - q3)
  
        A = np.matrix([[a00, a01, a02],[a10, a11, a12],[a20, a21, a22]])
        b = np.array([a00/2.0, a11/2.0, a22/2.0])

        if A.determinant() == 0:

            candidate_indices.pop()
            candidate_center, candidate_indices = mcsc(points, candidate_indices)

        else:
            l = b*A.inverted()
  
            l_sum = l[0] + l[1] + l[2]
            l3 = 1.0 - l_sum
        
            l_list = [l[0], l[1], l[2], l3]
  
            drop_index = None
            minimum = 0
  
            for index, number in enumerate(l_list):
                if number < minimum:
                    drop_index = index
                    minimum = number
     
            if drop_index != None:
                candidate_indices.pop(drop_index)
                candidate_center, candidate_indices = mcsc(points, candidate_indices)

            else:
                candidate_center = l[0]*q0 + l[1]*q1 + l[2]*q2 + l3*q3  
        
    return candidate_center, candidate_indices

#reduce covering sphere along line p*(t - center)
#new candidates have p values between 0 exclusive and 1 exclusive
#candidate with smallest value is choosen for next covering sphere
#if minimum_covering sphere_candidate returns 4 points then minimum found
#if no new candidate points are found then minimum is found
#if new_center == center minimum is found

def find_next_candidate(points, center, candidate_indices):

    t = 0.0
    t, candidate_indices = mcsc(points, candidate_indices)

    if len(candidate_indices) == 4:
        return True, t, candidate_indices

    p = []
    for point in points:
        p.append(1)
  
    for index, point in enumerate(points):
        if index not in candidate_indices:
            d = (t - center).dot(points[candidate_indices[0]] - point)
            if d > 0:
                p[index] = (((points[candidate_indices[0]] + point)/2.0 - center)).dot( \
                            ((points[candidate_indices[0]] - point)/d))
    minimum = 1
    min_index = None
    for index, number in enumerate(p):
        if number < minimum and number > 0:
            min_index = index
            minimum = number
 
    if minimum == 1:
        return True, t, candidate_indices
 
    new_center = center + p[min_index]*(t - center)
 
    candidate_indices.insert(0, min_index) 
 
    if new_center.all() == center.all():
        return True, new_center, candidate_indices

    return False, new_center, candidate_indices

#given set of points defined
#return center and radius of minimum covering sphere
#define starting center as first point in points
#define p_1 as point farthest away from center
#first set of candidate points contains only p_1
#run until no new candidate points are found or until
#four candidate points are found that lie on minimum covering sphere

def minimum_covering_sphere(points):
        
    point_0 = points[0]
    new_center = point_0
    
    max_d = 0.0
    point_1_index = 0
    for index, point in enumerate(points):
        d = math.sqrt((point_0 - point).dot(point_0 - point))
        if d > max_d:
            point_1_index = index
            max_d = d
   
    candidate_indices = [point_1_index]
    finished = False

    while not finished:
        center = new_center
        finished, new_center, candidate_indices = \
        	find_next_candidate(points, center, candidate_indices)

    diff = points[candidate_indices[0]] - new_center
    radius = math.sqrt(diff.dot(diff))

    return new_center, radius


#calculate the minimum enclosing sphere volume of all elements (in um^3)
#uses the end points of filaments
#minimum covering sphere based on www.mel.nist.gov/msidlibrary/doc/hopp95.pdf
#returns a list of the time and volume pairs
def minEnclosingVolume(SnapshotList):

	minEnclosingVols = []

	for i in xrange(2, len(SnapshotList), 2):
		Snapshot = SnapshotList[i]

		points = []

		for filamentID, FS in Snapshot.filaments.iteritems():

			points.append(FS.coords[0])
			points.append(FS.coords[len(FS.coords) - 1])

		center, radius = minimum_covering_sphere(points)

		#add to list
		minEnclosingVols.append((Snapshot.time, (4/3) * math.pi * pow(radius * 0.001,3)))

	return minEnclosingVols


#Calculate minimum enclosing sphere volumes from a number of trajectories (in um^3)
#uses the end points of filaments
#plot the data over time. If file is provided, save to that file
#minimum covering sphere based on www.mel.nist.gov/msidlibrary/doc/hopp95.pdf
def minEnclosingVolumes(SnapshotLists, saveFile=''):

	minEnclosingVolumesList = []

	#get all data
	for SnapshotList in SnapshotLists:

		minEnclosingVolumes = minEnclosingVolume(SnapshotList)
		minEnclosingVolumesList.append(minEnclosingVolumes)

	#now combine data
	finalMinEnclosingVolumes = []

	for i in xrange(0, len(minEnclosingVolumes)):

		finalTime = 0
		finalMinEnclosingVolume = 0

		minEnclosingVolumesSet = []

		for minEnclosingVolumes in minEnclosingVolumesList:

			finalTime += minEnclosingVolumes[i][0] 
			finalMinEnclosingVolume += minEnclosingVolumes[i][1]

			minEnclosingVolumesSet.append(minEnclosingVolumes[i][1])

		#average
		finalTime /= len(SnapshotLists)
		finalMinEnclosingVolume /= len(SnapshotLists)

		error = numpy.std(minEnclosingVolumesSet)

		finalMinEnclosingVolumes.append((finalTime, finalMinEnclosingVolume, error))

	#plot
	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	times  = [x[0] for x in finalMinEnclosingVolumes]
	vols   = [x[1] for x in finalMinEnclosingVolumes]
	#errors = [x[2] for x in finalMinEnclosingVolumes]

	plt.errorbar(times, vols)

	plt.title(r'$\mathrm{Minimum\/spherical\/enclosing\/volume}\/(um^3)$')
	plt.xlabel(r'$\mathrm{Time}\/(s)$')
	plt.ylabel(r'$\mathrm{Volume}\/(um^3)$')

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


#Calculate nematic order parameter of cylinders at each timestep
#from order parameter definition at http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
#takes the largest eigenvalue of Q, a second rank ordering tensor over cylinder pairs
#returns a list of time and value pairs
def orderParameter(SnapshotList, snapshot=1):

	orderParams = []

	Snapshot = SnapshotList[snapshot]

	Q = np.zeros((len(Snapshot.filaments), len(Snapshot.filaments)))

	i = 0
	for filIDA, FSA in Snapshot.filaments.iteritems():

		#calculate direction
		magA = sqrt((FSA.coords[len(FSA.coords) - 1] - FSA.coords[0]). \
					 dot(FSA.coords[len(FSA.coords) - 1] - FSA.coords[0]))
		directionA = (FSA.coords[len(FSA.coords) - 1] - FSA.coords[0]) / magA

		j = 0
		for filIDB, FSB in Snapshot.filaments.iteritems():

			#calculate direction
			magB = sqrt((FSB.coords[len(FSB.coords) - 1] - FSB.coords[0]). \
						 dot(FSB.coords[len(FSB.coords) - 1] - FSB.coords[0]))
			directionB = (FSB.coords[len(FSB.coords) - 1] - FSB.coords[0]) / magB

			#pairwise delta
			if filIDA == filIDB:
				delta = 1
			else:
				delta = 0

			#q value
			Q[i][j] = (1 / float(len(Snapshot.filaments))) * \
					  (1.5 * np.dot(directionA, directionB) - 0.5 * delta)

			j+=1

		i+=1

	#take largest eigenvalue
	largestEigenvalue, largestEigenVector = la.eigs(Q,k=1)
	#return
	return (Snapshot.time,largestEigenvalue[0])


##Calculate nematic order parameter of cylinders from a number of trajectories
#from order parameter definition at http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
#takes the largest eigenvalue of Q, a second rank ordering tensor over cylinder pairs
#plot the data over time. If file is provided, save to that file
def orderParameters(SnapshotLists, snapshot=1):

	orderParamsList = []

	#get all data
	for SnapshotList in SnapshotLists:

		orderParams = orderParameter(SnapshotList, snapshot)
		orderParamsList.append(orderParams)

	finalTime = 0
	finalOrderParam = 0

	orderParamSet = []

	for orderParam in orderParamsList:

		finalTime += orderParam[0]
		finalOrderParam += orderParam[1]

		orderParamSet.append(orderParam[1])

	#average
	finalTime /= len(orderParamsList)
	finalOrderParam /= len(orderParamsList)

	error = numpy.std(orderParamSet)

	return (finalTime,finalOrderParam,error)

#Calculate the radius of gyration over all cylinders 
#returns a time and Rg pair for the chosen snapshot
def radiusOfGyration(SnapshotList, snapshot=1):

	Snapshot = SnapshotList[snapshot]

	#get center of mass of Snapshot
	com = np.array([0.0,0.0,0.0])
	beadCount = 0

		#go through filaments
	for filamentID, FS in Snapshot.filaments.iteritems():
		for coord in FS.coords:

			com += coord
			beadCount += 1

	#normalize
	com /= beadCount

	RgSquare = 0
	numCylinders = 0

	#Loop through cylinders
	for cylinderID, CS in Snapshot.cylinders.iteritems():

		numCylinders += 1
		cam = (CS.coords[1] + CS.coords[0]) / 2

		RgSquare += (cam - com).dot(cam - com)

	RgSquare /= numCylinders
	Rg = sqrt(RgSquare * 0.00001)

	return (Snapshot.time, Rg)


#Calculate the radius of gyration over all cylinders from a number of trajectories
def radiusOfGyrations(SnapshotLists, snapshot=1):

	RgList = []

	#get all data
	for SnapshotList in SnapshotLists:

		Rg = radiusOfGyration(SnapshotList, snapshot)
		RgList.append(Rg)

	#now combine data
	finalTime = 0
	finalRg = 0

	for Rg in RgList:

		finalTime += Rg[0]
		finalRg += Rg[1]

	#average
	finalTime /= len(RgList)
	finalRg /= len(RgList)

	Rgs = [x[1] for x in RgList]

	perc_25 = np.percentile(Rgs, 25)
	perc_75 = np.percentile(Rgs, 75)

	return (finalTime,finalRg,perc_25,perc_75)

#Calculate contractile velocity towards the center of mass
#returns a list of time and velocity pairs
def contractileVelocity(SnapshotList):

	velocities = []

	for SnapshotIndex in xrange(0, len(SnapshotList) - 1):

		SnapshotI = SnapshotList[SnapshotIndex]
		SnapshotF = SnapshotList[SnapshotIndex + 1]

		deltaT = SnapshotF.time - SnapshotI.time

		#get center of mass of Snapshot i
		com = np.array([0.0,0.0,0.0])
		beadCount = 0

		#go through filaments
		for filamentID, FS in SnapshotI.filaments.iteritems():
			for coord in FS.coords:

				com += coord
				beadCount += 1

		#normalize
		com /= beadCount

		velocity = 0
		numCylinders = 0
		#Loop through cylinders
		for cylinderIDA, CSA in SnapshotF.cylinders.iteritems():

			#find in initial. if it doesnt exist, leave out
			try:
				CSB = SnapshotI.cylinders[cylinderIDA]
			except Exception, e:
				continue

			numCylinders += 1

			#calc velocity of midpoints of cylinders
			cam = (CSA.coords[1] + CSA.coords[0]) / 2
			cbm = (CSB.coords[1] + CSB.coords[0]) / 2

			vel = (cam - cbm)/ deltaT

			velocity += vel.dot((com - cbm)/((com - cbm).dot(com - cbm)))

		#normalize
		if numCylinders == 0: continue

		velocity /= numCylinders
		velocity *= 0.001 #units

		velocities.append((SnapshotF.time, velocity))

	#return velocities
	return velocities

#Calculate contractile velocity towards the center of mass from a number of trajectories
#plot the data over time. If file is provided, save to that file
def contractileVelocities(SnapshotLists, saveFile=''):

	velocitiesList = []

	#get all data
	for SnapshotList in SnapshotLists:

		velocities = contractileVelocity(SnapshotList)
		velocitiesList.append(velocities)

	#now combine data
	finalVelocities = []

	for i in xrange(0, len(velocities)):

		finalTime = 0
		finalVelocity = 0

		velocitySet = []

		for velocities in velocitiesList:

			finalTime += velocities[i][0]
			finalVelocity += velocities[i][1]

			velocitySet.append(velocities[i][1])

		#average
		finalTime /= len(SnapshotLists)
		finalVelocity /= len(SnapshotLists)

		error = numpy.std(velocitySet)

		finalVelocities.append((finalTime, finalVelocity, error))

	#plot
	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	times  = [x[0] for x in finalVelocities]
	vels   = [x[1] for x in finalVelocities]
	errors = [x[2] for x in finalVelocities]

	plt.errorbar(times, vels, yerr=errors)

	plt.title("Contractile velocity of filaments")
	plt.xlabel("Time (s)")
	plt.ylabel("Velocity (um / s)")

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)

#finds index of coordinate by linear mapping
def getIndex(coordinate, grid, compartment):

    index = 0
    i = 0

    for x in coordinate:
        if i == 0 :
            index += int(x / compartment[0])
        elif i == 1:
            index += int(x / compartment[1]) * grid[0];
        else:
            index += int(x / compartment[2]) * grid[0] * grid[1];
        i += 1

    return index


#Calculates density of filaments in each compartment for a Snapshot
#Uses number of cylinders in each compartment
#Param grid is the size of the grid (NX, NY, NZ)
#Param compartment is the size of compartment (nx, ny, nz)
#So, the entire grid is NX*nx * NY*ny * NZ*nz 
def density(Snapshot, grid, compartment):

	#create dict of compartment index and num cylinders by dimensions
	numCylinders = {}
	#create dict of compartment index and compartment
	compartmentMap = {}

	for i in xrange(0, grid[0]):

		for j in xrange(0, grid[1]):

			for k in xrange(0, grid[2]):

				coordinate = [i * compartment[0], \
							  j * compartment[1], \
							  k * compartment[2]]

				#add to map
				index = getIndex(coordinate, grid, compartment)
				compartmentMap[index] = coordinate

				#init num cylinders
				numCylinders[index] = 0

	#now, loop through cylinders, populate numCylinders
	for cylinderID, CS in Snapshot.cylinders.iteritems():

		#get coordinate, find index
		coordinate = (CS.coords[1] + CS.coords[0]) / 2
		index = getIndex(coordinate, grid, compartment)
	
		#add length of cylinder in um
		numCylinders[index] += 0.001 * sqrt((CS.coords[1] - CS.coords[0]). \
										 dot(CS.coords[1] - CS.coords[0]))

	#make final mapping of coordinates to num cylinders
	compartmentCoords = []
	densityValues = []

	for index, coordinate in compartmentMap.iteritems():

		compartmentCoords.append(compartmentMap[index])
		densityValues.append(numCylinders[index])

	#Plot
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	xVals = [x[0] + compartment[0] / 2 for x in compartmentCoords]
	yVals = [x[1] + compartment[1] / 2 for x in compartmentCoords]
	zVals = [x[2] + compartment[2] / 2 for x in compartmentCoords]

	#plt.xlim([0, compartment[]])

	plt.scatter(xVals, yVals, zs=zVals, c=densityValues)


###############################################
#
#
#
#				PLOS FIGURES
#
#
#
###############################################

def calculateRgs(snapshot=1):

	FrameLists = readTrajectories([])


	return radiusOfGyrations(FrameLists, snapshot)

def calculateRgVsT():

	FrameLists = readTrajectories([])

	Rgs1 = []

	for snap in xrange(2, len(FrameLists[0]), 2):
		Rgs1.append(radiusOfGyrations(FrameLists, snapshot=snap))

	del FrameLists

	FrameLists = readTrajectories([])

	Rgs2 = []

	for snap in xrange(2, len(FrameLists[0]), 2):
		Rgs2.append(radiusOfGyrations(FrameLists, snapshot=snap))

	del FrameLists

	FrameLists = readTrajectories([])

	Rgs3 = []

	for snap in xrange(2, len(FrameLists[0]), 2):
		Rgs3.append(radiusOfGyrations(FrameLists, snapshot=snap))

	del FrameLists

	FrameLists = readTrajectories([])

	Rgs4 = []

	for snap in xrange(2, len(FrameLists[0]), 2):
		Rgs4.append(radiusOfGyrations(FrameLists, snapshot=snap))

	del FrameLists

	FrameLists = readTrajectories([])

	Rgs5 = []

	for snap in xrange(2, len(FrameLists[0]), 2):
		Rgs5.append(radiusOfGyrations(FrameLists, snapshot=snap))

	del FrameLists


	FrameLists = readTrajectories([])

	Rgs6 = []

	for snap in xrange(2, len(FrameLists[0]), 2):
		Rgs6.append(radiusOfGyrations(FrameLists, snapshot=snap))

	del FrameLists

	return [Rgs1, Rgs2, Rgs3, Rgs4, Rgs5, Rgs6]


def plotRgHeatMap(saveFile=''):

	#read data
	x_labels = [0.01,0.02,0.05,0.1,0.2,0.5]
	y_labels = [0.02, 0.01, 0.005]

	data = np.array([[1.460,1.392,1.199,1.144,1.105,1.091],
				     [1.600,1.518,1.387,1.332,1.305,1.268],
			         [1.713,1.652,1.512,1.437,1.449,1.410]])

	fig, ax = plt.subplots()
	heatmap = ax.pcolor(data, cmap=plt.cm.YlGnBu_r)

	cbar = plt.colorbar(heatmap)

	cbar.ax.get_yaxis().set_ticks([])
	for j, lab in enumerate(['$1.0$','$1.2$','$1.4$','$1.6$']):
   		 cbar.ax.text(1.2, (2 * j + 1) / 8.0, lab, ha='left', va='center')
	cbar.ax.get_yaxis().labelpad = 15
	cbar.ax.set_ylabel(r'$\mathrm{R_g\/(\mu m)}$', fontsize=16)


	# put the major ticks at the middle of each cell
	ax.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
	ax.set_yticks(np.arange(data.shape[0])+0.5, minor=False)

	# want a more natural, table-like display
	ax.invert_yaxis()

	ax.set_xticklabels(x_labels, minor=False)
	ax.set_yticklabels(y_labels, minor=False)
	plt.show()

	plt.ylabel(r'$\mathrm{R_{m:a}}$', fontsize=20)
	plt.xlabel(r'$\mathrm{R_{\alpha:a}}$', fontsize=20)

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


def plotRgVsT(data, saveFile=''):

	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	#organize data
	data1 = data[0]
	data2 = data[1]
	data3 = data[2]
	data4 = data[3]
	data5 = data[4]
	data6 = data[5]

	time1 = [x[0] for x in data1]
	time2 = [x[0] for x in data2]
	time3 = [x[0] for x in data3]
	time4 = [x[0] for x in data4]
	time5 = [x[0] for x in data5]
	time6 = [x[0] for x in data6]

	Rg1_0 = data1[0][1]
	Rg2_0 = data2[0][1]
	Rg3_0 = data3[0][1]
	Rg4_0 = data4[0][1]
	Rg5_0 = data5[0][1]
	Rg6_0 = data6[0][1]

	Rg1 = [x[1]/Rg1_0 for x in data1]
	Rg2 = [x[1]/Rg2_0 for x in data2]
	Rg3 = [x[1]/Rg3_0 for x in data3]
	Rg4 = [x[1]/Rg4_0 for x in data4]
	Rg5 = [x[1]/Rg5_0 for x in data5]
	Rg6 = [x[1]/Rg6_0 for x in data6]
	
	perc_25_1 = [x[2]/Rg1_0 for x in data1]
	perc_25_2 = [x[2]/Rg2_0 for x in data2]
	perc_25_3 = [x[2]/Rg3_0 for x in data3]
	perc_25_4 = [x[2]/Rg4_0 for x in data4]
	perc_25_5 = [x[2]/Rg5_0 for x in data5]
	perc_25_6 = [x[2]/Rg6_0 for x in data6]

	perc_75_1 = [x[3]/Rg1_0 for x in data1]
	perc_75_2 = [x[3]/Rg2_0 for x in data2]
	perc_75_3 = [x[3]/Rg3_0 for x in data3]
	perc_75_4 = [x[3]/Rg4_0 for x in data4]
	perc_75_5 = [x[3]/Rg5_0 for x in data5]
	perc_75_6 = [x[3]/Rg6_0 for x in data6]

	plt.plot(time1, Rg1, 'r', label=r'$\mathrm{R_{\alpha:a}=0.01}$')
	plt.plot(time2, Rg2, 'g', label=r'$\mathrm{R_{\alpha:a}=0.02}$')
	plt.plot(time3, Rg3, 'b', label=r'$\mathrm{R_{\alpha:a}=0.05}$')
	plt.plot(time4, Rg4, 'y', label=r'$\mathrm{R_{\alpha:a}=0.1}$')
	plt.plot(time5, Rg5, 'm', label=r'$\mathrm{R_{\alpha:a}=0.2}$')
	plt.plot(time6, Rg6, 'c', label=r'$\mathrm{R_{\alpha:a}=0.5}$')

	fill_between(time1, perc_25_1, perc_75_1, alpha=0.25, linewidth=0, color='r')
	fill_between(time2, perc_25_2, perc_75_2, alpha=0.25, linewidth=0, color='g')
	fill_between(time3, perc_25_3, perc_75_3, alpha=0.25, linewidth=0, color='b')
	fill_between(time4, perc_25_4, perc_75_4, alpha=0.25, linewidth=0, color='y')
	fill_between(time5, perc_25_5, perc_75_5, alpha=0.25, linewidth=0, color='m')
	fill_between(time6, perc_25_6, perc_75_6, alpha=0.25, linewidth=0, color='c')

	plt.xlabel(r'$\mathrm{Time\/(s)}$', fontsize=20)
	plt.ylabel(r'$\mathrm{R_{g,f} / R_{g,i}\/(\mu m)}$', fontsize=20)
	plt.xlim(0,2000)

	plt.legend()

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)






