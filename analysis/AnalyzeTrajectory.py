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

#A complete frame
class Frame:
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

########################################################
#
#			         PARSING FUNCTIONS
#
########################################################


#Read all data from a trajectory file
#returns a list of frames with all data
def readTrajectory(filename):

	#Open the traj file
	traj_file=open(filename)

	FrameList=[]
	first_frame_line=True
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

		#get frame data

		if(first_frame_line):
			F=Frame()
			F.step, F.time, F.n_filaments, F.n_linkers, \
			F.n_motors, F.n_branchers = map(double,line.split())
			F.n_filaments = int(F.n_filaments)
			F.n_linkers = int(F.n_linkers)
			F.n_motors = int(F.n_motors)
			F.n_branchers = int(F.n_branchers)
			first_frame_line=False
			line_number+=1
			continue

		if(len(line)==0):
			first_frame_line=True
			first_filament_line=True
			FrameList.append(F)
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

	#return the final frame list
	return FrameList	

#Read a number of trajectories, return array of frame lists
def readTrajectories(fileNames):

	FrameLists = []

	for trajFile in fileNames:
		FrameLists.append(readTrajectory(trajFile))

	return FrameLists


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
def avgEndToEndDistance(FrameList):

	avgDists = []

	for frame in FrameList:

		totalDist = 0;

		for filamentID, FS in frame.filaments.iteritems():

			#get first and last bead
			b1 = FS.coords[0]
			b2 = FS.coords[len(FS.coords) - 1]

			#calculate distance
			dist = sqrt((b2 - b1).dot(b2 - b1))
			totalDist += dist

		#average
		totalDist /= len(frame.filaments)
		totalDist *= 0.001 #to um

		#add to list
		avgDists.append((frame.time, totalDist))

	return avgDists

#Calculate average filament end to end distance from a number of trajectories
#plot the data over time. If file is provided, save to that file
def avgEndToEndDistances(FrameLists, saveFile=''):

	avgDistsList = []

	#get all data
	for FrameList in FrameLists:

		avgDists = avgEndToEndDistance(FrameList)
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
		finalTime /= len(FrameLists)
		finalAvgDist /= len(FrameLists)
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
def avgLength(FrameList):

	avgLengths = []

	for frame in FrameList:

		totalLength = 0;

		for filamentID, FS in frame.filaments.iteritems():

			filLength = 0;

			for i in xrange(0, len(FS.coords) - 1):

				#get first and last bead
				b1 = FS.coords[i]
				b2 = FS.coords[i+1]

				#calculate distance
				filLength += sqrt((b2 - b1).dot(b2 - b1))
			
			totalLength += filLength

		#average
		totalLength /= len(frame.filaments)

		#add to list
		avgLengths.append((frame.time,totalLength))

	return avgLengths

#Calculate average lengths from a number of trajectories
#plot the data over time. If file is provided, save to that file
def avgLengths(FrameLists, saveFile=''):

	avgLengthsList = []

	#get all data
	for FrameList in FrameLists:

		avgLengths = avgLength(FrameList)
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
		finalTime /= len(FrameLists)
		finalAvgLength /= len(FrameLists)

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
def minEnclosingVolume(FrameList):

	minEnclosingVols = []

	for frame in FrameList:

		points = []

		for filamentID, FS in frame.filaments.iteritems():

			points.append(FS.coords[0])
			points.append(FS.coords[len(FS.coords) - 1])

		center, radius = minimum_covering_sphere(points)

		#add to list
		minEnclosingVols.append((frame.time, (4/3) * math.pi * pow(radius * 0.001,3)))

	return minEnclosingVols


#Calculate minimum enclosing sphere volumes from a number of trajectories (in um^3)
#uses the end points of filaments
#plot the data over time. If file is provided, save to that file
#minimum covering sphere based on www.mel.nist.gov/msidlibrary/doc/hopp95.pdf
def minEnclosingVolumes(FrameLists, saveFile=''):

	minEnclosingVolumesList = []

	#get all data
	for FrameList in FrameLists:

		minEnclosingVolumes = minEnclosingVolume(FrameList)
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
		finalTime /= len(FrameLists)
		finalMinEnclosingVolume /= len(FrameLists)

		error = numpy.std(minEnclosingVolumesSet)

		finalMinEnclosingVolumes.append((finalTime, finalMinEnclosingVolume, error))

	#plot
	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	times  = [x[0] for x in finalMinEnclosingVolumes]
	vols   = [x[1] for x in finalMinEnclosingVolumes]
	errors = [x[2] for x in finalMinEnclosingVolumes]

	plt.errorbar(times, vols, yerr=errors)

	plt.title("Minimum spherical enclosing volume (um^3)")
	plt.xlabel("Time (s)")
	plt.ylabel("Volume (um^3)")

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


#Calculate nematic order parameter of cylinders at each timestep
#from order parameter definition at http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
#takes the largest eigenvalue of Q, a second rank ordering tensor over cylinder pairs
#returns a list of time and value pairs
def orderParameter(FrameList):

	orderParams = []

	for frame in FrameList:

		if len(frame.cylinders) <= 1: 

			orderParams.append((frame.time, 1.0))
			continue

		Q = np.zeros((len(frame.cylinders), len(frame.cylinders)))

		i = 0
		for cylinderIDA, CSA in frame.cylinders.iteritems():

			#calculate direction
			magA = sqrt((CSA.coords[1] - CSA.coords[0]).dot(CSA.coords[1] - CSA.coords[0]))
			directionA = (CSA.coords[1] - CSA.coords[0]) / magA

			j = 0
			for cylinderIDB, CSB in frame.cylinders.iteritems():

				#calculate direction
				magB = sqrt((CSB.coords[1] - CSB.coords[0]).dot(CSB.coords[1] - CSB.coords[0]))
				directionB = (CSB.coords[1] - CSB.coords[0]) / magB

				#pairwise delta
				if cylinderIDA == cylinderIDB:
					delta = 1
				else:
					delta = 0

				#q value
				Q[i][j] = (1 / float(len(frame.cylinders))) * \
						  (1.5 * np.dot(directionA, directionB) - 0.5 * delta)

				j+=1

			i+=1

		#take largest eigenvalue
		largestEigenvalue, largestEigenVector = la.eigs(Q,k=1)
		#add to list
		orderParams.append((frame.time,largestEigenvalue[0]))

	return orderParams 


##Calculate nematic order parameter of cylinders from a number of trajectories
#from order parameter definition at http://cmt.dur.ac.uk/sjc/thesis_dlc/node65.html
#takes the largest eigenvalue of Q, a second rank ordering tensor over cylinder pairs
#plot the data over time. If file is provided, save to that file
def orderParameters(FrameLists, saveFile=''):

	orderParamsList = []

	#get all data
	for FrameList in FrameLists:

		orderParams = orderParameter(FrameList)
		orderParamsList.append(orderParams)

	#now combine data
	finalOrderParams = []

	for i in xrange(0, len(orderParams)):

		finalTime = 0
		finalOrderParam = 0

		orderParamSet = []

		for orderParams in orderParamsList:

			finalTime += orderParams[i][0]
			finalOrderParam += orderParams[i][1]

			orderParamSet.append(orderParams[i][1])

		#average
		finalTime /= len(FrameLists)
		finalOrderParam /= len(FrameLists)

		error = numpy.std(orderParamSet)

		finalOrderParams.append((finalTime,finalOrderParam,error))

	#plot
	fig = plt.figure(figsize=(10.0,10.0))
	host = fig.add_subplot(111)

	times  = [x[0] for x in finalOrderParams]
	vols   = [x[1] for x in finalOrderParams]
	errors = [x[2] for x in finalOrderParams]

	plt.errorbar(times, vols, yerr=errors)

	plt.title("Nematic order parameter of filaments")
	plt.xlabel("Time (s)")
	plt.ylabel("Order parameter")

	#if file provided, save
	if saveFile != '':
		fig.savefig(saveFile)


#Calculate contractile velocity towards the center of mass
#returns a list of time and velocity pairs
def contractileVelocity(FrameList):

	velocities = []

	for frameIndex in xrange(0, len(FrameList) - 1):

		frameI = FrameList[frameIndex]
		frameF = FrameList[frameIndex + 1]

		deltaT = frameF.time - frameI.time

		#get center of mass of frame i
		com = [0.0,0.0,0.0]
		beadCount = 0

		#go through filaments
		for filamentID, FS in frameI.filaments.iteritems():
			for coord in FS.coords:

				com += coord
				beadCount += 1

		#normalize
		com /= beadCount

		velocity = 0
		numCylinders = 0
		#Loop through cylinders
		for cylinderIDA, CSA in frameF.cylinders.iteritems():

			#find in initial. if it doesnt exist, leave out
			try:
				CSB = frameI.cylinders[cylinderIDA]
			except Exception, e:
				continue

			numCylinders += 1

			#calc velocity of midpoints of cylinders
			cam = (CSA.coords[1] + CSA.coords[0]) / 2
			cbm = (CSB.coords[1] + CSB.coords[0]) / 2

			vel = (cam - cbm)/ deltaT

			velocity += vel.dot((com - cbm)/((com - cbm).dot(com - cbm)))

		#normalize
		velocity /= numCylinders
		velocity *= 0.001 #units

		velocities.append((frameF.time, velocity))

	#return velocities
	return velocities

#Calculate contractile velocity towards the center of mass from a number of trajectories
#plot the data over time. If file is provided, save to that file
def contractileVelocities(FrameLists, saveFile=''):

	velocitiesList = []

	#get all data
	for FrameList in FrameLists:

		velocities = contractileVelocity(FrameList)
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
		finalTime /= len(FrameLists)
		finalVelocity /= len(FrameLists)

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


#Calculates density of filaments in each compartment for a frame
#Uses number of cylinders in each compartment
#Param grid is the size of the grid (NX, NY, NZ)
#Param compartment is the size of compartment (nx, ny, nz)
#So, the entire grid is NX*nx * NY*ny * NZ*nz 
def density(Frame, grid, compartment):

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
	for cylinderID, CS in Frame.cylinders.iteritems():

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

