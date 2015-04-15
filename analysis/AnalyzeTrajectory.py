from numpy import *

#A filament snapshot
class FilamentSnapshot:
	def __init__(self):
		self.id=None
		self.length=None
		self.delta_left=None	
		self.delta_right=None	
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
		self.linkers={}
		self.motors={}
		self.branchers={}


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
			FS.coords=FS.coords.reshape(N/3,3).T
			first_line=True
			reading_filament=False

		if(reading_linker):
			LS.coords=array(line.split(),'f')
			N=len(LS.coords)
			LS.coords=LS.coords.reshape(N/3,3).T
			first_line=True
			reading_linker=False

		if(reading_motor):
			MS.coords=array(line.split(),'f')
			N=len(MS.coords)
			MS.coords=MS.coords.reshape(N/3,3).T
			first_line=True
			reading_motor=False

		if(reading_brancher):
			BS.coords=array(line.split(),'f')
			N=len(BS.coords)
			BS.coords=BS.coords.reshape(N/3,3).T
			first_line=True
			reading_brancher=False

		line_number+=1

	#return the final frame list
	return FrameList	









