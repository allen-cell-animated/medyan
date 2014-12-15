from numpy import *
from mayavi import mlab

mlab.figure(1, size=(600, 600), bgcolor=(0, 0, 0))
mlab.show()

filename='/Users/jameskomianos/Code/M3SYM/M3SYM/filamentoutput.txt'
traj_file=open(filename)

class FilamentSnapshot:
	def __init__(self):
		self.id=None
		self.length=None
		self.delta_left=None	
		self.delta_right=None	
		self.coords=None
		self.connections=None

class LinkerSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.coords=None
		self.connections=None

class MotorSnapshot:
	def __init__(self):
		self.id=None
		self.type=None;
		self.coords=None
		self.connections=None

class Frame:
	def __init__(self):
		self.step=None
		self.time=None
		self.n_filaments=None
		self.n_linkers=None
		self.n_motors=None
		self.filaments={}
		self.linkers={}
		self.motors={}

FrameList=[]
first_frame_line=True
first_line=True
reading_filament=False
reading_linker=False
reading_motor=False
n_beads_filament=0
n_beads_linker=0
n_beads_motor=0

for line in traj_file:
	line = line.strip()
	#print "["+line+"]"
	if(first_frame_line):
		F=Frame()
		F.step, F.time, F.n_filaments, F.n_linkers, F.n_motors = map(double,line.split())
		F.n_filaments = int(F.n_filaments)
		F.n_linkers = int(F.n_linkers)
		F.n_motors = int(F.n_motors)
		first_frame_line=False
		continue

	if(len(line)==0):
		#print "Detected an empty line"
		first_frame_line=True
		first_filament_line=True
		FrameList.append(F)
		assert F.n_filaments == len(F.filaments)
		assert F.n_linkers == len(F.linkers)
		assert F.n_motors == len(F.motors)
		n_beads_filament=0
		n_beads_linker=0
		n_beads_motor=0
		continue
			
	if(first_line):
		line_split = line.split()
		if(line[0] == "F"):
			FS=FilamentSnapshot()
			FS.id, FS.length, FS.delta_left, FS.delta_right = map(int,line_split[1:])
			first_line=False
			F.filaments[FS.id]=FS
			FS.connections=vstack([arange(n_beads_filament, n_beads_filament + FS.length - 1.5), 
				arange(n_beads_filament + 1, n_beads_filament + FS.length - .5)]).T
			n_beads_filament+=FS.length
			reading_filament=True
			continue
		if(line[0] == "L"):
			LS = LinkerSnapshot()
			LS.id, LS.type = map(int, line_split[1:])
			first_line = False
			F.linkers[LS.id] = LS
			LS.connections=vstack([arange(n_beads_linker, n_beads_linker + 0.5), 
				arange(n_beads_linker + 1, n_beads_linker + 1.5)]).T
			n_beads_linker+=2
			reading_linker=True
			continue
		if(line[0] == "M"):
			MS = MotorSnapshot()
			MS.id, MS.type = map(int, line_split[1:])
			first_line = False
			F.motors[MS.id] = MS
			MS.connections=vstack([arange(n_beads_motor, n_beads_motor + 0.5), 
				arange(n_beads_motor + 1, n_beads_motor + 1.5)]).T
			n_beads_motor+=2
			reading_motor=True
			continue

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

@mlab.show
def show_frame(frame_number=-1):
	local_frame=FrameList[frame_number]
	mlab.figure(1, size=(600, 600), bgcolor=(0, 0, 0))
	mlab.clf()

	#DISPLAYING FILAMENTS
	if(len(local_frame.filaments) != 0):
		x=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for fid in sorted(local_frame.filaments.keys()):
				q.append(local_frame.filaments[fid].coords[i])
			x.append(hstack(q))

		for fid in sorted(local_frame.filaments.keys()):
				connections.append(local_frame.filaments[fid].connections)

		connections = vstack(connections)

		# Create the points
		src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])
		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24)

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=0.25)
		tube.filter.number_of_sides=12
		mlab.pipeline.surface(tube, color=(0.8, 0.3, 0.3))

	#DISPLAYING LINKERS
	if(len(local_frame.linkers) != 0):
		x=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for lid in sorted(local_frame.linkers.keys()):
				q.append(local_frame.linkers[lid].coords[i])
			x.append(hstack(q))

		for lid in sorted(local_frame.linkers.keys()):
				connections.append(local_frame.linkers[lid].connections)

		connections = vstack(connections)

		# Create the points
		src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])
		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24)

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=0.25)
		tube.filter.number_of_sides=12
		mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))

	#DISPLAYING MOTORS
	if(len(local_frame.motors) != 0):
		x=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for mid in sorted(local_frame.motors.keys()):
				q.append(local_frame.motors[mid].coords[i])
			x.append(hstack(q))

		for mid in sorted(local_frame.motors.keys()):
				connections.append(local_frame.motors[mid].connections)

		connections = vstack(connections)

		# Create the points
		src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])
		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24)

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=0.25)
		tube.filter.number_of_sides=12
		mlab.pipeline.surface(tube, color=(0.8, 0.4, 0.8))


@mlab.animate(delay=10, ui=True)
def anim():
	for i in range(0,len(FrameList), 1):
		show_frame(i)
		yield
		

