from numpy import *
from mayavi import mlab

#SPECIFY THE TRAJ FILE AND THE COLOR FILE
#If no color file is specified, the default coloring will be used	
traj_filename='/Users/jameskomianos/Code/M3SYM/M3SYM/Output/snapshot.traj'
color_filename = '/Users/jameskomianos/Code/M3SYM/M3SYM/Output/stresses.traj'
#color_filename = ''

#Open the traj file
traj_file=open(traj_filename)

#Do we have color?
color=False
if(color_filename != ''):
	color_file = open(color_filename)
	color = True

class FilamentSnapshot:
	def __init__(self):
		self.id=None
		self.length=None
		self.delta_left=None	
		self.delta_right=None	
		self.coords={}
		self.colors={}
		self.connections=None

class LinkerSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.coords={}
		self.colors={}
		self.connections=None

class MotorSnapshot:
	def __init__(self):
		self.id=None
		self.type=None
		self.coords={}
		self.colors={}
		self.connections=None

class BrancherSnapshot:
	def __init__(self):
		self.id=None
		self.type=None;
		self.coords={}

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

#read color file
if(color):
	color_lines = [line.strip() for line in color_file]

for line in traj_file:
	line = line.strip()

	if(color):
		color_line = color_lines[line_number]

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
		assert F.n_linkers == len(F.linkers)
		assert F.n_motors == len(F.motors)
		assert F.n_branchers == len(F.branchers)
		n_beads_filament=0
		n_beads_linker=0
		n_beads_motor=0
		n_beads_brancher=0;
		line_number+=1
		continue
			
	if(first_line):
		line_split = line.split()
		if(line[0] == "F"):
			FS=FilamentSnapshot()
			FS.id, FS.length, FS.delta_left, FS.delta_right = map(int,line_split[1:])
			first_line=False
			F.filaments[FS.id]=FS
			FS.connections=vstack(
				[arange(n_beads_filament, n_beads_filament + FS.length - 1.5), 
				arange(n_beads_filament + 1, n_beads_filament + FS.length - .5)]).T
			n_beads_filament+=FS.length
			reading_filament=True
			line_number+=1
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
			line_number+=1
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

	if(reading_filament):
		FS.coords=array(line.split(),'f')
		if(color):
			FS.colors=array(color_line.split(), 'f')
		N=len(FS.coords)
		FS.coords=FS.coords.reshape(N/3,3).T
		first_line=True
		reading_filament=False

	if(reading_linker):
		LS.coords=array(line.split(),'f')
		if(color):
			LS.colors=array(color_line.split(), 'f')
		N=len(LS.coords)
		LS.coords=LS.coords.reshape(N/3,3).T
		first_line=True
		reading_linker=False

	if(reading_motor):
		MS.coords=array(line.split(),'f')
		if(color):
			MS.colors=array(color_line.split(), 'f')
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

@mlab.show
def show_frame(frame_number=-1):

	#PARAMETERS TO SET FOR VISUAL
	#for color scaling
	MAXVAL = 2.00
	MINVAL = 0.00
	SCALETITLE = 'Stresses (pN/m)'
	#SCALETITLE = 'Birth Times (s)'
	COLORMAP = 'OrRd'

	#default color, in RGB
	DBEADCOLOR    = (1.0,1.0,1.0) 
	DFILCOLOR     = (1.0,0.0,0.0)
	DLINKERCOLOR  = (0.3,0.8,0.4)
	DMOTORCOLOR   = (0.2,0.4,0.8)

	local_frame=FrameList[frame_number]
	mlab.figure(1, size=(600, 600), bgcolor=(0, 0, 0))
	mlab.clf()

	#DISPLAYING FILAMENTS
	if(len(local_frame.filaments) != 0):
		x=[]
		c=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for fid in sorted(local_frame.filaments.keys()):
				q.append(local_frame.filaments[fid].coords[i])
			x.append(hstack(q))

		for fid in sorted(local_frame.filaments.keys()):
			for color in local_frame.filaments[fid].colors:
				c.append(color)

		for fid in sorted(local_frame.filaments.keys()):
				connections.append(local_frame.filaments[fid].connections)

		connections = vstack(connections)

		# Create the points
		if(len(c) != 0):
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2], c)
		else:
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])

		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
								    scale_mode='none', color=DBEADCOLOR)

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=0.25)
		tube.filter.number_of_sides=12

		if(len(c) != 0):
			surface = mlab.pipeline.surface(tube, colormap=COLORMAP,
										    vmax = MAXVAL, vmin = MINVAL)
			mlab.colorbar(object=surface, orientation='vertical', 
						  title=SCALETITLE, label_fmt='%.3f')

		else:
			surface = mlab.pipeline.surface(tube, color=DFILCOLOR)

		#display time
		time = 'Time = ' + str(round(local_frame.time,3)) + "s"
		mlab.text(0.6, 0.9, time)

	#DISPLAYING LINKERS
	if(False):
	#if(len(local_frame.linkers) != 0):
		x=[]
		c=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for lid in sorted(local_frame.linkers.keys()):
				q.append(local_frame.linkers[lid].coords[i])
			x.append(hstack(q))

		for lid in sorted(local_frame.linkers.keys()):
			for color in local_frame.linkers[lid].colors:
				c.append(color)

		for lid in sorted(local_frame.linkers.keys()):
				connections.append(local_frame.linkers[lid].connections)

		connections = vstack(connections)

		# Create the points
		if(len(c) != 0):
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2], c)
		else:
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])

		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
								    scale_mode='none', color=DBEADCOLOR)

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=0.25)
		tube.filter.number_of_sides=12

		if(len(c) != 0):
			surface = mlab.pipeline.surface(tube, colormap=COLORMAP, 
											vmax = MAXVAL, vmin = MINVAL)
		else:
			surface = mlab.pipeline.surface(tube, color=DLINKERCOLOR)

	#DISPLAYING MOTORS
	if(False):
	#if(len(local_frame.motors) != 0):
		x=[]
		c=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for mid in sorted(local_frame.motors.keys()):
				q.append(local_frame.motors[mid].coords[i])
			x.append(hstack(q))

		for mid in sorted(local_frame.motors.keys()):
			for color in local_frame.motors[mid].colors:
				c.append(color)

		for mid in sorted(local_frame.motors.keys()):
				connections.append(local_frame.motors[mid].connections)

		connections = vstack(connections)

		# Create the points
		if(len(c) != 0):
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2], c)
		else:
			src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])

		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
									scale_mode='none', color=DBEADCOLOR)

		# Connect them
		src.mlab_source.dataset.lines = connections

		# Finally, display the set of lines
		tube=mlab.pipeline.tube(src, tube_radius=0.25)
		tube.filter.number_of_sides=12

		if(len(c) != 0):
			surface = mlab.pipeline.surface(tube, colormap=COLORMAP, 
											vmax = MAXVAL, vmin = MINVAL)
		else:
			surface = mlab.pipeline.surface(tube, color=DMOTORCOLOR)

	#DISPLAYING BRANCHERS (NO COLOR)
	if(len(local_frame.branchers) != 0):
		x=[]
		connections=[]

		for i in 0, 1, 2:
			q=[]
			for bid in sorted(local_frame.branchers.keys()):
				q.append(local_frame.branchers[bid].coords[i])
			x.append(hstack(q))

		# Create the points
		src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])
		gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24, 
									scale_factor=2.0, color=DBEADCOLOR)


@mlab.animate(delay=10, ui=True)
def anim():
	for i in range(0,len(FrameList), 1):
		show_frame(i)
		yield
		

