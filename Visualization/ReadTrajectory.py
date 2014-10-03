from numpy import *
from mayavi import mlab

mlab.figure(1, size=(600, 600), bgcolor=(0, 0, 0))
mlab.show()

filename='/Users/jameskomianos/Code/CytoSim-Repo/Cyto/beadoutput.txt'
traj_file=open(filename)

class FilamentSnapshot:
	def __init__(self):
		self.id=None
		self.length=None
		self.delta_left=None	
		self.delta_right=None	
		self.coords=None
		self.connections=None

class Frame:
	def __init__(self):
		self.step=None
		self.time=None
		self.n_filaments=None
		self.filaments={}

FrameList=[]
first_frame_line=True
first_filament_line=True
n_beads=0
for line in traj_file:
	line = line.strip()
	#print "["+line+"]"
	if(first_frame_line):
		F=Frame()
		F.step, F.time, F.n_filaments = map(double,line.split())
		F.n_filaments = int(F.n_filaments)
		first_frame_line=False
		continue
	if(len(line)==0):
		#print "Detected an empty line"
		first_frame_line=True
		first_filament_line=True
		FrameList.append(F)
		assert F.n_filaments == len(F.filaments)
		n_beads=0
		continue
	if(first_filament_line):
		FS=FilamentSnapshot()
		FS.id, FS.length, FS.delta_left, FS.delta_right = map(int,line.split())
		first_filament_line=False
		F.filaments[FS.id]=FS
		FS.connections=vstack([arange(n_beads, n_beads + FS.length - 1.5), arange(n_beads + 1, n_beads + FS.length - .5)]).T
		n_beads+=FS.length
		continue
	FS.coords=array(line.split(),'f')
	N=len(FS.coords)
	FS.coords=FS.coords.reshape(N/3,3).T
	first_filament_line=True

@mlab.show
def show_frame(frame_number=-1):
	local_frame=FrameList[frame_number]
	x=[]
	connections=[]

	for i in 0, 1, 2:
	    #q=[f.coords[i] for f in local_frame.filaments]
	    q=[]
	    for fid in sorted(local_frame.filaments.keys()):
	        q.append(local_frame.filaments[fid].coords[i])
	    x.append(hstack(q))

	for fid in sorted(local_frame.filaments.keys()):
	        connections.append(local_frame.filaments[fid].connections)

	connections = vstack(connections)

	mlab.figure(1, size=(600, 600), bgcolor=(0, 0, 0))
	mlab.clf()

	# Create the points
	src = mlab.pipeline.scalar_scatter(x[0], x[1], x[2])
	gsphere=mlab.pipeline.glyph(src, mode="sphere", resolution=24)

	# Connect them
	src.mlab_source.dataset.lines = connections

	# Finally, display the set of lines
	tube=mlab.pipeline.tube(src, tube_radius=0.25)
	tube.filter.number_of_sides=12
	mlab.pipeline.surface(tube, color=(0.8, 0.3, 0.3))


@mlab.animate(delay=10, ui=True)
def anim():
	for i in range(0,len(FrameList)):
		show_frame(i)
		yield
		

