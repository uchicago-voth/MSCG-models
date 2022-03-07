import numpy as np
import networkx as nx
import mdtraj as md
import sys as sys
import math as math

# USAGE:
# analyze_assembly_and_extract_clusters.py trajbasename startindex endindex clusterMolID dtau_per_frame 
# assumes trajectory is LAMMPSTRJ format (xs, ys, zs)

##################### USER SPECIFIED INFO ##############################

trajname=sys.argv[1] # The name of the trajectory (base+index.lammpstrj)
startindex=int(sys.argv[2])
endindex=int(sys.argv[3])
mid=int(sys.argv[4]) # The molecule ID of the primary cluster
dt=int(sys.argv[5]) # Num timesteps between frames

##################### END USER SPECIFIED INFO ##############################

gcontact=18.5 #distance for contacts (gag-gag)

FSTATS=open("stats.dat","w")
FAST_FIND_BIN = 30.0
DEBUG_FLAG = 0

FSTATS.write("#time clustersize\n")

# write lammps stuff
def write_lammps_header(tfile, time, natoms, dims, fields) :
    tfile.write("ITEM: TIMESTEP\n")
    tfile.write("%d\n" % time)
    tfile.write("ITEM: NUMBER OF ATOMS\n")
    tfile.write("%d\n" % natoms)
    tfile.write("ITEM: BOX BOUNDS pp pp pp\n")
    tfile.write("0.0 %f\n" % dims[0])
    tfile.write("0.0 %f\n" % dims[1])
    tfile.write("0.0 %f\n" % dims[2])
    tfile.write("ITEM: ATOMS ")
    for f in fields :
        tfile.write("%s " % f)
    tfile.write("\n")

def write_lammps_atom(tfile, index, type, x, y, z, vx, vy, vz) :
    tfile.write("%d %d %f %f %f %f %f %f\n" % (index, type, x, y, z, vx, vy, vz))

#returns array ( nframes, natoms, (aid,type,x,y,z) )
def load_atom_traj(trajname) :
    FTRAJ=open(trajname,"r")
    all_lines = FTRAJ.readlines()
    natoms = int(all_lines[3].split()[0])
    nframes = len(all_lines) / (natoms+9)
    traj = np.zeros((nframes, natoms, 5))
    XLEN=0.0
    YLEN=0.0
    ZLEN=0.0
    for f in xrange(nframes) :
        start_line = f * (natoms+9)
        [xlo, xhi] = [float(i) for i in all_lines[start_line+5].split()] 
        [ylo, yhi] = [float(i) for i in all_lines[start_line+6].split()] 
        [zlo, zhi] = [float(i) for i in all_lines[start_line+7].split()] 
        XLEN = xhi - xlo
        YLEN = yhi - ylo
        ZLEN = zhi - zlo
        start_line += 9
        for a in xrange(natoms) :
            line = all_lines[start_line+a].split()
            aid = int(line[0])
            atype = int(line[1])
            xs = float(line[2])
            ys = float(line[3])
            zs = float(line[4])
            x = xs*XLEN + xlo
            y = ys*YLEN + ylo
            z = zs*ZLEN + zlo
            traj[int(f),int(a),0] = aid
            traj[int(f),int(a),1] = atype
            traj[int(f),int(a),2] = x
            traj[int(f),int(a),3] = y
            traj[int(f),int(a),4] = z

    return traj, [XLEN, YLEN, ZLEN], nframes, natoms

def calc_dist(r1, r2, dims) :
    dx = abs(r2[0] - r1[0])
    dy = abs(r2[1] - r1[1])
    dz = abs(r2[2] - r1[2])
    if(dx > dims[0]/2.0) :
        dx - dims[0]
    if(dy > dims[1]/2.0) :
        dy - dims[1]
    if(dz > dims[2]/2.0) :
        dz - dims[2]
    dist = (dx*dx + dy*dy + dz*dz)**0.5
    return dist

##################### COLLECT ALL TRAJECTORY DATA FIRST ######################################
nframes = 0
traj, dims, nframes, natoms = load_atom_traj("%s%d.lammpstrj" % (trajname, startindex))
for i in range(endindex - startindex) :
    tempname = "%s%d.lammpstrj" % (trajname, int(startindex+i+1))
    traj2, dims2, nframes2, null = load_atom_traj(tempname)    
    traj = np.concatenate((traj, traj2[1:,:,:]), axis=0)
    nframes = nframes + nframes2 - 1
print "Appended trajectory has shape:"
print traj.shape


##################### ANALYZE ALL TRAJECTORY DATA ######################################
# analyze cluster size

outfile = open("extracted.lammpstrj","w")

for i in xrange(nframes) :
    gag1_indices = set()
    gag2_indices = set()

    for ff, f in enumerate(traj[i,:,1]) :
        if(f == 26 or f == 61) :
            gag1_indices.add(ff)
        if(f == 31 or f == 66) :
            gag2_indices.add(ff)

    ngag = len(gag1_indices) #num mols

    # fast neighbor search
    nbinX = int( dims[0] / FAST_FIND_BIN )
    nbinY = int( dims[1] / FAST_FIND_BIN )
    nbinZ = int( dims[2] / FAST_FIND_BIN )
 
    fast_find_grid = np.zeros((nbinX, nbinY, nbinZ, natoms+1)) #0-index has length of list, rest is list of indices in grid
    for ff, (x,y,z) in enumerate(traj[i,:,2:5]) :
        ffx = int(x/FAST_FIND_BIN)
        ffy = int(y/FAST_FIND_BIN)
        ffz = int(z/FAST_FIND_BIN)
        if(ffx < 0) :
            ffx = ffx + nbinX
        if(ffx >= nbinX) :
            ffx = ffx - nbinX
        if(ffy < 0) :
            ffy = ffy + nbinY
        if(ffy >= nbinY) :
            ffy = ffy - nbinY
        if(ffz < 0) :
            ffz = ffz + nbinZ
        if(ffz >= nbinZ) :
            ffz = ffz - nbinZ
        ffcount = int(fast_find_grid[ffx,ffy,ffz,0])
        fast_find_grid[ffx,ffy,ffz,0] = ffcount + 1
        fast_find_grid[ffx,ffy,ffz,ffcount+1] = ff
    
    # process Gag connectivity graph
    print "Processing graph for frame %d" % i
    G = nx.Graph() #gag-gag
    for m in xrange(nbinX) :
        for n in xrange(nbinY) :
            for l in xrange(nbinZ) :
                nn1_num = int(fast_find_grid[m,n,l,0])
                for nn1_id in xrange(nn1_num) :
                    aid1 = int(fast_find_grid[m,n,l,nn1_id+1])
                    #check within +-1 of the current grid
                    for mm in xrange(3) :
                        for nn in xrange(3) :
                            for ll in xrange(3) :
                                mtemp = m - 1 + mm
                                ntemp = n - 1 + nn
                                ltemp = l - 1 + ll
                                if(mtemp < 0) :
                                    mtemp = mtemp + nbinX
                                if(mtemp >= nbinX) :
                                    mtemp = mtemp - nbinX
                                if(ntemp < 0) :
                                    ntemp = ntemp + nbinY
                                if(ntemp >= nbinY) :
                                    ntemp = ntemp - nbinY
                                if(ltemp < 0) :
                                    ltemp = ltemp + nbinZ
                                if(ltemp >= nbinZ) :
                                    ltemp = ltemp - nbinZ
                                nn2_num = int(fast_find_grid[mtemp,ntemp,ltemp,0])
                                for nn2_id in xrange(nn2_num) :
                                    aid2 = int(fast_find_grid[mtemp,ntemp,ltemp,nn2_id+1])
                                    # now check for contacts
                                    if(aid1 == aid2) :
                                        continue
                                    if(aid1 in gag1_indices and aid2 in gag1_indices) :
                                        r1 = traj[i, aid1, 2:5]
                                        r2 = traj[i, aid2, 2:5]
                                        m1 = int(math.ceil( float(traj[i, aid1, 0]) / 35.0 ) )
                                        m2 = int(math.ceil( float(traj[i, aid2, 0]) / 35.0 ) )
                                        dist = calc_dist(r1, r2, dims)
                                        if ( dist < gcontact ) :
                                            G.add_edges_from([(m1, m2)])
                                    elif(aid1 in gag2_indices and aid2 in gag2_indices) :
                                        r1 = traj[i, aid1, 2:5]
                                        r2 = traj[i, aid2, 2:5]
                                        m1 = int(math.ceil( float(traj[i, aid1, 0]) / 35.0 ) )
                                        m2 = int(math.ceil( float(traj[i, aid2, 0]) / 35.0 ) )
                                        dist = calc_dist(r1, r2, dims)
                                        if ( dist < gcontact ) :
                                            G.add_edges_from([(m1, m2)])

    cluster_size = 0.0 #equivalent to nc_gag
    total_size = 0.0
    nc_gag = 0.0
    G_max = set() #will be a set

    try :
        G_max = nx.node_connected_component(G, mid)
        cluster_size = len(G_max)
    except :
        print "No cluster found"

    time = dt*i

    # only continue if we have a cluster size > 1, otherwise, skip printing?
    if ( cluster_size > 0 ) :
        
        for nn in G_max :
            if(nn <= ngag) :
                nc_gag += 1.0

        # print padded trajectory with clustered Gag
        tempatoms = int(nc_gag*35)
        write_lammps_header(outfile, time, tempatoms, dims, ["id", "type", "x", "y", "z", "vx", "vy", "vz"])
        for ii, (aid, atype, xx, yy, zz) in enumerate(traj[i,:,:]) :
            
            # get molID then determine if its in the GRI_max cluster
            # if so, print it
            tempmid = 0
            if ( aid <= ngag*35) :
                tempmid = int(math.ceil( aid / 35.0 ))
            if ( tempmid in G_max ) :
                write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 0, 0)
                    

    FSTATS.write("%d %f\n" % (time, nc_gag))        

FSTATS.close()
outfile.close()
