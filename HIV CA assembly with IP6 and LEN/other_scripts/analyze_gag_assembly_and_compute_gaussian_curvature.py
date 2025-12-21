import numpy as np
import networkx as nx
import mdtraj as md
import sys as sys
import math as math
from scipy.spatial import Delaunay
import matplotlib as mpl
import matplotlib.pyplot as plt
from collections import defaultdict


trajname=sys.argv[1] # The name of the trajectory (.lammpstrj)
dt=int(sys.argv[2]) # Num timesteps between frames

gcontact = 20.0 #distance for contacts (gag-gag)
stats = ("%s_stats.dat" % (trajname))
FSTATS = open(stats,"w")
FAST_FIND_BIN = 40.0
DEBUG_FLAG = 0
out2_file = open("reduced_gaussian_curvature_%s.lammpstrj" % (trajname),"w")
xrange = range

def write_lammps_header(tfile, time, natoms, dims, fields) :
    tfile.write("ITEM: TIMESTEP\n")
    tfile.write("%d\n" % time)
    tfile.write("ITEM: NUMBER OF ATOMS\n")
    tfile.write("%d\n" % natoms)
    tfile.write("ITEM: BOX BOUNDS pp pp pp\n")
    tfile.write("0.0 900.0\n") #%f\n" % dims[0])
    tfile.write("0.0 900.0\n") #%f\n" % dims[1])
    tfile.write("0.0 900.0\n") #%f\n" % dims[2])
    tfile.write("ITEM: ATOMS ")
    for f in fields :
        tfile.write("%s " % f)
    tfile.write("\n")

def write_lammps_atom(tfile, index, type, x, y, z, vx, vy,vz) :
    tfile.write("%d %d %f %f %f %i %i %i\n" % (index, type, x, y, z, vx, vy, vz))

def load_atom_traj(trajname) :
    FTRAJ=open(trajname,"r")
    all_lines = FTRAJ.readlines()
    natoms = int(all_lines[3].split()[0])
    nframes = int(len(all_lines) / (natoms+9))
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

def sorted_triangle(triangle):
   return tuple(sorted(triangle))

def gaussian_curvature(cluster):

   #print(cluster)
   points = []
   gaussian_curvature = {}
   for arr in cluster :
      #print(int(arr[0]))
      if int(arr[0]) == 1 and int(arr[3]) == 55 :
         points.append(arr[-3:])
         gaussian_curvature[tuple(arr[-3:])] = [0,arr[1]]
      if int(arr[0]) == 2 and int(arr[3]) == 55 :
         points.append(arr[-3:])
         gaussian_curvature[tuple(arr[-3:])] = [0,arr[1]]
   points = np.array(points)
   print("points",len(points))
   if len(points) == 0:
      print("skip gc computation")
   else:
      tri = Delaunay(points)
      triangles = tri.simplices
      Nnodes = len(points)
      triangle_count = {}
    
      for simplex in triangles:
         faces = [
            sorted_triangle([simplex[0], simplex[1], simplex[2]]),
            sorted_triangle([simplex[0], simplex[1], simplex[3]]),
            sorted_triangle([simplex[0], simplex[2], simplex[3]]),
            sorted_triangle([simplex[1], simplex[2], simplex[3]])
         ]
    
    
         for face in faces:
            if face in triangle_count:
               triangle_count[face] += 1
            else:
               triangle_count[face] = 1
      surface_triangles = [face for face, count in triangle_count.items() if count == 1]
      surface_triangles_coordinates = np.array([[points[vertex] for vertex in triangle] for triangle in surface_triangles])
    
      vertex_to_triangles = defaultdict(list)
    
      for v,triface  in enumerate(surface_triangles_coordinates):
         for vertex in triface:
            vertex_to_triangles[tuple(vertex)].append(v)
    
    
      surface_coordinates = np.array([points[vertex] for triangle in surface_triangles for vertex in triangle])
    
      for vertex, triangle_indices in vertex_to_triangles.items():
         total_angle = 0
         for i in triangle_indices:
            facevertex = surface_triangles_coordinates[i,:,:]
            edges = []
            for t in facevertex:
               #print("t:",t)
               if not np.array_equal(t, np.array(vertex), equal_nan=True):
                  edges.append(t)
            v1, v2 = edges[0] - np.array(vertex), edges[1] - np.array(vertex)
            angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            total_angle += angle
         angle_deficit = 2 * np.pi - total_angle
         area = 0
    
         for j in triangle_indices:
            ftriangle = surface_triangles_coordinates[j,:,:]
            a, b, c = np.array(ftriangle[0]), np.array(ftriangle[1]), np.array(ftriangle[2])
            area += 0.5 * np.linalg.norm(np.cross(b-a, c-a))
    
            vv = gaussian_curvature[vertex]
            gaussian_curvature[tuple(vertex)] = (angle_deficit/area, vv[1])
      values = [v[0] for v in gaussian_curvature.values()]
      plt.hist(values, bins=10, edgecolor='black')
      plt.savefig("plot_gc_{}.png".format(trajname))
   return gaussian_curvature

def lammps_traj_file(coordinate,curvature):
    currTime = 0
    nGag = len(coordinate)

    out2_file.write("ITEM: TIMESTEP\n")
    out2_file.write("%d\n" % currTime)
    out2_file.write("ITEM: NUMBER OF ATOMS\n")
    out2_file.write("%d\n" % nGag)
    out2_file.write("ITEM: BOX BOUNDS pp pp pp\n")
    out2_file.write("0.0 900.0\n")
    out2_file.write("0.0 900.0\n")
    out2_file.write("0.0 900.0\n")
    out2_file.write("ITEM: ATOMS id type x y z vx vy vz\n")
    sr = 0
    clusterid = {}
    mollist = []
    for key, (value1, value2) in curvature.items():
        clusterid[int(value2)] = value1
        mollist.append(int(value2))

    for subarr in coordinate:
        vvy = subarr[0]
        gagid = int(subarr[1])
        cgtype = subarr[3]
        if gagid in mollist :
            vvz = clusterid[gagid]
            out2_file.write("%d %d %f %f %f %f %f %f\n" % (sr+1,cgtype,subarr[4],subarr[5],subarr[6],0,vvy,vvz))
        else:
            vvz = 0.0
            out2_file.write("%d %d %f %f %f %f %f %f\n" % (sr+1,cgtype,subarr[4],subarr[5],subarr[6],0,vvy,vvz))
        sr += 1
    print(sr)


traj, dims, nframes, natoms = load_atom_traj("%s.lammpstrj" % (trajname))

out_traj = ("final_%s_extracted.lammpstrj" % (trajname))
outfile = open(out_traj,"w")

for i in xrange(nframes) :
    gag1_indices = set()
    gag2_indices = set()
    ip6_indices = set()
    len_indices = set()
    for ff, f in enumerate(traj[i,:,1]) :
        if( f == 7) :
            gag1_indices.add(ff)
        if( f == 7) :
            gag2_indices.add(ff)
        if(f == 140) :
            ip6_indices.add(ff)
        if( f > 140) :
            len_indices.add(ff)

    ngag = len(gag2_indices) #num mols
    nip6 = len(ip6_indices)
    nlen = len(len_indices)
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
        ffcount = int(fast_find_grid[ffx,ffy,ffz,0]) # gives number of neighbor atoms for a partiulat atom
#        print("ffcount is:",ffcount)
        fast_find_grid[ffx,ffy,ffz,0] = ffcount + 1 # increases neighbor count by one
    #    print("fast_find_grid:", fast_find_grid[ffx,ffy,ffz,0])
        fast_find_grid[ffx,ffy,ffz,ffcount+1] = ff # atomid for neighboring atom is stored here
    #    print("ff: ", ff)
    # process Gag connectivity graph
    print("Processing graph for frame %d" % i)
    G = nx.Graph() #gag-gag
    Gcyc = nx.Graph()
#    I = nx.Graph()
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
                                        m1 = int(math.ceil( float(traj[i, aid1, 0]) / 134.0 ) )
                                        m2 = int(math.ceil( float(traj[i, aid2, 0]) / 134.0 ) )
                                        dist = calc_dist(r1, r2, dims)
                                        if ( dist < gcontact ) :
                                            G.add_edges_from([(m1, m2)])
                                            if(Gcyc.has_edge(m1,m2) == False) :##### subgraph analysis
                                                Gcyc.add_edges_from([(m1,m2)]) ##### subgraph
                                    if(aid1 in gag2_indices and aid2 in gag2_indices) :
                                        r1 = traj[i, aid1, 2:5]
                                        r2 = traj[i, aid2, 2:5]
                                        m1 = int(math.ceil( float(traj[i, aid1, 0]) / 134.0 ) )
                                        m2 = int(math.ceil( float(traj[i, aid2, 0]) / 134.0 ) )
                                        dist = calc_dist(r1, r2, dims)
                                        if ( dist < gcontact) :
                                            G.add_edges_from([(m1, m2)])


    cluster_size = 0.0 #equivalent to nc_gag
    nc_gag = 0.0
    G_max = set() #will be a set
    try :
        G_max = max(nx.connected_components(G), key=len)
        G_max = G
        cluster_size = len(G_max)
        print("Cluster size:", cluster_size)
    except :
        print("No cluster found")
    time = dt*i
    pent = []
    hexamer = []
    for ii in range(ngag*134) :
        try :
            cyc_edge = nx.find_cycle(Gcyc,source=ii,orientation="ignore")
            length = int(len(cyc_edge))
            if(length == 5) :
                for tup in cyc_edge :
                    for ele in tup :
                        if (ele not in pent) :
                            pent.append(ele)
            if(length == 6) :
                for t6 in cyc_edge :
                    for ele6 in t6 :
                        if (ele6 not in hexamer) :
                            hexamer.append(ele6)
        except :
            pass
    if ( cluster_size > 0 ) :
        for nn in G_max :
            if(nn <= ngag) :
                nc_gag += 1.0
        # print padded trajectory with clustered Gag
        tempatoms = int(ngag*134) + int(nip6) + int(nlen)
        write_lammps_header(outfile, time, natoms, dims, ["id", "type", "x", "y", "z", "vx", "vy", "vz"])
        ##############################################################
        cx = 0
        cy = 0
        cz = 0
        N = 0
        box = list(range(900))
        #cut_y = set()
        x_co = set() # stores int of x
        y_co = set()
        z_co = set()
        for ii, (aid, atype, xx, yy, zz) in enumerate(traj[i,:,:]) :
            tempmid = 0
            if ( aid <= ngag*134) :
                tempmid = int(math.ceil( aid / 134.0 ))
                if ( tempmid in G_max) :
                    if int(xx) not in x_co :
                        x_co.add(int(xx))
                    if int(yy) not in y_co :
                        y_co.add(int(yy))
                    if int(zz) not in z_co :
                        z_co.add(int(zz))
        cut_x = list(set(box) - set(sorted(x_co)))
        off_x = min(cut_x) if cut_x else 0
        cut_y = list(set(box) - set(sorted(y_co)))
        off_y = min(cut_y) if cut_y else 0
        cut_z = list(set(box) - set(sorted(z_co)))
        off_z = min(cut_z) if cut_z else 0
        # print(cut_y)
        for ii, (aid, atype, xx, yy, zz) in enumerate(traj[i,:,:]) :
            tempmid = 0
            if ( aid <= ngag*134) :
                tempmid = int(math.ceil( aid / 134.0 ))
                if ( tempmid in G_max) :
                    if(xx < off_x) :
                        xx = xx + 900.0
                        traj[i,ii,2] = xx
                    if(yy < off_y) :
                        yy = yy + 900.0
                        traj[i,ii,3] = yy
                    if(zz < off_z) :
                        zz = zz + 900.0
                        traj[i,ii,4] = zz
                    cx = cx + xx
                    cy = cy + yy
                    cz = cz + zz
                    N = N + 1
        ccx = cx/N
        ccy = cy/N
        ccz = cz/N
        print(ccx,ccy,ccz,N)
        for ii, (aid, atype, xx, yy, zz) in enumerate(traj[i,:,:]) :
            tempmid = 0
            if ( aid <= ngag*134) :
                tempmid = int(math.ceil( aid / 134.0 ))
                if ( tempmid in G_max) :
                    xx = (xx - ccx + 450)
                    traj[i,ii,2] = xx
                    yy = (yy - ccy + 450)
                    traj[i,ii,3] = yy
                    zz = (zz - ccz + 450)
                    traj[i,ii,4] = zz

        for ii, (aid, atype, xx, yy, zz) in enumerate(traj[i,:,:]) :

            tempmid = 0
            if ( aid <= ngag*134) :
                tempmid = int(math.ceil( aid / 134.0 ))
                if ( tempmid in G_max and tempmid in pent ) :
                    write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 1, 0)
                if ( tempmid in G_max and tempmid in hexamer and tempmid not in pent ) : 
                    write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 2, 0)
                if ( tempmid in G_max and tempmid not in hexamer and tempmid not in pent ) :
                    write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 3, 0)
                if ( tempmid not in G_max) :
                    write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 0, 0)
            if ( atype == 139) :
                write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 6, 0)
            if ( atype == 140) :
                write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 4, 0)
            if ( atype > 140) :
                write_lammps_atom(outfile, aid, atype, xx, yy, zz, 0, 5, 0)
    FSTATS.write("%d %d\n" % (time, nc_gag))

    capsid = np.zeros((ngag*134,7))
    k = 0
    
    xlen = 900.0
    ylen = 900.0
    zlen = 900.0
    
    for ii, (aid, atype, xx, yy, zz) in enumerate(traj[i,:,:]) :
        gagid = 0
        if (aid <= ngag*134) :
            gagid = int(math.ceil( aid / 134.0 ))
            if (gagid in G_max and gagid not in hexamer and gagid not in pent):
                nn = 3.0
                capsid[k,0],capsid[k,1],capsid[k,2],capsid[k,3],capsid[k,4],capsid[k,5],capsid[k,6] = nn, gagid,aid, atype, xx, yy, zz
            if (gagid in G_max and gagid in hexamer):
                nn = 2.0
                capsid[k,0],capsid[k,1],capsid[k,2],capsid[k,3],capsid[k,4],capsid[k,5],capsid[k,6] = nn, gagid, aid, atype,xx, yy, zz
            if (gagid in G_max and gagid in pent ):
                nn = 1.0
                capsid[k,0],capsid[k,1],capsid[k,2],capsid[k,3],capsid[k,4],capsid[k,5],capsid[k,6] = nn, gagid, aid, atype,xx, yy, zz

            if (gagid not in G_max) :
                nn = 0.0
                capsid[k,0],capsid[k,1],capsid[k,2],capsid[k,3],capsid[k,4],capsid[k,5],capsid[k,6] = nn, gagid,aid, atype, xx, yy, zz
            k += 1
    capsid = np.array([xx for xx in capsid if np.array_equal(xx, np.zeros(7)) == False ])

    print(capsid.shape) 
    
    curvature_data = gaussian_curvature(capsid)
    lammps_traj_file(capsid,curvature_data)


out2_file.close()

FSTATS.close()
outfile.close()

