#!/usr/bin/env python

import numpy as np
import math as math
import sys as sys
import networkx as nx
import scipy.stats as sp
import matplotlib.pyplot as plt


if(len(sys.argv) != 4) :
	print("ERROR dump.traj ATOMS_PER_GAG numTOD")
	exit(0)

name = sys.argv[1]
ATOMS_PER_GAG = int(sys.argv[2])
numTOD = int(sys.argv[3])

rMax = 30.0
rTol_Min = 0.01
rTol = 20.0
in_file = open("%s.lammpstrj" % (name),'r')
out1_file = open("direct_histo_TOD_%s.dat" % (name),'w')

out1_file.write("%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" % ("#Time(step)","num[1]","num[2]","num[3]","num[4]","num[5]","num[6]","num[6+]","numTOD","EVC_max","EVC_min","EVC_mean","EVC_stdev","EVC_skew","EVC_kurt", "CC_max", "CC_min", "CC_mean", "CC_stdev", "CC_skew", "CC_kurt", "over80", "below80"))

out2_file = open("reduced_eigen_centrality_%s.lammpstrj" % (name),"w")

histoBins = 22 # [1,2,3,4,5,6,6+] #1 is monomer, 6 is hexamer, 7 is TOD, 8/9/10/11 is eigenvector centrality max, mean, stdev, skew, num

FAST_FIND_BIN = 40.0



line = in_file.readline()
line = in_file.readline()
time0 = int(line.split()[0])
line = in_file.readline()
line = in_file.readline()
nAtoms = int(line.split()[0])
nLines = 4

nGag = 0
GagIdentifier = 7 #CTD dimer


# this block will read all lines and count the cluster size. cluster size is stored in nGag. 

counter = -5
for line in in_file :
	nLines=nLines+1
	if (counter < nAtoms and counter > -1) :
		aType = int(line.split()[1])
		if(aType == GagIdentifier) : #or any other representative type
			nGag = nGag + 1
	counter=counter+1

nGag = int(nGag)
print("nLines:",nLines)
nFrames = int(nLines/(nAtoms+9))

in_file.close()
in_file = open("%s.lammpstrj" % (name),'r')

print("Starting time is %d" % (time0))
print("Number of atoms is %d and number of frames is %d" % (nAtoms, nFrames))
print("Number of Gag monomers is %d" % (nGag))


numContacts = 1
contacts_set1 = [7] #start
contacts_set2 = [102] #end

gag_ghost_pos = np.zeros((nGag,numContacts,3))
gag_site_pos = np.zeros((nGag,numContacts,3))
gag_on_flag = np.zeros((nGag))
gag_mol_flag =np.zeros((nGag))

xLEN = 900
yLEN = 900
zLEN = 900 
#Gag-Gag distance based on Gag ID and contact ID (cType)

def calc_dist(id1, id2, cTypes, dflag) :

	distSum = 0
	numTypes = float(len(cTypes)) 
	for ic, cType in enumerate(cTypes) :
		if dflag == contacts_set1[0]:
			x1 = gag_ghost_pos[id1,cType,0]
			y1 = gag_ghost_pos[id1,cType,1]
			z1 = gag_ghost_pos[id1,cType,2]
			x2 = gag_ghost_pos[id2,cType,0]
			y2 = gag_ghost_pos[id2,cType,1]
			z2 = gag_ghost_pos[id2,cType,2]
		if dflag == contacts_set2[0]:
			x1 = gag_site_pos[id1,cType,0]
			y1 = gag_site_pos[id1,cType,1]
			z1 = gag_site_pos[id1,cType,2]
			x2 = gag_site_pos[id2,cType,0]
			y2 = gag_site_pos[id2,cType,1]
			z2 = gag_site_pos[id2,cType,2]

		xd = abs(x1 - x2)
		yd = abs(y1 - y2)
		zd = abs(z1 - z2)
		while(xd > xLEN/2.0) :
			 xd = xd - xLEN
		while(yd > yLEN/2.0) :
			yd = yd - yLEN
		while(zd > zLEN/2.0) :
			zd = zd - zLEN
		dist = (xd*xd + yd*yd + zd*zd)**0.5
		distSum = distSum + dist
	distSum /= numTypes
	return distSum

for i in range(nFrames) :
	print(nFrames)
	#extract chuck of lines for this frame
	chunk_lines = []
	for j in range(nAtoms+9) :
		chunk_lines.append(in_file.readline())

	currTime = int(chunk_lines[1].split()[0])
	print("analyzing currTime:",currTime)
	xlo = float(chunk_lines[5].split()[0])
	xhi = float(chunk_lines[5].split()[1])
	ylo = float(chunk_lines[6].split()[0])
	yhi = float(chunk_lines[6].split()[1])
	zlo = float(chunk_lines[7].split()[0])
	zhi = float(chunk_lines[7].split()[1])
	xLEN = xhi - xlo
	yLEN = yhi - ylo
	zLEN = zhi - zlo
	
	nbinX = int(xLEN/FAST_FIND_BIN)
	nbinY = int(yLEN/FAST_FIND_BIN)
	nbinZ = int(zLEN/FAST_FIND_BIN)
	print("For current time %d, xLEN = %f, yLEN = %f, and zLEN = %f" % (currTime, xLEN, yLEN, zLEN))
	
	histo = np.zeros(histoBins)
	
	gag_ghost_pos = np.zeros((nGag,numContacts,3))   #  gag index, ghost atom index, coordinate
	gag_site_pos = np.zeros((nGag,numContacts,3))    #  gag index, tether site atom index, coordinate
	
	gag_on_flag = np.zeros(nGag)    #  gag index (contains zero-based index of tethered site's gag)
	gag_mol_flag = np.zeros(nGag)
	fast_find_grid = np.zeros((nbinX,nbinY,nbinZ,nGag+1)) #to accelerate neighbor search, index 3 contains # of gag and list of indices
	
	G = nx.Graph()
	Gcyc = nx.Graph()
	
	for j in range(nGag) :
		gagLine = 8 + j*ATOMS_PER_GAG #gagLine determins which line in the FILE corresponds to this gag at this frame
		for k in range(len(contacts_set1)) :
			line = chunk_lines[gagLine + contacts_set1[k]].split()
			xs = float(line[2])
			ys = float(line[3])
			zs = float(line[4])
			x = xs
			y = ys
			z = zs
			aType = int(line[1])
			gag_ghost_pos[j][k][0] = x # ghost particle x
			gag_ghost_pos[j][k][1] = y # ghost particle y
			gag_ghost_pos[j][k][2] = z # ghost particle z
		for l in range(len(contacts_set2)) :
			line = chunk_lines[gagLine + contacts_set2[l]].split()
			xs = float(line[2])
			ys = float(line[3])
			zs = float(line[4])
			x = xs
			y = ys
			z = zs
			aType = int(line[1])
			#print "Populating site_pos with atom type %d from Gag %d for set %d" % (aType, j, k)
			gag_site_pos[j][l][0] = x # site particle x
			gag_site_pos[j][l][1] = y # site particle y
			gag_site_pos[j][l][2] = z # site particle z
		
		com_index = 66 #arbitrary, used to determine fast indexing bins (eg, fake cell list)
		line = chunk_lines[gagLine + com_index].split()
		x = float(line[2])
		y = float(line[3])
		z = float(line[4])
		
		vx = int(line[5])
		vy = int(line[6])
		if (vx == 0) :
			gag_on_flag[j] = 1
			gag_mol_flag[j] = vy
	
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
		fast_find_grid[ffx,ffy,ffz,ffcount+1] = j #gag index (starting from zero)
		#print(fast_find_grid[ffx,ffy,ffz,:]) #fast_find_grid[xgrid,ygrid,zgrid,:]=index 0:num of Gag and,rest:molID
	
	print("Performing neighbor search now")
	for m in range(nbinX) :
		for n in range(nbinY) :
			for l in range(nbinZ) :
				num_NN1 = int(fast_find_grid[m,n,l,0])
				#if num_NN1 > 0 :
				#	print("numNN1 is %d" % num_NN1)
				for nnrow in range(num_NN1) :
					gagID = int(fast_find_grid[m,n,l,nnrow+1])
					#print(gagID)
					#print("gag_on_flag:",gag_on_flag[gagID])
					if(gag_on_flag[gagID] != 1) :
						continue
				        #check within +-1 of the current grid
					for ii in range(3) :
						for jj in range(3) :
							for kk in range(3) :
								mtemp = m - 1 + ii
								ntemp = n - 1 + jj
								ltemp = l - 1 + kk
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
								num_NN2 = int(fast_find_grid[mtemp,ntemp,ltemp,0])
								#print("numNN2 ia %d" % num_NN2)


								for nn2_id in range(num_NN2) :
									gagID2 = int(fast_find_grid[mtemp,ntemp,ltemp,nn2_id+1])
									# now check for contacts
									if(gagID == gagID2) :
										continue
									if(gagID != gagID2):# in gag1_indices and gagID2 in gag1_indices) :
									dist_ntd = calc_dist(gagID, gagID2, [0], contacts_set1[0])
										dist_ctd = calc_dist(gagID, gagID2, [0], contacts_set2[0])
										if ( dist_ntd < rTol ) :
											G.add_edges_from([(gagID,gagID2)])
											if(Gcyc.has_edge(gagID,gagID2) == False) :##### subgraph analysis
												Gcyc.add_edges_from([(gagID,gagID2)]) ##### subgraph
										if ( dist_ctd < rTol ) :
											G.add_edges_from([(gagID,gagID2)])

	
	print("Performing connectivity search now")
	cFlag = True
	try : 
		node_list = G.nodes()
		G_max = max(nx.connected_components(G), key=len)
		G_keep = G.subgraph(G_max)
		G = G_keep
		print(len(G_max))
		centrality = nx.current_flow_betweenness_centrality(G)
		sort_eigen = np.asarray(sorted(centrality.values())) #sorted eigenvalues, max is last
		keep_eigen = sort_eigen[:]
		min_eigen = keep_eigen.min()
		max_eigen = keep_eigen.max()
		# Largest Eigen, Min Eigen, Mean Eigen, Stdev Eigen, Skew Eigen, # of nodes in graph
		# eigenvalues only mean anything as ratios
		histo[8] = keep_eigen[-1] / max_eigen
		histo[9] = keep_eigen[0] / max_eigen
	
		keep_eigen[:] = keep_eigen[:] / keep_eigen.sum()
		histo[10] = keep_eigen.mean() #/ max_eigen
		histo[11] = keep_eigen.std() #/ max_eigen
		histo[12] = sp.skew(keep_eigen[:])
		histo[13] = sp.kurtosis(keep_eigen[:]) #len(keep_eigen)
	
		centrality2 = nx.eccentricity(G) 
		sort_eigen = np.asarray(sorted(centrality2.values())) #sorted eigenvalues, max is last
	
		keep_eigen = sort_eigen[:] #/ max_eigen2
		max_eigen2 = float(keep_eigen.max())
	
		histo[14] = keep_eigen[-1] #/ min_eigen2
		histo[15] = keep_eigen[0] #/ min_eigen2
		histo[16] = keep_eigen.mean() #/ min_eigen2
		histo[17] = keep_eigen.std() #/ min_eigen2
		histo[18] = sp.skew(keep_eigen) 
		histo[19] = sp.kurtosis(keep_eigen) 
	
		norm_eigen = keep_eigen[:] / max_eigen2
	
		print("Maximum eigenvalue is %f" % max_eigen)
		print("Maximum eccentricity is %f" % max_eigen2)
		print("Num samples is %d" % len(norm_eigen))
		print("Num over 0.81 is %d" % len(norm_eigen[norm_eigen[:] > 0.81]))
		print("Num below 0.81 is %d" % len(norm_eigen[norm_eigen[:] < 0.81]))
	
		over80 = len(norm_eigen[norm_eigen[:] > 0.81])
		below80 = len(norm_eigen[norm_eigen[:] < 0.81])
	
		histo[20] = over80
		histo[21] = below80
	except :
		print("exception raised; assign eccentricity 0")
		histo[8:] = 0.0
		cFlag = False
	
	out2_file.write("ITEM: TIMESTEP\n")
	out2_file.write("%d\n" % currTime)
	out2_file.write("ITEM: NUMBER OF ATOMS\n")
	out2_file.write("%d\n" % nGag)
	out2_file.write("ITEM: BOX BOUNDS pp pp pp\n")
	out2_file.write("0.0 900.0\n")
	out2_file.write("0.0 900.0\n")
	out2_file.write("0.0 900.0\n")
	out2_file.write("ITEM: ATOMS id type x y z vx vy vz\n")
	#ec_count = -1
	for ii in range(nGag) :
		xx=gag_site_pos[ii,0,0]
		yy=gag_site_pos[ii,0,1]
		zz=gag_site_pos[ii,0,2]
		vvx=1
		vvy=0
		vvz=gag_mol_flag[ii]
		if(gag_on_flag[ii] == 1) :
			vvx=0
			if(vvz < numTOD) :
				histo[7] = histo[7] + 1.0
			
		if cFlag :
			if ii in centrality2 :
				vvy=float(centrality2[ii]) / max_eigen2
		out2_file.write("%d %d %f %f %f %f %f %f\n" % (ii+1,1,xx,yy,zz,vvx,vvy,vvz))
	
	
	for ii in range(nGag) :
		try :
			cyc_edge = nx.find_cycle(Gcyc,source=ii,orientation="ignore")

			length = int(len(cyc_edge))
			if(length <= 6) :
				histo[length-1] = histo[length-1] + 1.0
			else :
				histo[6] = histo[6] + 1.0
		except :
			histo[6] = histo[6] + 1.0


	extra_data = 14

	out1_file.write("%d    " % (currTime))
	for ii in range(histoBins-extra_data) :
		out1_file.write("   %d   " % math.ceil(histo[ii]))
	for ii in range(extra_data) :
		out1_file.write("   %f   " % (histo[ii+histoBins-extra_data]))
	out1_file.write("\n")

in_file.close()

out1_file.close()
out2_file.close()

