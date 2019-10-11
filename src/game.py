'''
Created on Sep 22, 2009

@author: bareno
'''
from numpy import *
import random
import pickle
import pprint

debug = False

class GameBoard(object):
    ''' container class for LiCoMn game simulation'''
     
    def __init__(self, bSize=1):
        '''
        Constructor
        Remember that not all sites will be occupied. I'm using a lattice oriented along alpha sites. Alpha sites occupy [3 0] 
        positions. Beta sites occupy [1 2] positions. [1 0], [2 0], [1 1] and [2 2] are empty allways 
        '''
        self.boardSize = 3*bSize # last site is equiv to origin, by trans simmetry
        self.board = zeros((self.boardSize, self.boardSize)) #All atomic sites
        self.pieces = ["..", "Li", "Mn", "Co", "H"]
        # Stores integer codes for the atoms that can occupy lattice sites
        self.totBoardSites = 9*bSize**2
        self.totAtomSites = 3*bSize**2 #total number of sites in board
        self.totAtoms = [0, 0, 0, 0] #Number of [occupied sites, Li, Mn, Co]
        self.targetComp =[1, 0.2, 0.4, 0.4] #Target [occupancy, Li-conc, Mn_conc, Co_conc]
        self.actualComp =[0, 0, 0, 0] #Actual values [same as above]
        self.neighbor_stats = zeros((3,9))
        #nearest neighbor stats table
        # 1st index is central atom: Li, Mn, Co
        # 2n index is stat: #of atoms, #Li neighbors, #Mn neigh, #Co neigh, #li-2nd neigh, #Mn 2nd n, #Co 2nd n, #bulk-like, #edge
        
        
        ''' Actually generates the board '''
         #generate list of empty alpha sites
        freeAlpha=[]
        for x in range(bSize):
            for y in range(bSize):
                freeAlpha.append((3*x, 3*y))
        
       #pick random empty alpha site
        #first is always isolated
        site_num = random.randrange(len(freeAlpha))
        site = freeAlpha[site_num]
        del(freeAlpha[site_num]) #no longer free
        #pick type
        type = random.randrange(1,4,2) #1=Li, 3=Co
        
        self.board[site[0], site[1]] = type
        for neigh in self.get_Alpha_Neighbors(site):
            self.fill_site(neigh, type, self.getNeighbors(neigh))
            if neigh in freeAlpha:
                freeAlpha.remove(neigh)
            '''else:
                if neigh[0] < 0: cpneigh = (neigh[0] + self.boardSize, neigh[1])
                if neigh[1] < 0: cpneigh = (neigh[0], neigh[1] + self.boardSize)
                freeAlpha.remove(cpneigh)'''
        
            
        #update board composition stats
        self.totAtoms[0] = 31
        if type == 1:
            self.totAtoms = [31, 7, 24, 0]
            self.actualComp = [31./self.totAtomSites, 7./self.totAtoms[0], 24./self.totAtoms[0], 0]
        else:
            self.totAtoms = [31, 0, 0, 31]
            self.actualComp= [31./self.totAtomSites,0, 0, 31./self.totAtoms[0]]
            self.neighbor_stats = zeros((3,9))
            
        #Now each site picked needs to decide whether is isolated or not...
        
        while len(freeAlpha) > 0:
        
            site_num = random.randrange(len(freeAlpha))
            site = freeAlpha[site_num]
            del(freeAlpha[site_num]) #no longer free
            
            if self.board[site[0], site[1]] !=0 : continue
            #if site was occupied, do not handle it. ALready removed from freeAlpha 
            
            #1st determine whether it is isolated, edge, or boundary. Look at Alpha Neighbours
            AlphaNeighbors = self.get_Alpha_Neighbors(site)
            AlphaKinds = [0,0,0,0] # Empty, Li, Mn, Co
            for AN in AlphaNeighbors:
                AlphaKinds[int(self.board[AN[0], AN[1]])] += 1
                            
            if AlphaKinds[0] == 6:
                #6 empty neighbors, isolated atom
                self.fill_isolated(site, self.getNeighbors(site), self.get_Alpha_Neighbors(site))
            elif AlphaKinds[1]*AlphaKinds[3] == 0:
                #Only has Li or Co neighbors (product of both kinds is zero)
                # It is an edge site => Grow grain
                if AlphaKinds[1] != 0:
                    type = 1
                else:
                    type = 3
                self.fill_edge(site, type,  self.getNeighbors(site))
            # Boundary sites are filled by random atoms kind. Likelyhood proportional to number of like neighbors
            elif random.random() <= AlphaKinds[1]/(6-AlphaKinds[0]):
                self.fill_boundary(site, 1, self.getNeighbors(site), self.get_Alpha_Neighbors(site))
            else:
                self.fill_boundary(site, 2, self.getNeighbors(site), self.get_Alpha_Neighbors(site))
            
        self.get_stats()
            
        if debug:
            for site in freeAlpha:
                if self.board[site[0], site[1]] ==0 :
                    self.board[site[0], site[1]] = 4
                
                
    def getNeighbors(self, site):
        '''
        Returns list of 1st neighbors to ALPHA site (x,y)
        Takes care of periodic boundary conditions
        Used to decide site type
        '''
        
        if site[0] > self.boardSize-5:
            a = site[0] - self.boardSize
        else:
            a=site[0]
        
        if site[1] > self.boardSize-5:
            b= site[1]-self.boardSize
        else:
            b = site[1]
        
        ret = [(a+1,b-1),(a+2,b+1),(a+1,b+2),(a-1,b+1),(a-2,b-1),(a-1,b-2)]
        
        return ret
    
    def get_Alpha_Neighbors(self, site):
        '''
        Returns list of 2nd neighbors to ALPHA site (x,y)
        Takes care of periodic boundary conditions
        Used to decide site type
        '''
        
                 
        neighbor_map = { \
                        "(1, 1)": lambda a,b : [(a+1,b-1),(a+2,b+1),(a+1,b+2),(a-1,b+1),(a-2,b-1),(a-1,b-2)],\
                        "(2, 1)": lambda a,b : [(0,b-1),(1,b+1),(0,b+2),(a-1,b+1),(a-2,b+1),(a-1,b-2)],\
                        "(1, 2)": lambda a,b : [(a+1,b-1),(a+2,0),(a+1,1),(a-1,0),(a-2,b-1),(a-1,b-2)],\
                        "(2, 2)": lambda a,b : [(0,b-1),(1,0),(0,1),(a-1,0),(a-2,b+1),(a-1,b-2)],\
        }
        
        if site[0] > self.boardSize-5:
            a = site[0] - self.boardSize
        else:
            a=site[0]
        
        if site[1] > self.boardSize-5:
            b= site[1]-self.boardSize
        else:
            b = site[1]
        
        ret = [(a+3,b),(a+3,b+3),(a,b+3),(a-3,b),(a-3,b-3),(a,b-3)]
        
        return ret

    def fill_isolated(self, site, neighbors, alpha_neighbors):
        ''' Fills an isolated Alpha site and all neighbors with atom drawn random to move close to target composition.
        '''
        if self.totAtoms[3] == 0:
            type = 3
        else:
            #What if I pick Li?
            Li_Li2Co = (self.totAtoms[1] + 7.)/ (self.totAtoms[3])
            Li_err = abs(Li_Li2Co - self.targetComp[1]/self.targetComp[3])
            #What if I pick Co?
            Mn_Li2Co = (self.totAtoms[1])/ (self.totAtoms[3]+31.)
            Mn_err = abs(Mn_Li2Co - self.targetComp[1]/self.targetComp[3])
            # Prob is inversely prop to error
            Tot_err = Mn_err + Li_err
            if (random.random() * Tot_err) < Li_err:
                #pick Co
                type = 3
            else:
                #pick Li
                type = 1
        self.board[site[0], site[1]] = type
        for neigh in alpha_neighbors:
            self.fill_site(neigh, type, self.getNeighbors(neigh))
         
    def fill_edge(self, site, type, neighbors):
        #only has one kind of neighbors
        self.fill_site(site, type, self.get_empty_neighbors(neighbors))
    
    def fill_boundary(self, site, type, neighbors, alpha ):
        self.fill_site(site, type, self.get_empty_neighbors(neighbors))
    
    def fill_site(self, site, type, empty_neighbors):
        ''' called by all fill_type
        returns number and kind of atoms written'''
        if type == 1: # Li
            self.board[site[0], site[1]] = 1
            for nsite in empty_neighbors:
                self.board[nsite[0], nsite[1]] = 2 #Mn surrounds Li
            ret = [len(empty_neighbors)+1, 1, len(empty_neighbors), 0  ] #[total, Li atoms, Mn, Co]
        else: # type = 3, Co
            self.board[site[0], site[1]] = 3
            for nsite in empty_neighbors:
                self.board[nsite[0], nsite[1]] = 3 #Co surrounds Co
            ret = [len(empty_neighbors)+1, 0, 0, len(empty_neighbors) +1  ] #[total, Li atoms, Mn, Co]
        return ret

    def get_empty_neighbors(self, neighbors):
        empty_neighbors =[]
        for neigh in neighbors:
            if self.board[neigh[0], neigh[1]] == 0:
                empty_neighbors.append(neigh)
        return empty_neighbors
        
    def to_Jmol(self, ofilename):
        ''' Provides .xyz representation of board, for debug purposes '''
        a = 1.6166 * array([1,0])
        b = 1.6166 * array([cos(radians(120)), sin(radians(120))])
        #Note: periodicity of lattive needs to be adjusted so nearest occupied neighbor dist = 2.8 angs
        ofile = open(ofilename, 'w')
        #ofile.write("%d\n\n" %(self.totAtoms[0] ))
        out_stream = []
        out_count = 0
        for j in range(self.boardSize):
            for i in range(self.boardSize):
                t =  int(self.board[i,j])
                if t !=0:
                    r = i * a + j * b
                    out_stream.append("%s\t" %(self.pieces[t]) +"%.3f\t%.3f\t" %tuple(r) + "0.0\n")
                    out_count +=1
                elif debug:
                    r = i * a + j * b
                    out_stream.append("H\t" +"%.3f\t%.3f\t" %tuple(r) + "0.0\n")
                    out_count +=1
        ofile.write("%d\n\n" %(out_count))
        for s in out_stream: ofile.write(s)        
        ofile.close()
    
    def to_hex_rep(self, ofilename, type=2):
        ''' Provides hexagonal representation of board by applying symmetry translations to unit cell.
        Origin of cell becomes center of hexagon.
        Allows export to xyz file (type = 1) and Jmol (type = 2, default)'''
        a = 1.6166 * array([1,0])
        b = 1.6166 * array([cos(radians(120)), sin(radians(120))]) 
        hex_coords = []
        hex_type = []
        num = self.boardSize
        for j in range(self.boardSize):
            for i in range(self.boardSize):
                t =  int(self.board[i,j])
                if t !=0:
                    r= i*a + j*b
                    hex_coords.append(r)
                    hex_type.append(self.pieces[t])
                    r = (i-num) * a + (j-num)*b
                    hex_coords.append(r)
                    hex_type.append(self.pieces[t])
                    if j>=i:
                        r = i*a + (j-num)*b
                        hex_coords.append(r)
                        hex_type.append(self.pieces[t])
                    if i>= j:
                        r=(i-num)*a + j*b
                        hex_coords.append(r)
                        hex_type.append(self.pieces[t])                                             
        ofile=open(ofilename, 'w')
        
        num = len(hex_coords)
        if type == 1:            
            ofile.write("%d\n\n" %(num))
            for j in range(num):
                ofile.write(hex_type[j] + "\t%.3f\t%.3f\t" %tuple(hex_coords[j]) + "0.0\n")
        else:
            for j in range(num):
                ofile.write("PutAtom_" + hex_type[j] + "(%.3f, %.3f, 0.0)\n" %tuple(hex_coords[j]))
        ofile.close()
        
    
    def get_stats(self):
        ''' saves game board to file 
        wipes the stats table first, for safety
        self.neighbor_stats[i,j] :
        i: Li, Mn, Co; central atom type
        j is stat: number, numb_Li-neigh, numb_Mn-neigh, numb_Co_neigh,
                   numb_Li_2ndneigh, numb_Mn_2nd_neigh, numb_co_2nd_neigh,
                   numb_bulk-like_atoms of i type, numb_edge_atoms of i type'''
        self.neighbor_stats = zeros((3,9))
        for x in range(self.boardSize):
            for y in range(self.boardSize):
                site = (x,y)
                kind = self.board[x,y]
                if kind == 0:
                    continue
                self.neighbor_stats[kind - 1, 0] += 1
                neigh_1st = self.getNeighbors(site)
                neigh_2nd = self.get_Alpha_Neighbors(site)
                bulk = True
                for nsite in neigh_1st:
                    nkind = int(self.board[nsite[0],nsite[1]])
                    self.neighbor_stats[kind - 1, nkind] +=1
                for nsite in neigh_2nd:
                    nkind = int(self.board[nsite[0],nsite[1]])
                    self.neighbor_stats[kind - 1, nkind + 3] +=1
                #also check whether edge or bulk atom:
                bulk = True
                if kind != 3:
                    #Li or Mn bulk if no Co arround
                    for nsite in neigh_1st:
                        if int(self.board[nsite[0],nsite[1]]) == 3:
                            #Co neighbor
                            bulk = False
                            break
                else:
                    #Co, bulk if only Co arround
                    for nsite in neigh_1st:
                        if int(self.board[nsite[0],nsite[1]]) != 3:
                            #Co neighbor
                            bulk = False
                            break
                if bulk:
                    self.neighbor_stats[kind-1,7] +=1
                else:
                    self.neighbor_stats[kind-1,8] +=1                    
                              
        for i in range(3):
            for j in range(8):
                self.neighbor_stats[i,j+1] /= self.neighbor_stats[i,0]
                #normalize stats (number of nieghbors of each knid) to number of central atoms of current kind
                #and bulk and edge atoms to fraction
            
        return self.neighbor_stats
    
    def harvest_grains(self, color=False, ofname = ""):
        '''Counts number of disconected grains in a board'''
        grains=[[],[]]
                
        #Each element of grains[0] is a list of coordinates of atoms in a LiMn grain, grains[1] same for Co grains
        grain_map = [[[3,0] for i in range(self.boardSize)] for j in range(self.boardSize)]
        # grain_map[x][y][t] gives info on atom at x,y. t=0 is grain type: 0 for LiMn, 1 for Co. 
        #t=1 is grain number of correct type; i.e. grains[ grain_map[x][y][t=0]][ grain_map[x][y][t=1] ] =
        # = list of [x,y] of similar atoms connected to it 
        
        current_grain=[4,0]
        # grain types 3,4 are shorthand for no type
        
        #auxiliary function to check for connectivity of current grain
        def neigh_grains(x,y,type):
            site = [x,y]
            if site[0] > self.boardSize-5:
                a = site[0] - self.boardSize
            else:
                a=site[0]
            if site[1] > self.boardSize-5:
                b= site[1]-self.boardSize
            else:
                b = site[1]
            low_neigh=[[a+1,b-1],[a-2,b-1],[a-1,b-2]]
            
            neigh_grainlist = []
            for neigh in low_neigh:
                if type == grain_map[neigh[0]][neigh[1]][0]:
                    neigh_grain = grain_map[neigh[0]][neigh[1]][1]
                    if not(neigh_grain in neigh_grainlist): 
                        neigh_grainlist.append(neigh_grain)
            #now I have a list of all grains of same type neighboring this grain
            return neigh_grainlist
        
        def fuse_grains(g1, g2, type):
            if g1 == g2:
                pass
            else:
                for atom in grains[type][g2]:
                    grains[type][g1].append(atom)
                    grain_map[atom[0]][atom[1]][1] = g1
                grains[type][g2]=[]
        
            
            
            
        #start walking the board
        for x in range(self.boardSize):
            for y in range(self.boardSize):
                t = self.board[x,y]
                if t == 1 or t == 2:
                    g = 0
                    #g is grain type, 0 = LiMn
                elif t == 3:
                    g = 1 #Co
                else:
                    continue
                
                #assign grain type
                grain_map[x][y][0]=g
                
                if g==current_grain[0]:
                    # Same grain as previous atom scanned
                    grain_map[x][y][1]=current_grain[1]
                    grains[g][current_grain[1]].append([x,y])
                    #check neighbors for connectivity
                    ng = neigh_grains(x, y, g)
                    while current_grain[1] in ng:
                        ng.remove(current_grain[1])
                    l=len(ng)
                    if l==0:
                        #no neighboring grains
                        pass
                    elif l == 1:
                        #one neighboring grain => fuse 'em
                        fuse_grains(current_grain[1], ng[0], g)
                    else:
                        fuse_grains(ng[0],ng[1],g)
                        fuse_grains(current_grain[1], ng[0], g)
                        #print "Something went wrong when scanning grains\n"
                        #print "An atom of current grain found 2 or more neighboring grains"
                        
                else:
                    # grain ended, other phase reached
                    current_grain[0]=g
                    #check neighbors
                    ng = neigh_grains(x, y, g)                    
                    l=len(ng)
                    if l==0:
                        #no neighboring grains
                         current_grain[1] = len(grains[g])
                         grain_map[x][y][1] =current_grain[1]
                         grains[g].append([[x,y]])
                    elif l ==1:
                        #neighboring existing grain: add to it
                        current_grain[1] = ng[0]
                        grain_map[x][y][1] =current_grain[1]
                        grains[g][ng[0]].append([x,y])
                    elif l ==2:
                        #bridging two grains => fusse
                        current_grain[1] = ng[0]
                        fuse_grains(ng[0], ng[1], g)
                        grain_map[x][y][1] =current_grain[1]
                        grains[g][ng[0]].append([x,y])
                    else:
                        print "Something went wrong when scanning grains\n"
                        print "got three neighboring atoms assigned to different grains of same type"
        
        #wrap up grain stats results
        LiMn_g = []
        Co_g = []
        for g_list in grains[0]:
            l = len(g_list)
            if  l > 0:
                LiMn_g.append(l)
                
        LiMn_g.sort()
        LiMn_hist =[]
        while len(LiMn_g)>0:
            size = LiMn_g[0]
            num = LiMn_g.count(size)
            LiMn_hist.append([size, num])            
            for i in range(num):
                LiMn_g.remove(size)
        for g_list in grains[1]:
            l = len(g_list)
            if  l > 0:
                Co_g.append(l)
        Co_g.sort()
        Co_hist =[]
        while len(Co_g)>0:
            size = Co_g[0]
            num = Co_g.count(size)
            Co_hist.append([size,num])
            for i in range(num):
                Co_g.remove(size)
                
        #if need to color grains, do so
        if color:
            a = 1.6166 * array([1,0])
            b = 1.6166 * array([cos(radians(120)), sin(radians(120))]) 
            color_file = open(ofname, 'w')
            for i in [0,1]:
                num=0
                for j in range(len(grains[i])):
                    l=len(grains[i][j])
                    if l>0:
                        num = num+1
                        for coord in grains[i][j]:
                            xycoord = coord[0]*a+coord[1]*b
                            color_str = "PutAtom_(%d, %d, " %(i,l)
                            color_str = color_str + "%.3f, %.3f, 0.0)\n" %tuple(xycoord)
                            color_file.write(color_str)
            color_file.close()
        
        return [LiMn_hist, Co_hist]
                     
                    
             
        
if __name__ == '__main__':
#    stats = open('stast.dat','w')
#    for i in range(4):
#        ofname = "test_" + str(i)+".xyz"
#        kk = GameBoard(20)
#        kk.to_Jmol(ofname)
#        stats.write('Stats for %s:\n' %(ofname))
#        stats.write(str(kk.totAtoms) +"\n")
#        stats.write(str(kk.actualComp) +"\n\n")
#    stats.close() 
#    kk = GameBoard(20)
#    file = open('pickle3.txt','w')
#    pickle.dump(kk, file)
#    file.close()
#    file=open('pickle3.txt','r')
#    board = pickle.load(file)
#    file.close()
#    board.to_Jmol('pickle3.xyz')
#    pprint.pprint(kk.neighbor_stats.transpose())
#    
#    print "That's all folks"

    '''self.neighbor_stats[i,j] :
        i: Li, Mn, Co; central atom type
        j is stat: number, numb_Li-neigh, numb_Mn-neigh, numb_Co_neigh,
                   numb_Li_2ndneigh, numb_Mn_2nd_neigh, numb_co_2nd_neigh,
                   numb_bulk-like_atoms of i type, numb_edge_atoms of i type'''

    result = []
    g_count = []
    threshold = 5
    #grain counted if > threshold atoms
    for size in range(5,26,2):
      for k in range(15):  
        for i in range(50):
            part_result = [size]
            kk = GameBoard(size)
            tot =0
            for j in range(3): tot+= kk.neighbor_stats[j,0]
            for j in range(3): part_result.append(kk.neighbor_stats[j,0] / tot)
            #part_result[1] +=1
            part_result.append(kk.neighbor_stats[1,2] + kk.neighbor_stats[1,3])
            part_result.append(kk.neighbor_stats[1,7])
            part_result.append(kk.neighbor_stats[2,2] + kk.neighbor_stats[2,3])
            part_result.append(kk.neighbor_stats[2,7])
            result.append(part_result)
            #now get grain count
#            if kk.neighbor_stats[0,0] == 0 or kk.neighbor_stats[1,0] == 0:
#                x_Mn = 0
#                x_Li = 0
#            else:                
#                r_Mn = kk.neighbor_stats[2,0] / 1.* kk.neighbor_stats[1,0]
#                r_Li = kk.neighbor_stats[2,0] / 1.* kk.neighbor_stats[0,0]
#                x_Mn = 3.0 / (3+2*r_Mn)
#                x_Li = 3.0 / (3+r_Li)
#                x = 0.5 * (x_Mn + x_Li)
#                Dx = 0.5* (x_Mn - x_Li)
            x = 1-part_result[3]
            Dx = (part_result[1]/part_result[2])
            
            harvest = kk.harvest_grains()
            Li_grains = 0
            Co_grains = 0
            g_sizes=[]
            for grain in harvest[0]:
                if grain[0] > 5:
                    Li_grains += grain[1]
                    for j in range(grain[1]):
                        g_sizes.append(grain[0])
            for grain in harvest[1]:
                if grain[0] > threshold:
                    Co_grains += grain[1]
                    for j in range(grain[1]):   
                        g_sizes.append(grain[0]) 
            g_count.append([x, Dx, size, Li_grains, Co_grains]+g_sizes)
            
            
            
        print "size " +str(size) + " 50 x " + str(k) + " done."
        print "size " + str(size) + " done"
    print "Calculation done, writing results"
    
    
#    opick = open("big_test.pic",'w')
#    pickle.dump(result, opick)
#    opick.close()
#    
    odata = open("Coordination.dat",'w')
    odata.write("LiCoMn game results file\n\n")
    odata.write('size\t[Li]\t[Mn]\t[Co]\tMn-M\t[bulk_Mn]\tCo-M\t[bulk_Co]\n\n')
    for line in result:
        for datum in line:
            odata.write("%f\t" %(datum))
        odata.write("\n")
    odata.close()
    
    odata = open("Grains.dat",'w')
    odata.write("LiCoMn game results file\n")
    odata.write("Grain count constrained to grain size >" +str(threshold) +" atoms\n\n")
    odata.write('x=1-[Co]\t[Li]/[Mn]\tsize\tLi_grains\tCo_grains\tGrain sizes (in order)\n\n')
    for line in g_count:
        for datum in line:
            odata.write("%.4f\t" %(datum))
        odata.write("\n")
    odata.close()
    
    print "Results files written. Time for coffee"


#    file=open('pickle3.txt','r')
#    board = pickle.load(file)
#    file.close()     
#    board.to_hex_rep("pickle_3hex.pov", 2)     
#    board.to_hex_rep("pickle_3hex.xyz", 1)
#       

#    statf = open('grain_stats.dat','w')
#    for i in range(5):
#        ofname = "test_" + str(i)
#        kk = GameBoard(20)
#        kk.to_Jmol(ofname+".xyz")
#        pickle_file = open(ofname+".pic",'w') 
#        pickle.dump(kk, pickle_file)
#        pickle_file.close()
#        
#        stat = kk.harvest_grains(True, ofname+".pov" )
#        statf.write("Grain stats on model #" + str(i) +"\n")
#        statf.write("LiMn6-like grains:\n")
#        for item in stat[0]:
#            statf.write("size = %d, number = %d\n" %tuple(item))
#        statf.write("Co-like grains:\n")
#        for item in stat[1]:
#            statf.write("size = %d, number = %d\n" %tuple(item))
#        statf.write("\n")
#    statf.close()
#    print "Done"   