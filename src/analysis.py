'''
Created on Oct 9, 2009

@author: bareno

Quick and dirty format parser to analyze LiCoMn data output

'''

import math

def parse_data(path):
    coord = open(path + "Coordination.dat", 'r')
    grains = open(path + "Grains.dat", 'r')
    data = []
    gdata = grains.readlines()
    cdata = coord.readlines()
    grains.close()
    coord.close()
    
    curr_size = 0
    curr_data = []
    for i in range(4, len(cdata)):
        size = float(cdata[i].split()[0])
        if size != curr_size:
            if len(curr_data) > 0:
                data.append([curr_size,curr_data])
            curr_data =[]
            curr_size = size                     
        
        tdata = cdata[i].split()[1:] + gdata[i+1].split()[0:2] + gdata[i+1].split()[3:5]
        curr_data.append(tdata)
    data.append([curr_size,curr_data])
    coord =[]
    grains=[]
    
    
    for dataset in data:
        coord_file = open(path+'coord_data_%2d.dat' %(dataset[0]),'w')
        coord_file_filtered = open(path+'coord_data_filt_%2d.dat'%(dataset[0]),'w')
        grain_set=[[0,0] for i in range(20)]    
        title_str = "[bulk_Mn] \t Mn-M \t[bulk_Co] \tCo-M\n"
        coord_file.write(title_str)
        coord_file_filtered.write(title_str)
       
        for data_line in dataset[1]:
            out_str = "\t".join(data_line[i] for i in [4,3,6,5]) + "\n"
            
            coord_file.write(out_str)
            if data_line[3] >= 0.4 and data_line[3] <= 0.6:
                coord_file_filtered.write(out_str)  
            
            j = int(round(float(data_line[2]) * 19))
            grain_set[j][0] +=1
            grain_set[j][1] += float(data_line[9])
              
        grains.append([dataset[0], grain_set])                       
        coord_file.close()
        coord_file_filtered.close()
    
    grain_file = open(path + "grain_stats.dat", 'w')
    grain_file.write("[Co]\tsize\tLiMn_graisn\n")
    for size in grains:
        for i in range(20):
            grain_file.write(str(size[0]) +"\t")
            grain_file.write(str((i+0.5)/20 ) + "\t")
            if size[1][i][0] > 0:
                grain_file.write(str(size[1][i][1] / size[1][i][0]) + "\n")
            else:
                grain_file.write("0\n")
    grain_file.close()
    
    grain_file = open(path + "grain_sampling.dat", 'w')
    grain_file.write("[Co]\tsize\tpoints\n")
    for size in grains:
        for i in range(20):
            grain_file.write(str(size[0]) +"\t")
            grain_file.write(str((i+0.5)/20 ) + "\t")
            grain_file.write(str(size[1][i][0]) + "\n")
            
    grain_file.close()   
 
def bin_x_y_data(data, num_bins, x_min_0 = 0, x_max_0 = 0):
    '''
    Bins data into num_bins equisized bins between x_min_0 and x_max_0
    If x valueues go out of this range, it gets expanded to actual dataset min and max
    data is a list of [x,y] values: [[x0,y0], [x1,y1], ...]
     
    output: x, <y>, Dy, num_points
    '''
    x_max = x_max_0
    x_min = x_min_0
    for data_point in data:
        if x_max < float(data_point[0]):
            x_max = float(data_point[0])
        if x_min > float(data_point[0]):
            x_min = float(data_point[0])
    
    bin_length = 1.0*(x_max - x_min)/num_bins
    data_out = [[(i+0.5)*bin_length, 0, 0, 0] for i in range(num_bins +1) ]
    
    for data_point in data:
        bin_num = int(math.floor((float(data_point[0])-x_min)/bin_length))
        data_out[bin_num][1] += float(data_point[1])
        data_out[bin_num][2] += float(data_point[1])**2
        data_out[bin_num][3] += 1
    #created extra bin for x_max, add to previous
    data_point = data_out[num_bins]
    data_out[-2][1] += data_out[-1][1]
    data_out[-2][2] += data_out[-1][2]
    data_out[-2][3] += data_out[-1][3]
    data_out.pop()
    
    for data_point in data_out:
        if data_point[3] != 0:
            data_point[1]/= data_point[3]
            data_point[2]/= data_point[3]
            data_point[2] = (data_point[2] - data_point[1]**2)** 0.5
        else:
            data_point[1] =0
            data_point[2] =0
            
    
    return data_out
    
def load_coord_data(fname):
    '''
    Loads coord data from file, returns in appropriate format for bin_x_y_data '''
    
    data_f = open(fname, 'r')
    data = data_f.readlines()
    data.pop(0)
    Mn_data =[]
    Co_data =[]
    for line in data:
        ln = line.split()
        Mn_data.append([float(ln[0]), float(ln[1])])
        Co_data.append([float(ln[2]), float(ln[3])])
    return [Mn_data, Co_data]
                
def grain_data(ifname, num_bins, x_min_0 = 0, x_max_0 = 0):
    '''
    Reads grain data and returns data binned in x
    for each size, [x, prob > 1grain, D_prob, tot_num models] '''
    ifile = open(ifname, 'r')
    if_data = ifile.readlines()
    ifile.close()
    
    gdata_bin=[]
    curr_size = 0
    curr_data = []
    for i in range(5, len(if_data)):
        data = if_data[i].split()
        size = float(data[2])
        if size != curr_size:
            if len(curr_data) > 0:
                gdata_bin.append([curr_size,curr_data])
            curr_data =[]
            curr_size = size  
        tdata = [data[0]]
        if float(data[3]) > 1:
            tdata.append(1.)
        else:
            tdata.append(0.)               
        curr_data.append(tdata)
    gdata_bin.append([curr_size,curr_data])
    
    g_data_out =[]
    for dataset in gdata_bin:
        sz = dataset[0]
        g_data_out.append([sz, bin_x_y_data(dataset[1], num_bins, x_min_0, x_max_0)])
    
    return g_data_out
    
def trace_frontier(gdata, p_threshold, min_sampling):
    ''' gdata is output from grain_data
    searches, for each size, x such that prob_>1grain > p_threshold
    maps x. Indicates if datapoint sampled > min_sampling times '''
    frontier = []
    for size in gdata:
        for i in range(len(size[1])):
            if size[1][i][3] < min_sampling:
                val = -1
                # point is undersampled
            else:
                val = size[1][i][1]
#            elif size[1][i][1] > p_threshold:
#                val = 2
#                # prob(>1 grain) > threshold
#            else:
#                val = 1
#                # prob(>1 grain) < threshold
            frontier.append([size[1][i][0], size[0], val ])
    return frontier
            
    
    
    
    
       
if __name__ == '__main__':
#     parse_data("C:\\Documents and Settings\\bareno\\Home-Local\\Eclipse\\LiCoMn_game\\src\\Game results 7 II\\")

    
#    data = load_coord_data("C:\\Documents and Settings\\bareno\\Home-Local\\Eclipse\\LiCoMn_game\\src\\Game results 7 II\\coord_data_13.dat")
#    Co_data = bin_x_y_data(data[1], 20, 0, 1)       
#    Mn_data = bin_x_y_data(data[0], 20, 0, 1)
#    
#    of = open("C:\\Documents and Settings\\bareno\\Home-Local\\Eclipse\\LiCoMn_game\\src\\Game results 7 II\\bin_Mn_coord_13.dat",'w')
#    of.write("X\tY\tDY\tNum\n")
#    for dp in Mn_data:
#        of.write("%.4f\t%.4f\t%.4f\t%d\n" %tuple(dp))
#    of.close()
#    of = open("C:\\Documents and Settings\\bareno\\Home-Local\\Eclipse\\LiCoMn_game\\src\\Game results 7 II\\bin_Co_coord_13.dat",'w')
#    of.write("X\tY\tDY\tNum\n")
#    for dp in Co_data:
#        of.write("%.4f\t%.4f\t%.4f\t%d\n" %tuple(dp))
#    of.close()

#    ifname = "C:\\Documents and Settings\\bareno\\Home-Local\\Eclipse\\LiCoMn_game\\src\\Game results 7 II\\Grains.dat"
#    ofname = "C:\\Documents and Settings\\bareno\\Home-Local\\Eclipse\\LiCoMn_game\\src\\Game results 7 II\\frontier.dat"
    ifname = "/Users/javierbareno/Documents/Eclipse/LiCoMn_game/src/Game results 7 II/Grains.dat"
    ofname = "/Users/javierbareno/Documents/Eclipse/LiCoMn_game/src/Game results 7 II/frontier.dat"
    
    data = grain_data(ifname, 20, 0, 1)
    front = trace_frontier(data, 0.95, 10)
    ofile = open(ofname,'w')
    ofile.write("x\tsize\tvalue\n")
    for datap in front:
        ofile.write("%.3f\t%.3f\t%.3f\n" %tuple(datap))
#    for datap in data:
#        for datapp in datap[1]:
#            ofile.write("%.3f\t" %(datap[0]))
#            ofile.write("%.3f\t%.3f\t%.3f\t%.3f\n" %tuple(datapp))
    ofile.close()
        
        