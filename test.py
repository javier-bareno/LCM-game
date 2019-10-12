'''
Created on Mar 4, 2014

@author: bareno
'''

if __name__ == '__main__':
    path = "/Users/bareno/Desktop/ZHC/new_tables/"
    fname = path +"file.txt"
    f=open(fname,'w')
    f.write("Filed opened and written to")
    f.close()