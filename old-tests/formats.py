'''
Created on Oct 6, 2009

@author: bareno
'''
from numpy import *
import random
import pickle
import pprint


def xyz_to_hex(ifname, ofname, type = 2):
    '''Reads xyz J_mol file and makes hex board out of it'''