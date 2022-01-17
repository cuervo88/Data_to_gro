#!/usr/bin/python3

#Code made for converting into a gro file from a dat file generated from h5: h5dump -d /particles/all/position/value

import pandas as pd
import numpy as np
import argparse
import sys
import os

parser = argparse.ArgumentParser(prog ="real_gross" ,description='Transform a data file of coordinates into a gro file f\
ormat for AZT simulations')
parser.add_argument('-f','--file',type=argparse.FileType('r',encoding='UTF-8'), required=True, help='File location')
parser.add_argument('-o','--output',type=str, help='Output file and location')
parser.add_argument('-az','--azt', type=int,required=True, help='Number of azobenzene-containing surfactant molecules')
parser.add_argument('-at','--atoms',type=int,required=True, help='Total number of atoms per time frame')
parser.add_argument('-tf','--tf', type=int, required=True, help='Number of time frames')
parser.add_argument('-b','--box',type=int, required=True, help='Box size of the simulation')
parser.add_argument('-tx', '--text', type=str, help='Description of the experiment')

args = parser.parse_args()
file_path = args.file
atom_num = int(args.atoms)
AZT_names= ["T", "T", "C", "B", "C", "S", "C", "B", "C", "O", "M", "N", "N", "N"]
AZT_molecules  = args.azt
ION = int(AZT_molecules)
AZT = len(AZT_names) * AZT_molecules
SOL = int(atom_num) - (AZT+ION)

text_name = "No file name supplied. All hail Ignacio"
if args.text is not None:
  text_name = args.text
time_frame = args.tf
box_size = args.box
file_name = file_path.name[:-3]+"gro"
if args.output is not None:
  file_name = args.output
if os.path.exists(file_name):
  raise Exception ("Output filename already exist. Please errase it or provide a new output name")
a = np.repeat(np.arange(1,AZT_molecules +1),len(AZT_names))
b = np.arange(AZT_molecules +1,AZT_molecules+SOL+ION+1)
c = list(a)+list(b)
d = AZT_names*AZT_molecules  + ['W']*SOL + ['F']*ION
e = ["AZT"]*AZT+ ["SOL"]*SOL + ['ION']*ION
f = np.arange(1,(AZT+SOL+ION)+1,1)
data1 = pd.read_csv(file_path, header=None, chunksize= (atom_num))
for i in range(time_frame):
  m = 1
  if i == 0:
    m=0
  skip_row = ((atom_num)*i)
  data = next(data1)
  data[['A',2]] = data[2].str.split(':', expand=True)
  data[2] = pd.to_numeric(data[2])
  data = data[[2,3,4]]
  data = data.applymap("{0:8.3f}".format)
  atoms = pd.DataFrame(index=range((atom_num+3)), columns=range(7))
  atoms[0][0] = "name"
  atoms[0][1] = atom_num
  atoms[0][2:(atom_num)+2] = ['{:>5}'.format(item)[-5:] for item in c]
  atoms[1][2:(atom_num)+2] = ['{: <5}'.format(item) for item in e]
  atoms[3][2:(atom_num)+2] = ['{:>5}'.format(item)[-5:] for item in f]
  atoms[2][2:(atom_num)+2] = ['{: >5}'.format(item) for item in d]
  atoms[4][2:(atom_num)+2] = data[2]
  atoms[5][2:(atom_num)+2] = data[3]
  atoms[6][2:(atom_num)+2] = data[4]
  atoms[0][atom_num+2] = str(box_size) + " " + str(box_size) + " " + str(box_size)
  atoms = atoms.fillna(" ")
  atoms2 = pd.DataFrame(index=range(atom_num+3), columns=range(1))
  atoms2[0]= atoms[0].apply(str)+ atoms[1].apply(str)+atoms[2].apply(str)+atoms[3].apply(str)+atoms[4].apply(str)+atoms[5].apply(str)+atoms[6].apply(str)
  atoms2.to_csv(file_name, header=None, mode='a', index=None)
  del data, atoms
