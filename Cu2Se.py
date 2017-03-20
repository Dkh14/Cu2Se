from ase import *
from ase import Atoms
from ase.build import bulk 
from ase.visualize import view
import numpy as np
import random
a=5.787
#construct ideal cell
Cu2Se=Atoms('Se4Cu8',
positions=[(0,0,0),(0,.5,.5),(.5,.5,0),(.5,0,.5),
(.25,.25,.25),(.25,.25,.75),(.25,.75,.25),(.25,.75,.75),(.75,.25,.25),(.75,.25,.75),(.75,.75,.25),(.75,.75,.75)],cell=(1,1,1),  pbc=True)
Cu2Se.set_cell((a,a,a),scale_atoms=True)
rep=np.array([2,2,2])
Cu2Se*=rep
#displace Cu atoms
m=12*rep[0]*rep[1]*rep[2]
i=4
n=np.array([0,0,0])
displace=np.array([0,0,0])
for i in range (4,m) :
    if i%12<4 :
        continue
    elif i%12==4 :
        n=(Cu2Se.positions[i]-a/4)/a #does this work?
        nx=n[0]
        ny=n[1]
        nz=n[2]
    r=random.random()
    temp=Cu2Se.positions[i]
    if r<.32 :
        displace=np.array([0,0,0])

    elif r>=.32 and r<=.87 :
        displace=np.array([.059*a,.059*a,.059*a])
        if r>=.32 and r<=.45 :
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
                displace*=(1,1,1)
            else:
                displace*=(-1,-1,-1)   
        elif r>=.46 and r<=.59 :
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
               displace*=(-1,-1,1)
            else : 
               displace*=(1,1,-1)
        elif r>=.60 and r <= .73 :
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
               displace*=(1,-1,-1)
            else:
               displace*=(-1,1,1)
        else:
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
               displace*=(-1,1,-1)
            else:
               displace*=(1,-1,1)

    else :
        displace = np.array([.134*a,.134*a,.134*a])
        if r>=.88 and r<=.90 :
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
               break
            else :
               displace*=(-1,-1,-1)
        elif r>=.91 and r<=.93 :
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
               displace*=(-1,-1,1)
            else :
               displace*=(1,1,-1)
        elif r>= .94 and r<=.96 :
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
               displace*=(1,-1,-1)
            else :
               displace*=(-1,1,1)
        else :
            if ((temp==[nx*a+a/4,ny*a+a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+3*a/4,nz*a+a/4]).all() or (temp==[nx*a+3*a/4,ny*a+a/4,nz*a+3*a/4]).all() or (temp==[nx*a+a/4,ny*a+3*a/4,nz*a+3*a/4]).all()):
               displace*=(-1,1,-1)
            else :
               displace*=(1,-1,1)
    
    Cu2Se.positions[i]+=displace
view(Cu2Se)
Cu2Se.write("Cu2Se.xyz") #save to.xyz so it can be converted to car
