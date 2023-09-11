import BiKi
import ast
import argparse

# parse input
parser = argparse.ArgumentParser(description='This script extract physical frames from your initial trajectory')
parser.add_argument('--traj','-t', type=str, required=True, help='xtc file with the trajectory')
parser.add_argument('--gro','-g', type=str, required=True, help='gro file of the system used to generate the trajectory')
args = parser.parse_args()

traj_file = args.traj
gro_file = args.gro

# extract frame list from physical_medoids.txt
f=open('physical_medoids.txt',"r")
lines=f.readlines()
frame_list=[]
for j in lines:
    frame_list.append(int(j.split()[0]))
f.close()
print(frame_list)

# extract physical frames from initial trajectory
for m in frame_list:
     print(m, frame_list.index(m)+1)
     new_str = BiKi.Structure()
     new_str.load(gro_file)
     trajLoader = BiKi.TrajectoryLoader(new_str)
     trajLoader.addTrajectory(traj_file)
     frame_new = trajLoader.getFrame(m)
     frame_new.save(str(frame_list.index(m)+1)+".gro")
