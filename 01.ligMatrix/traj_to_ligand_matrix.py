#################################################
try:
    import BiKi
except ImportError:
    print("The 'BiKi' module is not installed. Try:")
    print(". path-to-biki/BiKiLifeSciences/enableEssentials.source")
    print("or refererr to the BiKi manual")
    exit(1)
#################################################
import os
import argparse

# parse input
parser = argparse.ArgumentParser(description='This script align the protein part and export the aligned ligand trajectory from your initial trajectory. The selection is the same as in VMD')
parser.add_argument('--traj','-t', type=str, required=True, help='xtc file with the trajectory')
parser.add_argument('--ref','-r', type=str, required=True, help='reference file - gro or pdb')
parser.add_argument('--lig','-l', type=str, required=True, help='ligand name in reference file (e.g.: "AAA")')
parser.add_argument('--pro','-p', type=str, required=True, help='protein part to use for alignment (e.g.: "backbone", "resname AAA")')


args = parser.parse_args()

traj = args.traj
ref_file = args.ref
ligand_name = args.lig
pro_align = args.pro



# Load the box and trajectory
box = BiKi.Structure()
box.load(ref_file)
traj_loader = BiKi.TrajectoryLoader()
traj_loader.setStructure(box)
traj_loader.addTrajectory(traj)

print("Saving to file...")
with open("./ligandMatrix.txt", "w") as pw:
    N = traj_loader.getSize()

    temp = traj_loader.getFrame(0)
    ref = BiKi.Structure()
    ref.copy(temp)

    ms_ligand = BiKi.MultiStructure()
    for i in range(N):
        print("{}/{}".format(i, N))

        frame = traj_loader.getFrame(i)

        ad = BiKi.AlignmentData()

        frame_orig = BiKi.Structure()
        frame_orig.copy(frame)

        ligand = BiKi.Structure()
        ligand.copy(frame.select("resname " + ligand_name))

        ref.select(pro_align).align(frame.select(pro_align), ad)

        ligand.rotoTranslate(ad)
        ms_ligand.addStructure(ligand)

        NA = ligand.getSize()
        for j in range(NA):
            pos = ligand.getPosition(j)
            pw.write(str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " ")
        pw.write("\n")

# save traj of ligand
ms_ligand.saveMolecules("traj_ligand.pdb")
