# Requirement
import biotite
import biotite.structure as struc
import biotite.structure.io as strucio
from biotite.structure.io.xtc import XTCFile

import matplotlib.pyplot as plt
import numpy as np


def concatenate_parallel_simulation(trajectory_list):
    import mdtraj as md
    md.join(trajectory_list, discard_overlapping_frames=False, check_topology=True)


def rmsf_plot(topology, xtc_traj, start_frame=None, stop_frame=None):
    # Gromacs does not set the element symbol in its PDB files,
    # but Biotite guesses the element names from the atom names,
    # emitting a warning
    template = strucio.load_structure(topology)

    # The structure still has water and ions, that are not needed for our
    # calculations, we are only interested in the protein itself
    # These are removed for the sake of computational speed using a boolean
    # mask
    protein_mask = struc.filter_amino_acids(template)
    template = template[protein_mask]

    from mdtraj.formats.xtc import XTCTrajectoryFile
    # with XTCTrajectoryFile(xtc_traj) as f:

    xtc_file = XTCFile()
    xtc_file.read(xtc_traj, atom_i=np.where(protein_mask)[0], start=start_frame, stop=stop_frame)

    trajectory = xtc_file.get_structure(template)
    print(trajectory)
    time = xtc_file.get_time()  # Get simulation time for plotting purposes

    trajectory = struc.remove_pbc(trajectory)
    trajectory, transform = struc.superimpose(trajectory[0], trajectory)
    rmsd = struc.rmsd(trajectory[0], trajectory)

    figure = plt.figure(figsize=(6, 3))
    ax = figure.add_subplot(111)
    ax.plot(time, rmsd, color=biotite.colors["dimorange"])
    ax.set_xlim(time[0], time[-1])
    ax.set_ylim(0, 2)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("RMSD (Å)")
    figure.tight_layout()

    radius = struc.gyration_radius(trajectory)

    figure = plt.figure(figsize=(6, 3))
    ax = figure.add_subplot(111)
    ax.plot(time, radius, color=biotite.colors["dimorange"])
    ax.set_xlim(time[0], time[-1])
    ax.set_ylim(14.0, 14.5)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("Radius of gyration (Å)")
    figure.tight_layout()

    # In all models, mask the CA atoms
    ca_trajectory = trajectory[:, trajectory.atom_name == "CA"]
    rmsf = struc.rmsf(struc.average(ca_trajectory), ca_trajectory)

    figure = plt.figure(figsize=(6, 3))
    ax = figure.add_subplot(111)
    res_count = struc.get_residue_count(trajectory)
    ax.plot(np.arange(1, res_count + 1), rmsf, color=biotite.colors["dimorange"])
    ax.set_xlim(1, res_count)
    ax.set_ylim(0, 1.5)
    ax.set_xlabel("Residue")
    ax.set_ylabel("RMSF (Å)")
    figure.tight_layout()

    plt.show()


if __name__ == "__main__":
    # Convert_Trajectory_Format_From_Terminal_with_MDTraj : mdconvert <input_trajectory> -o <output_trajectory>
    # Put here the path of the files
    templ_file_path = 'C:\\Users\\HIbrahim\\Desktop\\OpenMM_PerTool\\examples\\lysozyme_md.pdb'
    traj_file_path = 'C:\\Users\\HIbrahim\\Desktop\\OpenMM_PerTool\\examples\\lysozyme_md.xtc'

    rmsf_plot(topology=templ_file_path, xtc_traj=traj_file_path)
