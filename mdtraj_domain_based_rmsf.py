"""
    Author: Halil ibrahim özdemir
    Loc: Marmara University / Bioengineering
"""


import mdtraj as md
import matplotlib.pyplot as plt


def get_rmsf_data(top, traj, start_frame, stop_frame, stride, selection, superpose, group_selections=None):
    """
        Make RMSF Calculation with Group Selections with mdtraj.

            Parameters
            ----------

            top: topology file

            traj : trajectory file

            start_frame : int
                Include after this snapshot for your analysis.

            stop_frame: int
                Include until this snapshot for your analysis

            stride: int
                It will take snapshots at intervals of the "stride" unit you specified for your analysis

            selection: str
                Normal Selection

            superpose: bool
                Will align Snapshots for best fitting

            group_selections: list
                list of domains selections

            Example
            ----------

            or_rmsf, domain_rmsf, time, residue_list = get_rmsf_data(top='test/protein.pdb', traj='test/50_frame.dcd',
                                                         start_frame=0, stop_frame=49, stride=1,
                                                         selection='backbone and name CA', name='aaa', superpose=True,
                                                         group_selections=["backbone and name CA and resid 0 to 20",
                                                                           "backbone and name CA and resid 21 to 59",
                                                                           "backbone and name CA and resid 60 to 100",
                                                                           "backbone and name CA and resid 101 to 115",
                                                                           "backbone and name CA and resid 116 to 142"])

        """

    global groupselections

    def group_selection(traj_universe, reference_universe, selection, group_selection):

        new_ref = reference_universe.atom_slice(topology.select(selection + ' and ' + group_selection))
        new_traj = traj_universe.atom_slice(topology.select(selection + ' and ' + group_selection))

        new_traj.superpose(reference=new_ref, parallel=True)
        print(new_traj)
        return md.rmsf(new_traj, new_ref, parallel=True) * 10

    try:
        traj_origin = md.load(traj, top=top, stride=stride)
        ref_origin = md.load(top)
        topology = traj_origin.topology

        traj = traj_origin.atom_slice(topology.select(selection))
        ref = ref_origin.atom_slice(topology.select(selection))

        if traj.n_atoms != ref.n_atoms:
            traj = traj.atom_slice(topology.select(selection))

        if group_selections is not None:
            groupselections = [group_selection(traj_origin, ref_origin, selection, s) for s in
                               group_selections]


        if superpose:
            traj.superpose(reference=ref, parallel=True)

        if (start_frame and stop_frame) is not None:
            traj = traj[start_frame:stop_frame]

        elif stop_frame is not None:
            traj = traj[:stop_frame]

        elif start_frame is not None:
            traj = traj[start_frame:]

        all_rmsf_data_struct = {'origin_RMSF': md.rmsf(traj, ref, parallel=True) * 10,
                                'groupSelection_RMSF': groupselections, 'time': traj.time,
                                'residues': list(range(0, traj.n_atoms))}

        ori_selection_rmsf = list(all_rmsf_data_struct['origin_RMSF'])
        selection_rmsf = all_rmsf_data_struct['origin_RMSF']
        domain_base_rmsf = all_rmsf_data_struct['groupSelection_RMSF']

        count = 0
        for i in range(len(groupselections)):
            selection_rmsf[count: len(domain_base_rmsf[i])+count] = domain_base_rmsf[i]
            count = len(domain_base_rmsf[i]) +count

        return ori_selection_rmsf, selection_rmsf, all_rmsf_data_struct['time'], all_rmsf_data_struct['residues']

    except Exception as Error:
        print(Error)
        print("problem in rmsf calculation")


# ---> For Plot
# fig = plt.figure(figsize=(4, 4))
# ax = fig.add_subplot(111)
# ax.plot(residue_list, domain_rmsf, 'r-', label="domain_based_rmsf")
#
# ax.plot(residue_list, or_rmsf, 'g-', label="Original (All)")
#
# ax.legend(loc="best")
# ax.set_xlabel("Residues")
# ax.set_ylabel(r"RMSF ($\AA$)")
# plt.show()
