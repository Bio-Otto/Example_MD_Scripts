
trj = '50_frame.dcd'
top = 'protein.pdb'

import MDAnalysis
import matplotlib.pyplot as plt
import numpy as np
from MDAnalysis.analysis.rms import rmsd

u = MDAnalysis.Universe(top, trj)

# one AtomGroup per domain
domains = {
    'Calpha RMSD': u.select_atoms("protein and (name C or name N or name CA)"),
    'All Atom RMSD (noH)': u.select_atoms("protein and not name H*"),
    }
colors = {'Calpha RMSD': 'green', 'All Atom RMSD (noH)': 'red'}

u.trajectory[0]   # rewind trajectory
xref0 = dict((name, g.positions - g.center_of_mass()) for name, g in domains.items())

nframes = len(u.trajectory)
results = dict((name, np.zeros((nframes, 2), dtype=np.float64)) for name in domains)

for iframe, ts in enumerate(u.trajectory):
    for name, g in domains.items():
        results[name][iframe, :] = (u.trajectory.frame,
                                    rmsd(g.positions, xref0[name],
                                         center=True, superposition=True))


# plot
fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
for name in "Calpha RMSD", "All Atom RMSD (noH)":
    data = results[name]
    ax.plot(data[:,0], data[:,1], linestyle="-", color=colors[name], lw=2, label=name)
ax.legend(loc="best")
ax.set_xlabel("Frame")
ax.set_ylabel(r"RMSD of Backbone ($\AA$)")

plt.show()
# for ext in ('svg', 'pdf', 'png'):
#     fig.savefig("AdK_domain_rigidity.{0}".format(ext))