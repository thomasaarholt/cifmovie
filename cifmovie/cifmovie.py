from pathlib import Path # for handling files

import math
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read # to read in crystallographic files
from ase.data import covalent_radii
from ase.data.colors import jmol_colors

import plato.draw.matplotlib

from scipy.spatial.transform import Rotation as R

def load(prefix="", directory="", extension='cif'):
    atoms = []
    for f in Path(directory).glob(prefix + "*" + extension):
        atoms.append(read(f))
    return atoms



def test_load(N=200):
    "Create a dummy test dataset of N cif files"
    import math
    
    test_download()
    # Duplicate the cif file and make lots of changes to it "over time"

    atoms_list = []
    for i in range(N):
        cell = read("sto3.cif").repeat(3) # using repeat to make it bigger - stack in xyz
        atoms_list.append(cell)

    a_length = cell.cell[0,0].copy() # original lattice a parameter
    x_positions = cell.positions[:,0].copy() # original atoms x position

    for i, atoms in enumerate(atoms_list):
        atoms[74].number = int(1 + i/4) # make one atom change its atomic number throughout series
        atoms.cell[0,0] = a_length * (1 + 0.3*math.sin(i/N * 2*math.pi)) # change lattice a parameter
        atoms.positions[:,0] = x_positions * (1 + 0.3*math.sin(i/N * 2*math.pi)) # change x-position of atoms

        # Make one atom spin around its original centre #justforfun
        r = 0.5
        eul = np.exp(i/200 * 2*np.pi * 1j)
        dx = eul.real
        dy = eul.imag

        atoms[2].position[0] += r*dx
        atoms[2].position[1] += r*dy
        
    return atom_list



def render(atoms_list, renderer='matplotlib', dpi=72, figsize='auto', pixel_scale=20):
    renderer = plato.draw.matplotlib
        
    p = Path("cif_images")
    p.mkdir(exist_ok=True)
    
    scene = renderer.Scene() # just used to get size
    real_size = scene.size * pixel_scale/dpi
    figure = plt.figure(figsize=real_size, dpi=dpi)
    
    for i, atom in enumerate(atoms_list):
        figure.clf()
        pos = atom.positions
        centre = (pos.min(0) + pos.max(0))/2
        pos -= centre
        numbers = atom.numbers
        radii = covalent_radii[numbers]
        colors = jmol_colors[numbers]
        colors = np.hstack((colors, np.ones((len(colors), 1)))) # format colors as matplotlib expects

        rotation = R.from_euler('z', (5*math.sin(i/100*(2*np.pi))), degrees=True).as_quat()
        scene = renderer.Scene(
            zoom=2, 
            antialiasing=True, 
            rotation=rotation, # comment this line out to remove rotation
        )
        sphere = renderer.Spheres(positions=pos, radii=radii/2, colors=colors, pixel_scale=20,)
        scene.add_primitive(sphere)
        scene.save('cif_images/vis{:03d}.png'.format(i), figure=figure);