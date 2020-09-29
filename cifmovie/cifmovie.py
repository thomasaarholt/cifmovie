from pathlib import Path # for handling files

import math
import numpy as np
import matplotlib.pyplot as plt

from ase.io import read # to read in crystallographic files
from ase.data import covalent_radii
from ase.data.colors import jmol_colors

import plato.draw.matplotlib

from scipy.spatial.transform import Rotation

import ffmpeg

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


def rotation_calc(axes='z', angles=(0,)):
    '''
    Multidimensional rotations
    
    Either:
    axes = 'z', angles = (90,)
    Or:
    axes = 'z', angles = (0,1,2,3...) of same length as number of frames
    Or:
    axes = 'zyx' angles = ([0,1,2, ...], [3,2,1, ...], [0,1,2,...]) # (zangles, yangles, xangles) or any combination

    See scipy.spatial.transform.Rotation for more info
    '''
    return Rotation.from_euler(axes, angles, degrees=True).as_quat()


def render(atoms_list, name='img', dpi=72, resolution=(1024, 1024), rotation=[1, 0, 0, 0], background_color='white', format='png'):

    renderer = plato.draw.matplotlib
    
    rotation = np.array(rotation)
    if np.ndim(rotation) == 1:
        rotation = np.stack(len(atoms_list)*(rotation,))

    p = Path("cif_images")
    p.mkdir(exist_ok=True)
    
    figsize = np.array(resolution) / dpi
    figure = plt.figure(figsize=figsize, dpi=dpi)
    
    for i, atom in enumerate(atoms_list):
        figure.clf()
        pos = atom.positions
        centre = (pos.min(0) + pos.max(0))/2
        pos -= centre
        numbers = atom.numbers
        radii = covalent_radii[numbers]
        colors = jmol_colors[numbers]
        colors = np.hstack((colors, np.ones((len(colors), 1)))) # format colors as matplotlib expects

        scene = renderer.Scene(
            zoom=1, 
            antialiasing=True, 
            rotation=rotation[i],
        )
        sphere = renderer.Spheres(positions=pos, radii=radii/2, colors=colors, )#pixel_scale=20,)
        scene.add_primitive(sphere)
        
        if background_color:
            ax = plt.gca()
            ax.set_facecolor(background_color)

        scene.save('cif_images/{}_{:03d}.{}'.format(name, i, format), figure=figure);

def movie(path='cif_images', name='img', format='png', duration=3, fps=None):
    if not fps:
        nImages = len(list(Path('cif_images').glob('{}*.{}'.format(name, format))))
        fps = nImages / duration

    ff_input = ffmpeg.input('cif_images/{}_%03d.{}'.format(name, format), framerate=fps)
    ff_output = ff_input.output('{}.mp4'.format(name), vcodec='copy')
    ff_output.run(overwrite_output=True)



    #ffmpeg -framerate 30 -i img_%03d.png -codec copy output.mkv