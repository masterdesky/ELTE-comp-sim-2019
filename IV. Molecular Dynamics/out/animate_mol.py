import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import os
import sys
from datetime import datetime
import imageio
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import mpl_toolkits.mplot3d.art3d as art3d

sns.set_style(style='whitegrid')
# Others
steps = 1
fps = 40


def ANIMATE_VIDEO(path, video_path, gen_video_title, packages):
    
    video_title = video_path + gen_video_title
    
    nrows=1
    ncols=1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*20,nrows*20), subplot_kw={'projection': '3d'})

    axislabelsize = 30
    labelpad = 20

    # Plot a trail of a planet position for the last trail_secs seconds.
    trail_secs = 5
    # This corresponds to max_trail time points.
    max_trail = int(trail_secs * fps)

    # Measuring time for progress bar
    starttime = datetime.now()
    # Width of progress bar, shown on sys.stdout
    progressbar_width = 20
    
    ##ANIMATION STUFF BEGINS HERE##
    # Plot and save an image of the twobody system for time point i
    def animation(i):

        # Axis labels
        axes.set_xlabel('Distance (x)', fontsize=axislabelsize, labelpad=labelpad)
        axes.set_ylabel('Distance (y)', fontsize=axislabelsize, labelpad=labelpad)
        axes.set_zlabel('Distance (z)', fontsize=axislabelsize, labelpad=labelpad)

        axes.tick_params(axis='both', which='major', labelsize=axislabelsize)
        
        # Centre the image on the fixed anchor point, and ensure the axes are equal
        axes.set_xlim(0,L)
        axes.set_ylim(0,L)
        axes.set_zlim(0,L)
        axes.set_aspect('equal', adjustable='box')

        # TRAILING EFFECT FROM: https://scipython.com/blog/the-double-pendulum/
        # The trail will be divided into ns segments and plotted as a fading line.
        ns = 20
        s = max_trail // ns

        for j in range(ns):
            imin = i - (ns-j)*s
            if imin < 0:
                imin = 0
            imax = imin + s + 1
            # The fading looks better if we square the fractional length along the trail
            alpha = (j/ns)**2

            for k in range(0, (data_set.shape[1]-1)//9):
                # Molecule trajectories
                if sys.argv[1] == '3':
                    axes.scatter(data_set[::steps,k*9][imin:imax],
                                 data_set[::steps,k*9+1][imin:imax],
                                 data_set[::steps,k*9+2][imin:imax], color=colors[k], alpha=alpha, s=10)

                else:
                    axes.scatter(data_set[::steps,k*9][imin:imax],
                                 data_set[::steps,k*9+1][imin:imax],
                                 data_set[::steps,k*9+2][imin:imax], color=colors[k], alpha=alpha, s=20)

                # Molecules itself
                axes.scatter(data_set[::steps,k*9][i],
                             data_set[::steps,k*9+1][i],
                             data_set[::steps,k*9+2][i], color='red', s=200)


        fig.tight_layout()

        #if sys.argv[1] == '1':
            # Don't show axes, only white background
            #axes.axis('off')

        plt.savefig(path + '_img{0:4d}.png'.format(i), dpi=80)    # Save next frame as png
        image = imageio.imread(path + '_img{0:4d}.png'.format(i)) # Load saved image
        writer.append_data(image)                                 # Append this image as the next frame to video

        # Clear the pyplot background for the next frame
        plt.cla()

        # Delete the now useless image from frames' folder
        os.unlink(path + '_img{0:4d}.png'.format(i))

    with imageio.get_writer(video_title + '.mp4', fps=fps) as writer:
        for i in range(MIN_LIM, MAX_LIM):
            animation(i)
            sys.stdout.write("Progress: [{0}{1}] {2:.3f}%\tElapsed: {3}\tRemaining: {4}\r".format(u'\u2588' * int((i - MIN_LIM +1)/(MAX_LIM-MIN_LIM) * progressbar_width),
                                                                                                  u'\u2591' * (progressbar_width - int((i - MIN_LIM +1)/(MAX_LIM-MIN_LIM) * progressbar_width)),
                                                                                                  (i - MIN_LIM +1)/(MAX_LIM-MIN_LIM) * 100,
                                                                                                  datetime.now()-starttime,
                                                                                                  (datetime.now()-starttime)/(i - MIN_LIM +1) * ((MAX_LIM-MIN_LIM) - (i - MIN_LIM +1))))
            sys.stdout.flush()


# =============== MAIN ===============
if sys.argv[1] == '1':
    data_set = np.genfromtxt('md1.dat')
if sys.argv[1] == '2':
    data_set = np.genfromtxt('md2.dat')
if sys.argv[1] == '3':
    data_set = np.genfromtxt('md3.dat')

# Some simulation parameters
if sys.argv[1] == '1':
    L = 10

if sys.argv[1] == '2':
    N=64
    rho=1.05
    L = pow(N / rho, 1.0/3)

if sys.argv[1] == '3':
    N=64
    rho=0.95
    L = pow(N / rho, 1.0/3)

# Limits of animated parts
MIN_LIM = int(sys.argv[2])
MAX_LIM = int(sys.argv[3])

if MIN_LIM < 0:
    MIN_LIM = 0
if MAX_LIM > data_set.shape[0]:
    MAX_LIM = data_set.shape[0]
if MIN_LIM > MAX_LIM:
    MIN_LIM_TEMP = MIN_LIM
    MIN_LIM = MAX_LIM
    MAX_LIM = MIN_LIM_TEMP


# Generate new colors
if sys.argv[5] == 'newcolors':
    # Define randomized colormaps for particles
    colors = np.empty(((data_set.shape[1]-1)//3,3))
    for i in range(0,len(colors)):
        colors[i] = np.array([random.random(), random.random(), random.random()])

    np.savetxt('colors.dat', colors, delimiter=' ')

else:
    colors = np.genfromtxt('colors.dat')

ANIMATE_VIDEO(path = '.\\frames\\', video_path='.\\videos\\', gen_video_title=sys.argv[4], packages=20)