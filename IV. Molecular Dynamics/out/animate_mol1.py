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
image_dpi = 150
image_format = 'pdf'

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
        axes.set_xlabel('Distance (x) [$m$]', fontsize=axislabelsize, labelpad=labelpad)
        axes.set_ylabel('Distance (y) [$m$]', fontsize=axislabelsize, labelpad=labelpad)
        axes.set_zlabel('Distance (z) [$m$]', fontsize=axislabelsize, labelpad=labelpad)

        #axes[1].set_xlabel('Velocity (x) [$\\frac{m}{s}$]', fontsize=axislabelsize, labelpad=labelpad)
        #axes[1].set_ylabel('Velocity (y) [$\\frac{m}{s}$]', fontsize=axislabelsize, labelpad=labelpad)
        #axes[1].set_zlabel('Velocity (z) [$\\frac{m}{s}$]', fontsize=axislabelsize, labelpad=labelpad)

        #axes[2].set_xlabel('Acceleration (x) [$\\frac{m}{s^{2}}$]', fontsize=axislabelsize, labelpad=labelpad)
        #axes[2].set_ylabel('Acceleration (y) [$\\frac{m}{s^{2}}$]', fontsize=axislabelsize, labelpad=labelpad)
        #axes[2].set_zlabel('Acceleration (z) [$\\frac{m}{s^{2}}$]', fontsize=axislabelsize, labelpad=labelpad)

        axes.tick_params(axis='both', which='major', labelsize=axislabelsize)
        #axes[1].tick_params(axis='both', which='major', labelsize=axislabelsize)
        #axes[2].tick_params(axis='both', which='major', labelsize=axislabelsize)
        
        # Centre the image on the fixed anchor point, and ensure the axes are equal
        axes.set_xlim(0,10)
        axes.set_ylim(0,10)
        axes.set_zlim(0,10)

        #axes[1].set_xlim(-3,3)
        #axes[1].set_ylim(-3,3)
        #axes[1].set_zlim(-3,3)

        #axes[2].set_xlim(-60,60)
        #axes[2].set_ylim(-60,60)
        #axes[2].set_zlim(-60,60)
        
        axes.set_aspect('equal', adjustable='box')
        #axes[1].set_aspect('equal', adjustable='box')
        #axes[2].set_aspect('equal', adjustable='box')

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

            for k in range(0, (data_set_1.shape[1]-1)//9):
                art3d.pathpatch_2d_to_3d(Circle((data_set_1[::steps,k*9][i],
                                                 data_set_1[::steps,k*9+1][i]), 0.1,
                                                 fc=colors[k], ec=colors[k], zorder=10),
                                                 z=data_set_1[::steps,k*9+2][i], zdir=i)
                
                axes.plot(data_set_1[::steps,k*9][imin:imax],
                          data_set_1[::steps,k*9+1][imin:imax],
                          data_set_1[::steps,k*9+2][imin:imax], lw=5, color=colors[k], alpha=alpha)
                #axes[1].plot(data_set_1[::steps,i*9+3][imin:imax], data_set_1[::steps,i*9+4][imin:imax], data_set_1[::steps,i*9+5][imin:imax], lw=3, alpha=alpha)
                #axes[2].plot(data_set_1[::steps,i*9+6][imin:imax], data_set_1[::steps,i*9+7][imin:imax], data_set_1[::steps,i*9+8][imin:imax], lw=3, alpha=alpha)


        fig.tight_layout()

        # Don't show axes, only white background
        #ax.axis('off')

        plt.savefig(path + '_img{0:4d}.png'.format(i), dpi=72)    # Save next frame as png
        image = imageio.imread(path + '_img{0:4d}.png'.format(i)) # Load saved image
        writer.append_data(image)                                 # Append this image as the next frame to video

        # Clear the pyplot background for the next frame
        plt.cla()

        # Delete the now useless image from frames' folder
        os.unlink(path + '_img{0:4d}.png'.format(i))

    with imageio.get_writer(video_title + '.mp4', fps=fps) as writer:
        for i in range(0, data_set_1.shape[0]):
            animation(i)
            sys.stdout.write("Progress: [{0}{1}] {2:.3f}%\tElapsed: {3}\tRemaining: {4}\r".format(u'\u2588' * int((i+1)/data_set_1.shape[0] * progressbar_width),
                                                                                                  u'\u2591' * (progressbar_width - int((i+1)/data_set_1.shape[0] * progressbar_width)),
                                                                                                  (i+1)/data_set_1.shape[0] * 100,
                                                                                                  datetime.now()-starttime,
                                                                                                  (datetime.now()-starttime)/(i+1) * (data_set_1.shape[0] - (i+1))))
            sys.stdout.flush()


# MAIN
data_set_1 = np.genfromtxt('md1.dat')

colors = np.empty(((data_set_1.shape[1]-1)//3,3))

for i in range(0,len(colors)):
    colors[i] = np.array([random.random(), random.random(), random.random()])

ANIMATE_VIDEO(path = '.\\frames\\', video_path='.\\videos\\', gen_video_title='threebody', packages=20)