import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import os
import sys
from datetime import datetime
import imageio

sns.set_style(style='whitegrid')
# Others
steps = 1
fps = 40


def ANIMATE_VIDEO(path, video_path, gen_video_title, packages):
    
    video_title = video_path + gen_video_title
    
    nrows=1
    ncols=1
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*20,nrows*20))

    # Measuring time for progress bar
    starttime = datetime.now()
    # Width of progress bar, shown on sys.stdout
    progressbar_width = 20
    
    ##ANIMATION STUFF BEGINS HERE##
    # Plot and save an image of the twobody system for time point i
    def animation(i):

        # Centre the image on the fixed anchor point, and ensure the axes are equal
        axes.set_xlim(0,width)
        axes.set_ylim(0,height)
        axes.set_aspect('equal', adjustable='box')

        # Show current image
        axes.imshow(data_cell_sp[height*i:height*(i+1)], cmap='hot')

        fig.tight_layout()

        # Don't show axes, only white background
        axes.axis('off')

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
            sys.stdout.write("Progress: [{0}{1}] {2:.3f}%\tElapsed: {3}\tRemaining: {4}\r".format(u'\u2588' * int((i - MIN_LIM + 1)/(MAX_LIM-MIN_LIM) * progressbar_width),
                                                                                                  u'\u2591' * (progressbar_width - int((i - MIN_LIM + 1)/(MAX_LIM-MIN_LIM) * progressbar_width)),
                                                                                                  (i - MIN_LIM + 1)/(MAX_LIM-MIN_LIM) * 100,
                                                                                                  datetime.now()-starttime,
                                                                                                  (datetime.now()-starttime)/(i - MIN_LIM + 1) * ((MAX_LIM-MIN_LIM) - (i - MIN_LIM + 1))))
            sys.stdout.flush()

# =============== MAIN ===============
# MAIN #
data_cell_sp = np.genfromtxt('sandpile.dat')

width = int(sys.argv[2])
height = int(sys.argv[3])

# Limits of animated parts
MIN_LIM = int(sys.argv[4])
MAX_LIM = int(sys.argv[5])

if MIN_LIM < 0:
    MIN_LIM = 0
if MAX_LIM > data_cell_sp.shape[0]/height - 1:
    MAX_LIM = data_cell_sp.shape[0]//height - 1
if MIN_LIM > MAX_LIM:
    MIN_LIM_TEMP = MIN_LIM
    MIN_LIM = MAX_LIM
    MAX_LIM = MIN_LIM_TEMP

ANIMATE_VIDEO(path = '.\\frames\\', video_path='.\\videos\\', gen_video_title=sys.argv[1], packages=20)