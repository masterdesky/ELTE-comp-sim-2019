import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class SandpileAnimation():

    def __init__(self, video_name, width, height, **kw):

        self.video_name = video_name

        sliced = 4
        self.width = width
        self.height = height
        self.scale_factor = 7
        self.game_arena = plt.figure(figsize=(self.width/sliced, self.height/sliced))

        self.test = True
        
        # Input data should be an array, so make sure, init_sandpile returns with that
        # Also save the current configuration
        self.starting_points = np.asarray(self.init_sandpile(width=self.width, height=self.height))
        self.current_state = self.starting_points

        # Create a full-zero background for the initialization
        self.Background = np.zeros_like(self.starting_points)

        # Set up a borderless frame for the sandpile
        self.ax = self.game_arena.add_axes([0, 0, 1, 1], xticks=[], yticks=[], frameon=False)

        # Defining the image frame, shown in the simulation
        self.Image = self.ax.imshow(self.starting_points, cmap='Greys', interpolation='nearest')


    #
    # Actual Abelian-sandpile mechanics
    #
    def init_sandpile(self, width, height):
    
        sandpile = np.zeros((height+2, width+2))
        
        for i in range(1, height+1):
            for j in range(1, width+1):
                sandpile[i][j] = self.scale_factor
        
        return(sandpile)

    def step_sandpile(self, sandpile, width, height):

        new_sandpile = sandpile

        for i in range(1, height+1):
            for j in range(1, width+1):
                if(sandpile[i][j] > 3):
                    new_sandpile[i][j] -= 4

                    # Add them to neighbours
                    new_sandpile[i+1][j] += 1
                    new_sandpile[i-1][j] += 1
                    new_sandpile[i][j+1] += 1
                    new_sandpile[i][j-1] += 1

        sandpile = new_sandpile
        
        return(sandpile)

    # Initializing the frame for the FuncAnimations
    def init_anim(self):
        self.Image.set_data(self.Background)
        return self.Image,

    # Stepping to the next data block
    def generation_update(self, x):
        self.Image.set_data(self.starting_points[1:self.width+1,1:self.height+1])
        self.starting_points = np.asarray(self.step_sandpile(sandpile=self.starting_points, width=self.width, height=self.height))

        #
        # Check if there are changes in the sandpile compared to the previous
        # state. If there are, then continue with the process, else stop.
        #
        # First set the checking to False, initially consider a static position.
        self.test = False
        for i in range(1, self.height+1):
            for j in range(1, self.width+1):
                # If a deviation is found, se the checking state to True to
                # indicate, there are deviations. In this case, jump out from
                # both of the loops.
                if(self.starting_points[i][j] != self.current_state[i][j]):
                    self.test = True
                    break
            if(self.test == True):
                break

        # If there are any deviation, continue to the next step, and save
        # the current configuration for the check in the next step.
        if(self.test == True):
            self.current_state = self.starting_points

        # If there aren't any deviation from the previous state, stop and save the video.
        else:
            self.anim.event_source.stop()

        return self.Image,

    def animate(self):

        # Set up movie formatter with ffmpeg
        self.Writer = animation.writers['ffmpeg']
        self.writer = self.Writer(fps=30, metadata=dict(artist='PB'), bitrate=1800)

        # Animate the function
        #self.starting_points = self.starting_points
        self.anim = animation.FuncAnimation(fig=self.game_arena, func=self.generation_update, frames=None, interval=50, init_func=self.init_anim, blit=True)
        self.anim.save('.\\videos\\' + self.video_name + '.mp4', writer=self.writer)

# MAIN
sandpileClass = SandpileAnimation(video_name='sandpile_fullpy_anim_cl', width=int(sys.argv[1]), height=int(sys.argv[2]))
sandpileClass.animate()