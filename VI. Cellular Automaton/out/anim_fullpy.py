import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import sys

class SandpileAnimation():

    def __init__(self, video_name, width, height, **kw):

        sliced = 4
        self.width = width
        self.height = height
        self.scale_factor = 7
        self.game_arena = plt.figure(figsize=(self.width/sliced, self.height/sliced))

        # Initially
        self.unstable = True
        
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

        # Set up movie formatter with ffmpeg
        self.Writer = animation.writers['ffmpeg']
        self.writer = self.Writer(fps=30, metadata=dict(artist='PB'), bitrate=1800)

        # Animate the function
        #self.starting_points = self.starting_points
        self.anim = animation.FuncAnimation(fig=self.game_arena, func=self.generation_update, frames=self.gen_frame(), interval=50, init_func=self.init_anim, blit=True)
        self.anim.save('.\\videos\\' + video_name + '.mp4', writer=self.writer)

    # Generating frames
    def gen_frame(self):
        i = 0
        while(self.unstable):
            i += 1
            yield i

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

        # Check if the stable state is reached (every pile is smaller than 4 grains)
        self.unstable = False
        for i in range(1, self.height+1):
            for j in range(1, self.width+1):
                # If there are unstable pile, then stop the checking
                if(self.starting_points[i][j] > 3):
                    self.unstable = True
                    break
            if(self.unstable == True):
                break

        # If there are any unstable pile found, continue to the next step
        if(self.unstable == True):
            pass

        # If all piles are stable, then stop and save the video.
        else:
            print('stopped at {0}'.format(x))
            #self.anim.event_source.stop()

        return self.Image,

# MAIN
sandpileClass = SandpileAnimation(video_name='sandpile_fullpy_anim_cl', width=int(sys.argv[1]), height=int(sys.argv[2]))