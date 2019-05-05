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
        
        # INPUT DATA SHOULD BE AN ARRAY, SO WE CONVERT INTO IT FROM A LIST
        self.starting_points = np.asarray(self.init_sandpile(width=self.width, height=self.height))
        self.current_state = self.starting_points

        # CREATE A FULL NULL-MATRIX IN THE SIZE OF THE FULL ARRAY
        self.Background = np.zeros_like(self.starting_points)

        # Set up a borderless frame
        self.ax = self.game_arena.add_axes([0, 0, 1, 1], xticks=[], yticks=[], frameon=False)

        # FROM STACKOVERFLOW, PLOTTING BINARY MATRIX WITH WHITE AND BLACK SQUARES
        # ON STACKOVERFLOW IT WAS AN EXAMPLE FOR ANIMATING A FRACTAL
        self.Image = self.ax.imshow(self.starting_points, cmap='Greys', interpolation='nearest')

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

    # WE SHOULD INITIALIZE THE PLOTTING ANIMATION WITH A BACKGROUND FIRST IN MATPLOTLIB'S FuncAnimation
    # THIS IS NECESSARY, BECAUSE THE ANIMATION WOULD ONLY SHOW ONLY ONE FRAME, AND WE DON'T WANT IT
    def init_anim(self):
        self.Image.set_data(self.Background)
        return self.Image,

    # Stepping to the next data block
    def generation_update(self, x):
        self.Image.set_data(self.starting_points[1:self.width+1,1:self.height+1])
        self.starting_points = np.asarray(self.step_sandpile(sandpile=self.starting_points, width=self.width, height=self.height))

        # Update current state, if there are changes in the sandpile, else stop
        for i in range(1, self.height+1):
            for j in range(1, self.width+1):
                if(self.starting_points[i][j] != self.current_state[i][j]):
                    self.test = True
                    break
            if(self.test == True):
                break

        if(self.test == True):
            self.current_state = self.starting_points

        return self.Image,

    def gen(self):
        i = 0
        while(self.test == True):
            i += 1
            yield i

    def animate(self):

        # Set up movie formatter with ffmpeg
        self.Writer = animation.writers['ffmpeg']
        self.writer = self.Writer(fps=30, metadata=dict(artist='PB'), bitrate=1800)

        # Animate the function
        #self.starting_points = self.starting_points
        self.anim = animation.FuncAnimation(fig=self.game_arena, func=self.generation_update, frames=self.gen, interval=50, init_func=self.init_anim, blit=True)
        self.anim.save('.\\videos\\' + self.video_name + '.mp4', writer=self.writer)

# MAIN
sandpileClass = SandpileAnimation(video_name='sandpile_fullpy_anim_cl', width=int(sys.argv[1]), height=int(sys.argv[2]))
sandpileClass.animate()