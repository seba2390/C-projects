import matplotlib.pyplot as plt
from matplotlib import animation, rc
rc('animation', html='jshtml')
#import statistics
from argparse import ArgumentParser
from copy import deepcopy
import numpy as np
plt.rc("font", family=["Helvetica", "Arial"])
plt.rc("text", usetex=True)



def load_data(nr_particles, nr_steps):
    if not type(nr_particles) is int:
        nr_particles = int(nr_particles)
    if not type(nr_steps) is int:
        nr_steps = int(nr_steps)
    positions  = np.loadtxt('coordinates.txt')
    velocities = np.loadtxt('velocities.txt')
    radia      = np.loadtxt('attributes.txt')[:,1]
    partitioned_positions  = [] 
    partitioned_velocities = []
    for i in range(nr_steps):
       partitioned_positions.append(positions[i*nr_particles:(i+1)*nr_particles])
       partitioned_velocities.append(velocities[i*nr_particles:(i+1)*nr_particles])
    return partitioned_positions,partitioned_velocities, radia


def animate_box(positions, velocities, radia, box_width, save_animation=False):
    fig, ax = plt.subplots(1,1,figsize=(12,12))
    if box_width == 10:
        my_lw = 3 
    if box_width == 20:
        my_lw = 0.3 
    if box_width == 40:
        my_lw = 0.1

    ax.vlines(0,0,box_width,color='black',lw=my_lw), ax.vlines(box_width,0,box_width,color='black',lw=my_lw)
    ax.hlines(0,0,box_width,color='black',lw=my_lw), ax.hlines(box_width,0,box_width,color='black',lw=my_lw)
    ax.set_xlabel("X"),ax.set_ylabel("Y")
    ax.axis('off')
    if np.all(radia == radia[0]):
        ax.set_title("Constant mass and radius of "+str(len(radia))+" particles \n (constant density)",size= 25)
    else:
        ax.set_title("Variating mass and radius of "+str(len(radia))+" particles \n (constant density)",size= 25)
    colors = deepcopy(radia) # Set this as c in ax.scatter for different colors depending on size
    if box_width == 10:
        radia = 4*radia**2*4e3
    if box_width == 20:
        radia = 4*radia**2.85*4e3
    if box_width == 40:
        radia = 4*radia**3.6*4e3

    scatter = ax.scatter(positions[0][: , 0], positions[0][: , 1], s = radia, c = colors,edgecolors='black')

    def update(frame_number):
        scatter.set_offsets(positions[frame_number])
        return scatter,

    anim = animation.FuncAnimation(fig,
                                   update,
                                   frames=kwargs['nr_steps'],
                                   interval=1,
                                   blit=True,
                                   repeat=False)
    if save_animation:
        writervideo = animation.FFMpegWriter(fps=60) 
        anim.save("ParticlesInBox.mp4", writer=writervideo)    
    else:
        plt.show()

def plot_histogram(velocities,save_plot=False):
    fig, ax = plt.subplots(1,2,figsize=(14,5))
    ax[0].set_title("Initial speed distribution ("+str(len(velocities[0]))+" particles)",size=25)
    ax[1].set_title("Final speed distribution ("+str(len(velocities[0]))+" particles)",size=25)
    ax[0].set_xlabel(r"$||\vec{v}||$",size=20),ax[0].set_ylabel("nr. of particles",size=20)
    ax[1].set_xlabel(r"$||\vec{v}||$",size=20),ax[1].set_ylabel("nr. of particles",size=20)
    speeds = []
    for frame in range(len(velocities)):
        temp = []
        for particle in range(len(velocities[frame])):
            speed = np.sqrt(velocities[frame][particle][0]**2+velocities[frame][particle][1]**2)
            temp.append(speed)
        speeds.append(temp)
    #speeds = np.array(speeds)
    ax[0].hist(speeds[0],edgecolor='k',bins=20)
    ax[1].hist(speeds[-1],edgecolor='k',bins=20)
    if save_plot:
        plt.savefig("Speed_histogram.jpg",dpi=150)
    else:
        plt.show()


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-particles',type=int, dest='nr_particles')
    parser.add_argument('-steps',    type=int, dest='nr_steps')
    parser.add_argument('-bw',       type=float, dest='box_width')

    args = parser.parse_args()
    assert args.nr_particles > 0, "Trying to plot zero particles\n"
    assert args.nr_steps >= 0,    "Trying to take a negative number of time steps\n"
    print('args', args)
    kwargs = {}
    if args.nr_particles > 0:
        kwargs['nr_particles'] = args.nr_particles
    if args.nr_steps >= 0:
        kwargs['nr_steps'] = args.nr_steps
    if args.box_width > 0:
        kwargs['box_width'] = args.box_width

    positions, velocities, radia = load_data(kwargs['nr_particles'],kwargs['nr_steps'])
    
    animate_box(positions, velocities, radia, kwargs['box_width'], save_animation=True)
    plot_histogram(velocities,save_plot=True)
