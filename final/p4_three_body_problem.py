# PROBLEM 4
#
# We have three stars in a three-dimensional space.  Use the symplectic
# Euler method to show how the position and velocity of these stars change
# with time.
#

from udacityplots import *
from mpl_toolkits.mplot3d import Axes3D
import itertools
import numpy

star_1_mass = 1e30 # kg
star_2_mass = 2e30 # kg
star_3_mass = 3e30 # kg
gravitational_constant = 6.67e-11 # m3 / kg s2
masses = numpy.array((star_1_mass,star_2_mass,star_3_mass))

end_time = 10. * 365.26 * 24. * 3600. # s
h = 2. * 3600. # s
num_steps = int(end_time / h)
times = h * numpy.array(range(num_steps + 1))

epsilon = 1e-6


def accelerations(positions, masses):
    '''Params:
    - positions: numpy array of size (n,3)
    - masses: numpy array of size (n,)
    '''
    n_bodies = len(masses)
    accelerations = numpy.zeros([n_bodies,3])  # n_bodies * (x,y,z)

    # vectors from mass(i) to mass(j)
    D = numpy.zeros([n_bodies, n_bodies, 3])  # n_bodies * n_bodies * (x,y,z)
    for i, j in itertools.product(range(n_bodies), range(n_bodies)):
        D[i][j] = positions[j]-positions[i]

    # Acceleration due to gravitational force between each pair of bodies
    A = numpy.zeros((n_bodies, n_bodies, 3))
    for i, j in itertools.product(range(n_bodies), range(n_bodies)):
        if numpy.linalg.norm(D[i][j]) > epsilon:
            A[i][j] = gravitational_constant * masses[j] * D[i][j] \
            / numpy.linalg.norm(D[i][j])**3

    # Calculate net accleration of each body
    accelerations = numpy.sum(A, axis=1) #  sum of accel vectors for each body

    return accelerations

def three_body_problem():

    # three indices: which time, which star, xyz
    positions = numpy.zeros([num_steps + 1, 3, 3]) # m
    velocities = numpy.zeros([num_steps + 1, 3, 3]) # m / s
    positions[0] = numpy.array([[1., 3., 2.], [6., -5., 4.], [7., 8., -7.]]) * 1e11
    velocities[0] = numpy.array([[-2., 0.5, 5.], [7., 0.5, 2.], [-4., -0.5, -3.]]) * 1e3

    for step in range(num_steps):
        # Task: Implement the symplectic Euler Method for the motion
        # of three stars with the data provided above.

        # your code here
        acc = accelerations(positions[step],masses)
        velocities[step+1] = velocities[step] + h * acc
        positions[step+1] = positions[step] + h * velocities[step+1]

    return positions, velocities

positions, velocities = three_body_problem()

@show_plot
def plot_stars():
    axes = matplotlib.pyplot.gca()
    axes.set_xlabel('x in m')
    axes.set_ylabel('z in m')
    axes.plot(positions[:, 0, 0], positions[:, 0, 2])
    axes.plot(positions[:, 1, 0], positions[:, 1, 2])
    axes.plot(positions[:, 2, 0], positions[:, 2, 2])
    matplotlib.pyplot.axis('equal')

plot_stars()
