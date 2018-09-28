# PROBLEM 6
#
# We will use an activator-inhibitor model to show how an initially only
# slightly random distribution of chemicals in a square grid evolves into
# a stable pattern.  The inhibitor chemical causes a reduction of both the
# activator and the inhibitor, and it diffuses quickly.  The activator,
# however, causes production of both the activator and the inhibitor, and it
# diffuses slowly.  First, insert periodic boundary conditions for the grid
# in both dimensions.  Then apply the explicit finite-difference scheme to
# the activator-inhibitor model.  Although the video shows this in 1D, please
# write your solution in 2D.
#

from udacityplots import *

diffusion_coefficient_a = 0.0005 # m2 / s
diffusion_coefficient_b = 0.004 # m2 / s

length = 1. # meters; domain extends from 0 to length
size = 20 # number of points per dimension of the grid
dx = length / size

# Pick a time step below the threshold of instability
h = 0.2 * dx ** 2 / max(diffusion_coefficient_a, diffusion_coefficient_b) # s
end_time = 15. # s

@show_plot
def pattern():
    numpy.random.seed(42)
    a_old = numpy.random.normal(0., 0.03, [size, size])
    b_old = numpy.random.normal(0., 0.03, [size, size])
    a_new = a_old.copy()
    b_new = b_old.copy()

    num_steps = int(end_time / h)
    for step in range(num_steps):
        for j in range(size):
            j_plus = j + 1
            j_minus = j - 1
            # Task 1: Treat j_plus and j_minus so that they are the coordinates of the next cell in the corresponding
            # direction -- taking a possible wrap-around into account (periodic boundary conditions).

            # your code here
            j_plus %= size
            j_minus = (j_minus+size)%size

            for i in range(size):
                i_plus = i + 1
                i_minus = i - 1
                # Task 2: Similar to task 1, but now for i_plus and i_minus

                # your code here
                i_plus %= size
                i_minus = (i_minus+size)%size

                a0 = a_old[j, i]
                b0 = b_old[j, i]

                # Task 3: Implement the explicit finite-difference scheme for the activator-inhibitor model.
                # da/dt = diffusion_coefficient_a * d^2a/dx^2 + (1-a**2) * a - b
                # db/dt = diffusion_coefficient_b * d^2b/dx^2 + a - 0.3 * b


                da_dt = diffusion_coefficient_a / dx**2 \
                    * (a_old[j_minus, i] + a_old[j_plus, i]
                        + a_old[j, i_minus] + a_old[j, i_plus]
                        - 4. * a_old[j, i]) \
                    + (1. - a_old[j, i]**2) * a_old[j, i] - b_old[j, i]

                db_dt = diffusion_coefficient_b / dx**2 \
                    * (b_old[j_minus, i] + b_old[j_plus, i]
                        + b_old[j, i_minus] + b_old[j, i_plus]
                        - 4. * b_old[j, i]) \
                    + a_old[j, i] - 0.7 * b_old[j, i]

                a_new[j, i] = a_old[j, i] + h * da_dt
                b_new[j, i] = b_old[j, i] + h * db_dt

        a_old, b_old, a_new, b_new = a_new, b_new, a_old, b_old

    dimensions = [0, length, 0, length]

    axes = matplotlib.pyplot.subplot(121)
    matplotlib.pyplot.imshow(a_old, cmap = matplotlib.cm.winter, origin = 'lower', extent = dimensions)
    matplotlib.pyplot.colorbar()
    axes.set_title('a')
    axes.set_xlabel('x in m')
    axes.set_ylabel('y in m')

    axes = matplotlib.pyplot.subplot(122)
    matplotlib.pyplot.imshow(b_old, cmap = matplotlib.cm.winter, origin = 'lower', extent = dimensions)
    matplotlib.pyplot.colorbar()
    axes.set_title('b')
    axes.set_xlabel('x in m')
    axes.set_ylabel('y in m')

    return a_old, b_old

pattern()
