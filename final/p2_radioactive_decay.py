# PROBLEM 2
#
# This problem models the decay chain of uranium-238.
# Use the Backward Euler Method to show how the
# populations of bismuth, polonium, and lead atoms
# change with time.
#

import math
from udacityplots import *

h = 10. # days
lifetime_bi = 5. / math.log(2.) # days
lifetime_po = 138. / math.log(2) # days

end_time = 5.0 * 365. # days
num_steps = int(end_time / h)
times = h * numpy.array(range(num_steps + 1))

def radioactive_decay():
    # numbers of atoms
    bi = numpy.zeros(num_steps + 1)
    po = numpy.zeros(num_steps + 1)
    pb = numpy.zeros(num_steps + 1)

    bi[0] = 1e24
    po[0] = 0.
    pb[0] = 0.


    for step in range(num_steps):
        # Task: Insert the Backward (!) Euler computation for the radioactive decay.
        # For exponential decay (or growth) we have x'(t) = c * x(t)
        # In particular, if x decays with average lifetime of lifetime:
        #   x(t) = A * exp(-t / lifetime) [where A is starting value]
        #   x'(t) = -1/lifetime * x(t)
        # We simulate this with step size h as:
        #   x(step+1) = x(step) + h * x'(step) [Forward Euler]
        #   x(step+1) = x(step) + h * x'(step+1) [Backward Euler]
        bi[step + 1] = bi[step] / (1.+h/lifetime_bi)
        po[step + 1] = (po[step] + h * (1./lifetime_bi * bi[step+1]))/(1.+h / lifetime_po)
        pb[step + 1] = pb[step] + h * (1./lifetime_po * po[step+1])

    return bi, po, pb

bi, po, pb = radioactive_decay()

@show_plot
def plot_decay():
    bi_plot = matplotlib.pyplot.plot(times, bi, label='$^{210}$Bi')
    po_plot = matplotlib.pyplot.plot(times, po, label='$^{210}$Po')
    pb_plot = matplotlib.pyplot.plot(times, pb, label='$^{206}$Pb')
    matplotlib.pyplot.legend(('$^{210}$Bi', '$^{210}$Po', '$^{206}$Pb'), loc='upper right')
    axes = matplotlib.pyplot.gca()
    axes.set_xlabel('Time in days')
    axes.set_ylabel('Number of atoms')
    matplotlib.pyplot.xlim(xmin = 0.)
    matplotlib.pyplot.ylim(ymin = 0.)

plot_decay()
