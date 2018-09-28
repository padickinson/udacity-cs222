# PROBLEM 1
#
# Use the Forward Euler Method to create an SEIR model that also
# takes into account the birth rate and death rate listed below.
# The total population is thus not constant.  People from all
# compartments can reproduce, children are born susceptible, and
# all compartments have the same mortality rate.
#

from udacityplots import *

h = 0.5 # days
contacts_per_day = 8. # 1 / day
transmission_probability = 0.8 # dimensionless
latency_time = 7. # days
infectious_time = 7. # days
birth_rate = 0.023 / 365. # 1 / days
mortality_rate = 0.013 / 365. # 1 / days

end_time = 5. * 365. # days
num_steps = int(end_time / h)
times = h * numpy.array(range(num_steps + 1))

def seir_model_births_deaths():

    s = numpy.zeros(num_steps + 1)
    e = numpy.zeros(num_steps + 1)
    i = numpy.zeros(num_steps + 1)
    r = numpy.zeros(num_steps + 1)

    s[0] = 1e8 - 1e6 - 1e5
    e[0] = 0.
    i[0] = 1e5
    r[0] = 1e6

    for step in range(num_steps):
        # Task 1: Compute s2e, the number of infections per day, taking into account contacts_per_day, transmission_probability, and the current total number of persons.
        pop = s[step]+e[step]+i[step]+r[step]
        s2e = h * contacts_per_day * s[step] / pop * transmission_probability * i[step]
        e2i = h / latency_time * e[step]
        i2r = h / infectious_time * i[step]
        # Task 2: Include births and deaths in the following four lines.
        # Modify the code below:
        births = h * birth_rate * pop
        s[step + 1] = s[step] - s2e + births - h * mortality_rate * s[step]
        e[step + 1] = e[step] + s2e - e2i - h * mortality_rate * e[step]
        i[step + 1] = i[step] + e2i - i2r - h * mortality_rate * i[step]
        r[step + 1] = r[step] + i2r - h * mortality_rate * r[step]

    return s, e, i, r

s, e, i, r = seir_model_births_deaths()

@show_plot
def plot_me():
    s_plot = matplotlib.pyplot.plot(times, s, label='S')
    e_plot = matplotlib.pyplot.plot(times, e, label='E')
    i_plot = matplotlib.pyplot.plot(times, i, label='I')
    r_plot = matplotlib.pyplot.plot(times, r, label='R')
    n_plot = matplotlib.pyplot.plot(times, s + e + i + r, label='N')
    matplotlib.pyplot.legend(('S', 'E', 'I', 'R', 'N'), 'upper right')

    axes = matplotlib.pyplot.gca()
    axes.set_xlabel('Time in days')
    axes.set_ylabel('Number of persons')
    matplotlib.pyplot.xlim(xmin = 0.)
    matplotlib.pyplot.ylim(ymin = 0.)

plot_me()
