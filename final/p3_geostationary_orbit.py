# PROBLEM 3
#
# A rocket orbits the earth at an altitude of 200 km, not firing its engine.  When
# it crosses the negative part of the x-axis for the first time, it turns on its
# engine to increase its speed by the amount given in the variable called boost and
# then releases a satellite.  This satellite will ascend to the radius of
# geostationary orbit.  Once that altitude is reached, the satellite briefly fires
# its own engine to to enter geostationary orbit.  First, find the radius and speed
# of the initial circular orbit.  Then make the the rocket fire its engine at the
# proper time.  Lastly, enter the value of boost that will send the satellite
# into geostationary orbit.
#

import math
from udacityplots import *

earth_mass = 5.97e24 # kg
earth_radius = 6.378e6 # m (at equator)
gravitational_constant = 6.67e-11 # m3 / kg s2

total_duration = 9. * 3600. # s
marker_time = 0.25 * 3600. # s
tolerance = 100000. # m

# Task 1: Use Section 2.2 and 2.3 to determine the speed of the inital circular orbit.
initial_radius = earth_radius + 200000
# orbital_speed is found by setting centripetal force F_c = v^2/r * m_s equal to
# gravitational force F_g = G * m_1 * m_2 / r^2
#     v^2 / r * m_s = G * m_s * m_e / r^2
# =>  v^2 = (G * m_e) / r
# =>  v = sqrt( G * m_e / r)
initial_speed = math.sqrt(gravitational_constant * earth_mass / initial_radius)
final_radius = 42164e3

# Task 3: Which is the appropriate value for the boost in velocity? 2.453, 24.53, 245.3 or 2453. m/s?
# Change boost to the correct value.
boost = 2453. # m / s

def acceleration(position):
    return -gravitational_constant * earth_mass / numpy.linalg.norm(position)**3 * position

@show_plot
def send_to_geostationary():
    position = numpy.array([initial_radius , 0.]) # m
    velocity = numpy.array([0., initial_speed]) # m/s
    current_time = 0.
    h = 0.1 # s, set to initial step size but will store current step size
    h_new = h # s, will store the adaptive step size of the next step
    previous_marker_number = -1;
    geostationary_time = 0.

    boost_done = False
    radius_reached = False
    while current_time < total_duration:
        acceleration0 = acceleration(position)
        velocityE = velocity + h * acceleration0
        positionE = position + h * velocity
        velocityH = velocity + h * 0.5 * (acceleration0
                                + acceleration(positionE))
        positionH = position + h * 0.5 * (velocity + velocityE)

        # Task 2: When the rocket crosses the negative part of the x axis,
        # fire its engine to increase the speed by the amount given in the variable boost.

        if not boost_done and positionH[1] < 0.:
            velocityH += boost * velocityH / numpy.linalg.norm(velocityH)
            boost_done = True

        velocity = velocityH
        position = positionH

        error = numpy.linalg.norm(positionE - positionH) + total_duration \
                * numpy.linalg.norm(velocityE - velocityH)
        h_new = h * math.sqrt(tolerance / error)

        h_new = min(0.5 * marker_time, max(0.1, h_new))

        if current_time >= marker_time * previous_marker_number:
            previous_marker_number += 1
            matplotlib.pyplot.scatter(position[0], position[1], s = 2., facecolor = 'r', edgecolor = 'none')

        current_time += h
        h = h_new

        if not radius_reached and abs(numpy.linalg.norm(position) - final_radius) <= 150e3:
            radius_reached = True
            geostationary_time = current_time

    axes = matplotlib.pyplot.gca()
    axes.set_xlabel('Longitudinal position in m')
    axes.set_ylabel('Lateral position in m')
    axes.add_patch(matplotlib.patches.Circle((0., 0.), earth_radius, facecolor = 'b', edgecolor = 'none'))
    axes.add_patch(matplotlib.patches.Circle((0., 0.), final_radius, facecolor = 'none', edgecolor = 'g'))
    matplotlib.pyplot.axis('equal')

    return radius_reached, geostationary_time
send_to_geostationary()
