import numpy as np 
import math as mh 
import mpmath as mpm


def geometrical_resolution(energy,angle,vertical_beam,detector_pitch,bending_radius):
    """
    Function for calculating the geometrical resolution derived from Mori.et.al 2012
    """
    delta_theta = (vertical_beam+detector_pitch)/bending_radius
    geometric_energy_resolution = delta_theta*mpm.cot(np.radians(angle))*energy
    return geometric_energy_resolution


def main():

    v_beam_size = 0.015
    theta = 79.1
    ev = 6380
    pixel_pitch = 0.055
    bending_radius = 400

    a = geometrical_resolution(ev,theta,v_beam_size,pixel_pitch,bending_radius)
    b = geometrical_resolution(ev,theta,v_beam_size,0.025,bending_radius)
    print(a-b)


if __name__ == "__main__":

    main()