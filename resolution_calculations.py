import numpy as np 
import mpmath as mpm
from matplotlib import pyplot as plt


def geometrical_resolution(energy,angle,vertical_beam,detector_pitch,bending_radius):
    """
    Function for calculating the geometrical resolution derived from Mori.et.al 2012
    """
    delta_theta = (vertical_beam+detector_pitch)/bending_radius
    geometric_energy_resolution = delta_theta*mpm.cot(np.radians(angle))*energy
    return geometric_energy_resolution

def calc_Etot(a):
    return ((a**2)+(0.065**2)+(0.15**2)+(0.20**2))**0.5


def anthea(energy,angle,vertical_beam,detector_pitch,bending_radius):
    numerator = mpm.cot(np.radians(angle))*((vertical_beam+detector_pitch)*0.001)*energy
    return numerator/0.5


def main():

    v_beam_size = 0.1
    theta = 84.2
    ev = 6490
    pixel_pitch = 0.05
    bending_radius = 500

    jung_reso = geometrical_resolution(ev,theta,v_beam_size,0.025,bending_radius)
    medi_reso = geometrical_resolution(ev,theta,v_beam_size,0.055,bending_radius)
    pilatus_reso = geometrical_resolution(ev,theta,v_beam_size,0.172,bending_radius)
    anthea_geo = anthea(ev,theta,v_beam_size,pixel_pitch,bending_radius)
    print(anthea_geo)

    jung_E = round(calc_Etot(jung_reso),2)
    medi_E = round(calc_Etot(medi_reso),2)
    pilatus_E = round(calc_Etot(pilatus_reso),2)

    print(jung_E)
    print(medi_E)
    print(pilatus_E)


    # input_pixel_pitchs = np.linspace(0.015,0.125,28)
    # a = geometrical_resolution(ev,theta,v_beam_size,pixel_pitch,bending_radius)
    # # b = geometrical_resolution(ev,theta,v_beam_size,0.025,bending_radius)
    # b = [float(geometrical_resolution(ev,theta,i,0.055,bending_radius)) for i in input_pixel_pitchs]
    # c = [calc_Etot(i) for i in b]

    # plt.plot(input_pixel_pitchs,c)
    # plt.xlabel("pixel pitch, μm")
    # plt.ylabel("Δ Etot")
    # plt.plot([0.055,0.055],[np.min(c),np.max(c)],linestyle='--',color='r',label=f"Medipix, Etot {medi_E}")
    # plt.plot([0.025,0.025],[np.min(c),np.max(c)],linestyle='--',color='b',label=f"Jungfrau, Etot {jung_E}")
    # plt.legend()
    # plt.show()
    


if __name__ == "__main__":

    main()