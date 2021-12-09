import numpy as np 
import Energy_range_graph as engy
from scipy import constants


def pix2e (Offset, pixel,R,PixelPitch):
    PixelPitch = 55e-6
    ll = (Offset + pixel*PixelPitch)
    dd = engy.dspace_Calc("Si",[2,2,0])
    e = constants.e
    c = constants.c
    h = constants.h
    EnergyCalc = (h*c/e)/(2*dd*np.sin((constants.pi/2) - np.arctan(ll/(2*R))))
    # EnergiesCalc =  (ll/(2*R))
    return EnergyCalc


def main():
    print("boggle")


if __name__ == "__main__":
    main()
