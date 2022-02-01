from functools import reduce
from operator import sub
import numpy as np 
import pickle
import copy
from matplotlib import pyplot as plt
from numpy.core.fromnumeric import std
from numpy.core.numeric import convolve
from numpy.lib import histograms
from numpy.lib.histograms import histogram
import scipy as sp
from xes_dataset import XesDataset as xs
from scipy.optimize import curve_fit
from scipy.ndimage.filters import uniform_filter1d
from scipy import constants
import scipy.interpolate
import math
import xes_spectral_decomposition as xsd
import energy_range_graph as en


def PixToE(Offset, pixel):
    e = constants.e  # electron charge
    c = constants.c  # speed of light	
    h = constants.h  # Planck constant	
    pi = constants.pi
    R = 0.50								# Hard coded that the R value is here, XX Would need to change this XX
    # Pixel Pitch = pixel height again hard coded for the epics detector XX need to change this XX
    PixelPitch = 55e-6
    # ll, apears to locate a given pixel on the detector.
    ll = (Offset + pixel*PixelPitch)
    # http://www.psi.ch/sls/superxas/johann-xes-spectrometer#Si_733_1_crystal

    dd = en.dspace_Calc("Si",[4,4,0])*10**-10		
    # d_Ge440 = 0.5*2.00040508e-10
    # dd = d_Ge440
    # print(dd)		
    # calculate the energy for that pixel according to robertos paper.
    EnergiesCalc = (h*c/e)/(2*dd*np.sin((pi/2) - np.arctan(ll/(2*R))))
    # EnergiesCalc =  (ll/(2*R))
    return EnergiesCalc		

##READ FUNCTIONS ###
def get_object_list(pickle_file_list):
    object_list = []
    for i in pickle_file_list:
        object_list.append(get_reduced_data(i))
    return object_list

def get_reduced_data(pickle_file):
    open_pickle = open(pickle_file,'rb')
    return pickle.load(open_pickle)

def calc_tiling(num):
    subplot_tiles = [0,0]
    square_root = math.sqrt(num)
    if square_root == 1:
        subplot_tiles = [1,1]
    elif square_root > 1 and square_root < 1.5:
        subplot_tiles = [2,1]
    return subplot_tiles
        
def get_spectrum_variability(object_list):
    reduced_spectra = object_list[0].get_reduced_subtracted()
    feature_std_red,variabiliy = calc_diff_std(reduced_spectra)
    reduced_spectrum = np.zeros((len(object_list),len(variabiliy)))
    for i in range(len(object_list)):
        reduced_spectra = object_list[i].get_reduced_subtracted()
        feature_std_red,variability = calc_diff_std(reduced_spectra)
        reduced_spectrum[i,:] = variability
    return reduced_spectrum

def get_reduced_spectrum_subtracted(object_list):
    reduced_spectra = object_list[0].get_reduced_subtracted()
    reduced_spectrum = np.zeros((len(object_list),len(reduced_spectra)))
    for i in range(len(object_list)):
        reduced_spectra = object_list[i].get_reduced_subtracted()
        reduced_spectrum[i,:] = reduced_spectra
    return reduced_spectrum

### STDDEV FUNCTION ###
def scale_spectra(spectral_line):
    return np.divide(np.subtract(spectral_line, min(spectral_line)),max(spectral_line)-min(spectral_line))

def calc_feature_std(spectral_line):
    scaled = scale_spectra(spectral_line)
    diff_scaled = np.diff(scaled)
    return np.std(diff_scaled),diff_scaled

def calc_diff_std(feature_line):
    diff = np.diff(feature_line)
    return np.std(diff),diff



### BACKGROUND ASSESSMENT FUNCTIONS ###

def calc_photon_accumulation_rate(detector_area,exposure_time):
    print(exposure_time)
    return np.mean(detector_area)/exposure_time

def calculate_histogram(data,binning):
    return np.histogram(data,bins=np.linspace(data.min(),data.max(),binning))

def calculate_interpolation(spectra,moving_windows,smoothing):
    spectral_length = len(spectra)
    reduced_spectra_convolved = uniform_filter1d(spectra, size=moving_windows)
    x = np.linspace(0, spectral_length, num=spectral_length, endpoint=True)

    # Calculate smoothing
    # Inplemented automatic smoothness suggested here https://github.com/scipy/scipy/issues/11916
    smoothing = 1/np.std(spectra)
    print(smoothing)
    reduced_spectra_interp1d = scipy.interpolate.splrep(x, spectra,s=smoothing)
    newx = np.linspace(1,spectral_length,spectral_length)
    newy = scipy.interpolate.splev(newx,reduced_spectra_interp1d,der=0)
    return reduced_spectra_convolved,newx,newy


### OUTPUT PLOTTING FUNCTIONS ###
def process_std_deviation(xes_dataset_object):
    summed_dataimage = xes_dataset_object.get_summed_dataimage()
    reduced_spectra = xes_dataset_object.get_reduced_subtracted()
    bckg1_average = xes_dataset_object.get_bckgA_avg()
    bckg2_average = xes_dataset_object.get_bckgB_avg()
    red_average = xes_dataset_object.get_reduced_avg()
    background_spectra = np.mean([bckg1_average,bckg2_average],axis=0)

    feature_std_red,difference_scaled1 = calc_feature_std(reduced_spectra)
    feature_std_raw,difference_scaled2 = calc_feature_std(red_average)
    bckg_spectra_std,background_variation = calc_diff_std(background_spectra)

    plt.figure(1112)
    plt.subplot(131)
    plt.plot(background_spectra,label=f"Background, sigma={bckg_spectra_std:.5}")
    plt.plot(red_average,label=f"Signal, sigma={bckg_spectra_std:.5}")
    plt.legend()

    plt.subplot(132)
    plt.plot(reduced_spectra,label="background subtracted")
    plt.legend()

    plt.subplot(133)
    plt.plot(difference_scaled1,label=f"$Std_feat$ = {feature_std_red:.3}")
    plt.legend()

    plt.show()


def histogram_plots(object_list,name_list):

    plt.figure(1116)
    for indx,i in enumerate(object_list):
        plt.subplot(len(object_list),1,indx+1)
        bckg1_average = i.get_bckgA_avg()
        bckg2_average = i.get_bckgB_avg()
        background_spectra = np.mean([bckg1_average,bckg2_average],axis=0)
        bckg_spectra_std,background_variation = calc_diff_std(background_spectra)
        counts,bins = calculate_histogram(background_variation,150)
        plt.hist(bins[:-1],bins,weights=counts,label=f"{name_list[indx]}, {bckg_spectra_std:.5}",alpha=0.8)
        plt.xlim(-700,700)
        plt.xticks(fontsize=6)
        plt.yticks(fontsize=6)
        plt.legend(fontsize=8)
    plt.show()

def plot_interpolation(spectra,convolved,spline_newx,spline_newy,ROI):
    plt.figure("Fittings methods")
    plt.subplot(211)
    plt.plot(spectra,label="Reduced spectra",linewidth=1.0,color='k')
    plt.plot(convolved,label="Moving average, window size 20",linewidth=2,color='r')
    plt.ylabel("Intensity (a.u.)",fontsize=18)
    plt.tick_params(labelsize=12)
    plt.legend(fontsize=16)

    plt.subplot(212)
    plt.plot(spectra,label="Reduced spectra",linewidth=1.0,color='k')
    plt.plot(spline_newx,spline_newy,label="Cublic spline",linewidth=2,color='r')
    plt.xlabel("Pixels",fontsize=18)
    plt.ylabel("Intensity (a.u.)",fontsize=18)
    plt.tick_params(labelsize=12)
    plt.legend(fontsize=16)
    plt.show()

    plt.figure("ROIs")
    plt.subplot(211)
    plt.plot(spectra,linewidth=1.0,color='k')
    plt.plot(convolved,label="Moving average, window size 20",linewidth=2,color='r')
    plt.ylabel("Intensity (a.u.)",fontsize=18)
    plt.xlabel("Pixels",fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim(ROI[0],ROI[1])
    plt.legend(fontsize=12)

    plt.subplot(212)
    plt.plot(spectra,linewidth=1.0,color='k')
    plt.plot(spline_newx,spline_newy,label="Cublic spline",linewidth=2,color='r')
    plt.xlabel("Pixels",fontsize=18)
    plt.tick_params(labelsize=12)
    plt.xlim(ROI[0],ROI[1])
    plt.legend(fontsize=12)
    plt.show()


def main():
    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Dark_test_01_1_reduced.pickle"
    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Anomolous_scattering_test_14_reduced.pickle"
    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Anomolous_scattering_test_4_reduced.pickle"
    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_26_reduced.pickle"
    # pickle_input2 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Energy_thresholding_test_1_reduced.pickle"

    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_1_reduced.pickle"
    # pickle_input2 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Detector_test_17_reduced.pickle"
    # pickle_input3 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Detector_test_18_reduced.pickle"
    # pickle_input4 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Detector_test_19_reduced.pickle"
    
    # file_pickle_list = [pickle_input1,pickle_input2,pickle_input3,pickle_input4]
    # object_list = get_object_list(file_pickle_list)  
    # name_list = ["Before","Votage change 1","Voltage change 2","Voltage change 3"]
    # histogram_plots(object_list,name_list)


    # exposure_time = 60
    # dataset_object = get_object_list([pickle_input1])[0]
    # process_std_deviation(dataset_object)
    # bckg1_average = dataset_object.get_bckgA_avg()
    # bckg2_average = dataset_object.get_bckgB_avg()
    # background_spectra = np.mean([bckg1_average,bckg2_average],axis=0)
    # photon_per_sec = calc_photon_accumulation_rate(background_spectra,[300,1033],exposure_time)
    # total_photons = np.mean(background_spectra[300:1033])
    # stddev,_ = calc_diff_std(background_spectra)
    # print(total_photons)
    # print(photon_per_sec)
    # print(stddev)
    # dataset_object.get_reduced_spectra()



    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_avrg_1_reduced.pickle"
    # pickle_input2 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_avrg_2_reduced.pickle"
    # pickle_input3 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_avrg_3_reduced.pickle"
    # pickle_input4 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_avrg_4_reduced.pickle"
    # pickle_input5 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_avrg_7_reduced.pickle"
    # pickle_input6 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_avrg_6_reduced.pickle"

    # Exposure time conclusion
    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_1_reduced.pickle"
    # pickle_input2 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_2_reduced.pickle"
    # pickle_input3 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_3_reduced.pickle"
    # pickle_input4 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_4_reduced.pickle"
    # pickle_input5 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_5_reduced.pickle"
    # pickle_input6 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_6_reduced.pickle"
    # pickle_input7 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_7_reduced.pickle"

    pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/MnCl2_1_reduced.pickle"
    pickle_input2 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/MnCl2_2_reduced.pickle"
    pickle_input3 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/MnCl2_3_reduced.pickle"
    pickle_input4 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/MnCl2_4_reduced.pickle"
    pickle_input5 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/MnCl2_5_reduced.pickle"
    pickle_input6 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/MnCl2_6_reduced.pickle"
    pickle_input7 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/MnCl2_7_reduced.pickle"

    #Calibrating the detector
    # pickle_input1 = "/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/VMXi_reduced_data/Manganese_foil_expl_7_reduced.pickle"



    # file_pickle_list = [pickle_input1,pickle_input2,pickle_input3,pickle_input4,pickle_input5,pickle_input6,pickle_input7]
    file_pickle_list = [pickle_input7]
    object_list = get_object_list(file_pickle_list)
    # name_list = ["100 sec","50 sec","25 sec","10 sec","1 sec","0.25 sec","0.04 sec"]
    name_list = ["0.05mM","0.5mM","5.0mM","50mM","100mM","250mM","500mM"]

    # reduced_spectra = object_list[-1].get_reduced_subtracted()

    # time = np.array([100,50,25,10,1,0.2,0.04])
    # vertROI = [500,750]
    # stddev_variation(object_list,name_list,time,vertROI,0.5,True,False)

    # Energy Calibration
    # convolved,newx,newy = calculate_interpolation(reduced_spectra,20,500)
    # plot_interpolation(reduced_spectra,convolved,newx,newy,[350,850])

    # energiesx =[ PixToE(0.0705,i) for i in newx ]

    # Slide 20 plots
    # plt.figure("Energy calibration")
    # plt.subplot(211)
    # plt.plot(reduced_spectra,label="Reduced spectra",linewidth=1.0,color='k')
    # plt.plot(newx,newy,label="Cublic spline",linewidth=2,color='r')
    # plt.xlabel("Pixels",fontsize=18)
    # plt.ylabel("Intensity (a.u.)",fontsize=18)
    # plt.tick_params(labelsize=12)
    # plt.legend(fontsize=16)

    # plt.subplot(212)
    # plt.plot(energiesx,reduced_spectra,label="Reduced spectra",linewidth=1.0,color='k')
    # plt.plot(energiesx,newy,label="Cublic spline",linewidth=2,color='r')
    # plt.xlabel("Energy",fontsize=18)
    # plt.ylabel("Intensity (a.u.)",fontsize=18)
    # plt.tick_params(labelsize=12)
    # plt.legend(fontsize=16)
    # plt.show()


    # plt.plot(energiesx,reduced_spectra,label="Reduced spectra",linewidth=1.0,color='k')
    # plt.plot(energiesx,newy,label="Mn Foil",linewidth=2,color='r')
    # plt.xlabel("Energy (eV)",fontsize=14)
    # plt.ylabel("Intensity (a.u.)",fontsize=14)
    # plt.tick_params(labelsize=12)
    # plt.legend(fontsize=16)
    # plt.xlim([6474,6501])
    # plt.show()    


    reduced_spectra = object_list[0].get_reduced_subtracted()
    
    # plot_interpolation(reduced_spectra,convolved,newx,newy,[350,850])

    for i in np.linspace(0,10000,10):
        convolved,newx,newy = calculate_interpolation(reduced_spectra,40,i)
        energiesx =[ PixToE(0.0705,i) for i in newx ]
        plt.plot(energiesx,newy,label=f"{i}",linewidth=2)
        plt.plot(energiesx,reduced_spectra,linewidth=0.5,color='k')

    plt.xlabel("Energy (eV)",fontsize=18)
    plt.ylabel("Intensity (a.u.)",fontsize=18)
    plt.tick_params(labelsize=12)
    plt.legend(fontsize=16)
    plt.xlim([6474,6501])
    plt.show()








    # plt.figure(1111)
    # max_rvec = np.max(R_mat[:,0])
    # max_PC1vec = np.max(reduced_spectrum[:,3])
    # scale = max_PC1vec/max_rvec
    # PC1_spectra = R_mat[:,0]*scale
    # scaled_reduced_spectrum = reduced_spectrum[:,3]
    # tmp = np.vstack([PC1_spectra,scaled_reduced_spectrum])
    # average_stuff = np.mean(tmp,axis=0)
    

        # reduced_spectra_adjusted[indx] = reduced_spectra_adjusted[indx]-PC1_spectra[indx]/2
        
    
    # plt.plot(PC1_spectra,label="PC1")
    # plt.plot(scaled_reduced_spectrum,label="100sec")
    # plt.plot(average_stuff,label="average")
    # plt.legend()
    # plt.show()

    # time = np.array([2,5,10,20,50,100])
    # vertROI = [500,750]
    # stddev_variation(object_list,name_list,time,vertROI)


    ### SVD ANALYSIS, WORK IN PROGRESS ###
    # reduced_spectrum = get_spectrum_variability(object_list)
    # reduced_spectrum = get_reduced_spectrum_subtracted(object_list)
    # # reduced_spectrum = reduced_spectrum[0:3,:]
    # reduced_spectrum = np.transpose(reduced_spectrum)
    # print(reduced_spectrum.shape)

    # xsd.plot_difference_matrix(reduced_spectrum,name_list)

    # U,S,Vh,S_mat,R_mat = xsd.calculate_SVD(reduced_spectrum)

    # xsd.plot_S_matrix(S)
    # xsd.plot_R_matrix(6,1033,50,R_mat*-1)
    # plt.xlabel("pixels",fontsize=14)
    # plt.ylabel("IArb",fontsize=14)
    # plt.show()




if __name__ == "__main__":

    main()