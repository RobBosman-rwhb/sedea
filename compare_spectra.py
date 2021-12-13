import os
import fnmatch
import argparse
import pickle
import numpy as np
from pathlib import Path

from numpy.lib.npyio import save
import xes_signal_analysis as xsa
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def process_args():
    """ 
    Function to process input arguments.

    """
    input_args = argparse.ArgumentParser()
    input_args.add_argument("-i","--input",type=str,help="Input regex for reuced data to compare")
    input_args.add_argument("-o","--output_folder",type=str,help="Output folder location for the created pickle file")
    input_args.add_argument("-p","--output_prefix",type=str,help="Set a new preference for the pickle file")
    
    input_args.add_argument("-sc","--do-scale",action='store_true',help="add to scale spectra")
    input_args.add_argument("-s","--save-comparison_set",action='store_true',help="add to save comparison set")

    args = input_args.parse_args()

    if not args.input:
        print("No input regex or pickle given exiting")
        exit()

    # Set outputs for the optional arguements.
    # output folder
    if not args.output_folder:
        args.output_folder = os.getcwd()

    # output prefix
    if not args.output_prefix:
        filelist,_ = get_xes_datasets_from_regex("*comparison_set*.pickle",args.output_folder)
        next_iterator_number = len(filelist)+1
        args.output_prefix = f"Comparison_set_{next_iterator_number}"
        print(f"No output prefix set using {args.output_prefix}")
        


    return args


def get_xes_datasets_from_regex(regex,folder):
    """ 
    Function to return list of data object given an input regex
    """
    input_filelist =  [ file for file in os.listdir(folder) if fnmatch.fnmatch(file,regex) ]
    
    if folder == None:
        folder = os.path.dirname(regex)
        local_regex = os.path.basename(regex)
        print(folder)
        input_filelist =  [ os.path.join(folder,file) for file in os.listdir(folder) if fnmatch.fnmatch(file,local_regex) ]

    print(input_filelist)

    if len(input_filelist) == 0:
        print("No files found, please check regex")
        exit()

    xes_dataset_list = xsa.get_object_list(input_filelist)
    return xes_dataset_list,input_filelist


def plot_variation_comparison(object_list,name_list,time,vertROI,displacement_I,scaled,absolute):
    
    stddev_all = np.zeros((len(object_list)))
    stddev_ROI = np.zeros((len(object_list)))

    plt.figure(1113)
    for indx,i in enumerate(object_list):
        reduced_spectra = i.get_reduced_subtracted()
        if scaled == False:
            feature_std_red,difference_scaled1 = xsa.calc_diff_std(reduced_spectra)
        elif scaled == True:
            feature_std_red,difference_scaled1 = xsa.calc_feature_std(reduced_spectra)
        stddev_all[indx] = feature_std_red

        reduced_spectra_verticalROI = reduced_spectra[vertROI[0]:vertROI[1]]
        if scaled == False:
            feature_std_redROI,difference_scaledROI = xsa.calc_diff_std(reduced_spectra_verticalROI)
        elif scaled == True:
            feature_std_redROI,difference_scaledROI = xsa.calc_feature_std(reduced_spectra_verticalROI)
        stddev_ROI[indx] = feature_std_redROI
        
        if absolute == True:
            plt.plot(xsa.scale_spectra(reduced_spectra)+(displacement_I*(indx+1)),label=f"{name_list[indx]}, σ{feature_std_red:.3}")
            plt.text(x=0,y=displacement_I*(indx+1-displacement_I/2),s=f"{name_list[indx]}")
        elif absolute == False:
            plt.plot(difference_scaled1+(displacement_I*(indx+1)),label=f"{name_list[indx]}, σ{feature_std_red:.3}")
            plt.text(x=0,y=displacement_I*(indx+1-displacement_I/2),s=f"{name_list[indx]}")            
        # plt.legend()
    plt.axvline(vertROI[0],linestyle='--',color='k')
    plt.axvline(vertROI[1],linestyle='--',color='k')
    plt.xlabel("Pixels",fontsize=12)
    plt.ylabel("I arb",fontsize=12)

    print(stddev_all)
    plt.figure(1114)
    pars1, cov = curve_fit(f=xsa.func_powerlaw, xdata=time, ydata=stddev_all, p0=[0, 0], bounds=(-np.inf, np.inf))
    fitted_curve_all = [ xsa.func_powerlaw(i,pars1[0],pars1[1]) for i in np.linspace(1,100,100) ]
    pars2, cov = curve_fit(f=xsa.func_powerlaw, xdata=time, ydata=stddev_ROI, p0=[0, 0], bounds=(-np.inf, np.inf))
    fitted_curve_ROI = [ xsa.func_powerlaw(i,pars2[0],pars2[1]) for i in np.linspace(1,100,100) ]
    residuals_all = stddev_all - xsa.func_powerlaw(time, *pars1)
    residuals_ROI = stddev_ROI - xsa.func_powerlaw(time, *pars2)

    plt.plot(fitted_curve_all)
    plt.plot(fitted_curve_ROI)
    plt.scatter(time,stddev_all,label="σAll")
    plt.scatter(time,stddev_ROI,label="σROI")
    plt.xlabel("Exposure Time (sec)",fontsize=12)
    plt.ylabel("Pixel variation stddev σ",fontsize=12)
    plt.legend()

    plt.figure(1115)
    pixel_values = np.linspace(vertROI[0],vertROI[1],vertROI[1]-vertROI[0])
    reduced_spectra = object_list[0].get_reduced_avg()
    if scaled == True:
        plt.plot(pixel_values,xsa.scale_spectra(reduced_spectra[vertROI[0]:vertROI[1]]),label=f"{name_list[0]}, scaled")
    elif scaled == False: 
        plt.plot(pixel_values,xsa.scale_spectra(reduced_spectra[vertROI[0]:vertROI[1]]),label=f"{name_list[0]}")
    reduced_spectra = object_list[-1].get_reduced_avg()
    if scaled == True:
        plt.plot(pixel_values,xsa.scale_spectra(reduced_spectra[vertROI[0]:vertROI[1]]),label=f"{name_list[-1]}, scaled")
    elif scaled == False:
        plt.plot(pixel_values,xsa.scale_spectra(reduced_spectra[vertROI[0]:vertROI[1]]),label=f"{name_list[-1]}")
    plt.ylabel("Intensity",fontsize=12)
    plt.xlabel("Pixels",fontsize=12)
    plt.legend()
    plt.show()


def save_comparison_set(args,stored_files):
    """
    Saves out a comparison set as a pickle file can be used for reloading comparison.
    """
    complete_file = args.output_prefix+"_comparison_set.pickle"
    full_output = os.path.join(args.output_folder,complete_file)
    print(f"Saving output as {full_output}")
    with open(full_output,'wb') as output_handle:
        pickle.dump(stored_files,output_handle)
    output_handle.close()



def main():
    input_args = process_args()
    object_list,filelist = get_xes_datasets_from_regex(input_args.input,folder=None)
    name_list = [ i.get_dataset_name()  for i in object_list ]
    time_list = [ i.get_single_exposure_time() for i in object_list ]


    plot_variation_comparison(object_list,name_list=name_list,
                                time=time_list,
                                scaled=input_args.do_scale,absolute=True,
                                vertROI=[500,750],displacement_I=0.45)

    if input_args.save_comparison_set == True:               
        save_comparison_set(input_args,filelist)

if __name__ == "__main__":
    main()
