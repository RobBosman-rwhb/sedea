import os
import fnmatch
import argparse
import pickle
import numpy as np
from pathlib import Path

from numpy.lib.npyio import save
from xes_dataset import XesDataset
import xes_signal_analysis as xsa
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit


def process_args():
    """ 
    Function to process input arguments.
    Returns processed args objects for use in main.

    """
    input_args = argparse.ArgumentParser()
    input_args.add_argument("-i","--input",type=str,help="Input regex for reuced data to compare")
    input_args.add_argument("-o","--output_folder",type=str,help="Output folder location for the created pickle file")
    input_args.add_argument("-p","--output_prefix",type=str,help="Set a new preference for the pickle file")
    input_args.add_argument("-d","--difference_dark",type=str,help="dataset to subtract from the others")
    
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

    if not args.difference_dark:
        print(f"No dark given subtraction will not run")
        args.difference_dark == None

    return args


def calculate_difference_spectra(
        objects:list,
        dark:XesDataset
    ):
    """
    Function that takes a list of XesDataset objects and
    calculates difference spectra for each one.
    """
    dark_reduced = dark.get_reduced_subtracted()
    reduced_spectra = np.empty((len(objects),len(dark_reduced)))

    [ np.vstack((reduced_spectra,i.get_reduced_subtracted())) for i in objects ]
    difference_matrix = np.tile(dark_reduced,(len(objects),1))
    subtracted_difference_matrix = reduced_spectra-difference_matrix
    displacement = np.max(subtracted_difference_matrix)
    movement = np.max(np.diff(subtracted_difference_matrix))

    for i in range(len(objects)):
        plt.plot(subtracted_difference_matrix[i,:]+(displacement*(i)+movement))
        # plt.text(x=0,y=displacement*(indx+1),s=f"{name_list[indx]}")
    # [ plot_spectra_with_placement(i,displacement,name_shift=movement) for i in iter(difference_matrix)]

    ## Still trying to get this plotted out ##

    plt.show()


def get_xes_datasets_from_regex(regex,folder):
    """ 
    Function to return list of data objects from an input regex
    The regex follows the simple bash based pattern recognition.
    A folder variable is given for a situation where the files are
    looked for in another directory.
    """
    print(regex)
    input_filelist =  [ print(file) for file in os.listdir(folder) if fnmatch.fnmatch(file,regex) ]
    
    if folder == None:
        folder = os.path.dirname(regex)
        local_regex = os.path.basename(regex)
        input_filelist =  [ os.path.join(folder,file) for file in \
            os.listdir(folder) if fnmatch.fnmatch(file,local_regex) ]

    if len(input_filelist) == 0:
        print("No files found, please check regex")
        exit()

    xes_dataset_list = xsa.get_object_list(input_filelist)

    return xes_dataset_list,input_filelist


def plot_absolute_spectra_comparison(
        normalise_strategy:str,
        object_list:list,
        vertical_roi:tuple
        ):
    """
        Plots each spectra with an offset to see each clearly
        while still making visual comparison possible. Returns
        plot object that can be plotted as required.
    """
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

        total_variation[indx] = np.max(reduced_spectra)-np.min(reduced_spectra)
        variation_index[indx] = np.max(np.diff(reduced_spectra))
        
    displacement = np.mean(total_variation)
    movement = np.max(variation_index)

    for indx,i in enumerate(object_list):
        plt.plot(i.get_reduced_subtracted()+(displacement*(indx)+movement),label=f"{name_list[indx]}, σ{feature_std_red:.3}")
        plt.text(x=0,y=displacement*(indx+1),s=f"{name_list[indx]}")
           
    # plt.legend()
    plt.axvline(vertROI[0],linestyle='--',color='k')
    plt.axvline(vertROI[1],linestyle='--',color='k')
    plt.xlabel("Pixels",fontsize=12)
    plt.ylabel("I arb",fontsize=12)


def plot_comparison(object_list,name_list,time,vertROI,scaled):
    
    stddev_all = np.zeros((len(object_list)))   # stddev should be calculated for the datasets at the spectral reduction stage?
    stddev_ROI = np.zeros((len(object_list)))   # same as above, should create a xesDataset class with series of class methods on creation or recipet of data.
    total_variation = np.zeros((len(object_list)))
    variation_index = np.zeros((len(object_list)))

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

        total_variation[indx] = np.max(reduced_spectra)-np.min(reduced_spectra)
        variation_index[indx] = np.max(np.diff(reduced_spectra))
        
    displacement = np.mean(total_variation)
    movement = np.max(variation_index)

    for indx,i in enumerate(object_list):
        plt.plot(i.get_reduced_subtracted()+(displacement*(indx)+movement),label=f"{name_list[indx]}, σ{feature_std_red:.3}")
        plt.text(x=0,y=displacement*(indx+1),s=f"{name_list[indx]}")
           
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

    if input_args.difference_dark is not None:
        reduced_dark_data = xsa.get_reduced_data(input_args.difference_dark)
        calculate_difference_spectra(object_list,reduced_dark_data)


    plot_variation_comparison(object_list,name_list=name_list,
                                time=time_list,
                                scaled=input_args.do_scale,
                                vertROI=[500,750])

    if input_args.save_comparison_set == True:               
        save_comparison_set(input_args,filelist)

if __name__ == "__main__":
    main()
