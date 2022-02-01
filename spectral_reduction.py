#!/usr/bin/python3

from matplotlib import pyplot as plt
from matplotlib import colors
from numpy.core.fromnumeric import std
from scipy import interpolate
from pathlib import Path
from xes_dataset import XesDataset as xs
from xes_signal_analysis import calc_photon_accumulation_rate as photon_rate
import xes_signal_analysis as xsa
import numpy as np
import warnings
import argparse
import pickle
import copy
import h5py
import time
import sys
import os

## Loading code functions

def process_args():
    """
    Handles the argparse workflow for spectral resolution 
    returns the args object for use in the code later
    """

    #   Arguements related to inputs.
    input_args = argparse.ArgumentParser()
    input_args.add_argument(
        "-i","--input_datafile",
        type=str,
        help="Input hdf5 or nexus file with XES data"
        )

    input_args.add_argument("-o","--output_folder",
        type=str,
        help="""Output folder for output files, will use the pwd 
        if not specified"""
        )

    input_args.add_argument("-p","--output_prefix",
        type=str,
        help="""Prefix for output files, will use Input datafile if 
        not specified"""
        )

    input_args.add_argument("-fo","--focal_orientation",
        type=int,
        help="""Set whether the orientation the 
        first or second array dimension"""
        )

    input_args.add_argument("-rmin","--roi_min",
        type=int,
        help="Set True if you want us to plot, default is 0"
        )
    
    input_args.add_argument("-rmax","--roi_max",
        type=int,
        help="Set True if you want us to plot, default is 0"
        )

    input_args.add_argument("-ff","--flatfield",
        type=str,
        help="File to apply flatfield correction."
        )

    input_args.add_argument("-bgg","--background_gap",
        type=int,
        help="Seperation between ROI and background area, by deafult 0."
        )

    input_args.add_argument("-bgw","--background_width",
        type=int,
        help="Background section width, by default same as ROI width."
        )

    #    boolen arguements asking the program to do something else.
    input_args.add_argument("-x","--plot_reduced",
        action='store_true',
        help="Set, if you want us to plot")

    input_args.add_argument("-s","--save_output",
        action='store_true',
        help="""Set true to save the reduced data, produces pickle 
        file with dataset object"""
        )

    input_args.add_argument("-ab","--average_borders",
        action='store_true',
        help="include to  average the borders of the image"
        )

    input_args.add_argument("-pc","--plot_colourmap",
        action='store_true',
        help="include to  average the borders of the image"
        )

    input_args.add_argument("-th","--thresholding",
        action='store_true',
        help="include to perform averaging"
        )

    args = input_args.parse_args()

    if not args.input_datafile:
        print("No input file given exiting")
        exit()

    # Set outputs for the optional arguements.
    # output folder
    if not args.output_folder:
        args.output_folder = os.getcwd()
    # output prefix
    if not args.output_prefix:
        args.output_prefix = Path(args.input_datafile).stem
        print(f"No output prefix set using {args.output_prefix}")
    # focus orientation
    if not args.focal_orientation:
        args.focal_orientation = None

    if not args.plot_reduced:
        args.plot_reduced = False

    if not args.roi_min and not args.roi_max:
        args.roi_min = False
        args.roi_max = False

    if not args.background_gap:
        args.background_gap = 0

    if not args.background_width:
        args.background_width = args.roi_max - args.roi_min

    if not args.flatfield:
        args.flatfield = False
        print("No flatfield given flat field not applied")

    return args

def load_data(input_args):
    datafile = input_args.input_datafile
    data = h5py.File(datafile,'r') #open raw data .hdf file						# open data file
    datainfo = dict()
    dataset = data["entry/data/data"][()]
    dataset_shape = np.array(dataset.shape)
    datainfo['dataset_shape'] = dataset_shape
    datainfo['exp_time'] = np.array(data["/entry/instrument/NDAttributes/Shutter Time"])
    datainfo['mean_values'] = np.array(data["/entry/instrument/NDAttributes/StatsMean"])
    datainfo['frame_no'] = len(np.array(data["/entry/instrument/NDAttributes/Frame Number"]))
    datainfo['frame_index'] = int(list(np.where(datainfo['frame_no']==datainfo['dataset_shape']))[0])
    
    image_dimensions = np.delete(dataset_shape,datainfo['frame_index'])
    datainfo['image_dimensions'] = image_dimensions

    if input_args.focal_orientation == None:
        focal_orientation = np.where(np.max(image_dimensions)==image_dimensions)[0]
        input_args.focal_orientation = focal_orientation
        datainfo['focal_oritenation'] = focal_orientation
    else:
        focal_orientation = input_args.focal_orientation

    if input_args.roi_min is False:
        roi_selection_dim = np.delete(image_dimensions,focal_orientation)
        roi_selection_dim = np.where(np.max(roi_selection_dim))
        input_args.roi_min = (image_dimensions[roi_selection_dim]/2)-10
        input_args.roi_max = (image_dimensions[roi_selection_dim]/2)+10
        input_args.background_width = input_args.roi_max - input_args.roi_min

    return dataset,datainfo,input_args

def load_flatfield(input_flatfield):
    """
    Applying flatfield to the dataset, take input flatfield file
    """
    
    if not os.path.isfile(input_flatfield):
        print ("Flatfield file not found, skipping flatfield")
        return False

    with open(input_flatfield) as f:
        FF_data = f.readlines()
        rows = len(FF_data)
    f.close()
    columns = len(FF_data[0].split(" "))
    FlatField = np.zeros((rows,columns),dtype=float)

    for i in range(rows):
        FF_data_line = FF_data[i].split(" ")
        for indx,a in enumerate(FF_data_line):
            FlatField[i,indx] = float(a.strip('/n'))
    
    return FlatField

## Determine value functions


## Image correction and processing functions

def apply_thresholds(datread_filt):
    datread_filt[datread_filt<threshold] = 0								            # datread_filt set everything below threshold to 0
    datread_filt[((datread_filt>threshold_upper) * (datread_filt<threshold_2ph))] = 0	# datread_filt thresholding both timescale
    datread_filt[datread_filt>threshold_upper_2ph] = 0							    # datread_filt removes data above the upper theshold for 2photons
    return datread

def average_over_gap(dataimage,sensor_edge):
    lower_bound = sensor_edge
    upper_bound = lower_bound+5
    steps = upper_bound-lower_bound
    for i in range(steps):
        preside = np.average(dataimage[:,lower_bound-steps-1+i:lower_bound-1+i],axis=1)
        afterside = np.average(dataimage[:,upper_bound+i+1:upper_bound+steps+i+1],axis=1)
        to_average_vectors = np.stack([preside,afterside],axis=1)
        new_value = np.average(to_average_vectors,axis=1)
        dataimage[:,lower_bound+i] = new_value
    return dataimage

def average_image_gaps(dataimage):
    sensor_gaps = 3
    for a in range(sensor_gaps):
        sensor_points = [255,514,773]
        dataimage = average_over_gap(dataimage,sensor_points[a])
    return dataimage

def average_borders(dataset):
    return average_image_gaps(dataset)
    
def apply_flatfield(image,FlatField):
    ff_size = FlatField.size
    data_size = image.size
    if ff_size != data_size:
        print("Flatfield array is inconsistent image size, not applying flat field")
        return image
    ff_corrected = np.multiply(image,FlatField)
    return ff_corrected

## Main code functions

def reduce_data(dataset,info,args):

    # Calculate bckglimits
    ROImin, ROImax = args.roi_min,args.roi_max
    ROI_bg1_min = ROImin - (args.background_gap + args.background_width)
    ROI_bg1_max = ROImin - args.background_gap
    ROI_bg2_max = ROImax + (args.background_gap + args.background_width)
    ROI_bg2_min = ROImax + args.background_gap

    dataset_images = len(dataset)
    print('Number of pulses in scan: ' + str(dataset_images))							# print number of pulse

    # Create output variables reduced data
    summed_image = np.zeros(np.shape(dataset[0]))							    # create an empty data table type thing
    summed_image_filt = np.zeros(np.shape(dataset[0]))						    # create an empty dataset for the filtered
    reduced_spectra_matrix = np.zeros((len(dataset),np.max(info['image_dimensions'])))					# create the empty for output
    bg1_reduced = np.zeros((len(dataset),np.max(info['image_dimensions'])))					# create empty for output
    bg2_reduced = np.zeros((len(dataset),np.max(info['image_dimensions'])))

    if args.flatfield is not False:
        ff_matrix = load_flatfield(args.flatfield)
        if ff_matrix is False:
            args.flatfield = False

    # Reduce data in loop
    for i in range(len(dataset)):
        single_shot_data = np.array(dataset[i])
        single_shot_data_filt = copy.deepcopy(np.array(dataset[i]))

        if args.flatfield is not False:
            single_shot_data = apply_flatfield(single_shot_data,ff_matrix)
            single_shot_data_filt = apply_flatfield(single_shot_data,ff_matrix)

        # If set, apply energy thresholds
        if args.thresholding == True:
            single_shot_data_filt = apply_thresholds(single_shot_data_filt)
            summed_image_filt = summed_image_filt + single_shot_data_filt

        # Apply zero averaging
            # Find the index for all the zeros
            # cycle through all the zeros and create an average of those around them.

        # Creating the ROI spectra, and associated background spectra.
        if args.focal_orientation == 0:
            reduced_spectra_matrix[i,:] = np.mean(single_shot_data_filt[ROImin:ROImax,:],axis=0)    
            bg1_reduced[i,:] = np.mean(single_shot_data_filt[ROI_bg1_min:ROI_bg1_max,:],axis=0) 
            bg2_reduced[i,:] = np.mean(single_shot_data_filt[ROI_bg2_min:ROI_bg2_max,:],axis=0) 
        if args.focal_orientation == 1:
            reduced_spectra_matrix[i,:] = np.mean(single_shot_data_filt[ROImin:ROImax,:],axis=0)      
            bg1_reduced[i,:] = np.mean(single_shot_data_filt[ROI_bg1_min:ROI_bg1_max,:],axis=0)   
            bg2_reduced[i,:] = np.mean(single_shot_data_filt[ROI_bg2_min:ROI_bg2_max,:],axis=0) 

        summed_image = summed_image + single_shot_data

    if args.average_borders == True:
        summed_image = average_borders(summed_image)
        reduced_spectra_matrix = average_borders(reduced_spectra_matrix)
        bg1_reduced = average_borders(bg1_reduced)
        bg1_reduced = average_borders(bg2_reduced)   

    # Subtraction on a shot-by-shot basis
    shot_by_shot_subtracted = reduced_spectra_matrix - np.mean([bg1_reduced,bg2_reduced],axis=0)
    reduced_subtracted = np.mean(shot_by_shot_subtracted,axis=0)


    to_return_xes_dataset_object = xs(reduced_spectra_matrix,
                                    bg1_reduced,
                                    bg2_reduced,
                                    summed_image,
                                    [ROImin,ROImax,ROI_bg1_min,ROI_bg1_max,ROI_bg2_min,ROI_bg2_max],
                                    shot_by_shot_subtracted=shot_by_shot_subtracted,
                                    reduced_subtracted=reduced_subtracted,
                                    dataset_name=args.output_prefix,total_images=dataset_images)
    
    return to_return_xes_dataset_object

def calculate_timing_information(dataset_info):
    """
    Function extract the timing information from the detector
    """
    shutter_time_array = dataset_info["exp_time"]
    assert np.all(shutter_time_array == shutter_time_array[0])
    single_exposure = float(shutter_time_array[0])
    total_exposure = np.sum(shutter_time_array)
    return total_exposure,single_exposure    

## process after the image

def interpolate_line(line,multiplier,smooth):
    interpolation_sampling = (len(line))*multiplier
    x_axis_evaluation = np.linspace(1,len(line),len(line))
    interpolation_points = np.linspace(1,len(line),interpolation_sampling)
    cubic_interp_coefficients = interpolate.splrep(x_axis_evaluation, line, s=smooth)
    return interpolate.splev(interpolation_points, cubic_interp_coefficients, der=0),x_axis_evaluation,interpolation_points

def interpolate_peak(reduced_data,signal_ROImin,signal_ROImax):
    reduced_subtracted_spectra = reduced_data.get_reduced_subtracted()
    peak_spectra = reduced_subtracted_spectra[signal_ROImin:signal_ROImax]
    interpolated_line,old_xaxis,new_xaxis = interpolate_line(peak_spectra,5,30)
    reduced_data.set_interpolated_peak(interpolated_line)
    return reduced_data

def focally_interpolate_spectra(dataset_obj,image,ROI_min,ROI_max,focal_dimension):
    focally_interpolated_spectra = np.zeros((focal_dimension))
    for i in range(focal_dimension):
        line = image[ROI_min:ROI_max,i]
        interpolated_line,old_xaxis,new_xaxis=interpolate_line(line,3,0)
        focally_interpolated_spectra[i] = max(interpolated_line)/(len(line)*5)
    dataset_obj.set_focal_interpolation(focally_interpolated_spectra)
    return 

## Plotting functions
# def function_fits(x_data,y_data):
#     power_pars, cov = curve_fit(f=xsa.func_powerlaw, xdata=x_data, ydata=y_data,
#                              p0=[0.05, -0.1], bounds=(-1, 1),method='dogbox')
#     power_curve_fit = [ xsa.func_powerlaw(i,power_pars[0],power_pars[1]) for i in image_sampling ]
#     poly_fit_1 = 



def plot_selection_histogram(dataset,threshold,threshold_upper,threshold_2ph,threshold_upper_2ph):
    histodata1, histoedges = np.histogram(dataset[:50], bins=np.arange(2000,6000,10))				# create a histogram with bin set at 0.1 intervals between -5 to 30
    plt.figure(1113)												# create a plt with name 1113 (#matlab)
    plt.plot(histoedges[:-1], histodata1)									# plot the histogram that you just created
    plt.xlim([2000, 6000])												# set limits to what we set the ranges to
    plt.yscale('log')												# give it a log scale 
    plt.axvline(x=threshold, color = 'r')									# set sections defined by the threshold as 'r'
    plt.axvline(x=threshold_upper, color = 'r')									# See above
    plt.axvline(x=threshold_2ph, color = 'r')									# See above
    plt.axvline(x=threshold_upper_2ph, color = 'r')								# See above
    # plt.savefig(Figsavepath + 'adva_histo_Run' +runstr + '.png', dpi = 600)					# See above
    plt.show()

def plot_detector_image(dataset_obj,args):
    image = dataset_obj.get_summed_dataimage()
    reduced_subtracted = dataset_obj.get_reduced_subtracted()
    bg1 = dataset_obj.get_bg1_spectra()
    bg2 = dataset_obj.get_bg2_spectra()
    bckg_mean = (np.mean(bg1,axis=0)+np.mean(bg2,axis=0))/2
    roi_selections = dataset_obj.get_ROI()

    single_exp_time = dataset_obj.get_single_exposure_time()
    photons_per_sec = photon_rate(bckg_mean,single_exp_time)

    # Calculate the commit times for wonderful things!!!

    power_fit = dataset_obj.get_power_fit_dict()
    running_diff_std = power_fit["running_std"]
    number_photons = power_fit["tot_photons_array"]
    est_double_std = power_fit["est_double_std"]
    percentage_improvement = power_fit["expected_improv_percent"]
    power_curve_fit = power_fit["fitted_std_curve"]
    cutoff_photons = power_fit["cutoff_photon"]

    # poly_fit_1 = 
    printable_final_std = round(running_diff_std[-1],8)
    printable_est_std = round(est_double_std,8)
    printable_percent_imp = round(percentage_improvement,2)


    # plt.figure(1113)
    plt.subplot(324)
    plt.scatter(number_photons,running_diff_std,marker="1",s=15,label=f"Final std = {printable_final_std}")
    plt.plot(number_photons,power_curve_fit,color='red')
    plt.plot([0,number_photons[-1]],[est_double_std,est_double_std],
            label=f"2x std = {printable_est_std}, % improvement={printable_percent_imp})",color='k')
    # plt.plot([cutoff_photons,cutoff_photons],[est_double_std,running_diff_std[0]],lineStyle='--',color='k')
    # Needs reworking probably with the algo
    plt.ylabel("pixel-to-pixle σ")
    plt.xlabel("Total number photons")
    plt.legend()


    # plt.subplot(212)
    # plt.plot(np.diff(running_diff_std))

    moving_average,interp_x,interp_y=xsa.calculate_interpolation(reduced_subtracted,10,2)
    max_interp = round(np.max(interp_y),4)
    max_average = round(np.max(moving_average),4)

    linear_interp = interpolate.UnivariateSpline(interp_x,interp_y-max_interp/2,s=0)
    
    r1,r2 = linear_interp.roots()
    fwhm1 = r2-r1


    plt.subplot(322)
    plt.plot(reduced_subtracted,label=f"Reduced spectra",color='k',linewidth=0.1)
    plt.plot(moving_average,label=f"Moving average, peak max={max_average}",linestyle='--',color='r')
    plt.plot(interp_x,interp_y,label=f"Cubic spline, peak max={max_interp}",linestyle='--',color='b')
    plt.axvspan(r1,r2,label=f"FWHM = {fwhm1}",alpha=0.5)
    plt.legend()


    # plt.figure(1114)
    plt.subplot(321)

    if args.plot_colourmap == True:
        plt.imshow(image,cmap='jet',norm=colors.LogNorm())
    else:
        plt.imshow(image, vmin=np.percentile(image,1), vmax=np.percentile(image,99))
    # central selection
    plt.plot([0,1032],[args.roi_min,args.roi_min],color='r',linestyle='dashed')
    plt.plot([0,1032],[args.roi_max,args.roi_max],color='r',linestyle='dashed')
    # upper roi selection
    plt.plot([0,1032],[roi_selections[2],roi_selections[2]],color='k',linestyle='dotted')
    plt.plot([0,1032],[roi_selections[3],roi_selections[3]],color='k',linestyle='dotted')
    # lower roi selection
    plt.plot([0,1032],[roi_selections[4],roi_selections[4]],color='k',linestyle='dotted')
    plt.plot([0,1032],[roi_selections[5],roi_selections[5]],color='k',linestyle='dotted')

    plt.subplot(323)
    plt.plot(reduced_subtracted,label="reduced spectra")
    plt.legend()

    plt.subplot(325)
    plt.plot(bckg_mean,label=f"{photons_per_sec} ph/s")
    plt.legend()
    plt.show()


## Output function the dark lord of sauron

def save_output(dataset_obj,args):
    complete_file = args.output_prefix+"_reduced.pickle"
    full_output = os.path.join(args.output_folder,complete_file)
    print(f"Saving output as {full_output}")
    with open(full_output,'wb') as output_handle:
        pickle.dump(dataset_obj,output_handle)
    output_handle.close()


def main():
    start_timer = time.process_time()
    args = process_args()
    data,dataset_info,args = load_data(args)
    # focal_dimension, focal_orientation = set_focal_orientation(data, args.focal_orientation)
    reduced_data = reduce_data(data, dataset_info,args)
    # calculate_stats = (data,dataset_info,reduced_data,args)
    # focally_interpolate_spectra(reduced_data,reduced_data.get_summed_dataimage(),
                                        # reduced_data.get_ROI()[0],
                                        # reduced_data.get_ROI()[1],
                                        # focal_dimension)
    # reduced_data = interpolate_peak(reduced_data,525,760)
    total_exposure,single_exposure = calculate_timing_information(dataset_info)
    reduced_data.set_total_exposure_time(total_exposure)
    reduced_data.set_single_exposure_time(single_exposure)

    # assess_noise_suppression(reduced_data)

    if args.plot_reduced == True:
        plot_detector_image(reduced_data,args)

    if args.save_output == True:
        save_output(reduced_data,args)

    end_timer = time.process_time()
    print(f"Total run time : {end_timer-start_timer}")



if __name__ == "__main__":

    main()

