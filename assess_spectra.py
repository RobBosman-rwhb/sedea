import numpy as np
import argparse
import xes_signal_analysis as xsa
from matplotlib import pyplot as plt
from scipy import interpolate



def process_args():
    input_args = argparse.ArgumentParser()
    input_args.add_argument(
        "-i","--input_datasetfile",
        type=str,
        help="""Input dataset file with reduced XES data"""
        )
    
    input_args.add_argument(
        "-s","--smoothing",
        type=float,
        help="""input for numpy.splrep for subic interpolation see 
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.splrep.html for details"""
        )

    input_args.add_argument(
        "-a","--average_window",
        type=int,
        help="""set the window size for moveing average"""
        )

    args = input_args.parse_args()

    if args.input_datasetfile == None:
        print("Input file not recognised")
        exit()

    if args.input_datasetfile.endswith(".pickle"):
        print(f"Input data file == {args.input_datasetfile}")
    else:
        print("Input file does not appear to be a pickle file")
        exit()
    
    if args.smoothing is None:
        print("Smoothing value is not set")

    return args


def plot_interpolations(dataset_obj,args):

    error_codes = [-1,-2]
    reduced_subtracted = dataset_obj.get_reduced_subtracted()
    power_fit = dataset_obj.get_power_fit_dict()

    if np.isin(power_fit,error_codes) == False:
        running_diff_std = power_fit["running_std"]
        number_photons = power_fit["tot_photons_array"]
        est_double_std = power_fit["est_double_std"]
        percentage_improvement = power_fit["expected_improv_percent"]
        power_curve_fit = power_fit["fitted_std_curve"]
        cutoff_photons = power_fit["cutoff_photon"]

        printable_final_std = round(running_diff_std[-1],8)
        printable_est_std = round(est_double_std,8)
        printable_percent_imp = round(percentage_improvement,2)
    elif power_fit == -2:
        print("Single image skipping stddev calc")
    elif power_fit == -1:
        print("Couldn't fit stddev")


    moving_average,max_average,y_average = dataset_obj.do_moving_average_for_plots()

    fwhm_dict = dataset_obj.do_fwhm_assesment_for_plots()

    cubic_dict = dataset_obj.do_cubic_interpolation_for_plots()

    if fwhm_dict is not -1:
        stddevs_sampling = fwhm_dict["images_used"]
        stddevs = fwhm_dict["fwhm_stddevs"]

    if cubic_dict is not -1:
        max_interp = cubic_dict["max"]
        fwhm1 = cubic_dict["fwhm"]
        rvals = cubic_dict["rvals"]
        interp_x = cubic_dict["x_fit"]
        interp_y = cubic_dict["y_fit"]

    plt.figure(1111)

    
    plt.subplot(224)
    plt.scatter(y_average,reduced_subtracted,label=f"Reduced spectra",color='grey',s=10,marker='.')
    plt.plot(moving_average,label=f"Moving average, peak max={max_average}",linestyle='--',color='r')

    if cubic_dict is not -1:
        plt.plot(interp_x,interp_y,label=f"Cubic spline, peak max={max_interp}",linestyle='--',color='b')
        plt.axvspan(rvals[0],rvals[1],label=f"FWHM = {fwhm1}",alpha=0.5)

    plt.legend()

    if np.isin(power_fit,error_codes) == False:
        plt.subplot(222)
        plt.scatter(number_photons,running_diff_std,marker="1",s=15,label=f"Final std = {printable_final_std}")
        plt.plot(number_photons,power_curve_fit,color='red')
        plt.plot([0,number_photons[-1]],[est_double_std,est_double_std],
                label=f"2x std = {printable_est_std}, % improvement={printable_percent_imp})",color='k')
        # plt.plot([cutoff_photons,cutoff_photons],[est_double_std,running_diff_std[0]],lineStyle='--',color='k')
        # Needs reworking probably with the algo
        plt.ylabel("pixel-to-pixle σ")
        plt.xlabel("Total number photons")
        plt.legend()

    if fwhm_dict is not -1:
        plt.subplot(121)
        plt.scatter(stddevs_sampling,stddevs)

    # plt.subplot(212)
    # plt.plot(np.diff(running_diff_std))



    plt.show()


def main():
    args = process_args()
    dataset_obj = xsa.get_reduced_data(args.input_datasetfile)

    plot_interpolations(dataset_obj,args)

if __name__ == "__main__":

    main()