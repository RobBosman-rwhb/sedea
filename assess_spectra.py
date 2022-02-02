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
    
    if type(args.smoothing) is not float:
        print("smoothing value should be a number, a cuddly number")

    return args
    

def plot_interpolations(dataset_obj,args):

    reduced_subtracted = dataset_obj.get_reduced_subtracted()
    power_fit = dataset_obj.get_power_fit_dict()
    running_diff_std = power_fit["running_std"]
    number_photons = power_fit["tot_photons_array"]
    est_double_std = power_fit["est_double_std"]
    percentage_improvement = power_fit["expected_improv_percent"]
    power_curve_fit = power_fit["fitted_std_curve"]
    cutoff_photons = power_fit["cutoff_photon"]

    printable_final_std = round(running_diff_std[-1],8)
    printable_est_std = round(est_double_std,8)
    printable_percent_imp = round(percentage_improvement,2)

    moving_average,interp_x,interp_y=xsa.calculate_interpolation(reduced_subtracted,args.average_window,args.smoothing)
    max_interp = round(np.max(interp_y),4)
    max_average = round(np.max(moving_average),4)

    linear_interp = interpolate.UnivariateSpline(interp_x,interp_y-max_interp/2,s=0)
    
    r1,r2 = linear_interp.roots()
    fwhm1 = r2-r1

    plt.figure(1111)

    plt.subplot(211)
    plt.plot(reduced_subtracted,label=f"Reduced spectra",color='k',linewidth=0.1)
    plt.plot(moving_average,label=f"Moving average, peak max={max_average}",linestyle='--',color='r')
    plt.plot(interp_x,interp_y,label=f"Cubic spline, peak max={max_interp}",linestyle='--',color='b')
    plt.axvspan(r1,r2,label=f"FWHM = {fwhm1}",alpha=0.5)
    plt.legend()

    plt.subplot(212)
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



    plt.show()


def main():
    args = process_args()
    dataset_obj = xsa.get_reduced_data(args.input_datasetfile)

    plot_interpolations(dataset_obj,args)

if __name__ == "__main__":

    main()