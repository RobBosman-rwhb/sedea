from operator import index
import numpy as np
import random
import xes_signal_analysis as xsa
from numpy.core.numeric import rollaxis
from scipy.optimize import curve_fit
from scipy import interpolate
from matplotlib import pyplot as plt


class XesDataset:

    def __init__(self,reduced_spectra,bg1_spectra=None,
                    bg2_spectra=None,summed_dataimage=None,
                    ROI=None,reduced_spectra_filtered=None,
                    shot_by_shot_subtracted=None,
                    focal_interpolation=None,reduced_subtracted=None,
                    interpolated_peak=None,dataset_name=None,
                    total_exposure_time=None,single_exposure_time=None,
                    total_images=None,power_fit_dict=None):

        self.reduced_spectra = reduced_spectra
        self.bg1_spectra = bg1_spectra
        self.bg2_spectra = bg2_spectra
        self.summed_dataimage = summed_dataimage
        self.ROI = ROI
        self.reduced_spectra_filtered = reduced_spectra_filtered
        self.shot_by_shot_subtracted = shot_by_shot_subtracted
        self.focal_interpolation = focal_interpolation
        self.reduced_subtracted = reduced_subtracted
        self.interpolated_peak = interpolated_peak
        self.dataset_name = dataset_name
        self.total_exposure_time = total_exposure_time
        self.single_exposure_time = single_exposure_time
        self.total_images = total_images
        self.power_fit_dict = power_fit_dict

    ## CLASS FUNCTIONS ##


    ## SET FUNCTIONS ##
    def set_reduced_spectra(self,reduced_spectra):
        self.reduced_spectra=reduced_spectra

    def set_bg1_spectra(self,bg1_spectra):
        self.bg1_spectra=bg1_spectra

    def set_bg2_spectra(self,bg2_spectra):
        self.bg2_spectra=bg2_spectra

    def set_summed_dataimage(self,summed_dataimage):
        self.summed_dataimage=summed_dataimage

    def set_shot_by_shot(self,shot_by_shot_subtracted):
        self.shot_by_shot_subtracted=shot_by_shot_subtracted
    

    ## Information variables ##

    def set_dataset_name(self,dataset_name):
        self.dataset_name=dataset_name

    def set_ROI(self,ROI):
        self.ROI=ROI

    def set_total_exposure_time(self,total_exposure_time):
        self.total_exposure_time=total_exposure_time

    def set_single_exposure_time(self,single_exposure_time):
        self.single_exposure_time=single_exposure_time
    
    def set_total_images(self,total_images):
        self.total_images = total_images

    ## Post processing variables ##

    def set_focal_interpolation(self,focal_interpolation):
        self.focal_interpolation=focal_interpolation

    def set_reduced_subtracted(self,reduced_subtracted):
        self.reduced_subtracted=reduced_subtracted

    def set_interpolated_peak(self,interpolated_peak):
        self.interpolated_peak=interpolated_peak



    ## GET FUNCTIONS ##
    def get_reduced_spectra(self):
        return self.reduced_spectra

    def get_reduced_spectra_filtered(self):
        return self.reduced_spectra_filtered

    def get_bg1_spectra(self):
        return self.bg1_spectra

    def get_bg2_spectra(self):
        return self.bg2_spectra

    def get_summed_dataimage(self):
        return self.summed_dataimage

    def get_shot_by_shot(self):
        return self.shot_by_shot_subtracted

    ## Information variables ##

    def get_dataset_name(self):
        return self.dataset_name

    def get_ROI(self):
        return self.ROI

    def get_total_exposure_time(self):
        return self.total_exposure_time

    def get_single_exposure_time(self):
        return self.single_exposure_time

    def get_total_images(self):
        return self.total_images

    ## Post processing analysis variables ##

    def get_focal_interpolation(self):
        return self.focal_interpolation

    def get_reduced_subtracted(self):
        return self.reduced_subtracted
    
    def get_interpolated_peak(self):
        return self.interpolated_peak

    def get_power_fit_dict(self):
        if self.power_fit_dict == None:
            return self.do_std_calc_for_plots()
        elif type(self.power_fit_dict) == "dict": # Needs to have an fitting failure check
            return self.power_fit_dict


    # Functions that should really be public? or in xsa?

    def func_powerlaw(self,x, a, b, c):
        return (a*np.power(x, -b))+c


    def calc_percentage_improvement(self,x_data,asymptote_level):
        return ((x_data[-1]-asymptote_level) / (x_data[0]-x_data[-1])) * 100


    def power_func_fit(self,x,y):
        try:
            print(x,y)
            pars1, cov = curve_fit(f=self.func_powerlaw, xdata=x, ydata=y,
                            p0=[0, 0, 0], bounds=(-1000, 1000),method='dogbox')
        except RuntimeError:
            if len(x) or len(y) <= 10 :
                return -1,-1,-1
            pars1,_,_ = self.power_func_fit(x[:-10],y[:-10])

        return pars1,len(x),len(y)


    def generate_image_averaging(self,sampling_freq,images):
        """ Function to calculate a sampling of images and handle non-divisible inputs."""
        rounded_images = images - (images%sampling_freq)
        sample_no = int(rounded_images/sampling_freq)
        image_sampling = np.linspace(sampling_freq,rounded_images,sample_no)
        return image_sampling,sample_no


    def generate_recombined_spectra(self,sbs_matrix,num):
        index_numbers = random.sample(range(0,len(sbs_matrix[:,0])),int(num))
        return np.mean(sbs_matrix[index_numbers,:],axis=0)

    # Private functions that must not be touched?


    def set_sampling_frequency(self):
        images = self.total_images
        if images < 1 and images > 10:
            std_sampling_frequency = images
        if images > 10:
            std_sampling_frequency = round(images/100,-1)
        return std_sampling_frequency


    def fit_cubic_spline(self,spectra):
        """Function to fit cublic spline without setting the s variable
        does a grid scan over the variable. Principally used with FWHM std calc"""
        # print(np.linspace(1,len(spectra),int(len(spectra)*10)))
        a = np.logspace(0,3,30)-1
        for i in a:
            try:
                interp_x,interp_y=xsa.calculate_cubic_interpolation(spectra,i)
                linear_interp = interpolate.UnivariateSpline(interp_x,interp_y-round(np.max(interp_y),4)/2,s=0)
                root1,root2 = linear_interp.roots()
                print("Assessed cubic spline and fitted roots")
                cubic_dict = { "x_fit":interp_x,
                                "y_fit":interp_y,
                                "rvals":[root1,root2]}
                return cubic_dict
            except ValueError:
                continue
        print("Found no cubic spline")
        return -1


    def do_moving_average_for_plots(self):
        spectra = self.reduced_subtracted
        spectral_length = len(spectra)
        moving_average = xsa.calculate_moving_average(spectra,10)
        max_average = round(np.max(moving_average),4)
        y = np.linspace(0, spectral_length, num=spectral_length, endpoint=True)
        return moving_average,max_average,y


    def do_fwhm_assesment_for_plots(self):
        ## Calculate the std as a function of something
        shot_by_shot = self.shot_by_shot_subtracted
        images = self.total_images

        if images < 100:
            return -1

        repeats = 5
        sample_freq = self.set_sampling_frequency()
        sampling_array,samples = self.generate_image_averaging(sample_freq,images)
        # sampling_array = np.insert(sampling_array,0,1,axis=0)
        # samples = samples +1
        fwhm_stddevs = np.zeros((samples))

        for indx,i in enumerate(sampling_array):
            recombined_spectra = self.generate_recombined_spectra(shot_by_shot,i)
            fittings = [ self.fit_cubic_spline(recombined_spectra) for _ in range(repeats) ]

            if np.isin(-1,fittings):
                break

            roots = [ i["rvals"] for i in fittings ]
            fwhms = [ roots[d][1] - roots[d][0] for d in range(repeats) ]
            fwhm_stddevs[indx] = np.mean(fwhms)
        fwhm_assessment = {"images_used":sampling_array,
                            "fwhm_stddevs":fwhm_stddevs}
        return fwhm_assessment


    def do_cubic_interpolation_for_plots(self):
        spectra = self.reduced_subtracted
        cubic_dict = self.fit_cubic_spline(spectra)

        if cubic_dict is -1:
            return cubic_dict

        y_fit = cubic_dict["y_fit"]
        roots = cubic_dict["rvals"]
        max_interp = round(np.max(y_fit),4)
        fwhm = roots[1]-roots[0]
        cubic_dict["fwhm"] = fwhm
        cubic_dict["max"] = max_interp

        return cubic_dict


    def do_std_calc_for_plots(self):
        ## Grabs some data that we need    
        shot_by_shot_matrix = self.shot_by_shot_subtracted
        images = self.total_images

        if images == 1:
            return -2

        std_sampling_frequency = self.set_sampling_frequency()

        minimal_percentage = 0.05
        print(std_sampling_frequency)
        image_sampling_vector,samples = self.generate_image_averaging(std_sampling_frequency,images)
        print(image_sampling_vector,samples)
        number_photons = np.zeros((samples))
        running_pixel_std = np.zeros((samples))
        for indx,i in enumerate(image_sampling_vector):
            a = int(i)
            reduced_spectra = np.mean(shot_by_shot_matrix[0:a,:],axis=0)
            running_pixel_std[indx] = np.std(reduced_spectra)
            number_photons[indx] = np.sum(shot_by_shot_matrix[0:a,:])

        pars1, _, _ = self.power_func_fit(number_photons,running_pixel_std)
        
        if pars1[0] == -1:
            return -1

        fitted_std_curve = [ self.func_powerlaw(i,pars1[0],pars1[1],pars1[2]) for i in number_photons ]

        # back calculate the improvement in 
        estimated_double_stdsig = self.func_powerlaw(number_photons[-1]*2,pars1[0],pars1[1],pars1[2])
        percentage_improvement = self.calc_percentage_improvement(running_pixel_std,estimated_double_stdsig)
        cutoff_std = (minimal_percentage*running_pixel_std[0])+estimated_double_stdsig/(1+minimal_percentage)
        cutoff_photons = ((cutoff_std/pars1[0])**(1/pars1[1]))
    
        power_fit_dict = {"opt_params":pars1,
                            "running_std":running_pixel_std,
                            "tot_photons_array":number_photons,
                            "fitted_std_curve":fitted_std_curve,
                            "est_double_std":estimated_double_stdsig,
                            "expected_improv_percent":percentage_improvement,
                            "cutoff_std":cutoff_std,
                            "cutoff_photon":cutoff_photons}

        self.power_fit_dict = power_fit_dict
        return self.power_fit_dict