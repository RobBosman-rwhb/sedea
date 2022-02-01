import numpy as np
from numpy.core.numeric import rollaxis
from scipy.optimize import curve_fit

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
            return self.do_std_fitting_for_plots()
        elif type(self.power_fit_dict) == "dict": # Needs to have an fitting failure check
            return self.power_fit_dict


    # Active Functions

    def func_powerlaw(self,x, a, b, c):
        return (a*np.power(x, -b))+c

    def calc_percentage_improvement(self,x_data,asymptote_level):
        return ((x_data[-1]-asymptote_level) / (x_data[0]-x_data[-1])) * 100


    def do_std_fitting_for_plots(self):
        ## Grabs some data that we need
    
        shot_by_shot_matrix = self.get_shot_by_shot()
        images = self.get_total_images()
        
        ## Creates some empty space that would be nice to have
        std_sampling_frequency = 10
        minimal_percentage = 0.05

        std_sample_no = images - (images%std_sampling_frequency)
        rounded_images = images - (images%std_sample_no)
        number_photons = np.zeros((std_sample_no))
        running_pixel_std = np.zeros((std_sample_no))
        image_sampling = np.linspace(std_sampling_frequency,rounded_images,std_sample_no)

        ## Calculate the std as a function of somethin

        for indx,i in enumerate(image_sampling):
            a = int(i)
            reduced_spectra = np.mean(shot_by_shot_matrix[0:a,:],axis=0)
            running_pixel_std[indx] = np.std(reduced_spectra)
            number_photons[indx] = np.sum(shot_by_shot_matrix[0:a,:])

        pars1, cov = curve_fit(f=self.func_powerlaw, xdata=number_photons[:], ydata=running_pixel_std[:],
                                p0=[0, 0, 0], bounds=(-1000, 1000),method='dogbox')

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