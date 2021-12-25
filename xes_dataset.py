import numpy as np
from numpy.core.numeric import rollaxis

class XesDataset:

    def __init__(self,reduced_spectra,bg1_spectra=None,
                    bg2_spectra=None,summed_dataimage=None,
                    ROI=None,reduced_spectra_filtered=None,
                    shot_by_shot_subtracted=None,
                    focal_interpolation=None,reduced_subtracted=None,
                    interpolated_peak=None,dataset_name=None,
                    total_exposure_time=None,single_exposure_time=None):

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

    ## CLASS FUNCTIONS ##


    ## SET FUNCTIONS ##
    def set_reduced_spectra(self,reduced_spectra):
        self.reduced_spectra=reduced_spectra

    def get_reduced_spectra_filtered(self):
        return self.reduced_spectra_filtered

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

    ## Post processing analysis variables ##

    def get_focal_interpolation(self):
        return self.focal_interpolation

    def get_reduced_subtracted(self):
        return self.reduced_subtracted
    
    def get_interpolated_peak(self):
        return self.interpolated_peak