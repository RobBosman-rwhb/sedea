import numpy as np
from numpy.core.numeric import rollaxis

class XesDataset:

    def __init__(self,reduced_spectra,bg1_spectra=None,
                    bg2_spectra=None,summed_dataimage=None,
                    ROI=None,reduced_spectra_filtered=None,
                    reduced_avg=None,bckgA_avg=None,bckgB_avg=None,
                    focal_interpolation=None,reduced_subtracted=None,
                    interpolated_peak=None,dataset_name=None):

        self.reduced_spectra = reduced_spectra
        self.bcA_spectra = bg1_spectra
        self.bg2_spectra = bg2_spectra
        self.summed_dataimage = summed_dataimage
        self.ROI = ROI
        self.reduced_spectra_filtered = reduced_spectra_filtered
        self.reduced_avg = reduced_avg
        self.bckgA_avg = bckgA_avg
        self.bckgB_avg = bckgB_avg
        self.focal_interpolation = focal_interpolation
        self.reduced_subtracted = reduced_subtracted
        self.interpolated_peak = interpolated_peak
        self.dataset_name = dataset_name

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
    
    def set_reduced_avg(self,reduced_avg):
        self.reduced_avg=reduced_avg

    def set_bckgA_avg(self,bckgA_avg):
        self.bckgA_avg=bckgA_avg

    def set_bckgB_avg(self,bckgB_avg):
        self.bckg_avg=bckgB_avg

    ## Information variables ##

    def set_dataset_name(self,dataset_name):
        self.dataset_name=dataset_name

    def set_ROI(self,ROI):
        self.ROI=ROI

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

    def get_reduced_avg(self):
        return self.reduced_avg

    def get_bckgA_avg(self):
        return self.bckgA_avg

    def get_bckgB_avg(self):
        return self.bckgB_avg

    ## Information variables ##

    def dataset_name(self):
        return self.dataset_name

    def get_ROI(self):
        return self.ROI

    ## Post processing analysis variables ##

    def get_focal_interpolation(self):
        return self.focal_interpolation

    def get_reduced_subtracted(self):
        return self.reduced_subtracted
    
    def get_interpolated_peak(self):
        return self.interpolated_peak

