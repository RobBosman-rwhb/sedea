import numpy as np
from numpy.core.numeric import rollaxis

class XesDataset:

    def __init__(self,reduced_spectra,bckgA_spectra=None,
                    bckgB_spectra=None,summed_dataimage=None,
                    ROI=None,reduced_spectra_filtered=None,
                    reduced_avg=None,bckgA_avg=None,bckgB_avg=None,
                    focal_interpolation=None,reduced_subtracted=None,
                    interpolated_peak=None,dataset_name=None):

        self.reduced_spectra = reduced_spectra
        self.bckgA_spectra = bckgA_spectra
        self.bckgB_spectra = bckgB_spectra
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

    # Set functions
    def set_reduced_spectra(self,reduced_spectra):
        self.reduced_spectra=reduced_spectra

    def set_bckgA_spectra(self,bckgA_spectra):
        self.bckgA_spectra=bckgA_spectra

    def set_bckgB_spectra(self,bckgB_spectra):
        self.bckgB_spectra=bckgB_spectra

    def set_summed_dataimage(self,summed_dataimage):
        self.summed_dataimage=summed_dataimage

    def set_ROI(self,ROI):
        self.ROI=ROI

    def set_reduced_spectra_filtered(self,reduced_spectra_filtered):
        self.reduced_spectra_filtered=reduced_spectra_filtered
    
    def set_reduced_avg(self,reduced_avg):
        self.reduced_avg=reduced_avg

    def set_bckgA_avg(self,bckgA_avg):
        self.bckgA_avg=bckgA_avg

    def set_bckgB_avg(self,bckgB_avg):
        self.bckg_avg=bckgB_avg

    def set_focal_interpolation(self,focal_interpolation):
        self.focal_interpolation=focal_interpolation

    def set_reduced_subtracted(self,reduced_subtracted):
        self.reduced_subtracted=reduced_subtracted

    def set_interpolated_peak(self,interpolated_peak):
        self.interpolated_peak=interpolated_peak

    def set_dataset_name(self,dataset_name):
        self.dataset_name=dataset_name

    # Get functions
    def get_reduced_spectra(self):
        return self.reduced_spectra

    def get_bckgA_spectra(self):
        return self.bckgA_spectra

    def get_bckgB_spectra(self):
        return self.bckgB_spectra

    def get_summed_dataimage(self):
        return self.summed_dataimage

    def get_ROI(self):
        return self.ROI

    def get_reduced_spectra_filtered(self):
        return self.reduced_spectra_filtered

    def get_reduced_avg(self):
        return self.reduced_avg

    def get_bckgA_avg(self):
        return self.bckgA_avg

    def get_bckgB_avg(self):
        return self.bckgB_avg

    def get_focal_interpolation(self):
        return self.focal_interpolation

    def get_reduced_subtracted(self):
        return self.reduced_subtracted
    
    def get_interpolated_peak(self):
        return self.interpolated_peak

    def dataset_name(self):
        return self.dataset_name