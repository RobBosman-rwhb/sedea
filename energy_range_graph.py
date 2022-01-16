from typing import Callable
from analyser_crystal import analyser_crystal
import numpy as np 
import math
import analyser_array_position as an
import scipy.constants as cons
import matplotlib.pyplot as plt
import mpmath as mpm

def signal_above_midpoint(bragg,bending_radius):
    return math.tan(np.radians(90)-np.radians(bragg))*bending_radius

def calculate_bragg_from_height(height,bending_radius):
    return 90-math.degrees(math.atan(height/bending_radius))

def E_Calc (bragg_angle,dspace):
    h = cons.h/cons.electron_volt          # calculations from RT4XES
    energy_factor = h*cons.c*(1e10) 
    bragg_lambda = 2*dspace*math.sin(np.radians(bragg_angle))
    energy = energy_factor/bragg_lambda
    return energy

def Calc_BraggFromE(E,dspace):
    h = cons.h/cons.electron_volt          # calculations from RT4XES
    energy_factor = h*cons.c*(1e10)
    bragg_lambda = energy_factor/E
    bragg_angle = np.degrees(math.asin(bragg_lambda/(2*dspace)))
    return bragg_angle

def dspace_Calc (material,hkl):
    if material == 'Si':
        a = 5.4309
    elif material == 'Ge':
        a = 5.6574
    else:
        print ("Material not recognised returning Si")
        a = 5.4309
    denomintor_term = math.sqrt( pow(hkl[0],2) +pow(hkl[1],2) +pow(hkl[2],2) )
    return a/denomintor_term

def crystal_E_range(analyser_height,bending_radius,bragg_angle,analyser_crystal):
    # function that given analyser crystal height dimension, angle range and indicies calculates max/min energy
    crystal_centre = an.calculate_midpoint(bragg_angle,bending_radius)
    top_crystal_vertice = crystal_centre+(analyser_height/2)
    bottom_crystal_vertice = crystal_centre-(analyser_height/2)

    bragg_top = np.degrees(float(mpm.acot(top_crystal_vertice/bending_radius)))
    bragg_bottom = np.degrees(float(mpm.acot(bottom_crystal_vertice/bending_radius)))
    E_top = E_Calc(bragg_top,dspace_Calc(analyser_crystal[0],analyser_crystal[1]))
    E_bottom = E_Calc(bragg_bottom,dspace_Calc(analyser_crystal[0],analyser_crystal[1]))
    E_range = E_top-E_bottom
    return E_top, E_bottom, E_range

def crystal_E_range_over_braggRange(bragg_angle_range,Rs,crystal_height,analyser_crystal):
    # function that given Analyser crystal height dimension, angle_range and indices calculates max/min energy over a bragg angle range
    analyser_Si440_dspace = dspace_Calc(analyser_crystal[0],analyser_crystal[1])
    # highest_energy,_,_ = crystal_E_range(crystal_height,Rs,bragg_angle_range[0],analyser_crystal)
    # lowest_energy,_,_ = crystal_E_range(crystal_height,Rs,bragg_angle_range[1],analyser_crystal)
    higher_E_midpoint = an.calculate_midpoint(bragg_angle_range[0],Rs)
    lower_E_midpoint = an.calculate_midpoint(bragg_angle_range[1],Rs)
    higher_E_min = higher_E_midpoint+(crystal_height/2)
    lower_E_max = lower_E_midpoint-(crystal_height/2)
    highest_bragg = calculate_bragg_from_height(lower_E_max,Rs) 
    lowest_bragg = calculate_bragg_from_height(higher_E_min,Rs)
    lowest_energy = E_Calc(highest_bragg,analyser_Si440_dspace)
    highest_energy = E_Calc(lowest_bragg,analyser_Si440_dspace)
    return [lowest_energy,highest_energy]


def calc_crystal_signal_range_above_midpoint(energy1,energy2,bending_radius,analyser_crystal):
    bragg1 = Calc_BraggFromE(energy1,dspace_Calc(analyser_crystal[0],analyser_crystal[1]))
    bragg2 = Calc_BraggFromE(energy2,dspace_Calc(analyser_crystal[0],analyser_crystal[1]))
    height1 = signal_above_midpoint(bragg1,bending_radius)
    height2 = signal_above_midpoint(bragg2,bending_radius)
    signal_range_above_midpoint = height1-height2
    return height1,height2,signal_range_above_midpoint,bragg1,bragg2,(bragg1+bragg2)/2

def for_plot_bragg_Erange(bragg_range,analyser_height,bending_radius,analyser_crystal):
        all_Etop = []
        all_Ebottom = []
        for i in bragg_range:
            E_top,E_bottom,E_range = crystal_E_range(analyser_height,bending_radius,i,analyser_crystal)
            all_Etop.append(E_top)
            all_Ebottom.append(E_bottom)
        return all_Etop,all_Ebottom,E_range

def plot_hub_spectrometer_operating_energies():
    bragg_range = [70,90]
    angle_ranges = np.linspace(bragg_range[0],bragg_range[1],1000)
    Ge620_Erange = [ E_Calc(i,dspace_Calc('Ge',[6,2,0])) for i in np.nditer(angle_ranges)]
    Ge440_Erange = [ E_Calc(i,dspace_Calc('Ge',[4,4,0])) for i in np.nditer(angle_ranges)]
    # Ge555_Erange = [ E_Calc(i,dspace_Calc('Ge',[5,5,5])) for i in np.nditer(angle_ranges)]
    Si440_Erange = [ E_Calc(i,dspace_Calc('Si',[4,4,0])) for i in np.nditer(angle_ranges)]
    Si444_Erange = [ E_Calc(i,dspace_Calc('Si',[4,4,4])) for i in np.nditer(angle_ranges)]
    Si553_Erange = [ E_Calc(i,dspace_Calc('Si',[5,5,3])) for i in np.nditer(angle_ranges)]
    # Si642_Erange = [ E_Calc(i,dspace_Calc('Si',[6,4,2])) for i in np.nditer(angle_ranges)]

    plt.plot(Ge620_Erange,angle_ranges, label='Ge620')
    plt.plot(Ge440_Erange,angle_ranges, label='Ge440')
    # plt.plot(Ge555_Erange,angle_ranges, label='Ge555')
    plt.plot(Si440_Erange,angle_ranges, label='Si440')
    plt.plot(Si444_Erange,angle_ranges, label='Si444')
    plt.plot(Si553_Erange,angle_ranges, label='Si553')
    # plt.plot(Si642_Erange,angle_ranges, label='Si642')

    plt.axhspan(75,85,alpha=0.2,color='green')
    plt.plot(6485,Calc_BraggFromE(6485,dspace_Calc('Si',[4,4,0])),'ro')
    plt.plot(7058,Calc_BraggFromE(7058,dspace_Calc('Ge',[6,2,0])),'ro')
    plt.plot(8905,Calc_BraggFromE(8905,dspace_Calc('Si',[5,5,3])),'ro')
    plt.plot(6404,Calc_BraggFromE(6404,dspace_Calc('Ge',[4,4,0])),'go')
    plt.plot(6391,Calc_BraggFromE(6391,dspace_Calc('Ge',[4,4,0])),'go')
    plt.plot(8028,Calc_BraggFromE(8028,dspace_Calc('Si',[4,4,4])),'go')
    plt.plot(8048,Calc_BraggFromE(8048,dspace_Calc('Si',[4,4,4])),'go')

    plt.text(6485,Calc_BraggFromE(6485,dspace_Calc('Si',[4,4,0])),'MnKß1,3',horizontalalignment='right')
    plt.text(7058,Calc_BraggFromE(7058,dspace_Calc('Ge',[6,2,0])),'FeKß1,3',horizontalalignment='right')
    plt.text(8905,Calc_BraggFromE(8905,dspace_Calc('Si',[5,5,3])),'CuKß1,3',horizontalalignment='right')
    plt.text(6404,Calc_BraggFromE(6404,dspace_Calc('Ge',[4,4,0])),'FeKα1',horizontalalignment='right')
    plt.text(6391,Calc_BraggFromE(6391,dspace_Calc('Ge',[4,4,0])),'FeKα2',horizontalalignment='right')
    plt.text(8028,Calc_BraggFromE(8028,dspace_Calc('Si',[4,4,4])),'CuKα1',horizontalalignment='right')
    plt.text(8048,Calc_BraggFromE(8048,dspace_Calc('Si',[4,4,4])),'CuKα2',horizontalalignment='right')

    plt.legend()
    plt.xlim([6000,10000])
    plt.ylim(bragg_range)
    plt.xlabel("Energy range (eV)")
    plt.ylabel("Bragg Angle $^o$")
    plt.show()

def plot_crystal_E_range_in_bragg():
    analyser_height = 25
    bending_radius = 400

    bragg_angles = np.linspace(75,85,21)
    analyser_crystals = [["Ge",[4,4,0]],["Ge",[6,2,0]],["Si",[4,4,0]],["Si",[4,4,4]],["Si",[5,5,3]],["Si",[6,2,0]],["Si",[5,5,1]]]
    Analyser_crystalGe440 = for_plot_bragg_Erange(bragg_angles,analyser_height,bending_radius,analyser_crystals[0])
    Analyser_crystalGe620 = for_plot_bragg_Erange(bragg_angles,analyser_height,bending_radius,analyser_crystals[1])
    Analyser_crystalSi440 = for_plot_bragg_Erange(bragg_angles,analyser_height,bending_radius,analyser_crystals[2])
    Analyser_crystalSi444 = for_plot_bragg_Erange(bragg_angles,analyser_height,bending_radius,analyser_crystals[3])
    Analyser_crystalSi553 = for_plot_bragg_Erange(bragg_angles,analyser_height,bending_radius,analyser_crystals[4])
    Analyser_crystalSi620 = for_plot_bragg_Erange(bragg_angles,analyser_height,bending_radius,analyser_crystals[5])
    Analyser_crystalSi551 = for_plot_bragg_Erange(bragg_angles,analyser_height,bending_radius,analyser_crystals[6])

    Fekalpha12 = [6385,6403,6410]                           # https://sci-hub.se/https://doi.org/10.1021/acs.inorgchem.0c01620
    Fekbeta13 = [7035,7058,7070]                            # https://sci-hub.se/https://doi.org/10.1021/acs.inorgchem.0c01620
    MnKBeta13 = [6476,6492,6500]                            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3422323/
    MnKa1 = [5894,5899,5906]                                # https://pubs.acs.org/doi/10.1021/acs.inorgchem.8b01674
    CuKa12 = [8015,8048,8060]                               # K alpha and K beta x-ray emission spectra of copper
    CuKBeta13 = [8890,8905,8920]                            # K alpha and K beta x-ray emission spectra of copper
    NiKa12 = [7450,7478,7485]                               # High-resolution studies of the K emission spectra of nickel
    NiKBeta13 = [8250,8263,8275]                            # Sorum, H., & Bremer, J. (1982). High-resolution studies of the K emission spectra of nickel. Journal of Physics F: Metal Physics, 12(11), 2721–2728. doi:10.1088/0305-4608/12/11/029 

    print(Analyser_crystalGe440)
    
    ## Ploting for Feka 
    plt.subplot(241)
    plt.title("FeKα (Ge440)")
    plt.plot(bragg_angles,Analyser_crystalGe440[0],color='k')
    plt.plot(bragg_angles,Analyser_crystalGe440[1],color='k')
    plt.plot([75.5,75.5],[Fekalpha12[0],Fekalpha12[2]],color='g')
    plt.text(75.7,6383,"Δ 25eV",horizontalalignment='left')
    plt.fill_between(bragg_angles,Analyser_crystalGe440[0],Analyser_crystalGe440[1],alpha=0.3,color='blue')
    plt.xlim([75,85])
    plt.ylim([6200,6500])
    
    ## Ploting for Fekbeta13
    plt.subplot(242)
    plt.title("FeKß1,3 (Ge620)")
    plt.plot(bragg_angles,Analyser_crystalGe620[0],'k')
    plt.plot(bragg_angles,Analyser_crystalGe620[1],'k')
    plt.plot([79.1,79.1],[Fekbeta13[0],Fekbeta13[2]],color='g')
    plt.text(79.3,7037,"Δ 35eV",horizontalalignment='left')
    plt.fill_between(bragg_angles,Analyser_crystalGe620[0],Analyser_crystalGe620[1],alpha=0.3,color='blue')
    plt.xlim([75,85])

    ## Ploting for MnKbeta13
    plt.subplot(243)
    plt.title("MnKß1,3 (Si440)")
    plt.plot(bragg_angles,Analyser_crystalSi440[0],'k')
    plt.plot(bragg_angles,Analyser_crystalSi440[1],'k')
    plt.plot([83.9,83.9],[MnKBeta13[0],MnKBeta13[2]],color='g')
    plt.text(83.7,6500,"Δ 25eV",horizontalalignment='right')
    plt.fill_between(bragg_angles,Analyser_crystalSi440[0],Analyser_crystalSi440[1],alpha=0.3,color='blue')
    plt.xlim([75,85])
    plt.ylim([6450,6750])

    ## Ploting for CuKa12
    plt.subplot(244)
    plt.title("CuKα (Si444)")
    plt.plot(bragg_angles,Analyser_crystalSi444[0],'k')
    plt.plot(bragg_angles,Analyser_crystalSi444[1],'k')
    plt.plot([80.1,80.1],[CuKa12[0],CuKa12[2]],color='g')
    plt.text(80.3,8011,"Δ 45eV",horizontalalignment='left')
    plt.fill_between(bragg_angles,Analyser_crystalSi444[0],Analyser_crystalSi444[1],alpha=0.3,color='blue')
    plt.xlim([75,85])
    plt.ylim([7900,8300])

    ## Ploting for CuKBeta13
    plt.subplot(245)
    plt.title("CuKß1,3 (Si553)")
    plt.plot(bragg_angles,Analyser_crystalSi553[0],'k')
    plt.plot(bragg_angles,Analyser_crystalSi553[1],'k')
    plt.plot([80.0,80.0],[CuKBeta13[0],CuKBeta13[2]],color='g')
    plt.text(80.2,8880,"Δ 30eV",horizontalalignment='left')
    plt.fill_between(bragg_angles,Analyser_crystalSi553[0],Analyser_crystalSi553[1],alpha=0.3,color='blue')
    plt.xlim([75,85])
    plt.ylim([8750,9200])

    ## Plotting for NiKa12
    plt.subplot(246)
    plt.title("NiKα (Si620)")
    plt.plot(bragg_angles,Analyser_crystalSi620[0],'k')
    plt.plot(bragg_angles,Analyser_crystalSi620[1],'k')
    plt.plot([75.7,75.7],[NiKa12[0],NiKa12[2]],color='g')
    plt.text(75.5,7430,"Δ 35eV",horizontalalignment='left')
    plt.fill_between(bragg_angles,Analyser_crystalSi620[0],Analyser_crystalSi620[1],alpha=0.3,color='blue')
    plt.xlim([75,85])
    plt.ylim([7200,7600])

    ## Plotting for NiKBeta13
    plt.subplot(247)
    plt.title("NiKBeta13 (Si551)")
    plt.plot(bragg_angles,Analyser_crystalSi551[0],'k')
    plt.plot(bragg_angles,Analyser_crystalSi551[1],'k')
    plt.plot([80.0,80.0],[NiKBeta13[0],NiKBeta13[2]],color='g')
    plt.text(80.2,8250,"Δ 25eV",horizontalalignment='left')
    plt.fill_between(bragg_angles,Analyser_crystalSi551[0],Analyser_crystalSi551[1],alpha=0.3,color='blue')
    plt.xlim([75,85]) 
    plt.ylim([8150,8550])

    plt.subplot(248)
    plt.plot(bragg_angles,Analyser_crystalGe440[0],color='b')
    plt.plot(bragg_angles,Analyser_crystalGe440[1],color='b')   
    plt.fill_between(bragg_angles,Analyser_crystalGe440[0],Analyser_crystalGe440[1],alpha=0.3,color='blue')
    plt.text(83,6320,"Ge440",horizontalalignment='left')

    plt.plot(bragg_angles,Analyser_crystalGe620[0],color='g')
    plt.plot(bragg_angles,Analyser_crystalGe620[1],color='g')
    plt.fill_between(bragg_angles,Analyser_crystalGe620[0],Analyser_crystalGe620[1],alpha=0.3,color='green')
    plt.text(83,7050,"Ge620",horizontalalignment='left') 

    plt.plot(bragg_angles,Analyser_crystalSi440[0],color='r')
    plt.plot(bragg_angles,Analyser_crystalSi440[1],color='r')  
    plt.fill_between(bragg_angles,Analyser_crystalSi440[0],Analyser_crystalSi440[1],alpha=0.3,color='red')
    plt.text(83,6620,"Si440",horizontalalignment='left') 

    plt.plot(bragg_angles,Analyser_crystalSi444[0],color='c')
    plt.plot(bragg_angles,Analyser_crystalSi444[1],color='c') 
    plt.fill_between(bragg_angles,Analyser_crystalSi444[0],Analyser_crystalSi444[1],alpha=0.3,color='cyan')
    plt.text(83,8050,"Si444",horizontalalignment='left') 

    plt.plot(bragg_angles,Analyser_crystalSi553[0],color='m')
    plt.plot(bragg_angles,Analyser_crystalSi553[1],color='m')
    plt.fill_between(bragg_angles,Analyser_crystalSi553[0],Analyser_crystalSi553[1],alpha=0.3,color='purple')
    plt.text(83,8940,"Si553",horizontalalignment='left') 

    plt.plot(bragg_angles,Analyser_crystalSi620[0],color='y')
    plt.plot(bragg_angles,Analyser_crystalSi620[1],color='y')
    plt.fill_between(bragg_angles,Analyser_crystalSi620[0],Analyser_crystalSi620[1],alpha=0.3,color='yellow')
    plt.text(83,7400,"Si620",horizontalalignment='left') 

    plt.plot(bragg_angles,Analyser_crystalSi551[0],color='k')
    plt.plot(bragg_angles,Analyser_crystalSi551[1],color='k')
    plt.fill_between(bragg_angles,Analyser_crystalSi551[0],Analyser_crystalSi551[1],alpha=0.3,color='k')
    plt.text(83,8300,"Si551",horizontalalignment='left') 

    plt.ylim([6200,9200])
    plt.xlim([75,85])
    plt.show()
 
def plot_detector_height():
    Rs = 400
    half_crystal = 12.5

    Fekalpha12 = [6385,6403,6410]                           # https://sci-hub.se/https://doi.org/10.1021/acs.inorgchem.0c01620
    Fekbeta13 = [7035,7058,7070]                            # https://sci-hub.se/https://doi.org/10.1021/acs.inorgchem.0c01620
    Fekbeta25 = [7100,00,7115]
    MnKBeta13 = [6475,6492,6500]                            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3422323/
    MnKa1 = [5894,5899,5906]                                # https://pubs.acs.org/doi/10.1021/acs.inorgchem.8b01674
    CuKa12 = [8015,8048,8060]                               # K alpha and K beta x-ray emission spectra of copper
    CuKBeta13 = [8890,8905,8920]                            # K alpha and K beta x-ray emission spectra of copper
    NiKa12 = [7450,7478,7485]                               # High-resolution studies of the K emission spectra of nickel
    NiKBeta13 = [8250,8263,8275]                            # Sorum, H., & Bremer, J. (1982). High-resolution studies of the K emission spectra of nickel. Journal of Physics F: Metal Physics, 12(11), 2721–2728. doi:10.1088/0305-4608/12/11/029 
    # ZnKa12 = [7450 7478 7490]                               # Kα1,2 and Kβ1,3 x-ray emission lines of the 3d transition metals
    # ZnKBeta13 = [8240 8265 8280]                            # Kα1,2 and Kβ1,3 x-ray emission lines of the 3d transition metals

    Fekalpha12_Erange = calc_crystal_signal_range_above_midpoint(Fekalpha12[2],Fekalpha12[0],Rs,["Ge",[4,4,0]])
    Fekbeta13_Erange = calc_crystal_signal_range_above_midpoint(Fekbeta13[2],Fekbeta13[0],Rs,["Ge",[6,2,0]])
    Fekbeta25_Erange = calc_crystal_signal_range_above_midpoint(Fekbeta25[2],Fekbeta25[0],Rs,["Ge",[6,2,0]])
    MnKBeta13_Erange = calc_crystal_signal_range_above_midpoint(MnKBeta13[2],MnKBeta13[0],Rs,["Si",[4,4,0]])
    CuKa12_Erange = calc_crystal_signal_range_above_midpoint(CuKa12[2],CuKa12[0],Rs,["Si",[4,4,4]])
    CuKBeta13_Erange = calc_crystal_signal_range_above_midpoint(CuKBeta13[2],CuKBeta13[0],Rs,["Si",[5,5,3]])
    NiKa12_Erange = calc_crystal_signal_range_above_midpoint(NiKa12[2],NiKa12[0],Rs,["Si",[6,2,0]])
    NiKBeta13_Erange = calc_crystal_signal_range_above_midpoint(NiKBeta13[2],NiKBeta13[0],Rs,["Si",[5,5,1]])

    print(Fekalpha12_Erange)

    Fekalpha12_mid = an.calculate_midpoint(Fekalpha12_Erange[3],Rs)
    Fekbeta13_mid = an.calculate_midpoint(Fekbeta13_Erange[3],Rs)
    Fekbeta25_mid = an.calculate_midpoint(Fekbeta25_Erange[3],Rs)
    MnKBeta13_mid = an.calculate_midpoint(MnKBeta13_Erange[3],Rs)
    CuKa12_mid = an.calculate_midpoint(CuKa12_Erange[3],Rs)
    CuKBeta13_mid = an.calculate_midpoint(CuKBeta13_Erange[3],Rs)
    NiKa12_mid = an.calculate_midpoint(NiKa12_Erange[3],Rs)
    NiKBeta13_mid = an.calculate_midpoint(NiKBeta13_Erange[3],Rs)
    # an.calculate_midpoint(Fekalpha12_Erange[3],["Ge",[4,4,0]])


    Fekbeta13_peak = an.calculate_midpoint(Calc_BraggFromE(Fekbeta13[1],dspace_Calc("Ge",[6,2,0])),Rs)*2
    Fekbeta13_perfect_array_height = an.calculate_midpoint(Calc_BraggFromE(Fekbeta13[1],dspace_Calc("Ge",[6,2,0])),Rs)
    print(f"The Fekbeta1,3 peak should be {Fekbeta13_peak}mm from the sample point")
    print(f"The Fekbeta1,3 midpoint {Fekbeta13_perfect_array_height}mm from the sample point")

    variable_names = ["$FeKα$","$FeKβ_{1,3}$","$FeKβ_{2,5}$","$MnKβ_{1,3}$","$CuKα$","$CuKβ_{1,3}$","$NiKα$","$NiKβ_{1,3}$"]
    # plt.axhspan(MnKBeta13_Erange[1]+MnKBeta13_mid,MnKBeta13_Erange[1]+MnKBeta13_mid+75,alpha=0.2,color='green')
    plt.axhspan(MnKBeta13_Erange[1]+MnKBeta13_mid,MnKBeta13_Erange[1]+MnKBeta13_mid+150,alpha=0.2,color='blue')
    plt.plot([-1,9],[MnKBeta13_Erange[1]+MnKBeta13_mid+75,MnKBeta13_Erange[1]+MnKBeta13_mid+75],color='red',lineWidth=.3)
    # plt.axhspan(CuKa12_Erange[1]+CuKa12_mid,CuKa12_Erange[1]+CuKa12_mid+75,alpha=0.2,color='red')
    plt.plot([0,0],[Fekalpha12_Erange[0]+Fekalpha12_mid+half_crystal,Fekalpha12_Erange[1]+Fekalpha12_mid-half_crystal])
    plt.plot([1,1],[Fekbeta13_Erange[0]+Fekbeta13_mid+half_crystal,Fekbeta13_Erange[1]+Fekbeta13_mid-half_crystal])
    plt.plot(1,Fekbeta13_peak,marker='o')
    plt.plot([2,2],[Fekbeta25_Erange[0]+Fekbeta25_mid+half_crystal,Fekbeta25_Erange[1]+Fekbeta25_mid-half_crystal])
    plt.plot([3,3],[MnKBeta13_Erange[0]+MnKBeta13_mid+half_crystal,MnKBeta13_Erange[1]+MnKBeta13_mid-half_crystal])
    plt.plot([4,4],[CuKa12_Erange[0]+CuKa12_mid+half_crystal,CuKa12_Erange[1]+CuKa12_mid-half_crystal])
    plt.plot([5,5],[CuKBeta13_Erange[0]+CuKBeta13_mid+half_crystal,CuKBeta13_Erange[1]+CuKBeta13_mid-half_crystal])
    plt.plot([6,6],[NiKa12_Erange[0]+NiKa12_mid+half_crystal,NiKa12_Erange[1]+NiKa12_mid-half_crystal])
    plt.plot([7,7],[NiKBeta13_Erange[0]+NiKBeta13_mid+half_crystal,NiKBeta13_Erange[1]+NiKBeta13_mid-half_crystal])

    plt.xticks(np.arange(7),variable_names)
    plt.xlabel("Spectra")
    plt.xlim([-1,8])
    plt.ylabel("Distance from sampe point (mm)")
    plt.show()

def energy_references():
    energy_dict = {"Cu":{"Ka":[8018,8068],"Kb13":[8880,8930],"Kb25":[8952,9003]},\
                 "Fe":{"Ka":[6374,6424],"Kb13":[7033,7083],"Kb25":[7085,7135]},\
                 "Mn":{"Ka":[5869,5919],"Kb13":[6465,6515],"Kb25":[6512,6563]},\
                 "Co":{"Ka":[6900,6950],"Kb13":[7624,7674]},\
                 "Ni":{"Ka":[7448,7498],"Kb13":[8240,8290]},\
                 "Zn":{"Ka":[8609,8659],"Kb13":[8952,9003]},\
                 "Mo":{"Ka":[17325,17550],"Kb13":[19550,19650]},\
                 "Rd":{"Ka":[20074,20216],"Kb13":[22704,22744]}}
    return energy_dict

if __name__ == '__main__':


    plot_detector_height()


    Rs = 500
    MnKbeta13 = [6475,6500]
    crystal_height = 25
    Rs = 500
    bragg_angle_range = [75,85]
    crystal_height = 25

    print(Calc_BraggFromE(8977,dspace_Calc("Si",[5,5,3])))
    # print(calc_crystal_signal_range_above_midpoint(7110.59,7110.61,500,["Si",[5,3,1]])[0]+an.calculate_midpoint(71.75,500))

    

    top_position = an.calculate_midpoint(85,Rs)
    bottom_position = an.calculate_midpoint(75,Rs)

    top_bound = top_position-((25+10)*1.5)
    bottom_bound = bottom_position+((25+10)*1.5)

    print(an.calculate_midpoint(77.9,500))

    # print(top_position)
    # print(bottom_position)
    # print(top_position-bottom_position)
    # print(top_bound)
    # print(bottom_bound)
    # print(bottom_bound-top_bound)

    # lowest_energy, highest_energy = crystal_E_range_over_braggRange(bragg_angle_range,Rs,crystal_height,["Si",[4,4,0]])

    # energy_book = energy_references()
    # MoKbeta13 = energy_book['Mo']['Kb13']
    # Max_bragg = Calc_BraggFromE(MoKbeta13[0],dspace_Calc("Si",[12,12,0]))
    # Min_bragg = Calc_BraggFromE(MoKbeta13[1],dspace_Calc("Si",[12,12,0]))

    # print((Max_bragg+Min_bragg)/2)
    # print(Min_bragg)

    # function that given Analyser crystal height dimension, angle_range and indices calculates max/min energy

    # midpoint_VMXi_tests = calc_crystal_signal_range_above_midpoint(MnKbeta13[0],MnKbeta13[1],Rs,('Si',[4,4,0]))
    # the_midpoint = an.calculate_midpoint(83.9,Rs)

    # upper_level = the_midpoint+(crystal_height/2)+midpoint_VMXi_tests[1]
    # lower_level = the_midpoint-(crystal_height/2)+midpoint_VMXi_tests[0]

    # upper_level = the_midpoint+midpoint_VMXi_tests[1]
    # lower_level = the_midpoint+midpoint_VMXi_tests[0]

    # print(midpoint_VMXi_tests)
    # print(f'The midpoint:{the_midpoint}')
    # print(upper_level)
    # print(lower_level)

    # print(upper_level-lower_level)

    # NiKa12 = [7450,7478,7485]                           # High-resolution studies of the K emission spectra of nickel
    # NiKBeta13 = [8250,8263,8275]                        # Sorum, H., & Bremer, J. (1982). High-resolution studies of the K emission spectra of nickel. Journal of Physics F: Metal Physics, 12(11), 2721–2728. doi:10.1088/0305-4608/12/11/029 

    # braggAPS = Calc_BraggFromE(7050,dspace_Calc('Si',[5,3,1]))
    # print("APS thing")
    # print(braggAPS)

    # Ge620_crystal_for_use_min = min(list(crystal_E_range(25,400,85,["Ge",[18,6,0]])[0:1]))
    # Ge620_crystal_for_use_max = max(list(crystal_E_range(25,400,55,["Ge",[18,6,0]])[0:1]))
    # print(Ge620_crystal_for_use_max,Ge620_crystal_for_use_min)

    # Ge620_crystal_for_use = crystal_E_range(25,400,85,["Si",[6,6,0]])
    # print(Ge620_crystal_for_use)

    