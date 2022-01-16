# Radius of curvature demonstration
from re import M
from mpmath.functions.functions import cot
from numpy.lib.function_base import angle
import scipy.optimize as opt
import energy_range_graph as Ecalcs
import matplotlib.pyplot as plt
import numpy as np 
import copy
import math
import mpmath as mpm
import numpy as np

def X_position(theta,bending_radius,crystal_offset_in_X):
    # X_displacement_theta = math.asin(abs(crystal_offset_in_X)/bending_radius)
    # crystal_X_pos = bending_radius*math.cos(X_displacement_theta)
    c_square = math.pow(bending_radius,2)
    b_square = math.pow(crystal_offset_in_X,2)
    a_square = c_square - b_square
    crystal_X_pos = math.sqrt(a_square)
    return crystal_X_pos

def calc_yaw(opp,adj):
    yaw = math.atan(opp/adj)
    return yaw

def calc_crystal_yaws(array):
    cent_yaw = calc_yaw(array[2,0,2],array[2,0,0])
    out_yaw = calc_yaw(array[3,0,2],array[3,0,0])
    crystal_yawes = np.array([-out_yaw, -cent_yaw, cent_yaw, out_yaw])
    return crystal_yawes

def calc_pitch(crystal_height,central_theta_degrees,bending_radius):
    # Calculates the pitch (XRT) or roll (DLS-standard) required from a given 
    a = float(mpm.acot(crystal_height/bending_radius))      # calculate the bragg angle at the crystal 
    b = np.radians(central_theta_degrees)                   # convert theta to radians
    pitch = a-b                                             # Take the difference (roll required to make crystals meet.)
    pitch = np.radians(90)+pitch                            # add to a 90 degree offset, (this is for XRT as 90 is the flat plane)      
    print("Crystal Height = {0}".format(crystal_height))
    print("Crystal Bragg = {0}".format(np.degrees(a)))
    print("Central Bragg = {0}".format(np.degrees(b)))
    print("Pitch = {0}".format(np.degrees(pitch)))
    return pitch

def calc_crystal_pitches(Y_opp_array,theta,bending_radius):
    pitch_vector = [calc_pitch(a,theta,bending_radius) for a in np.nditer(Y_opp_array)]
    pitch_vector = np.array(pitch_vector)
    return pitch_vector

def calculate_midpoint(theta,bending_radius):
    theta_radians = np.radians(theta)
    return bending_radius*(1/math.tan(theta_radians))

def calculate_detector_height(theta,bending_radius):
    theta_radians = np.radians(theta)
    return bending_radius*(1/math.tan(theta_radians))*2

def generate_crystal_stack(bottom_crystal_height,Analyser_height,height_dimension_gap,No_crystals):
    height_gap_unit = Analyser_height+height_dimension_gap
    crystal_X_array = [ (height_gap_unit*i)+bottom_crystal_height for i in range(No_crystals)]
    crystal_X_array = np.array(crystal_X_array)
    print(crystal_X_array)
    return crystal_X_array

def generate_crystal_Z_positions(No_crystals_width,Analyser_width,Z_dimension_gap):

    Analyzer_Z_positions = np.zeros((1,No_crystals_width))
    width_and_Zgap = Analyser_width+Z_dimension_gap
    Half_crystals = No_crystals_width/2

    if No_crystals_width%2 == 0:
        for i in range(int(Half_crystals)):
            positive_index_position = int(Half_crystals+i)
            negative_index_position = int(math.floor(Half_crystals-0.1)-i)
            Analyzer_Z_positions[0,positive_index_position] = np.multiply((i+1/2),width_and_Zgap)
            Analyzer_Z_positions[0,negative_index_position] = -np.multiply((i+1/2),width_and_Zgap)
    elif No_crystals_width%2 == 1:
        for i in range(math.ceil(Half_crystals)):
            positive_index_position = int(math.floor(Half_crystals)+i)
            negative_index_position = int(math.floor(Half_crystals)-i)
            Analyzer_Z_positions[0,positive_index_position] = np.multiply(i,width_and_Zgap)
            Analyzer_Z_positions[0,negative_index_position] = -np.multiply(i,width_and_Zgap)
    
    return Analyzer_Z_positions

def calculate_crystal_X_positions(theta,bending_radius,crystals_Z_array):
    tmp_x_array = copy.deepcopy(crystals_Z_array)
    crystal_X_array = [ X_position(theta,bending_radius,i) for i in np.nditer(tmp_x_array) ]
    return crystal_X_array

def calculate_crystal_Y_positions(theta,bending_radius,Analyser_height,Y_dimension_gap,array_height_shift,No_crystals_height):
    midpoint = calculate_midpoint(theta,bending_radius)                     # Calculate the perpendicular offset
    height_gap_unit = Analyser_height+Y_dimension_gap                       # Combine the analyser height and planned gap between crystals (cryst pos multiple of this)
    if No_crystals_height == 0:                                                
        return midpoint
    elif No_crystals_height%2 == 0:
        offset_multiplier = (No_crystals_height/2)-0.5
        bottom_crystal = (midpoint-(offset_multiplier*height_gap_unit))+array_height_shift
        crystals_Y_array = generate_crystal_stack(bottom_crystal,Analyser_height,Y_dimension_gap,No_crystals_height)
    elif No_crystals_height%2 != 0:
        offset_multiplier = (No_crystals_height/2)
        bottom_crystal = (midpoint-(offset_multiplier*height_gap_unit))+array_height_shift
        crystals_Y_array = generate_crystal_stack(bottom_crystal,Analyser_height,Y_dimension_gap,No_crystals_height)
    else:
        print("Desired crystals stack appears to be negative ({0}), please reset".format(No_crystals_height))
    return crystals_Y_array

def generate_array_coordinates(theta,bending_radius,Array_dimensions_XY,Analyser_dimensions_XY,Y_offset,ZY_dimension_gaps):
    crystal_Z_positions = generate_crystal_Z_positions(Array_dimensions_XY[0], Analyser_dimensions_XY[0], ZY_dimension_gaps[0])
    crystal_X_positions = calculate_crystal_X_positions(theta, bending_radius, crystal_Z_positions)
    crystal_Y_positions = calculate_crystal_Y_positions(theta, bending_radius, Analyser_dimensions_XY[1], ZY_dimension_gaps[1], Y_offset, Array_dimensions_XY[1])
    
    crystals_pitch = calc_crystal_pitches(crystal_Y_positions,theta,bending_radius)
    array_crystal_positions = np.zeros((Array_dimensions_XY[1],Array_dimensions_XY[0],3))

    for i in range(Array_dimensions_XY[1]):
        array_crystal_positions[:,i,2] = crystal_Z_positions
        array_crystal_positions[i,:,1] = crystal_Y_positions
        array_crystal_positions[:,i,0] = crystal_X_positions
    
    return array_crystal_positions,crystals_pitch

def cot_term(bending_radius:float,theta:float):
    return 0.5*bending_radius*(1+(mpm.cot(mpm.radians(theta))))

def calculate_distance_adjustment(x,theta,bending_radius,height_offset):
    # calculation provided by John Sutter optics department phi = (arctan(R+yc)tan(thetab)) - arctan
    a = cot_term(bending_radius,theta)
    boggle = (x+a)**2 + (height_offset**2) - (a**2)
    return np.array(boggle,dtype='float64')

def calculate_angle_adjustment(Xc,Yc,bending_radius,theta):
    """
    Calculate the rotation adjustment for crystals positioned 
    above and below on the array.
    """
    Xc = float(Xc)
    a = bending_radius+Xc
    tangent_term = mpm.tan(np.radians(theta))
    first_term = float(((a)*tangent_term)/(bending_radius-Yc*tangent_term))
    second_term = float((a)*tangent_term/(bending_radius+Yc*tangent_term))
    return 0.5*(np.arctan(first_term)-np.arctan(second_term))

## Cryst E-Offset calculations.

def calculate_ray_height(theta_degrees,Rcryst,dsep,Rcentre):
    array_midpoint = calculate_midpoint(theta_degrees,Rcentre)
    h1_crystal = array_midpoint+dsep
    crystal_theta = np.degrees((math.atan(Rcryst/(h1_crystal))))
    crystal_angle_from_flat = (theta_degrees-crystal_theta)+theta_degrees
    h2_crystal = Rcryst/math.tan(np.radians(crystal_angle_from_flat))
    return h1_crystal+h2_crystal

def nearest_value_index_by_other(array_index_by,array_lookup_in,value):
    array = np.asarray(array_index_by)
    indx = (np.abs(array-value)).argmin()
    return array_lookup_in[indx]

def crystal_Rcryst_shift_approxiamation(theta_degrees,Rcryst,dsep):
    adjustment_to_Rcryst_movement = np.arange(-100,100,0.0001)
    forward_movement_offset = [calculate_ray_height(theta_degrees,Rcryst-a,dsep,Rcryst)-(2*calculate_midpoint(theta_degrees,Rcryst)) for a in adjustment_to_Rcryst_movement]
    forward_movement_offset = np.array(forward_movement_offset)
    offset = nearest_value_index_by_other(forward_movement_offset,adjustment_to_Rcryst_movement,0)
    # offset_accuracy = calculate_ray_height(theta,Rcryst-offset,dsep,Rcryst)-(2*calculate_midpoint(theta,Rcryst))
    # plt.plot(adjustment_to_Rcryst_movement,forward_movement_offset)
    # plt.xlabel("adjustment to R_cryst")
    # plt.ylabel("Ray height")
    # plt.show()
    return offset

def apply_Rcryst_X_shifts(current_crystal_position,bending_radius,offset):
    new_X_position = math.sqrt(((bending_radius-offset)**2)-(current_crystal_position[0]**2))
    new_yaw = -math.atan(current_crystal_position[0]/new_X_position)
    return new_X_position,new_yaw     

## Calculate crystal heights, 

##PLOTS for display purposes can be called as individual functions.

def analytical_cone_dispersion():
    
    shallowest_midpoint = calculate_midpoint(85,400)
    analyser_array_half_height = 25+25+5+2.5
    length_of_perpendicular_ray = math.sqrt((400**2)-(34.5**2))
    length_of_upper_most_ray = math.sqrt((length_of_perpendicular_ray**2)+(23**2))
    upward_angle_from_ninety = np.degrees(math.atan(23/length_of_perpendicular_ray))
    gradient_of_upward_ray = -(23/length_of_perpendicular_ray)
    reciprocal_gradient = -(1/gradient_of_upward_ray)
    b = 60/(math.sqrt((reciprocal_gradient**2)+1))
    a = b*reciprocal_gradient
    angle_from_downwards = np.degrees(math.sin(b/60))

    print(shallowest_midpoint)
    print(analyser_array_half_height)
    print(length_of_perpendicular_ray)
    print(gradient_of_upward_ray)
    print(a)
    print(b)
    print(angle_from_downwards)
    print(reciprocal_gradient)
    print("This is the length of the upper most ray {0}".format(length_of_upper_most_ray))
    print(upward_angle_from_ninety)
    print(shallowest_midpoint)

def plot_array_movements():
    theta= 75
    bending_radius = 400
    Array_dimensions_XY = [4,4]
    Analyser_dimensions_XY = [110,25]
    Y_offset = 25
    ZY_dimension_gaps = [5,5]

    mid_a = calculate_midpoint(75,400)
    mid_b = calculate_midpoint(85,400)

    print(mid_a)
    print(mid_b)

    a,a_pitch = generate_array_coordinates(75,bending_radius,Array_dimensions_XY,Analyser_dimensions_XY,Y_offset,ZY_dimension_gaps)
    b,b_pitch = generate_array_coordinates(85,bending_radius,Array_dimensions_XY,Analyser_dimensions_XY,Y_offset,ZY_dimension_gaps)
    full_array = crystal_bank(a)

    print(np.subtract(np.degrees(a_pitch),75))
    print(np.subtract(np.degrees(b_pitch),85))

    x = [1,2]
    plt.plot(0,0)
    plt.plot(2.1,82.5,'go')
    plt.plot(1.1,mid_a,'go')

    print(mid_a-82.5)

    plt.plot(x[0],a[0,0,1]-Analyser_dimensions_XY[1]/2,'gx')
    plt.plot(x[1],b[0,0,1]-Analyser_dimensions_XY[1]/2,'gx')
    plt.plot([0,3],[82.5,82.5],color='k',linestyle='--')
    plt.plot([0,3],[mid_a,mid_a],color='k',linestyle='--')
    plt.plot(x[0],a[0,0,1],'ro')
    plt.plot(x[0],a[0,1,1],'ro')
    plt.plot(x[0],a[0,2,1],'ro')
    plt.plot(x[0],a[0,3,1],'ro')
    plt.plot(x[1],b[0,0,1],'ro')
    plt.plot(x[1],b[0,1,1],'ro')
    plt.plot(x[1],b[0,2,1],'ro')
    plt.plot(x[1],b[0,3,1],'ro')

    plt.xlim([0,2.5])
    plt.ylim([0,160])

    plt.ylabel("Height above sample (mm)")

    plt.show()

def plot_detector_heights():
    signal_MnK13 = calculate_detector_height(Ecalcs.Calc_BraggFromE(6492.70,Ecalcs.dspace_Calc('Si',(4,4,0))),400)
    signal_FeKa1 = calculate_detector_height(Ecalcs.Calc_BraggFromE(6389.51,Ecalcs.dspace_Calc('Ge',(4,4,0))),400)
    signal_FeKa2 = calculate_detector_height(Ecalcs.Calc_BraggFromE(6403.13,Ecalcs.dspace_Calc('Ge',(4,4,0))),400)
    signal_FeK13 = calculate_detector_height(Ecalcs.Calc_BraggFromE(7059.9,Ecalcs.dspace_Calc('Ge',(6,2,0))),400)
    signal_CuKa1 = calculate_detector_height(Ecalcs.Calc_BraggFromE(8028.0,Ecalcs.dspace_Calc('Si',(4,4,4))),400)
    signal_CuKa2 = calculate_detector_height(Ecalcs.Calc_BraggFromE(8048.0,Ecalcs.dspace_Calc('Si',(4,4,4))),400)
    signal_CuK13 = calculate_detector_height(Ecalcs.Calc_BraggFromE(8903.0,Ecalcs.dspace_Calc('Si',(5,5,3))),400)
    
    print("Signal position for signal_MnK13 = {0}".format(signal_MnK13))
    print("Signal position for signal_FeKa1 = {0}".format(signal_FeKa1))
    print("Signal position for signal_FeKa2 = {0}".format(signal_FeKa2))
    print("Signal position for signal_FeK13 = {0}".format(signal_FeK13))
    print("Signal position for signal_CuKa1 = {0}".format(signal_CuKa1))
    print("Signal position for signal_CuKa2 = {0}".format(signal_CuKa2))
    print("Signal position for signal_CuK13 = {0}".format(signal_CuK13))

def plot_ray_height_vs_dsep():
    theta= 79.1
    Rcryst = 400
    plotting_dsep = np.linspace(-55,55,201)
    print(plotting_dsep)
    norm_height = calculate_detector_height(theta,Rs)
    central_array_point = calculate_detector_height(theta,Rs)
    crystal_heights = [ calculate_ray_height(theta,Rcryst,a,Rcryst) for a in plotting_dsep]

    plt.title("Central ray height agaisnt dsep")
    plt.xlabel("dsep (mm)")
    plt.ylabel("Ray height (mm)")
    plt.plot(0,central_array_point,"bo")
    plt.plot(plotting_dsep,crystal_heights)
    plt.show()

class crystal_bank:
    def __init__(self,crystal_array):
        self.crystal_array = crystal_array
        # Bank1
        self.cryst0 = (crystal_array[0,0,2],crystal_array[0,0,0],crystal_array[0,0,1])
        self.cryst1 = (crystal_array[0,1,2],crystal_array[0,1,0],crystal_array[0,1,1])
        self.cryst2 = (crystal_array[0,2,2],crystal_array[0,2,0],crystal_array[0,2,1])
        self.cryst3 = (crystal_array[0,3,2],crystal_array[0,3,0],crystal_array[0,3,1])
        # Bank2
        self.cryst4 = (crystal_array[1,0,2],crystal_array[1,0,0],crystal_array[1,0,1])
        self.cryst5 = (crystal_array[1,1,2],crystal_array[1,1,0],crystal_array[1,1,1])
        self.cryst6 = (crystal_array[1,2,2],crystal_array[1,2,0],crystal_array[1,2,1])
        self.cryst7 = (crystal_array[1,3,2],crystal_array[1,3,0],crystal_array[1,3,1])
        # Bank3   
        self.cryst8 = (crystal_array[2,0,2],crystal_array[2,0,0],crystal_array[2,0,1])
        self.cryst9 = (crystal_array[2,1,2],crystal_array[2,1,0],crystal_array[2,1,1])
        self.cryst10 = (crystal_array[2,2,2],crystal_array[2,2,0],crystal_array[2,2,1])
        self.cryst11 = (crystal_array[2,3,2],crystal_array[2,3,0],crystal_array[2,3,1])  
        # Bank4
        self.cryst12 = (crystal_array[3,0,2],crystal_array[3,0,0],crystal_array[3,0,1])
        self.cryst13 = (crystal_array[3,1,2],crystal_array[3,1,0],crystal_array[3,1,1])
        self.cryst14 = (crystal_array[3,2,2],crystal_array[3,2,0],crystal_array[3,2,1])
        self.cryst15 = (crystal_array[3,3,2],crystal_array[3,3,0],crystal_array[3,3,1])
        # Bank information 
        self.yaw_vector = calc_crystal_yaws(crystal_array)


if __name__ == '__main__':

    theta= 79.7
    bending_radius = 400
    cryst_width = 25
    sep = 6
    cryst_height = 1.5*(sep+cryst_width)

    Array_dimensions_XY = [4,4]
    Analyser_dimensions_XY = [110,25]
    Y_offset = 0
    ZY_dimension_gaps = [5,5]
    theta_radians = np.radians(theta)
    dsep = -15
    Rs = 400
    x = np.array(1,dtype='float64')



    b = opt.root(calculate_distance_adjustment,x0=x,args=(theta,bending_radius,cryst_height),method='hybr')

    # print(b)
    if b.success:
        print(f"Crystal heigh is {cryst_height}")
        print(f"Xc determined as {b.x}, final value {b.fun}")
        distance_adjustment = b.x
        angle_adjustment = calculate_angle_adjustment(distance_adjustment,cryst_height,bending_radius,theta)
        angle_adjustment = np.degrees(angle_adjustment)
        print(f"φ is {angle_adjustment}")

    

    # print(calculate_distance_adjustment(-2.8174,theta,bending_radius,37.5))

    # x = 2.817
    # a = -x*bending_radius*(mpm.csc(np.radians(theta)**2))-(x**2)+0
    # print(a)



    # array_height1 = calculate_midpoint(75,500)
    # array_height2 = calculate_midpoint(85,400)

    # print(array_height1)
    # print(array_height2)



    # # analytical_cone_dispersion()
    # opposite = math.cos(np.radians(85))*400
    # # opposite_upper = 
    # a = calculate_midpoint(85,bending_radius)
    # mid_ray_theta = math.atan(np.radians(85))*bending_radius
    # upper_ray_height = (a-Analyser_dimensions_XY[1])
    # # upper_ray_length = upper_ray_height/math.asin(np.radians(85))
    # lower_ray_height = a+Analyser_dimensions_XY[1]
    # # lower_ray_length = lower_ray_height/math.asin(np.radians(85))

    # print(opposite)
    # print(lower_ray_length)

    # offset,offset_accuracy = crystal_Rcryst_shift_approxiamation(theta,bending_radius,dsep)
    # print(offset)
    # print(offset_accuracy)

    # a = generate_array_coordinates(theta,bending_radius,Array_dimensions_XY,Analyser_dimensions_XY,Y_offset,ZY_dimension_gaps)


    # height = Rs*mpm.cot(np.radians(theta))
    # print(height)

    # p = Rs/(math.sin(theta)-(dsep/math.cos(theta)))
    # print(p)
    # p = 50.86
    # a = ((Rs/math.sin(theta_radians))-p)/math.cos(theta_radians)
    # print(a)

