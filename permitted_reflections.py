from operator import index
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import analyser_array_position as an
import energy_range_graph as er
from analyser_crystal import analyser_crystal as Cryst
from itertools import compress
from itertools import combinations_with_replacement
import matplotlib.patches as mpatches
import copy


def determine_permitted(index):
    all_even_test = np.sum([1 for i in index if i % 2 == 0])
    all_odd_test = np.sum([1 for i in index if i % 2 == 1])
    multiples_test = np.sum(index) % 4
    if all_even_test == 3 and multiples_test == 0:
        return True
    elif all_odd_test == 3:
        return True
    else:
        return False

def determine_harmonics(max_index,index_list):
    return_list = []
    harmonic_list = []
    
    a=-1
    while len(index_list) != 0:
        a=a+1
        # return_list.append(list(copy.deepcopy(i)))
        base_index = np.array(copy.deepcopy(index_list[0]))
        return_list.append(copy.deepcopy(harmonic_list))
        harmonic_list = []
        new_index = [0,0,0]
        i = 1
        while np.sum(np.array(new_index)>max_index) == 0: 
            print(np.sum(np.array(new_index)>max_index))
            new_index = list(np.array(base_index)*i)
            i=i+1
                
            if new_index in index_list:
                harmonic_list.append(copy.deepcopy(new_index))
                index_list.remove(new_index)
                print(return_list)
    return_list.pop(0)
    return return_list

# def generate_cubic_reflection_table(max_indicies):
#     a = [np.array([1,0,0])]
#     while np.sum(a[-1][-1]) <= max_indicies:
#         working_array = copy.deepcopy(a[-1])
#         min_index = np.argmin(working_array)
#         if working_array[1] < working_array[0]:
#             working_array[1] = working_array[1]+1
#             a.append(working_array)
#             continue
#         elif working_array[2] < working_array[1]:
#             working_array[2] = working_array[2]+1
#             a.append(working_array)
#             continue
#         elif working_array[2] == working_array[-1]:
#             no_max = working_array[2]+1
#             for i in range(no_max):
#                 a.append(np.array([no_max, no_max-i-1, no_max-i-1]))
#             # a.append(np.array([working_array[0]+1,0,0]))
#     return a

def generate_reflection_table(max_indicies):
    a = combinations_with_replacement(range(max_indicies),3)
    b = []
    for i in a:
        b.append(sorted(i,reverse=True))
    b = sorted(b, key=lambda x: x[0])
    return b

# def determine_harmonics(max_index,indicies_list):
    return_list = []
    i = 1
    while True:
        harmonic = indicies_list[0]*i
        i = i+1
        a = determine_permitted(harmonic)
        if a == True:
            return_list.append(harmonic)
        is_over_max_indicies = np.greater_equal(harmonic,max_index)
        if np.any(is_over_max_indicies):
            break
    
    for a in return_list:
        for indx,b in enumerate(indicies_list):
            if np.array_equal(a,b):
                indicies_list.pop(indx)
                break

    return return_list,indicies_list

def generate_permitted_reflections(highest_indicies,lattice_type):
    permitted_indicies_boolean = []
    cubic_indicies = generate_reflection_table(highest_indicies)

    print(cubic_indicies)
    for i in range(len(cubic_indicies)):
        if lattice_type == "FCC":
            permitted_indicies_boolean.append(determine_permitted(cubic_indicies[i]))
        elif lattice_type == "P":
            permitted_indicies_boolean.append(cubic_indicies[i])
        
    diamond_cubic_indicies = list(compress(cubic_indicies,permitted_indicies_boolean))
    diamond_cubic_indicies.pop(0)

    list_of_harmonic_lists = []
    ac_list = []
    # for i in range(highest_indicies):
        # permitted_harmonic_list,diamond_cubic_indicies = determine_harmonics(highest_indicies,diamond_cubic_indicies)
        # list_of_harmonic_lists.append(permitted_harmonic_list)
    list_of_harmonic_lists = determine_harmonics(highest_indicies,diamond_cubic_indicies)
    
    for i in list_of_harmonic_lists:
        a = Cryst(list(i[0]))
        a.set_harmonic_list(i)
        ac_list.append(a)
    return ac_list


def main():

    # analyser_Si440_dspace = er.dspace_Calc("Si",[4,4,0])
    Rs = 400
    angle_range = [75,85]
    crystal_height = 30
    indicies_depth = 15
    ac_objects_list = generate_permitted_reflections(indicies_depth,"P")
    a = ac_objects_list[1]
    
    # Silicon_energies = [ er.crystal_E_range_over_braggRange(angle_range,Rs,crystal_height,["Si",i]) for i in harmonics_list ]
    # Germanium_energies = [ er.crystal_E_range_over_braggRange(angle_range,Rs,crystal_height,["Ge",i]) for i in harmonics_list ]
    # print(Silicon_energies)
    # print(Germanium_energies)

    plt.figure(1)
    base_harmonic = []
    for indx,a in enumerate(ac_objects_list):
        harmonics_list = a.get_harmonic_list()
        Si_energies = [ er.crystal_E_range_over_braggRange(angle_range,Rs,crystal_height,["Si",i]) for i in harmonics_list ]
        Ge_energies = [ er.crystal_E_range_over_braggRange(angle_range,Rs,crystal_height,["Ge",i]) for i in harmonics_list ]

        base_harmonic.append(harmonics_list[0])
        for indx_2,i in enumerate(Si_energies):
            plt.hlines(y=(indx*2)-1, xmin=i[0], xmax=i[1], linewidth=2, color='r')
            print(harmonics_list[indx_2])
            plt.text(i[0],(indx*2)-1.5,"Si "+str(harmonics_list[indx_2]),fontsize=6)
            plt.hlines(y=indx*2, xmin=Ge_energies[indx_2][0], xmax=Ge_energies[indx_2][1], linewidth=2, color='g')
            plt.text(Ge_energies[indx_2][0],(indx*2)-0.5,"Ge "+str(harmonics_list[indx_2]),fontsize=6)
    # plt.yticks(list(np.linspace(-1,21,12)),base_harmonic)

    reference_edict = er.energy_references()
    metals = reference_edict.keys()
    for metal in metals:
        metals_spectra = reference_edict[metal].keys()
        for spectra in metals_spectra:
            energy_range = reference_edict[metal][spectra]
            if metal == "Cu":
                plt.axvspan(xmin=energy_range[0],xmax=energy_range[1],facecolor='b', alpha=0.5)
            elif metal == "Mo":
                plt.axvspan(xmin=energy_range[0],xmax=energy_range[1],facecolor='purple', alpha=0.5)
            elif metal == "Fe":
                plt.axvspan(xmin=energy_range[0],xmax=energy_range[1],facecolor='orange', alpha=0.5)
            elif metal == "Ni":
                plt.axvspan(xmin=energy_range[0],xmax=energy_range[1],facecolor='yellow', alpha=0.5)
            elif metal == "Mn":
                plt.axvspan(xmin=energy_range[0],xmax=energy_range[1],facecolor='magenta', alpha=0.5)
            elif metal == "Co":
                plt.axvspan(xmin=energy_range[0],xmax=energy_range[1],facecolor='lightblue', alpha=0.5)
            else:
                plt.axvspan(xmin=energy_range[0],xmax=energy_range[1],facecolor='k', alpha=0.5)

    Cu_patch = mpatches.Patch(color='blue', label='Cu')
    Mo_patch = mpatches.Patch(color='purple', label='Mo')
    Fe_patch = mpatches.Patch(color='orange', label='Fe')
    Ni_patch = mpatches.Patch(color='yellow', label='Ni')
    Mn_patch = mpatches.Patch(color='magenta', label='Mn')
    Co_patch = mpatches.Patch(color='lightblue', label='Co')

    plt.tick_params(left = False, labelleft = False)
    # plt.ylim(-2,indicies_depth*2)
    # plt.xlim(1500,20000)
    plt.xlabel("eV")

    plt.legend(handles=[Cu_patch,Mo_patch,Fe_patch,Ni_patch,Mn_patch,Co_patch],handlelength=1,fontsize=8)



    plt.show()

    # for indx,a in enumerate(ac_objects_list):
    #     harmonics_list = a.get_harmonic_list()
    #     Ge_energies = [ er.crystal_E_range_over_braggRange(angle_range,Rs,crystal_height,["Ge",i]) for i in harmonics_list ]
    #     for i in Ge_energies:
    #         plt.hlines(y=indicies_depth+indx+1, xmin=i[0], xmax=i[1], linewidth=2, color='g')
    # plt.ylim(0,10)
    # plt.show()




if __name__ == "__main__":

    main()