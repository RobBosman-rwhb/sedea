import numpy as np 
import matplotlib.pyplot as plt
import scipy.linalg as liny
from numpy.lib.function_base import diff

def load_spectra(filename):
    fobj = open(filename,'r')
    data = fobj.readlines()
    data =  [float(i.strip("\n")) for i in data ]
    data_array = np.array(data)
    return data_array

def calculate_SVD(data):
    dimension = data.shape[1]
    U,S,Vh = liny.svd(data,full_matrices=True)
    S_mat = np.zeros((data.shape[0],data.shape[1]))
    S_mat[:dimension,:dimension] = np.diag(S)
    R_matrix = np.dot(U,S_mat)
    return U,S,Vh,S_mat,R_matrix

def plot_spectra(IPNS_FeII_norm,IPNS_FeIII_norm,IPNS_1sec_norm,IPNS_Anaerobic_norm):
    plt.plot(IPNS_FeII_norm,label="FeII")
    plt.plot(IPNS_FeIII_norm,label="FeIII")
    plt.plot(IPNS_1sec_norm,label="1 sec")
    plt.plot(IPNS_Anaerobic_norm,label="Anaerobic")
    plt.legend()
    plt.show()

def plot_difference_matrix(matrix,matrix_label):
    plt.figure(1)
    spectra_no = matrix.shape
    spectra_length = spectra_no[0]
    spectra_no = spectra_no[1]
    arb_disp = 0.005
    print(spectra_no)
    for i in range(spectra_no):
        plt.plot(matrix[:,i]+(i*arb_disp),linewidth=1.5,color='k')
        plt.plot(np.zeros((spectra_length))+(i*arb_disp),'--',linewidth=0.5,color='r')
        plt.text(0,(i*arb_disp+0.0005),matrix_label[i])
    plt.show()

def plot_S_matrix(S):
    plt.figure(2)
    plt.imshow(np.diag(S))

def plot_R_matrix(spectra_no,spectra_length,arb_disp,R_matrix):
    plt.figure(3)
    for i in range(spectra_no):
        plt.plot(R_matrix[:,i]+(i*arb_disp),linewidth=1.5,color='k')
        plt.plot(np.zeros((spectra_length))+(i*arb_disp),'--',linewidth=0.5,color='r')


def main():
    filename1="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/fe_II_NORM.txt"
    filename2="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/fe_III_NORM.txt"
    filename3="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/Anaerobic_NORM.txt"
    filename4="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/1_sec_NORM.txt"
    filename5="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/2_sec_NORM.txt"
    filename6="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/3_sec_NORM.txt"
    filename7="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/4_sec_NORM.txt"
    filename8="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/6_sec_NORM.txt"
    filename9="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/8_sec_NORM.txt"
    filename10="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/10_sec_NORM.txt"
    filename11="/home/kyot/Documents/work/Diamond/Crystallina_project/XES_processing/AcdV_normalized_smoothed/12_sec_NORM.txt"

    IPNS_FeII_norm = load_spectra(filename1)
    IPNS_FeIII_norm = load_spectra(filename2)
    IPNS_Anaerobic_norm = load_spectra(filename3)
    IPNS_1sec_norm = load_spectra(filename4)
    IPNS_2sec_norm = load_spectra(filename5)
    IPNS_3sec_norm = load_spectra(filename6)
    IPNS_4sec_norm = load_spectra(filename7)
    IPNS_6sec_norm = load_spectra(filename8)
    IPNS_8sec_norm = load_spectra(filename9)
    IPNS_10sec_norm = load_spectra(filename10)
    IPNS_12sec_norm = load_spectra(filename11)

    data = np.array([IPNS_FeII_norm,IPNS_FeIII_norm,IPNS_Anaerobic_norm,IPNS_1sec_norm,IPNS_2sec_norm,IPNS_3sec_norm,IPNS_4sec_norm,IPNS_6sec_norm,IPNS_8sec_norm,IPNS_10sec_norm,IPNS_12sec_norm])
    matrix_label = ['1sec','2sec','3sec','4sec','6sec','8sec','10sec','12sec']
    m_dmin = IPNS_12sec_norm.shape[0]
    n_dmin = 8

    diff_data = np.zeros((8,IPNS_12sec_norm.shape[0]))
    for i in range(8):
        diff_data[i] = np.subtract(data[i+3,:],data[2,:])
        print(i+3)
    diff_data = np.transpose(diff_data)
    print(diff_data.shape)

    U,S,Vh,S_mat,R_matrix = calculate_SVD(diff_data)

    plot_difference_matrix(diff_data,matrix_label)
    plot_S_matrix(S)
    plot_R_matrix(8,diff_data.shape[0],0.0035,R_matrix)

    # print(U.shape,S_mat.shape, Vh.shape)
    

    # print(U.shape,S.shape, Vh.shape, R_matrix.shape)

    # IE=np.zeros((8,1))
    # for indx,i in enumerate(S):
    #     IE[indx] = (indx*np.sum(S[indx+1:]))/((m_dmin*n_dmin)*(n_dmin-indx))

    # REV=np.zeros((8,1))
    # REV_denoms=np.zeros((8,1))
    # for indx,eigen_raw in enumerate(S):
    #     REV_denoms[indx] = ((m_dmin-indx+1)*(n_dmin-indx+1))
    #     REV[indx] = eigen_raw/REV_denoms[indx]
    
    # F=np.zeros((8,1))
    # for indx in range(8):
    #     denominator_sum = np.sum(np.array(REV[indx+1:]))
    #     REVi_sum = np.sum(REV_denoms[indx-1:])
    #     F[indx] = (REV[indx]/denominator_sum)*(REVi_sum)
        
    # print(REV_denoms)

    # plt.figure(1)
    # plt.plot(IE)

    # plt.figure(2)
    # plt.plot(REV)

    # plt.figure(3)
    # plt.plot(F)
    # plt.show()






if __name__ == "__main__":
    main()
