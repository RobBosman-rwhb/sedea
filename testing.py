import numpy as np
from matplotlib import pyplot as plt
import itertools as it

def flatfield_test():
    input_file = "Mn_FFcoefficientMatrix_20211013.txt"
    with open(input_file) as f:
        FF_data = f.readlines()
        rows = len(FF_data)
    f.close()
    columns = len(FF_data[0].split(" "))
    FlatField = np.zeros((rows,columns),dtype=float)

    for i in range(rows):
        FF_data_line = FF_data[i].split(" ")
        for indx,a in enumerate(FF_data_line):
            FlatField[i,indx] = float(a.strip('/n'))
    
    plt.imshow(FlatField, vmax=np.percentile(FlatField,99))
    plt.show()


def main():
    a = it.combinations_with_replacement(range(15),3)
    b = []
    for i in a:
        b.append(sorted(i,reverse=True))
    b = sorted(b, key=lambda x: x[0])
    print(b)



if __name__ == "__main__":
    main()