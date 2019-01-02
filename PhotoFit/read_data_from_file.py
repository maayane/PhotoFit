import numpy as np
import pandas as pd


def read_data_into_numpy_array(path,delimiter=',',header=False,skiprows=None,no_repeat_rows=True):
    """Description: numpy genfrom txt but with mixed types
        Input  :- a path to an ascii file
                - delimiter
                - header: if false, then the output is simply a numpy array with mized types
                         if true, then there are 3 outputs: the numpy array, a numpy array with the fields names, and a dictionnary of which the keys
                         are the fields names and the values are the columns
        Output :- if header=True: an array with the data, an array with the header,
        and a dictionnary where the keys are the header values and the values are the corresponding columns
                - if header=False: the data array only
        a numpy array with mixed types
        with the corresponding data
        Tested : ?
             By : Maayane T. Soumagnac Nov 2016
            URL :
        Example: Swift=read_data_from_file.read_data_into_numpy_array('Swift.txt',delimiter=' ',header=True)[0]
                header_array=read_data_from_file.read_data_into_numpy_array('Swift.txt',delimiter=' ',header=True)[1]

        Reliable:  """
    if header==True:
        pd_array = pd.read_csv(path, header=0,
                                     index_col=False, delimiter=delimiter,skiprows=skiprows)
        pd_numpy_array=pd_array.values
        #print(pd_numpy_array)
        #print(pd_numpy_array[:, :-2])
        if no_repeat_rows==True:
            unique, index = np.unique(pd_numpy_array.astype("<U22"),axis=0,return_index=True)
            #print(unique)
            pd_numpy_array_norepeat=pd_numpy_array[index]

        pd_keys=np.array(pd_array.keys())
        pd_dict=dict()
        for i,k in enumerate(pd_keys):
            if no_repeat_rows==True:
                pd_dict[k]=pd_numpy_array_norepeat[:,i]
            else:
                pd_dict[k]=pd_numpy_array[:,i]
        return pd_numpy_array,pd_keys,pd_dict



    else:
        pd_array = pd.read_csv(path, index_col=False, delimiter=delimiter,skiprows=skiprows)
        pd_numpy_array = pd_array.values
        return pd_numpy_array
