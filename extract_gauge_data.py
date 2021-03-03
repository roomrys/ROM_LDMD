'''
Module to extract gauge data using gauges.data and _outputs/gauge000X.txt to form matrices used in Lagrangian DMD.
'''

import numpy as np
from scipy.io import savemat

def handle_init_data(reader, words, init_data, line_num, ngauges, eog):
    '''
    Handles the data for the gauges.data file
        Parameters:
            reader ():
            words (list of strings): words in line
            init_data (numpy array): matrix containing data from file
            line_num (int): current line number being read
            ngauges (int): number of gauges
            eog (int): expected end of gauges line
        Returns:
            init_data (numpy array): init data for gauges
            ngauges (int): number of gauges
            eog (int): expected end of gauges line
    '''
    if 'ngauges' in words:
        try:
            ngauges = int(words[0])
            eog = ngauges + line_num # line number where we expect gauge initial conditions to end
            init_data = np.inf * np.ones((ngauges, 5))
            return init_data, ngauges, eog
        except:
            print('Error: Could not find number of gauges.')
            raise
            
    if line_num <= eog:  # only occurs after finding ngauges
        words = [float(num) for num in words]
        init_data[line_num - eog - 1, :] = words
        
    return init_data, ngauges, eog

def handle_gauge_data(reader, words, gauge_data, line_num, ngauges, eog):
    '''
    Handles the data for the gauge00001.txt files
        Parameters:
            reader ():
            words (list of strings): words in line
            gauge_data (numpy array): matrix containing data from file
            line_num (int): current line number being read
            ngauges (int): N/A
            eog (int): N/A
        Returns:
            gauge_data (numpy array): data for specific gauge
            ngauges (int): N/A
            eog (int): N/A
    '''
    if line_num < 4:  # skip over commented lines
        return gauge_data, ngauges, eog
    
    words = [float(num) for num in words]
    if line_num == 4:
        gauge_data = np.array(words)
    else:
        gauge_data = np.vstack((gauge_data.reshape(-1, 6), words))
    
    return gauge_data, ngauges, eog


def get_file_data(file_path, handle_case=handle_init_data):
    '''
    Returns the data from file depending on the handle_case function.
        Parameters:
            file_path (str): path to file
            handle_case (function): specific task to perform for specific files
        Returns:
            file_data (numpy array): data for gauges
    '''
    reader = open(file_path, 'r')  # open the gauges.data file

    eof = False  # set our variable to detect end of file to false
    empty_line = False  # variable to store if previous line was empty
    max_iter = 500000000
    line_num = 0
    ngauges = -1
    eog = -1
    file_data = np.empty((1, 1))
    try:
        # READING LOOP 1
        while not (eof or line_num > max_iter):  # while we have not yet reached the eof...
            line_num += 1

            line = reader.readline()  # read file line by line
            if not line:  # current line is empty
                if empty_line:  # if prev line was also empty, then eof
                    eof = True
                empty_line = True  # set empty_line variable true
                continue
            
            empty_line = False
            words = line.split()  # split line into words based on whitespaces
            file_data, ngauges, eog = handle_case(reader, words, file_data, line_num, ngauges, eog)

    finally:
        reader.close()

    return file_data, ngauges

if __name__ == '__main__':
    g_all = 'gauges.data'  # gauge data file should be in current working directory
    (init_data, ngauges) = get_file_data(g_all)  # get inital data of all gauges
    print(f"init_data = {init_data}")
    print(f"init_data.shape = {init_data.shape}")
    print(f"ngauges = {ngauges}")
    
    path_to_output = '/home/jovyan/clawpack-v5.7.1/geoclaw/rom/_output/'
    base_fname = 'gauge'
    
    gauges_dict = {}
    for n in range(ngauges):
        gnum_str = str(n).rjust(5, '0')
        gauge_fname = path_to_output + base_fname + gnum_str + '.txt'
        print(f"gauge_fname = {gauge_fname}")

        (gauge_data, _) = get_file_data(gauge_fname, handle_gauge_data)  # get data of specific gauge
        gauges_dict['g_' + gnum_str] = gauge_data
        print(f"guage_data = {gauge_data}")
        print(f"gauge_data.shape = {gauge_data.shape}")
    
    savemat("matlab_matrix.mat", gauges_dict)  # export gauges0000X.txt data to matlab matrix file