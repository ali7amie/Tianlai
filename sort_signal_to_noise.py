import numpy as np
import numpy.ma as ma


# the signal to noise ratio here is the ratio of the average of central kernel and the std of the peripheric kernel, so this function compute the std and later the ratio


def compute_signal_to_noise_results(detection_tools,detected_src_coor_upper_pixcorner,map):
    
    #computing signal to noise ratio
    signal_to_noise_map = detection_tools[1][0]/detection_tools[1][2][2]
    signal_to_noise_list = signal_to_noise_map[detected_src_coor_upper_pixcorner]
    signal_list = map[detected_src_coor_upper_pixcorner]
    signal_to_noise_and_coor_list = np.transpose(( np.arange(0,len(detected_src_coor_upper_pixcorner[0])),signal_list,signal_to_noise_list,detected_src_coor_upper_pixcorner[1],detected_src_coor_upper_pixcorner[0] ))

    argsorted_signal_to_noise_index_list = np.flip(np.argsort(signal_to_noise_list))

    sorted_signal_to_noise_and_coor_list = np.transpose(( argsorted_signal_to_noise_index_list, signal_list[argsorted_signal_to_noise_index_list],signal_to_noise_list[argsorted_signal_to_noise_index_list], detected_src_coor_upper_pixcorner[1][argsorted_signal_to_noise_index_list], detected_src_coor_upper_pixcorner[0][argsorted_signal_to_noise_index_list], np.zeros_like(signal_to_noise_list) ))

    return (signal_to_noise_map,signal_to_noise_and_coor_list,sorted_signal_to_noise_and_coor_list)

    