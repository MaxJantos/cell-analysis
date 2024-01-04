# Outputting data to csv
# Max Jantos


import csv
import tifffile as tf
from os import path

# want to have an empty space flag
def get_data(well, file_d, nucpts_d, count_threshold, ceiling_threshold):
    img = tf.imread(file_d[well])
    max = img.max()
    count = len(nucpts_d[well])
    count_flag = False
    ceil_flag = False
    if count < count_threshold: count_flag = True
    if max > ceiling_threshold: ceil_flag = True
    
    return (max, count, count_flag, ceil_flag)
    
def get_flags(count_flag, ceil_flag):
    if count_flag:
        if ceil_flag:
            return "Low cell count, bright spots"
        return "Low cell count"
    if ceil_flag:
        return "Bright spots"
    return "No flags"

def export_data(dirname, filename, file_d, nucpts_d, count_threshold, 
                    ceiling_threshold):
    if filename.endswith('.csv') == False:
            filename = filename + ".csv"
    filepath = path.join(dirname, filename)
    with open(filepath, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        header = ["Well", "Count", "Max value", "Flag(s)", "Filename"]
        count_flags = 0
        ceil_flags = 0
        blank_flags = 0
        writer.writerow(header)
        for well in file_d.keys():
            (max, count, count_flag, ceil_flag) = get_data(well, file_d, nucpts_d, 
                                                            count_threshold, 
                                                            ceiling_threshold)
            flag = get_flags(count_flag, ceil_flag)
            if count_flag: count_flags += 1
            if ceil_flag: ceil_flags += 1
            row = [well, str(count), max, flag, file_d[well]]
            writer.writerow(row)
        
        summary = ["Files Analyzed", "Total flags", "Cell count flags", 
                    "Bright spot flags", "Blank spot flags"]
        writer.writerow(summary)
        summary_data = [str(len(file_d)), str(count_flags + ceil_flags), 
                        str(count_flags), str(ceil_flags), str(blank_flags)]
        writer.writerow(summary_data)
    return