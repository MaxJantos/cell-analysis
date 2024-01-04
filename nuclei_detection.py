# Nuclei Counting for 10x Images
# Max Jantos

from os import path, listdir

import numpy as np

from scipy.ndimage import convolve, maximum_filter, minimum_filter
from skimage.measure import block_reduce

import tifffile as tf

from tkinter import filedialog as fd
from tkinter.messagebox import showinfo


# ******************************************************************************
# File Management
# ******************************************************************************

# Open file dialog to select one file
def select_file():
    # only grab .tif files
    # filetypes can also be a list, or a tuple of one tuple with a trailing comma as below
    filetypes = (("tif file", "*.tif"),)
    filename = fd.askopenfilename(title='Open a file', initialdir='/', filetypes=filetypes)
    # returns empty string "" if cancelled or is using a thumb image
    if "_thumb" in filename: return ""
    return filename

# Get path of a directory
def select_directory():
    dirname = fd.askdirectory(title='Open directory', initialdir='/')
    return dirname

# TODO: break into scale image function?
def get_img(filename, maxthresh=17500):
    raw_img = tf.imread(filename)
    # eliminate brightspots
    raw_img[raw_img > maxthresh] = 0
    # rescale DAPI image
    u, v = np.min(raw_img), np.max(raw_img)
    raw_img = 255.0 * (raw_img - u) / (v - u)
    return raw_img

# determine if given file lines us with a desired image
def valid_file(f, rows, cols, site):
    if ("_thumb" in f) or (".tif" not in f): return (False, "None")
    # pattern 0: every well
    #         1: every other well
    #         2: 1/4 wells
    for r in rows:
        if (f'overview_{r}' not in f): continue
        for c in cols:
            # only retrieve the image at site 5
            if (f"_{r}{c}_s{site}" in f):
                return (True, f"{r}{c} s{site}") 
    return (False, "None")

# Get dict containing well_site : image path pairs from the given directory
def get_filenames(cur_path, pattern, site):
    # pattern 0: every well
    #         1: every other well?
    #         2: 1/4 wells
    rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
    cols = [n for n in range(25) if n > 0]
    if pattern == 1:
        rows = rows[::2]
    elif pattern == 2:
        rows = rows[::2]
        cols = cols[::2]
        
    file_dict = {valid_file(f, rows, cols, site)[1]:path.join(cur_path,f) for f in listdir(cur_path) if (path.isfile(path.join(cur_path,f)) and valid_file(f, rows, cols, site)[0]) }
    return file_dict

# ******************************************************************************
# Get Data
# ******************************************************************************

# Get nuclei centroids as coords
def detect_nuclei(img_in, minthresh=25, searchlen=21, mincellsize=2, minpeak=0.2):
    img = np.copy(img_in)
    img[img < minthresh] = 0

    # downsample the image to the factor based on accepted mincellsize
    multfactor = 2**(mincellsize - 1)
    # block_reduce downsamples, it does not change the images shape
    img = block_reduce(img, block_size=(multfactor, multfactor), func=np.mean)

    # expand mins and maxes
    ########### mess with the SIZE param, will impact accuracy but also speed
    max_img = maximum_filter(img, size=3) # expand nuc centers
    min_img = minimum_filter(img, size=3) # concentrate nuc centers

    # nuclei with "hollow" centers
    avoid = ((max_img + minthresh) - min_img) / (min_img + minthresh)
    # very small centers of nuclei
    cellness = (min_img) / (max_img + minthresh)

    x_grid, y_grid = np.mgrid[0:(2 * searchlen + 1), 0:(2 * searchlen + 1)]
    dist = np.sqrt(np.square(x_grid - searchlen) + np.square(y_grid - searchlen))
    krnl = searchlen / (1 + dist)
    cellness = convolve(cellness, krnl, mode='constant')
    avoid = convolve(avoid, krnl, mode='constant')
    optfcn = cellness/(avoid+0.1)

    img = block_reduce(img, block_size=(2,2), func=np.max)
    nuc_pts = []
    selected = np.zeros_like(optfcn)
    visited = np.zeros_like(optfcn)
    dr = [-1,-1,-1, 0, 0, 1, 1, 1]
    dc = [-1, 0, 1,-1, 1,-1, 0, 1]
    for r in range(2,img.shape[0]-2):
        for c in range(2,img.shape[1]-2):
            if img[r,c] > minthresh:
                # cur_nucpt is a coord
                cur_nucpt = np.array([2*r+1,2*c+1])
                # update the visited array to account for seeing this nuc point
                visited[cur_nucpt[0], cur_nucpt[1]] = 1
                skipped = False
                switched = True

                while switched:
                    switched = False
                    # check for a better adjacent maximum
                    peak = optfcn[cur_nucpt[0],cur_nucpt[1]]
                    adj_tiles = np.array([optfcn[cur_nucpt[0] + dr[i], cur_nucpt[1] + dc[i]] for i in range(len(dr))])
                    adj_max_indx = np.argmax(adj_tiles)
                    # if there is a better local peak, move to it
                    if adj_tiles[adj_max_indx]>peak:
                        # update coords
                        cur_nucpt += np.array([dr[adj_max_indx],dc[adj_max_indx]])
                        # check "visited" array for early termination
                        if visited[cur_nucpt[0], cur_nucpt[1]]:
                            skipped = True
                            break
                        visited[cur_nucpt[0], cur_nucpt[1]] = 1
                        # continue ascent while within array bounds
                        if cur_nucpt[0]>1 and cur_nucpt[1]>1 and cur_nucpt[0]<(optfcn.shape[0]-2) and cur_nucpt[1]<(optfcn.shape[1]-2):
                            switched = True

                # add found peak to list if large enough and not already seen
                if not skipped and selected[cur_nucpt[0],cur_nucpt[1]]==0 and peak>minpeak:
                    nuc_pts.append(multfactor*cur_nucpt)
                    selected[cur_nucpt[0],cur_nucpt[1]]=1

    return np.array(nuc_pts)

# Wrapper for detect nuclei on 10x images
def get_nuc_centers(img, params=None):
    if not params:
        # Default 10x parameters
        params = {
            'minthresh': 25, # threshold to eliminate noise from actual nuclei (smaller = more sensitive to noise / nuclei)
            #'maxthresh': 17500, # threshold to eliminate bright spots from actual nuclei (smaller = more sensitive to bright spots in scaling)
            'searchlen': 3,  # controls the filter that highlights nuceli centers (bigger = wider gradients to make one peak/nucleus, but more likelihood of merged nuclei)
            'mincellsize': 2,  # controls number of downsamples - will tend to remove smaller "nuclei"
            'minpeak': 0.02  # minimum allowed value of optimization function - higher values reject noise but may lose low intensity nuclei
        }
    # Find nuclei centroids
    try:
        nucpts = detect_nuclei(img, **params)
    except:
        return None
    return nucpts

# returns two dicts:
#   1) well : well's nuclei points as a list
#   2) well : number of detected nuclei
def get_well_nuc_pairs(well_file_dict, maxthresh):
    # dict mapping the name of the well and site of the image (key) 
    # to its list of nuc centers (value)
    images = dict(map(lambda x: (x[0], get_img(x[1], maxthresh)), well_file_dict.items() ))
    nuc_list_d = dict(map(lambda x: (x[0], get_nuc_centers(x[1])), images.items() ))
    nuc_count_d = dict(map(lambda x: (x[0], len(x[1])), nuc_list_d.items() ))
    return (nuc_list_d, nuc_count_d)

# seperates the well tag from the well site string
def get_well(s):
    (w, *_) = s.split(' ', 1)
    return w

# Not currently being used
# returns well estimate dictionary and an overall well average
def calc_well_data(nucCount_d, pattern, site):
    # TODO: fix calculations for the other pattern/site options
    # currently assuming only site 5
    if site > 0: # single site selection
        well_est_d = dict( map(lambda x: (get_well(x[0]), 9*x[1]), nucCount_d.items()) )
        well_avg = sum(well_est_d.values()) / len(well_est_d)
        return (well_est_d, well_avg)
    else:
        return None

# ******************************************************************************
# Misc
# ******************************************************************************

# given an image's nuc list, get data regarding the maxes at those points
def image_data_summary(img, nucpts):
    f = lambda x: img[x[0], x[1]]
    img_copy = np.copy(img)
    peaks = f(nucpts)
    max_peak = peaks.max()
    min_peak = peaks.min()
    avg_peak = np.sum(peaks) / peaks.size
    median_peak = np.sort(peaks)[peaks.size//2]
    return (peaks.max(), peaks.min(), avg_peak, median_peak, img.max(), img.min())


# ******************************************************************************
# Main Functions
# ******************************************************************************
# TODO: - threads?
#           - current concurrency through map is enough after investigating thread options
#       - parallelism?

def single_file_analysis(filename):
    if filename == "":
        showinfo(title='Error', message="No file selected")
        return None

    img = get_img(filename)
    nucpts = get_nuc_centers(img)
    nuc_count = len(nucpts)

    return (filename, nucpts, nuc_count)


def multi_file_analysis(dirname, pattern=0, maxthresh=17500, site=5):
    if dirname == "": return None

    total_files = len([f for f in listdir(dirname) if path.isfile(path.join(dirname,f))])
    if total_files == 0: 
        showinfo(title='Error', message="No files in selected directory")
        return None
    file_d = get_filenames(dirname, pattern, site)
    if len(file_d) == 0:
        showinfo(title='Error', message="No files that match the selected search parameters")
        return None

    (nucpts_d, nucCounts_d) = get_well_nuc_pairs(file_d, maxthresh)
    #well_est_d, well_avg = calc_well_data(nucCounts_d, pattern, site)

    return (file_d, nucpts_d, nucCounts_d, total_files)
