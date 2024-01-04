# Full GUI
# Max Jantos


import time
import platform

from tkinter import *
from tkinter import ttk
from tkinter import filedialog as fd
from tkinter.messagebox import showinfo

import matplotlib.pyplot as plt

from scrollable_window import *
from nuclei_detection import *
from csv_write import *


# ******************************************************************************
# Non-tkinter GUI tools
# ******************************************************************************

# Display image with the found nuclei centroids "highlighted"
def visualize_nucpts(scaled_img, nucpts):
    img = np.copy(scaled_img)
    for cur_nuc in nucpts:
        img[cur_nuc[0], cur_nuc[1]] = 500 # arbitrary number to make a peak
    
    img = maximum_filter(img, 5) # further highlight the detected points
    img_show = plt.imshow(img)
    plt.show()
    return

# Given a tif img as a np array
# Returns a scaled np array
def scale_img(raw_img):
    u, v = np.min(raw_img), np.max(raw_img)
    raw_img = 255.0 * (raw_img - u) / (v - u)
    return raw_img


# ******************************************************************************
# GUI Object
# ******************************************************************************
# TODO:
#       - color legend for buttons (yellow, red, green, etc)
#       - button to modify initial directory
#       - add option to select if you want to look for bright spots or empty spots
#           - cell sparsity / empty space issue
#               - investigate voroni diagram using openCV distanceTransform
#       - improve speed of post collection data analysis
class App(Tk):
    def __init__(self):
        super().__init__()

        # internal data
        self.well_dim = (16,24) # (row, col)
        self.total_files = 0
        self.files_analyzed = 0
        self.data = False
        self.dirname = None

        # setting initial pattern and site options
        self.pattern = -1 # (1 = every well, 2 = every other well, 3 = 1/4 wells), -1 means not selected
        self.site = -1 # 0 = all available imaged sites, 1-9 = only that specific site's image analyzed, -1 means not selected
        self.pattern_store = IntVar(value = -1)
        self.site_store = IntVar(value = -1)
        # setting initial thresholds
        self.count_thresh = 2000
        self.ceiling_thresh = 17500
        self.count_thresh_store = StringVar()
        self.ceiling_thresh_store = StringVar()

        # timer
        self.start_time = 0
        self.end_time = 0

        # Declare various data dictionaries for operating on a whole folder
        # well tag with site = {row tag}{col tag} s{site number} (i.e. C17 s5)
        self.file_d = None          # well tag with site : image's filepath (only for wells that have been analyzed)
        self.nucpts_d = None        # well tag with site : list of nucpt coordinates (only for wells that have been analyzed)
        self.nucCounts_d = None     # well tag with site : number of nuclei detected (only for wells that have been analyzed)
        self.sum_ips_d = None       # well tag with site : sum intensity projection (only for wells that have been analyzed)

        self.well_buttons = None    # well tag : corresponding button on wellplate
        self.data_labels = None

        # set up main window
        self.title('Well Plate GUI')
        self.resizable(True, True)
        self.geometry('1400x600')

        # Set up main gui
        altius_canvas = self.create_altius_canvas(self)
        altius_canvas.pack(side='top', anchor='n')
        # setup the buttons
        button_frame = self.create_button_frame(self) 
        button_frame.pack(side='top', anchor='w')

        # setup the output data
        output_frame = Frame(self)
        # get the data visual
        data_frame, self.data_labels = self.create_data_frame(output_frame)
        data_frame.grid(column=0, row=1, sticky='w')
        # get the settings visual
        settings_frame, self.settings_labels = self.create_settings_frame(output_frame)
        settings_frame.grid(column=0, row=0, sticky='w')
        output_frame.pack(side="right")

        # setup the wellplate
        wellplate_scrollFrame, self.well_buttons = self.create_wellplate_scrollFrame(self)
        wellplate_scrollFrame.pack(side='top', fill='both', expand=True)

    # copies given string to clipboard
    def copy(self, s):
        self.clipboard_clear()
        self.clipboard_append(s)
    
    # Command to exit any window
    def close(self, top):
        top.destroy()
        top.update()

    # Resets all internal/stored data regarding image analysis and parameters
    def reset_data(self):
        # Reset the internal data
        self.data = False
        self.total_files = 0
        self.files_analyzed = 0

        self.pattern = -1 # (1 = every well, 2 = every other well, 3 = 1/4 wells)
        self.site = -1 # 0 = all available imaged sites, 1-9 = only that specific site's image analyzed
        self.pattern_store.set(-1)
        self.site_store.set(-1)
        self.count_thresh = 2000
        self.ceiling_thresh = 17500
        self.count_thresh_store.set("")
        self.ceiling_thresh_store.set("")

        self.dirname = None

        # Reset operation timer
        self.start_time = 0
        self.end_time = 0

        # Reset dictionaries
        self.file_d = None
        self.nucpts_d = None
        self.nucCounts_d = None
        self.sum_ips_d = None

        self.update_gui()

    # Creates 2 popup windows: a data popup and a visualization popup
    def summary_popup(self, filename, nucpts, nuc_count, flag):
        newWindow = Toplevel(self)
        newWindow.title("Image Summary")
        newWindow.resizable(True, True)
        newWindow.geometry("500x350")
        raw_img = tf.imread(filename)
        (max_peak, min_peak, avg_peak, median_peak, max, min) = image_data_summary(raw_img, nucpts)

        data_title = Label(newWindow, text="Nuclei Peak Data", font="Helvetica 16 bold")
        count = Label(newWindow, text=f"Nuclei Count: {nuc_count}")
        sum_ip = Label(newWindow, text=f"Total Intensity: {np.sum(raw_img)}")
        max_val = Label(newWindow, text=f"Maximum value: {max}")
        min_val = Label(newWindow, text=f"Minimum value: {min}")
        max_peak_int = Label(newWindow, text=f"Maximum peak intensity: {max_peak}")
        min_peak_int = Label(newWindow, text=f"Minimum peak intensity: {min_peak}")
        avg_peak_int = Label(newWindow, text=f"Average peak intensity: {avg_peak}")
        median_peak_int = Label(newWindow, text=f"Median peak intensity: {median_peak}")
        flag = Label(newWindow, text=f"Flags: {flag}", font="Helvetica 12 bold")
        file_path = Label(newWindow, text=f"File Path: {filename}", justify="left", wraplength=300)

        #flag = check_for_flag(filename)
        #flag_label = Label(newWindow, text=f"Flagged for: {flag}")

        file_path.pack(side="top", anchor="nw")
        data_title.pack(side="top", anchor="nw")
        count.pack(side="top", anchor="nw")
        sum_ip.pack(side="top", anchor="nw")
        max_val.pack(side="top", anchor="nw")
        min_val.pack(side="top", anchor="nw")
        max_peak_int.pack(side="top", anchor="nw")
        min_peak_int.pack(side="top", anchor="nw")
        avg_peak_int.pack(side="top", anchor="nw")
        median_peak_int.pack(side="top", anchor="nw")
        flag.pack(side="top", anchor="nw")

        button_frame = Frame(newWindow)
        copy_button = Button(button_frame, text="Copy file path", command= lambda: self.copy(filename))
        quit_button = Button(button_frame, text="Quit", command=newWindow.destroy) #.pack(side="top", anchor="nw")

        copy_button.grid(column = 0, row = 0)
        quit_button.grid(column = 1, row = 0)
        
        button_frame.pack(side="top", anchor="nw")
        
        scaled_img = scale_img(raw_img)
        visualize_nucpts(scaled_img, nucpts)

    # Command run when you choose one file
    def single_file(self):
        filename = select_file()
        self.start_time = time.perf_counter()
        data = single_file_analysis(filename)
        self.end_time = time.perf_counter()
        if data != None:
            self.summary_popup(data[0], data[1], data[2], "None")
            return
        #showinfo(title='Error', message="Something went wrong")
        return

    # Command run when you choose a directory
    def multiple_files(self):
        self.dirname = select_directory()

        self.start_time = time.perf_counter()
        data = multi_file_analysis(self.dirname, self.pattern, self.ceiling_thresh, self.site)
        self.end_time = time.perf_counter()

        if data != None:
            self.data = True
            (self.file_d, self.nucpts_d, self.nucCounts_d, self.total_files) = data
            self.files_analyzed = len(self.nucpts_d)
            # update wellplate buttons according to newly collected data
            self.update_gui()
        return

    # Command to update/set the search and well parameters
    def update_params(self, top):
        if self.site_store.get() == 0:
            showinfo(title='Error', message="All sites option is not setup yet. Select a different site")
            return
        if self.pattern_store.get() == -1 or self.site_store.get() == -1:
            showinfo(title='Error', message="Must select a pattern and a site to continue")
            return
        self.pattern = self.pattern_store.get()
        self.site = self.site_store.get()
        self.pattern_store.set(-1)
        self.site_store.set(-1)

        top.destroy()
        top.update()

        self.multiple_files()
        return

    # Popup to choose site and pattern parameters
    def search_param_selection(self):
        newWindow = Toplevel(self)

        pattern_label = Label(newWindow, text="Pattern", font="Helvetica 18 bold")
        site_label = Label(newWindow, text="Site", font="Helvetica 18 bold")
        allFiles = Radiobutton(newWindow, text="Analyze every well", variable=self.pattern_store, value = 0)
        everyOther = Radiobutton(newWindow, text="Analyze every other well", variable=self.pattern_store, value = 1)
        everyFour = Radiobutton(newWindow, text="Analyze every four wells", variable=self.pattern_store, value = 2)

        pattern_label.grid(column = 0, row = 0)
        site_label.grid(column = 1, row = 0)
        allFiles.grid(column = 0, row = 1, sticky="w")
        everyOther.grid(column = 0, row = 2, sticky="w")
        everyFour.grid(column = 0, row = 3, sticky="w")
        
        for i in range(9):
            siteButton = Radiobutton(newWindow, text=f"Site {i + 1}", variable=self.site_store, value = i+1)
            siteButton.grid(column = 1, row = i + 1, sticky="w")

        confirm = Button(newWindow, text="Confirm", command= lambda: self.update_params(newWindow))
        cancel = Button(newWindow, text="Cancel", command= lambda: self.close(newWindow))

        confirm.grid(column = 0, row = 11)
        cancel.grid(column = 1, row = 11)

    def create_altius_canvas(self, master):
        canvas = Canvas(master, bg='black', height=52/800*600, width=216/1000*800)
        canvas.create_text(110/1000*800, 36/800*600, text="ALTIUS\nINSTITUTE\nFOR BIOMEDICAL SCIENCES\n", fill="white", font=('Helvetica 10 bold'))
        return canvas

    # Updates the thresholds based on the inputs in the threshold popup
    def update_thresholds(self, top):
        temp_count = self.count_thresh_store.get()
        temp_ceiling = self.ceiling_thresh_store.get()

        if (temp_count == "" or temp_ceiling == "" or 
            temp_count.isnumeric() == False or temp_ceiling.isnumeric() == False):
            showinfo(title='Error', message="Entries must be numbers to confirm")
            return
        self.count_thresh = int(temp_count)
        self.ceiling_thresh = int(temp_ceiling)
        self.count_thresh_store.set("")
        self.ceiling_thresh_store.set("")

        top.destroy()
        top.update()
        # updates the buttons and settings frames according to new thresholds
        self.update_gui()
        return

    # Creates popup to enter new threshold values
    def select_thresholds(self):
        newWindow = Toplevel(self)

        count = Label(newWindow, text="Nuclei Count Threshold:")
        ceiling = Label(newWindow, text="Maximum Intensity Threshold:")
        count_entry = Entry(newWindow, textvariable= self.count_thresh_store)
        ceiling_entry = Entry(newWindow, textvariable= self.ceiling_thresh_store)
        confirm = Button(newWindow, text="Confirm", command= lambda: self.update_thresholds(newWindow))
        cancel = Button(newWindow, text="Cancel", command= lambda: self.close(newWindow))

        count.grid(column=0, row=0, sticky='w')
        ceiling.grid(column=0, row=1, sticky='w')
        count_entry.grid(column=1, row=0, sticky='w')
        ceiling_entry.grid(column=1, row=1, sticky='w')
        confirm.grid(column=0, row=2, sticky='e')
        cancel.grid(column=1, row=2, sticky='w')

    def get_save_directory(self, dirname, dir_label):
        selected = fd.askdirectory(title = "Select save directory", initialdir='/')
        if selected == "":
            showinfo(title="No directory selected", message="No directory selected")
            return
        dirname.set(value=selected)
        dir_label.configure(text=selected)
        
    def write_wrapper(self, top, dirname_var, filename_var):
        dir = dirname_var.get()
        file = filename_var.get()
        if dir == '' or file == '':
            showinfo(title='Error', message='Must input valid filename and select valid directory')
            return
        top.destroy()
        export_data(dir, file, self.file_d, self.nucpts_d, self.count_thresh, 
                    self.ceiling_thresh)

    def export_data_popup(self, master):
        if self.data == False:
            showinfo(title="Error", message="No data to export")
            return
        newWindow = Toplevel(master)
        filename = StringVar()
        dirname = StringVar()

        file_label = Label(newWindow, text="Filename: ")
        dir_label = Label(newWindow, text="Directory: ")
        filename_input = Entry(newWindow, textvariable = filename)
        chosen_dir_label = Label(newWindow, text="None selected")
        dir_button = Button(newWindow, text="Select Directory", command= lambda: self.get_save_directory(dirname, chosen_dir_label))
        confirm_button = Button(newWindow, text="Export", command= lambda: self.write_wrapper(newWindow, dirname, filename))
        quit_button = Button(newWindow, text="Quit", command= newWindow.destroy)

        file_label.grid(column = 0, row = 0)
        dir_label.grid(column = 0, row = 1)
        filename_input.grid(column = 1, row = 0)
        chosen_dir_label.grid(column = 1, row = 1)
        dir_button.grid(column = 2, row = 1)
        confirm_button.grid(column = 0, row = 2)
        return

    def create_button_frame(self, master):
        frame = Frame(master)

        openFile_button = Button(frame, text='Select a file', command=self.single_file)
        openDir_button = Button(frame, text='Select a folder', command= self.search_param_selection)
        threshold_button = Button(frame, text='Change thresholds', command= self.select_thresholds)
        self.export_button = Button(frame, text='Export data', state='disabled', command= lambda: self.export_data_popup(self))
        reset_button = Button(frame, text='Reset', command=self.reset_data)
        quit_button = Button(frame, text='Quit', command= lambda: self.close(self))

        openFile_button.grid(column=0, row=0, padx=10, pady=5, sticky='w')
        openDir_button.grid(column=1, row=0, padx=10, pady=5, sticky='w')
        threshold_button.grid(column=2, row=0, padx=10, pady=5, sticky='w')
        self.export_button.grid(column=3, row=0, padx=10, pady=5, sticky='w')
        reset_button.grid(column=4, row=0, padx=10, pady=5, sticky='w')
        quit_button.grid(column=5, row=0, padx=10, pady=5, sticky='w')

        return frame

    def create_wellplate_scrollFrame(self, master):
        scrollFrame = ScrollFrame(master)

        rows = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P']
        cols = [n for n in range(25) if n > 0]
        #cols = ["01", "02", "03", "04"]

        well_buttons = {}

        # place labels and canvases for wells
        for r_i in range(self.well_dim[0] + 1):
            for c_i in range(self.well_dim[1] + 1):
                #print(str(r_i)+str(c_i))
                if c_i == 0:
                    if r_i != 0:
                        Label(scrollFrame.viewPort, text=rows[r_i - 1]).grid(column = c_i, row = r_i)
                    continue # early exit because c_i == 0 and r_i == 0, or the label was already placed
                elif r_i == 0: # c_i != 0 implied
                    Label(scrollFrame.viewPort, text=str(cols[c_i - 1])).grid(column = c_i, row = r_i)
                    continue
                else: # c_i != 0 and r_i != 0 implied
                    well_tag = rows[r_i - 1] + str(c_i)
                    # well_tag = well_tag creates a closure on this lambda call so that it is called with the correct well tag for that specific button
                    button = Button(scrollFrame.viewPort, text=well_tag, state="disabled", command= lambda well_tag=well_tag: self.get_well_summary(well_tag))
                    button.grid(column = c_i, row = r_i)
                    well_buttons[well_tag] = button

        return scrollFrame, well_buttons

    def pattern_as_str(self):
        if self.pattern == 0:
            return "All wells"
        elif self.pattern == 1:
            return "Every other well"
        elif self.pattern == 2:
            return "Every fourth well"
        return "None"

    def site_as_str(self):
        if self.site < 0:
            return "None"
        if self.site == 0:
            return "All available sites"
        return f"Site {self.site}"

    def update_settings_frame(self):
        self.settings_labels["intensity threshold"].configure(text=f"Maximum intensity threshold: {self.ceiling_thresh}")
        self.settings_labels["count threshold"].configure(text=f"Cell count threshold: {self.count_thresh}")
        self.settings_labels["pattern"].configure(text=f"Analysis pattern: {self.pattern_as_str()}")
        self.settings_labels["site"].configure(text=f"Sites selected: {self.site_as_str()}")
        return

    def create_settings_frame(self, master):
        frame = Frame(master)
        settings_labels = {}

        settings_title = Label(frame, text="Settings", font="Helvetica 16 bold")
        intensity_thresh = Label(frame, text=f"Maximum intensity threshold: {self.ceiling_thresh}")
        count_thresh = Label(frame, text=f"Cell count threshold: {self.count_thresh}")
        pattern = Label(frame, text=f"Analysis pattern: {self.pattern_as_str()}")
        site = Label(frame, text=f"Sites selected: {self.site_as_str()}")

        settings_labels["intensity threshold"] = intensity_thresh
        settings_labels["count threshold"] = count_thresh
        settings_labels["pattern"] = pattern
        settings_labels["site"] = site

        settings_title.grid(column = 0, row = 0, padx=10, pady=5, sticky='w')
        intensity_thresh.grid(column = 0, row = 1, padx=10, pady=5, sticky='w')
        count_thresh.grid(column = 0, row = 2, padx=10, pady=5, sticky='w')
        pattern.grid(column = 0, row = 3, padx=10, pady=5, sticky='w')
        site.grid(column = 0, row = 4, padx=10, pady=5, sticky='w')

        return frame, settings_labels

    # Frame containing label for overall data on selected folder
    def create_data_frame(self, master):
        frame = Frame(master)
        labels = {}
        frame.columnconfigure(0, weight=1)

        data_title = Label(frame, text="Data", font="Helvetica 16 bold")
        dirname = Label(frame, text=f"Given Directory: {self.dirname}",wraplength=200, justify="left")
        imgs_analyzed = Label(frame, text=f"{self.files_analyzed} out of {self.total_files} images analyzed")
        time_elapsed = Label(frame, text=f"Analysis took {self.end_time - self.start_time} seconds")

        #if self.files_analyzed != 0:
        #    nucCount_avg = sum(self.nucCounts_d.values()) / self.files_analyzed
        img_s5_avg = Label(frame, text=f"Average cells per analyzed image: NO DATA")

        #min_well, min_count = min(self.nucCounts_d, key=lambda x: x[1])
        #max_well, max_count = max(self.nucCounts_d., key=lambda x: x[1])
        img_min = Label(frame, text=f"Minimum nuclei count of NO DATA at NO DATA")
        img_max = Label(frame, text=f"Maximum nuclei count of NO DATA at NO DATA")

        labels["image count"] = imgs_analyzed
        labels["time elapsed"] = time_elapsed
        labels["img avg"] = img_s5_avg
        labels["min"] = img_min
        labels["max"] = img_max
        labels["directory"] = dirname

        data_title.grid(column = 0, row = 0, padx=10, pady=5, sticky='w')
        dirname.grid(column = 0, row = 1, padx=10, pady=5, sticky='w')
        imgs_analyzed.grid(column = 0, row = 2, padx=10, pady=5, sticky='w')
        time_elapsed.grid(column = 0, row = 3, padx=10, pady=5, sticky='w')
        img_s5_avg.grid(column = 0, row = 4, padx=10, pady=5, sticky='w')
        img_min.grid(column = 0, row = 5, padx=10, pady=5, sticky='w')
        img_max.grid(column = 0, row = 6, padx=10, pady=5, sticky='w')

        return frame, labels
    
    # Update the overall data based on most recently collected data
    def update_data_frame(self):
        self.data_labels["directory"].configure(text=f"Given Directory: {self.dirname}")

        self.data_labels["image count"].configure(text = f"{self.files_analyzed} out of {self.total_files} images analyzed")
        self.data_labels["time elapsed"].configure(text = f"Analysis took {self.end_time - self.start_time} seconds")

        min_well, min_count, max_well, max_count, nucCount_avg = ("NO DATA", "NO DATA", "NO DATA", "NO DATA", "NO DATA")
        if self.nucCounts_d != None:
            min_well, min_count = min(self.nucCounts_d.items(), key=lambda x: x[1])
            max_well, max_count = max(self.nucCounts_d.items(), key=lambda x: x[1])
            nucCount_avg = sum(self.nucCounts_d.values()) / self.files_analyzed

    
        self.data_labels["img avg"].configure(text=f"Average cells per analyzed image: {nucCount_avg}")
        self.data_labels["min"].configure(text=f"Minimum nuclei count of {min_count} at {min_well}")
        self.data_labels["max"].configure(text=f"Maximum nuclei count of {max_count} at {max_well}")
        return

    # Updates the wellplate buttons based on most recently collected data
    def update_buttons(self):
        if self.nucpts_d == None:
            self.export_button.configure(state='disabled')
            for button in self.well_buttons.values():
                button.configure(state = "disabled", fg = "black")
            return
        
        
        self.export_button.configure(state='normal')
        for well, button in self.well_buttons.items():
            new_state = "disabled"
            color = "black"
            adjusted_well = well + f" s{self.site}"
            if adjusted_well in self.nucCounts_d:
                color = "green"
                if self.nucCounts_d[adjusted_well] < self.count_thresh:
                    color = "red"
                filename = self.file_d[adjusted_well]
                img = tf.imread(filename)
                if img.max() > self.ceiling_thresh:

                    color = "yellow"
                #color = "green" if (self.nucCounts_d[adjusted_well] >= self.count_thresh) else "red"
                new_state = "normal" 

            button.configure(state = new_state, fg = color)
        return

    # Wrapper thats the main window gui
    def update_gui(self):
        self.update_settings_frame()
        self.update_data_frame()
        self.update_buttons()
        return

    # Wellplate button command to get a data summary of a specified well
    def get_well_summary(self, well_tag):
        if self.site != 0:
            key = well_tag + " s" + str(self.site)
            color = self.well_buttons[well_tag]['fg']
            flag = "None"
            if color == 'red':
                flag = "Insufficient cell count"
            elif color == 'yellow':
                flag = "Bright spots detected"

            self.summary_popup(self.file_d[key], self.nucpts_d[key], self.nucCounts_d[key], flag)
            return
        return


# ******************************************************************************
# Main
# ******************************************************************************
if __name__ == "__main__":
    app = App()
    app.mainloop()