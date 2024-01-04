# Scrollable Window
# Notes: 
#   - When adding widgets, use scrollBar.viewport as the parent as that is the frame
#   - Adjusted found github code to include horizontal scrollbar
#   - Does not have mousewheel functionality
# Max Jantos
# https://gist.github.com/mp035/9f2027c3ef9172264532fcd6262f3b01

import tkinter as tk
import platform

# ************************
# Scrollable Frame Class
# ************************
class ScrollFrame(tk.Frame):
    def __init__(self, parent):
        super().__init__(parent) # create a frame (self)

        self.canvas = tk.Canvas(self, borderwidth=0, background="#ffffff")              #place canvas on self
        self.viewPort = tk.Frame(self.canvas, background="#ffffff")                     #place a frame on the canvas, this frame will hold the child widgets 
        
        self.vsb = tk.Scrollbar(self, orient="vertical", command=self.canvas.yview)     #place a vertical scrollbar on self 
        self.canvas.configure(yscrollcommand=self.vsb.set)                              #attach vertical scrollbar action to scroll of canvas
        
        self.hsb = tk.Scrollbar(self, orient="horizontal", command=self.canvas.xview)   #place a horizontal scrollbar on self 
        self.canvas.configure(xscrollcommand=self.hsb.set)                              #attach horizontal scrollbar action to scroll of canvas

        # pack version (original) - does expand with given space
        self.vsb.pack(side="right", fill="y")                                       #pack vertical scrollbar to right of self
        self.hsb.pack(side="bottom", fill="x")                                      #pack horizontal scrollbar to bottom of self
        self.canvas.pack(side="left", fill="both", expand=True)                     #pack canvas to left of self and expand to fil
        self.canvas_window = self.canvas.create_window((4,4), window=self.viewPort, anchor="nw",            #add view port frame to canvas
                                  tags="self.viewPort")

        self.viewPort.bind("<Configure>", self.onFrameConfigure)                       #bind an event whenever the size of the viewPort frame changes.
        #self.canvas.bind("<Configure>", self.onCanvasConfigure)                       #bind an event whenever the size of the canvas frame changes.
            
        self.viewPort.bind('<Enter>', self.onEnter)                                 # bind wheel events when the cursor enters the control
        self.viewPort.bind('<Leave>', self.onLeave)                                 # unbind wheel events when the cursorl leaves the control

        self.onFrameConfigure(None)                                                 #perform an initial stretch on render, otherwise the scroll region has a tiny border until the first resize

    def onFrameConfigure(self, event):                                              
        '''Reset the scroll region to encompass the inner frame'''
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))                 #whenever the size of the frame changes, alter the scroll region respectively.

    def onCanvasConfigure(self, event):
        '''Reset the canvas window to encompass inner frame when required'''
        canvas_width = event.width
        canvas_height = event.height
        self.canvas.itemconfig(self.canvas_window, width = canvas_width)            #whenever the size of the canvas changes alter the window region respectively.
        self.canvas.itemconfig(self.canvas_window, height = canvas_height)

    def onMouseWheel(self, event):                                                  # cross platform scroll wheel event
        if platform.system() == 'Windows':
            self.canvas.yview_scroll(int(-1* (event.delta/120)), "units")
        elif platform.system() == 'Darwin':
            self.canvas.yview_scroll(int(-1 * event.delta), "units")
        else:
            if event.num == 4:
                self.canvas.yview_scroll( -1, "units" )
            elif event.num == 5:
                self.canvas.yview_scroll( 1, "units" )
    
    def onEnter(self, event):                                                       # bind wheel events when the cursor enters the control
        if platform.system() == 'Linux':
            self.canvas.bind_all("<Button-4>", self.onMouseWheel)
            self.canvas.bind_all("<Button-5>", self.onMouseWheel)
        else:
            self.canvas.bind_all("<MouseWheel>", self.onMouseWheel)

    def onLeave(self, event):                                                       # unbind wheel events when the cursorl leaves the control
        if platform.system() == 'Linux':
            self.canvas.unbind_all("<Button-4>")
            self.canvas.unbind_all("<Button-5>")
        else:
            self.canvas.unbind_all("<MouseWheel>")
