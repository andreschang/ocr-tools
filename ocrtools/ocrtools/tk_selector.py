#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCR TOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

# Google Maps polygon coordinates tool for future implementation
# https://codepen.io/jhawes/pen/ujdgK
# http://www.gis.osu.edu/misc/map-projections/
from tkinter import *
from PIL import Image, ImageTk
import os
dirname = os.path.dirname(__file__)


def get_dims():
    import numpy as np

    class CanvasEventsDemo:
        def __init__(self, parent=None):
            self.w, self.h = 611, 320
            self.merid = 180
            canvas = Canvas(width=self.w, height=self.h)
            canvas.pack()
            canvas.bind('<ButtonPress-1>', self.onStart)
            canvas.bind('<B1-Motion>', self.onGrow)
            canvas.bind('<Double-1>', self.onClear)
            canvas.bind('<ButtonRelease-1>', self.onRelease)
            self.canvas = canvas
            self.drawn = None
            self.kind = canvas.create_rectangle
            self.coords = []

            f_lmap = dirname + '/images/lambert_cylindrical.gif'
            self.tk_im = PhotoImage(file=f_lmap)
            self.canvas.create_image(0, 0, anchor="nw", image=self.tk_im)

        def onStart(self, event):
            self.shape = self.kind
            self.start = event
            self.drawn = None

        def onGrow(self, event):
            canvas = event.widget
            if self.drawn: canvas.delete(self.drawn)
            objectId = self.shape(self.start.x, self.start.y, event.x, event.y)
            self.drawn = objectId

        def onClear(self, event):
            event.widget.delete('all')

        def onRelease(self, event):
            box_xy = (self.start.x, self.start.y, event.x, event.y)
            self.coords.append(self.xy_to_coords(box_xy))

        def xy_to_coords(self, box_xy):
            x0, x1 = box_xy[0]/self.w - 0.5, box_xy[2]/self.w - 0.5
            y0, y1 = 2 * (0.5 - box_xy[1]/self.h), 2 * (0.5 - box_xy[3]/self.h)
            lon0, lon1 = x0 * 360 + self.merid, x1 * 360 + self.merid
            lat0, lat1 = np.rad2deg(np.arcsin(y0)), np.rad2deg(np.arcsin(y1))
            # note that dimensions are output in coordinate order (lat, lon)
            return((lat0, lon0, lat1, lon1))

    tk_selector = CanvasEventsDemo()
    mainloop()

    return(tk_selector.coords)
