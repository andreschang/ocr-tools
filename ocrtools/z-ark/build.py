#!/usr/bin/env python
# -*- coding: utf-8 -*-

########################################
#
#  OCRTOOLS - Making climate data malleable
#  Copyright (C) 2018 Andres Chang
#
########################################

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import ocrtools.stage as st
ndivs = st.ndivs


class build(object):

    def __init__(self, query, base_list, verbose=True, **kwargs):
        """
        Initializes a build object that can be used to modify and/or
        interpolate between existing climate simulation data products

        Args:
        * query (object): query object with attributes that are shared by build
         object
        * base_list (list): data product (usually timeseries) that is used as a
         basis for newly developed products

        Kwargs:
        * verbose (Boolean): print extra info for debugging etc.
        * dt (str): specify whether input data is monthly or daily
        * stage (str) [optional]: specify a stage for class to use. Otherwise,
         it will be set based on that of input query object
        * combine_steps (int) [optional]: how many consecutive datapoints
         should be used for the analysis (ex. January 1-5, combine_steps = 5)
        """

        # Modify the input query class for naming
        self.query = query
        self.query.src = 'build'

        self.base_list = base_list
        self.verbose = verbose
        try:
            self.dt = kwargs['dt']
        except KeyError:
            try:
                self.dt = query.dt
            except AttributeError:
                raise "Could not autofill dt. Please specify dt in build init \
                        with kwarg dt"
        self.ndiv = ndivs[self.dt]
        try:
            self.stage = kwargs['stage']
        except KeyError:
            self.stage = self.query.stage
        try:
            self.directories = self.stage.directories
        except AttributeError:
            self.directories = {"cesm-raw": "", "cesm-reformatted": "",
                                "other-reformatted": "", "scratch": "",
                                "plot": "", "csv": ""}
        try:
            self.subfolders = self.stage.subfolders
        except AttributeError:
            self.subfolders = {"cesm-raw": "", "cesm-reformatted": "",
                               "other-reformatted": "",
                               "scratch": "", "plot": "", "csv": ""}
        try:
            self.naming = self.stage.naming
        except AttributeError:
            self.naming = {"cesm-reformatted": ["var_name", "dt", "mem",
                           "time_slice"], "other-reformatted":
                           ["var_name", "dt", "time_slice"], "spatial-average":
                           ["var_name", "dt", "yr_range", "mean", "scratchId"],
                           "cesm-report": ["var_name", "dt", "mem"],
                           "other-report": ["var_name", "dt"]}
        try:
            self.time_as = self.stage.time_as
        except AttributeError:
            self.time_as = "date"
        try:
            self.base_yr = self.stage.base_yr
        except AttributeError:
            self.base_yr = 1920
        try:
            self.combine_steps = kwargs['combine_steps']
        except KeyError:
            pass

    def annual_cycle(self, **kwargs):
        """
        Returns average annual cycle (i.e. climatology) of input list

        Kwargs:
        list (str) [optional]: which list to use for calculating annual cycle
        combine_steps (int) [optional]: how many consecutive datapoints should
        be used for the analysis (ex. January 1-5, combine_steps = 5)
        savgol_window (float) [optional]: coefficient for Savitzky-Golay
        filter. Higher number leads to more smoothing
        """

        try:
            list0 = kwargs['list']
        except KeyError:
            list0 = self.base_list
        try:
            combine_steps = kwargs['combine_steps']
        except KeyError:
            try:
                combine_steps = self.combine_steps
            except AttributeError:
                combine_steps = 1
        try:
            savgol_window = kwargs['savgol_window']
        except KeyError:
            savgol_window = 0

        # List of each starting data index (ind0) to be used for annual cycle
        # If combine_steps = 1, then this is every div. If = 2, every other div
        all_step0 = np.arange(0, self.ndiv, combine_steps)
        pad_list = np.concatenate((list0, list0[-self.ndiv:]))
        all_vals, average_steps, div_average_steps = [], [], []

        for i in all_step0:
            # ind0: list of starting segment indices (ex. january year 0 OR
            # january+february year 0)
            # indf: list of ending segment indices (ex. february year f OR
            # march+april year f)
            ind0, indf, all_ind0, nstep = self.get_build_indices(
                list0, i, combine_steps)
            all_vals.append([pad_list[int(k)] for k in all_ind0])

        for j in range(len(all_vals)-1):
            average_steps.append(np.mean(all_vals[j+1])-np.mean(all_vals[j]))
        average_steps.append(np.mean(all_vals[0])-np.mean(all_vals[-1]))

        for i in range(len(all_vals)):
            for j in range(combine_steps):
                div_average_steps.append(average_steps[i])

        new_list = [np.mean(list0[0:combine_steps])]
        for i in range(len(list0)-1):
            year, div = divmod(i, self.ndiv)
            which_step = div
            step_i = div_average_steps[which_step]/(combine_steps)+new_list[i]
            new_list.append(step_i)

        # Smooth with Savitzky-Golay filter
        if savgol_window > 0:
            new_list = signal.savgol_filter(new_list, savgol_window, 2)

        if self.verbose is True:
            print('Annual cycle of length '+str(len(new_list)) +
                  ' calculated from list with combine_steps = ' +
                  str(combine_steps))

        return new_list

    def get_variance(self, **kwargs):
        """
        Returns step_var - variance of each equivalent div (ex. each January) -
        & blur_var - variance within each div
        (only non-zero if combine_steps > 0)

        Kwargs:
        list (str) [optional]: which list to use for calculating annual cycle
        combine_steps (int) [optional]: how many consecutive datapoints should
        be used for the analysis (ex. January 1-5, combine_steps = 5)
        """

        try:
            list0 = kwargs['list']
        except KeyError:
            list0 = self.base_list
        try:
            combine_steps = kwargs['combine_steps']
        except KeyError:
            try:
                combine_steps = self.combine_steps
            except AttributeError:
                combine_steps = 1

        self.ndiv = ndivs[self.dt]
        all_step0 = np.arange(0, self.ndiv, combine_steps)
        pad_list = np.concatenate((list0, list0[-self.ndiv:]))
        all_blur_vars, all_step_vars = [], []

        for i in all_step0:
            ind0, indf, all_ind0, nstep = self.get_build_indices(
                list0, i, combine_steps)
            div_vals = [pad_list[int(k)] for k in all_ind0]

            # get blur variances for each div
            all_blur_var = []
            all_blur_means = []

            for j in range(nstep+1):
                select_blur_range = np.arange(
                    j*combine_steps, (j+1)*combine_steps)
                # print(select_blur_range)
                # print([div_vals[z] for z in select_blur_range])
                step_vals = [div_vals[z] for z in select_blur_range]
                all_blur_var.append(np.var(step_vals))
                all_blur_means.append(np.mean(step_vals))

            all_blur_vars.append(all_blur_var)
            all_step_vars.append(np.var(signal.detrend(all_blur_means)))
            # print(all_blur_vars)
            # print(np.var(signal.detrend(all_blur_means)))

        if self.verbose:
            print(
                'Blur var of length '+str(len(all_blur_vars)) +
                ', with sublist of length '+str(len(all_blur_vars[0])) +
                ', and step var of length '+str(len(all_step_vars)) +
                ' calculated')

        return all_step_vars, all_blur_vars

    def new(self, list2, save_rand=False, **kwargs):
        """
        Makes new climate data based on existing modeled data and user options
        using a step-wise approach (i.e. OCR Tools calculates each timestep in
        sequence.)

        Args:
        * list2 (list): A list of climate data that is used as the
        'experimental' or end-member run compared to the list used to
        initialize class, which is similar to a 'control' run.
        Must be the same length
        * save_rand (bool): Also return list of random values calculated during
        build loop

        Kwargs:
        * a (float) [optional]: amplitude of change, i.e. how much the returned
        data aims to deviate from list1 into (or beyond) list2. a = 1 means
        that the new data should start equivalent to list1 and end equivalent
        to list2
        * step_var_a (float) [optional]: how much random variance based on
        detrended step variance to incorporate into new modeled data
        * blur_var_a (float) [optional]: how much random variance based on
        variance within each combine_step to incorporate into new modeled data
        * var_min (float or Bool) [optional]: set lower limit to calculate
        values
        * var_max (float or Bool) [optional]: set upper limit to calculate
        values
        * snap (float) [optional]: corrective amount for if data strays far
        from source
        * combine_steps (int) [optional]: how many consecutive datapoints
        should be used
        for the analysis (ex. January 1-5, combine_steps = 5)
        * head (int) [optional]: how many first years to use to set starting
        point of data
        * tail (int) [optional]: how many last years to use to set projected
        end of data
        * print_report (Bool) [optional]: saves a plot with list1 and list2
        plots
        * load_rand (list) [optional]: loads a list of random values instead
        of calculating random variance at each timestep
        * snap_atten (float) [optional]: how much to attenuate dynamic snapping
        * savgol_window (float) [optional]: coefficient for Savitzky-Golay
        filter used in annual_cycle(). Higher number leads to more smoothing
        * hist_stretch (bool) [optional]: enhances contrast of data series
        based on input value and max/min values of input data.
        * hist_dist (function) [optional]: distribution of histogram
        stretching. Should be 0-1.
        """

        try:
            combine_steps = kwargs['combine_steps']
        except KeyError:
            try:
                combine_steps = self.combine_steps
            except AttributeError:
                combine_steps = 1
        try:
            print_report = kwargs["print_report"]
        except KeyError:
            print_report = False
        try:
            step_var_a = kwargs["step_var_a"]
        except KeyError:
            step_var_a = 10.
        try:
            blur_var_a = kwargs["blur_var_a"]
        except KeyError:
            blur_var_a = 10.
        try:
            var_min = kwargs["var_min"]
        except KeyError:
            var_min = None
        try:
            var_max = kwargs["var_max"]
        except KeyError:
            var_max = None
        try:
            a = kwargs['a']
        except KeyError:
            a = 1
        try:
            snap = kwargs["snap"]
            if snap > 15:
                snap = 15
        except KeyError:
            snap = 10.
        try:
            snap_atten = kwargs["snap_atten"]
            if snap_atten > 50:
                snap_atten = 50
        except KeyError:
            snap_atten = 20
        try:
            head = kwargs["head"]
        except KeyError:
            head = 1
        try:
            tail = kwargs["tail"]
        except KeyError:
            tail = 1
        try:
            load_rand = kwargs["load_rand"]
        except KeyError:
            load_rand = None
        try:
            savgol_window = kwargs["savgol_window"]
        except KeyError:
            savgol_window = 0
        try:
            hist_stretch = kwargs["hist_stretch"]
        except KeyError:
            hist_stretch = False
        try:
            if kwargs["hist_dist"]:
                c_dist = kwargs["hist_dist"]
            else:
                def c_dist(x):
                    return(x)
        except KeyError:
            def c_dist(x):
                return(x)

        # First calculate the mean cycle (climatology) of both input lists
        list1 = self.base_list
        list1_base = self.annual_cycle(list=list1, combine_steps=combine_steps,
                                       savgol_window=savgol_window)
        list2_base = self.annual_cycle(list=list2, combine_steps=combine_steps,
                                       savgol_window=savgol_window)

        # calculate variances:
        # step var is detrended variance of each div
        # (ex. january or mean(january, february))
        # blur var is variance of each "combine_steps" segment
        # (ex. january (0) or var(january, february) year 0)
        list1_step_var, list1_blur_var = self.get_variance(
            list=list1, combine_steps=combine_steps)
        list2_step_var, list2_blur_var = self.get_variance(
            list=list2, combine_steps=combine_steps)

        # Make an average of the list1 and list2 bases
        list0_base = np.mean([list1_base, list2_base], axis=0)
        max0 = np.amax(list1 + list2)
        min0 = np.amin(list1 + list2)
        data_range = max0-min0

        # Some print-outs
        if self.verbose:
            print("List 1 step_var: ")
            print(list1_step_var)
            print("List 1 blur_var: ")
            print(list2_blur_var)
            print("List 2 step_var: ")
            print(list2_step_var)
            print("List 2 blur_var: ")
            print(list2_blur_var)

        if print_report:
            x = np.arange(0, len(list1_base))
            plt.plot(x, list1_base, label='list 1')
            plt.plot(x, list2_base, label='list 2')
            plt.legend(loc='best')
            f1 = self.query.outfile(mode='spatial_average-plot',
                                    custom_tags=['annual_cycle'])
            plt.savefig(f1, dpi=200)
            plt.close()

        # Array of each starting step in year (ex. if combine_steps = 3
        # all_step0 = [0, 3, 6])
        all_step0 = np.arange(0, self.ndiv, combine_steps)
        pad_list1 = np.concatenate((list1, list1[-self.ndiv:]))
        pad_list2 = np.concatenate((list2, list2[-self.ndiv:]))

        # Absolute 'full_step' is calculated, which is the change between
        # div(i~), yr0 and div(i~+1), yrf. Where i~ here means the full
        # 'combine_steps' div. Ex. January 1-5 1995, January 6-10 2000.
        # Opt steps is the full step divided by number of years
        full_steps, opt_steps = [], []

        for i in all_step0:
            for j in range(combine_steps):
                # ind0: list of starting segment indices
                # (ex. january year 0 OR january+february year 0)
                # indf: list of ending segment indices
                # (ex. february year f OR march+april year f)
                ind0, indf, all_ind0, nstep = self.get_build_indices(
                    list1, i, combine_steps)

                full_step = a*(np.mean([pad_list2[k-g*self.ndiv] for k
                                       in indf for g in range(tail)]) -
                               np.mean([pad_list1[m+n*self.ndiv]for m
                                       in ind0 for n in range(head)]))
                opt_step = full_step/(len(list1)/self.ndiv - 1)

                full_steps.append(full_step)
                opt_steps.append(opt_step)
                if (self.verbose is True and j == 0):
                    print('ind0, indf, nstep')
                    print(ind0, indf, nstep)
                    print('Build indices:')
                    print(all_ind0)
                    print('Full step')
                    print(full_step)
                    print('opt_step')
                    print(opt_step)

        # List with snap amounts
        self.snap_list = [0]

        # Generate new climate data starting at list1
        new_list = [list1[0]]
        list2_blur_var.append(list2_blur_var[0])
        new_rand_list = []

        for i in range(len(list1)-1):
            which_year, which_div = divmod(i, self.ndiv)
            which_blur = int(which_div/combine_steps)

            # Value at step_i+1 calculated from 'opt_step'
            step_i1 = new_list[i] + opt_steps[which_blur]

            # Add some randomness based on observed deviation in input
            # data and user args.
            # First, get average standard deviation in step i and i+1
            which_year_next, which_div_next = divmod(i+1, self.ndiv)
            which_blur_next = int(which_div_next/combine_steps)
            step_dev = ((list2_step_var[which_blur] +
                        list2_step_var[which_blur_next])/2)**0.5
            blur_dev = (
                (list2_blur_var[which_blur][which_year] +
                 list2_blur_var[which_blur_next][which_year_next])/2)**0.5

            if load_rand is None:
                r1, r2 = np.random.normal(), np.random.normal()
            else:
                r1, r2 = load_rand[i][0], load_rand[i][1]

            step_i1 = step_i1 + r1 * step_dev * (step_var_a/50.) + r2 * \
                blur_dev * (blur_var_a/50.)
            new_rand_list.append([r1, r2])

            # Normalize snapping envelope to 5 standard deviations
            dev_from_base = np.abs(
                step_i1 - (list0_base[i+1]))/(step_dev * (15-snap))
            if dev_from_base > 1:
                dev_from_base = 1

            snap_amt = dev_from_base**((50-snap_atten)/5)
            self.snap_list.append(snap_amt)
            step_i1 = (snap_amt * list0_base[i+1]) + (1-snap_amt)*step_i1

            # If var_min or var_max exceeded, set to min/max allowed value
            if var_min is not None:
                if step_i1 < var_min:
                    step_i1 = var_min
            if var_max is not None:
                if step_i1 > var_max:
                    step_i1 = var_max

            # add it to the list!
            new_list.append(step_i1)

        # enhance contrast
        if hist_stretch:
            in_mean = np.mean(new_list)
            in_dev = np.var(new_list)**0.5
            contrast_lims = [in_mean-in_dev, in_mean+in_dev]
            c_range = contrast_lims[1]-contrast_lims[0]

            if var_min is not None:
                out_min = var_min
            else:
                out_min = np.amin(list1+list2)
            if var_max is not None:
                out_max = var_max
            else:
                out_max = np.amax(list1+list2)

            for div in range(len(new_list)):
                if (
                      new_list[div] > contrast_lims[0] and
                      new_list[div] < contrast_lims[1]):
                    new_list[div] = (out_min +
                                     (c_dist((new_list[div] -
                                      contrast_lims[0]) /
                                      c_range)*(out_max-out_min)))

        if save_rand is True:
            return new_list, new_rand_list
        else:
            return new_list

    def magnet(self, list0, listM, **kwargs):
        try:
            magnet_a = kwargs["magnet_a"]
        except KeyError:
            magnet_a = 50.
        try:
            which_divs = kwargs["divs"]
        except KeyError:
            which_divs = np.arange(0, self.ndiv)

        new_list = list0

        for i in range(len(list0)-1):
            which_year, which_div = divmod(i, self.ndiv)
            if which_div in which_divs:

                new_list[i] = (list0[i]*(100.-magnet_a)/100.) + \
                              (listM[i]*magnet_a/100.)

        return new_list

    def get_build_indices(self, list0, div, combine_steps):
        """
        Returns index lists - ind0, indf, and all_ind0 - and nstep.
        ind0 is the first div under consideration and indf is the next
        regular step, taken from the last year. all_ind0 includes all the
        "ind0s" from every year in list (ex. if ind0 is Jan 1990,
        indf is February 1995, and all_ind0 is Jan 1990-1995)

        Args:
        * list0 (list): Input list (only used for length assessments)
        * div (int): Specify which div in the year under consideration
         (0-364 or 0-11)
        * combine_steps (int): Specify number of steps to include in each
        returned list, ind0, indf, and all_ind0. Introduces some distortion

        """

        # get num of full years and remainder divs in full list
        yrf, yrf_div = divmod(len(list0), self.ndiv)

        # get set of indexes for data to be averaged
        ind0 = np.arange(div, div+combine_steps)

        divf_ind0_plus_one = (yrf)*self.ndiv+(div+combine_steps)

        # make sure that all divf units (based on combine_step)
        # would be captured by the list
        if divf_ind0_plus_one+(combine_steps-1)+1 <= len(list0):
            divf = (yrf)*self.ndiv+(div+combine_steps)
        else:
            divf = (yrf-1)*self.ndiv+(div+combine_steps)

        # if yrf_div >= (div+2*combine_steps):
        #   divf = divf+self.ndiv
        #   # print(divf)

        indf = np.arange(divf, divf+combine_steps)
        nstep = int((divf-div-1)/self.ndiv)
        # print(nstep)
        # print(len(ind0))
        all_ind0 = []
        for j in range(nstep+1):
            for g in range(combine_steps):
                all_ind0.append(div+self.ndiv*j+g)
        return ind0, indf, all_ind0, nstep
