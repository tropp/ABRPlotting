#!/usr/bin/python
from __future__ import print_function
import sys
import os
import re
import os.path
import glob
from collections import OrderedDict
import pandas as pd
import numpy as np
import scipy.signal
import matplotlib.pyplot as mpl
import itertools
import seaborn as sns  # makes plot background light grey with grid, no splines. Remove for publication plots
import peakdetect  # from Brad Buran's project, but cloned and modified here
from matplotlib.backends.backend_pdf import PdfPages


from getcomputer import getcomputer # stub to return the computer and base directory
from ABR_Datasets import ABR_Datasets # just the dict describing the datasets

basedir, computer_name = getcomputer()

class ABR():
    """
    Read an ABR data set from the matlab program
    Plot the traces.
    
    Parameters
    ----------
    datapath : str
        Path to the datasets. The data sets are expected to be collected under this
        path into individual directories, each with the results of the ABR runs for
        one subject.
    mode : string ('clicks' or 'tones')
        Specify type of data in the data set
    info : dict
        Dictionary with the following keys:
            invert : boolean (default : False)
                Changes sign of waveform if True
            minlat : float (default 0.75)
                minimum latency for event detection
            term : float (default '\r')
                line terminator (may change if the data has been read by an editor or is
                on a windows vs. mac vs linux system)
    """
    
    def __init__(self, datapath, mode='clicks',  info={'invert': False, 'minlat': 0.75, 'term': '\r'}):

        # Set some default parameters for the data and analyses
        
        self.sample_rate = 1e-5 # standard interpolated sample rate for matlab abr program: 10 microsecnds
        self.sample_freq = 1./self.sample_rate
        self.hpf = 300.
        self.lpf = 1500.  # filter frequencies, Hz
        self.clickdata = {}
        self.tonemapdata = {}
        # build color map where each SPL is a color (cycles over 12 levels)
        bounds = np.linspace(0, 120, 25)  # 0 to 120 db inclusive, 5 db steps
        color_labels = np.unique(bounds)
        rgb_values = sns.color_palette("Set2", 25)
        # Map label to RGB
        self.color_map = dict(zip(color_labels, rgb_values))
        
        self.datapath = datapath
        self.term = info['term']
        self.minlat = info['minlat']
        self.invert = info['invert']  # data flip... depends on "active" lead connection to vertex (false) or ear (true).
        
        self.characterizeDataset()

        # build color map where each SPL is a color (cycles over 12 levels)
        bounds = np.linspace(0, 120, 25)  # 0 to 120 db inclusive, 5 db steps
        color_labels = np.unique(bounds)
        self.color_map = self.makeColorMap(25, color_labels)
        color_labels2 = range(25)
        self.summaryClick_color_map = self.makeColorMap(25, range(25))  # make a default map, but overwrite for true number of datasets
        self.psdIOPlot = False

    def characterizeDataset(self):
        """
        Look at the directory in datapath, and determine what datasets are present.
        A dataset consists of at least 3 files:
        yyyymmdd-HHMM-SPL.txt : the SPL levels in a given run, with year, month day hour and minute
        yyyymmdd-HHMM-{n,p}-[freq].txt : hold waveform data for a series of SPL levels
            For clicks, just click polarity - holds waveforms from the run
            typically there will be both an n and a p file because we use alternating polarities.
            For tones, there is also a number indicating the tone pip frequency used.
        
        The runs are stored in separate dictionaries for the click and tone map runs, including
        their SPL levels (and for tone maps, frequencies).
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        
        files = [f for f in os.listdir(self.datapath) if os.path.isfile(os.path.join(self.datapath, f))]        
        self.spls = self.getSPLs(files)
        self.freqs = self.getFreqs(files)
        # A click run will consist of an SPL, n and p file, but NO additional files.
        self.clicks = {}
        allfiles = glob.glob(self.datapath)
        for s in self.spls.keys():
            if len(glob.glob(os.path.join(self.datapath, s + '*.txt'))) == 3:
                self.clicks[s] = self.spls[s]

        # inspect the directory and get a listing of all the tone and click maps
        # that can be found.
        self.tonemaps = {}
        for i, f in enumerate(self.freqs.keys()):
            self.tonemaps[f] = {'stimtype': 'tonepip', 'Freqs': self.freqs[f], 'SPLs': self.spls[f[:13]]}
            if mode == 'tones':
                if i == 0:
                    print( '  Frequency maps: ')
                print ('    ', f, self.tonemaps[f])

        self.clickmaps = {}
        for i, s in enumerate(self.clicks.keys()):
            self.clickmaps[s] = {'stimtype': 'click', 'SPLs': self.spls[s]}
            if mode == 'clicks':
                if i == 0:
                    print ('\n  Click Intensity Runs: ')
                print( '    Run: {:s}'.format(s))
                print( '        {:s}'.format(self.clickmaps[s]))

    def makeColorMap(self, N, labels):
        """
        create a color map of N levels using the specified labels
        
        Parameters
        ----------
        N : int
            Number of color levels
        labels : tuples
            list of tuples of rgb values corresponding to color levels
        
        Returns
        -------
        dict for colormap
        """
        rgb_values = sns.color_palette("Set2", N)
        return dict(zip(labels, rgb_values))

    def getClickData(self, select):
        """
        Gets the click data for the current selection
        The resulting data is held in a dictionary structured as
        {mapidentity: dict of {waves, time, spls and optional marker}
        
        Parameters
        ----------
        select : which data to select (see main section at end of code)
        
        """
        select, freqs = self.adjustSelection(select)
        # get data for clicks and plot all on one plot
        self.clickdata = {}
        for i, s in enumerate(self.clickmaps.keys()):
            if select is not None:
                if s[9:] not in select:
                    continue
            if s[0:8] == '20170419':  # these should not be here... should be in an excel table 
                smarker = 'go-'
            elif s[0:8] in ['20170608', '20170609']:
                smarker = 'bs-'
            else:
                smarker = 'kx-'
            waves = self.get_combineddata(s, self.clickmaps[s])
            if waves is None:
                print('Malformed data set for run %s. Continuing' % s)
                continue
            waves = waves[::-1]  # reverse order to match spls
            t = np.linspace(0, waves.shape[1]*self.sample_rate*1000., waves.shape[1])
            spls = self.clickmaps[s]['SPLs']  # get spls
            self.clickdata[s] = {'waves': waves, 'timebase': t, 'spls': spls, 'marker': smarker}

    def getToneData(self, select):
        """
        Gets the tone map data for the current selection
        The resulting data is held in a dictionary structured as
        {mapidentity: dict of frequencies}
        Where each dictoffrequencies holds a dict of {waves, time, spls and optional marker}
        
        """
        select, freqs = self.adjustSelection(select)
        self.tonemapdata = {}
        # convert select to make life easier 
        # select list should have lists of strings ['0124', '0244'] or Nones...
        select, freqs = self.adjustSelection(select, tone=True)
        for i, s in enumerate(self.tonemaps.keys()):
            freqs = []
            if select is not None:
                if s[9:] not in select:
                    continue
            for f in self.tonemaps[s]['Freqs']:
                if f not in freqs:
                    freqs.append(f)
            freqs.sort()
            if len(freqs) == 0:  # check of no tone pip ABR data in this directory
                continue
            # now we can build the tonemapdata
            self.tonemapdata[s] = OrderedDict()
            for fr in self.tonemaps[s]['Freqs']:
                waves = self.get_combineddata(s, self.tonemaps[s], freq=fr)
                if waves is None:
                    print('Malformed data set for run %s. Continuing' % s)
                    continue
                t = np.linspace(0, waves.shape[1]*self.sample_rate*1000., waves.shape[1])
                spls = self.tonemaps[s]['SPLs']
                self.tonemapdata[s][fr] = {'waves': waves, 'timebase': t, 'spls': spls, 'marker': 'ko-'}

    def plotClicks(self, select=None, plottarget=None, IOplot=None, PSDplot=None, superIOPlot=None, colorindex=0):
        """
        Plot the click ABR intensity series, one column per subject, for one subject
        
        Parameters
        ----------
        select : list of str (default : None)
            A list of the times for the datasets to plot for each subject.
            If None, then all detected click ABR runs for that subject are superimposed.
            The list looks like [['0115'], None, None]  - for each directory in order. 
        
        plottarget : matplotlib axis object
            The axis to plot the data into.
        
        IOPlot : Matplotlib Axes object
            Input output plot target. If not None, then use the specified Mabplotlib Axes for the IO plot
        
        PSDPlot : Matplotlib Axes object
            Power spectral density plot target. If not None, then use the specified Mabplotlib Axes for the plot
        
        superIOPlot : Matplotlib Axes object
            Input output plot target. If not None, then use the specified Mabplotlib Axes for the IO plot

        """
        # get data for clicks and plot all on one plot
        A = Analyzer()
        thrs = {}
        icol = colorindex
        for s in self.clickdata.keys():
            datatitle = os.path.basename(os.path.normpath(self.datapath)) + '\n' + s
            #datatitle = datatitle.replace('_', '\_')  # if TeX is enabled, will need to escape the underscores
            waves = self.clickdata[s]['waves']
            t = self.clickdata[s]['timebase']
            spls = self.clickdata[s]['spls']
            
            A.analyze(t, waves)
            r, n = A.p1n1
            halfspl = np.max(spls)/2.
            latmap = []
            spllat = []
            for j in range(len(spls)):
                if spls[j] > halfspl:
                    latmap.append(t[r[j][0]]) # get latency for first value
                    spllat.append(spls[j])
            latp = np.polyfit(spllat, latmap, 1)
            fitline = np.polyval(latp, spls)
            thr_spl = A.threshold_spec(waves, spls, SD=3.5)
            thrs[s] = thr_spl
            linewidth=1.0
            IO = np.zeros(len(spls))
            for j in range(len(spls)):
                if spls[j] == thr_spl:
                    plottarget.plot(t, 0*waves[j]*1e6+spls[j], color=[0.5, 0.5, 0.5, 0.4], linewidth=5)
                try:
                    c = self.color_map[spls[j]]
                    plottarget.plot(t, 2*waves[j]*1e6+spls[j], color=c, linewidth=linewidth)
                except:
                    plottarget.plot(t, 2*waves[j]*1e6+spls[j], color='k', linewidth=linewidth)
                for p in r[j]:
                    plottarget.plot(t[p], 2*waves[j][p]*1e6+spls[j], 'ro', markersize=3.5)
                for p in n[j]:
                    plottarget.plot(t[p], 2*waves[j][p]*1e6+spls[j], 'bo', markersize=3.5)
                plottarget.plot(fitline, spls, 'g-', linewidth=0.7)
                if spls[j] >= thr_spl:
                    IO[j] = 1.e6*(waves[j][r[j][0]] - waves[j][n[j][0]])
                else:
                    ti = int(fitline[j]/(self.sample_rate*1000.))
                    IO[j] = 1.e6*waves[j][ti]
            
            if superIOPlot is not None:  # superimposed IO plots
                datatitle_short = os.path.basename(os.path.normpath(self.datapath))
                superIOPlot.plot(spls, 1e6*A.ppio, self.clickdata[s]['marker'], color=self.summaryClick_color_map[icol], label=datatitle_short)
                
                # print out the data for import into another plotting program, such as Prism or Igor
                print( '*'*20)
                print('dataset: ', s)
                print('t\tV')
                for i in range(len(spls)):
                    print('%.1f\t%.3f' % (spls[i], IO[i]))
                print('*'*20)
                
            if IOplot is not None:  # generic io plot for cell
                IOplot.set_title(datatitle, {'fontsize': 7})  # directory plus file
                IOplot.plot(spls, 1e6*A.ppio, marker=A.ppioMarker, color=self.summaryClick_color_map[icol], label='P-P')
                IOplot.plot(spls, 1e6*A.rms_response, marker=A.rmsMarker, color=self.summaryClick_color_map[icol], label='RMS signal')
                IOplot.plot(spls, 1e6*A.rms_baseline, marker=A.baselineMarker, color='k', label='RMS baseline')
                if self.psdIOPlot:
                    ax2 = IOplot.twinx()
                    ax2.plot(spls, A.psdwindow, marker=A.psdMarker, color='r', label='PSD signal')
                    ax2.tick_params('y', colors='r')
                handles, labels = IOplot.get_legend_handles_labels()
                legend = IOplot.legend(loc='upper left')
                for label in legend.get_texts():
                    label.set_fontsize(6)
                
            if PSDplot is not None:  # power spectral density
                for j in range(len(spls)):
                    PSDplot.semilogy(np.array(A.fr), np.array(A.psd[j])) # color=self.color_map[spls])
                PSDplot.set_ylim(1e-6, 0.01)
                PSDplot.set_xlim(100., 2000.)
        plottarget.set_xlim(0, 8.)
        plottarget.set_ylim(10., 115.)
        plottarget.title.set_text(datatitle)
        plottarget.title.set_size(7)
        if superIOPlot is not None:
            legend = superIOPlot.legend(loc='upper left')
            for label in legend.get_texts():
                label.set_fontsize(7)
            
        print("")
        for s in thrs.keys():
            print ("dataset: %s  thr=%.0f" % (s, thrs[s]))
        print ("")
        
    
    def plotTones(self, select=None, pdf=None):
        """
        Plot the tone ABR intensity series, one column per frequency, for one subject
        
        Parameters
        ----------
        select : list of str (default : None)
            A list of the times for the datasets to plot for each subject.
            If None, then all detected tone ABR runs for that subject/frequency
            comtination are superimposed.
            The list looks like [['0115'], None, None]  - for each directory in order. 
        
        pdf : pdfPages object
            The pdfPages object that the plot will be appended to. Results in a multipage
            pdf, one page per subject. 

        """
        self.thrs = {}  # holds thresholds for this dataset
        freqs = []
        for run in self.tonemapdata.keys():  # just passing one dict... 
            freqs.extend(self.tonemapdata[run].keys())
        freqs.sort()
        if len(freqs) == 0:  # check of no tone pip ABR data in this directory
            return

        f, axarr = mpl.subplots(1, len(freqs), figsize=(12,6), num="Tones")
        datatitle = os.path.basename(os.path.normpath(self.datapath))[15:]
        datatitle = datatitle.replace('_', '\_')
        f.suptitle(datatitle)
        A = Analyzer()
        for i, s in enumerate(self.tonemapdata.keys()):
            thr_spls = np.zeros(len(self.tonemaps[s]['Freqs']))
            for k, fr in enumerate(self.tonemaps[s]['Freqs']):  # next key is frequency
                waves = self.tonemapdata[s][fr]['waves']
                t = self.tonemapdata[s][fr]['timebase']
                spls = self.tonemapdata[s][fr]['spls']
                A.analyze(t, waves)
                thr_spl = A.threshold_spec(waves, spls, SD=3.5)
                thr_spls[k] = thr_spl
                plottarget = axarr[freqs.index(fr)]
                for j in range(len(spls))[::-1]:
                    if spls[j] == thr_spl:
                        plottarget.plot(t, 0*waves[j]*1e6+spls[j], color=[0.5, 0.5, 0.5, 0.4], linewidth=5)
                    plottarget.plot(t, 2*waves[j]*1e6+spls[j], color=self.color_map[spls[j]])
                plottarget.set_xlim(0, 8.)
                plottarget.set_ylim(10., 110.)
                frtitle = '%.1f kHz' % (float(fr)/1000.)
                plottarget.title.set_text(frtitle)
                plottarget.title.set_size(9)
            self.thrs[s] = [self.tonemaps[s]['Freqs'], thr_spls] 
        self.cleanAxes(axarr)
        if pdf is not None:
            pdf.savefig()
            mpl.close()

    def plotToneThresholds(self, allthrs, num):
        """
        Make a nice plot of the tone thresholds for all of the datasets
        Data are plotted against a log frequency scale (2-64kHz)
        Data is plotted into the current figure.
        
        Parameters
        ----------
        allthrs : dict
        A dictionary holding all the threshold information. The following
        structure is required:
            Keys: filenames for each dataset
            Values a dict of thresholds. The keys are the names of the tone maps
            (because more than one tone map may be combined)
            The values are tuples of (frequency, threshold)
        
        Returns
        -------
        Nothing
        """
        f, ax = mpl.subplots(nrows=1, ncols=1, num=num)
        ax.set_xscale('log', nonposx='clip', base=2)
        n_datasets = len(allthrs.keys())
        c_map = self.makeColorMap(n_datasets, allthrs.keys())
        
        thrfrs = {}
        for i, d in enumerate(allthrs):  # for all the datasets
            for m in allthrs[d]:  # for all the maps in the dataset combined
                ax.scatter(np.array(allthrs[d][m][0])/1000., allthrs[d][m][1], color=c_map[d], s=12)
                for j, f in enumerate(allthrs[d][m][0]):
                    if f not in thrfrs.keys():
                        thrfrs[f] = [allthrs[d][m][1][j]]
                    else:
                        thrfrs[f].append(allthrs[d][m][1][j])
        # sort the threshold list
        thrs_sorted = OrderedDict(sorted(thrfrs.items(), key=lambda t: t[0]))
        frmean = np.zeros(len(thrs_sorted.keys()))
        frstd = np.zeros(len(thrs_sorted.keys()))

        for i, f in enumerate(thrs_sorted):
            print ('i, f: ', i, f)
            print (thrs_sorted[f])
            frmean[i] = np.nanmean(thrs_sorted[f])
            frstd[i] = np.nanstd(thrs_sorted[f])
        ax.errorbar(np.array(thrs_sorted.keys())/1000., frmean, yerr=frstd, fmt='o')
        ax.set_xlim(1.8, 65.)
        xt = [2., 4., 8., 16., 32., 64.]
        mpl.xticks(xt, [str(x) for x in xt])

    def adjustSelection(self, select, tone=False):
        freqs = []
        if select is None:
            return select, freqs
        for i, s in enumerate(select):
            if s is None:
                continue
            if isinstance(s, int):
                select[i] = '%04d' % s  # convert to string
            if isinstance(s, str) and len(s) < 4:  # make sure size is ok
                select[i] = '%04d' % int(s)
            if tone:
                base = self.tonemaps.keys()[0][:9]
                for f in self.tonemaps[base + select[i]]['Freqs']:
                    if f not in freqs:
                        freqs.append(f)
        return select, freqs

    def cleanAxes(self, axl):
        """
        Remove axis splines on top and right
        
        Parameters
        ----------
        axl : list of axes objects
        
        """
        if type(axl) is not list:
            axl = [axl]
        axl = list(itertools.chain(*axl))
        for ax in axl:
            if ax is None:
                continue
            for loc, spine in ax.spines.iteritems():
                if loc in ['left', 'bottom']:
                    pass
                elif loc in ['right', 'top']:
                    spine.set_color('none')
                    # do not draw the spine
                else:
                    raise ValueError('Unknown spine location: %s' % loc)
                # turn off ticks when there is no spine
                ax.xaxis.set_ticks_position('bottom')
                ax.yaxis.set_ticks_position('left')  # stopped working in matplotlib 1.10

    def getSPLs(self, dir):
        """
        Return all the spl files in the directory. There is one spl file
        per intensity run.
        """
        spls = [f[:-8] for f in os.listdir(self.datapath) if f[14:17] == 'SPL']
        rundict = {}
        for run in spls:
             SPLfile = open(os.path.join(self.datapath, run+'-SPL.txt'), 'r')
             rundict[run] = [float(spl) for spl in SPLfile]
        return rundict
 
    def getFreqs(self, dir):
        """
        Return all the tonepip files in the directory. There is one tonepip file
        per frequency for its intensity run.
        
        """
        # return all the tone response files in the directory.
        freqs = [f[:-8] for f in os.listdir(self.datapath) if f[14:17] == 'kHz']
        rundict = {}
        for run in freqs:
          Freqfile = open(os.path.join(self.datapath, run+'-kHz.txt'), 'r')
          
          rundict[run] = [float(khz) for khz in Freqfile if khz[0] != '\t'] # handle old data with blank line at end
        return rundict

    def get_combineddata(self, datasetname, dataset, freq=None):
        """
        Read the data sets and combine the p (condensation) and n
        (rarefaction) data sets for alternating polarity stimuli.
        
        Parameters
        ----------
        datasetname : str
            yyyymmdddd-time format for start of dataset name
        dataset : dict
            dictionary for the dataset that identifies the stimulus
            type, the SPLs in the dataset, and the frequency if the dataset
            is a tone pip run
        freq : float (default: None)
            for tone maps, the specific frequency intensity series to return
        
        Returns
        -------
        waves : numpy array
            Waveforms, as a nxm array, where n is the number of intensities,
            and m is the length of each waveform
        
        """
#        print dataset
        if dataset['stimtype'] == 'click':
            fnamepos = datasetname + '-p.txt'
            fnameneg = datasetname + '-n.txt'
            waves = self.read_dataset(fnamepos, fnameneg)
            return waves
        if dataset['stimtype'] == 'tonepip':
            fnamepos = datasetname + '-p-%.3f.txt' % freq
            fnameneg = datasetname + '-n-%.3f.txt' % freq
            waves = self.read_dataset(fnamepos, fnameneg)
            return waves

    def read_dataset(self, fnamepos, fnameneg):
        """
        Read a dataset, combining the positive and negative recordings,
        which are stored in separate files on disk. The waveforms are averaged
        which helps to minimize the CAP contribution. 
        The waveforms are then bandpass filtered to remove the low frequency
        "rumble" and excess high-frequency noise.
        
        Parameters
        ----------
        fnamepos : str
             name of the positive (condensation) file
        fnameneg : str
            name of the negative (rarefaction) file
        
        Returns
        -------
        waveform
            Waveform, as a nxm array, where n is the number of intensities,
            and m is the length of each waveform
         
        """

        posf = pd.io.parsers.read_csv(os.path.join(self.datapath, fnamepos), delim_whitespace=True, lineterminator='\r',
            skip_blank_lines=True, header=0)
        negf = pd.io.parsers.read_csv(os.path.join(self.datapath, fnameneg), delim_whitespace=True, lineterminator='\r',
            skip_blank_lines=True, header=0)
        negseries=[]
        posseries=[]
        for col in negf.columns:
            for i in range(len(negf)):
                negseries.append(negf[col][i])

        for col in posf.columns:
            for i in range(len(posf)):
                posseries.append(posf[col][i])
        
        wvfmdata = [(x + y)/2 for x, y in zip(negseries, posseries)]
#        wvfmdata = negseries  # just get one polarity
        d1 = len(wvfmdata)/len(posf[col])
        waves = np.reshape(wvfmdata,(d1,len(posf[col])))
        for i, w in enumerate(waves):
            waves[i, -1] = waves[i, -2]  # remove nan from end of waveform...
            waves[i,:] = self.filter(waves[i,:], 4, self.lpf, self.hpf, samplefreq=self.sample_freq)
            if self.invert:
                waves[i,:] = -waves[i,:]
        return waves


    def filter(self, data, order, lowpass, highpass, samplefreq, ftype='butter'):
        """
        Returns waveform filtered using filter paramters specified. Since
        forward and reverse filtering is used to avoid introducing phase delay,
        the filter order is essentially doubled.
        
        Parameters
        ----------
        data : m numpy array of floats
            the data waveforms, as a 1-d array
        
        order : filter order.
        
        lowpass : float
            Low pass filter setting, in Hz
        
        highpass : float
            High pass filter setting, in Hz
        
        samplefreq : float
            sample frequency, in Hz (inverse of sample rate)
        
        ftype : str (default: 'butter')
            Type of filter to implement. Default is a Butterworth filter.
        
        Returns
        -------
        signal : numpy array of floats
            the filtered signal waveform.
        """

        Wn = highpass/(samplefreq/2.), lowpass/(samplefreq/2.)
        kwargs = dict(N=order, Wn=Wn, btype='band', ftype=ftype)
        b, a = scipy.signal.iirfilter(output='ba', **kwargs)
        zpk = scipy.signal.iirfilter(output='zpk', **kwargs)
        try:
            self._zpk.append(zpk)
        except:
            self._zpk = [zpk]
        self.signal = scipy.signal.filtfilt(b, a, data, padlen=int(len(data)/10))
        return self.signal
        
    def SignalFilter(self, signal, LPF, HPF, samplefreq, debugFlag=True):
        """Filter signal within a bandpass with elliptical filter.

        Digitally filter a signal with an elliptical filter; handles
        bandpass filtering between two frequencies. 

        Parameters
        ----------
        signal : array
            The signal to be filtered.
        LPF : float
            The low-pass frequency of the filter (Hz)
        HPF : float
            The high-pass frequency of the filter (Hz)
        samplefreq : float
            The uniform sampling rate for the signal (in seconds)

        Returns
        -------
        w : array
            filtered version of the input signal
        """
        print ('nans: ', np.argwhere(np.isnan(signal)))
        if debugFlag:
            print ("sfreq: %f LPF: %f HPF: %f" % (samplefreq, LPF, HPF))
        flpf = float(LPF)
        fhpf = float(HPF)
        sf = float(samplefreq)
        sf2 = sf/2.
        wp = [fhpf/sf2, flpf/sf2]
        ws = [0.5*fhpf/sf2, 2*flpf/sf2]
        if debugFlag:
            print( "signalfilter: samplef: %f  wp: %f, %f  ws: %f, %f lpf: %f  hpf: %f" % (
               sf, wp[0], wp[1], ws[0], ws[1], flpf, fhpf))
        filter_b, filter_a = scipy.signal.iirdesign(wp, ws,
                gpass=1.0,
                gstop=60.0,
                ftype="ellip")
        msig = np.nanmean(signal)
        signal = signal - msig
        w = scipy.signal.lfilter(filter_b, filter_a, signal) # filter the incoming signal
        signal = signal + msig
        if debugFlag:
            print ("sig: %f-%f w: %f-%f" % (np.amin(signal), np.amax(signal), np.amin(w), np.amax(w)))
        return(w)


class Analyzer(object):
    """
    Provide analysis functions for ABRs. 
    """
    def __init__(self):
        self.ppioMarker = 's'
        self.rmsMarker = 'o'
        self.psdMarker = '*'
        self.baselineMarker = '+'
    
    def analyze(self, timebase, waves):
        
        self.sample_rate = 0.001*(timebase[1]-timebase[0])
        self.sample_freq = 1./self.sample_rate
        self.waves = waves
        self.timebase = timebase
        
        response_range = [2.5, 7] # window, in msec... 
        baseline = [0., 2.0]

        self.P1N1()
        self.ppio = self.peaktopeak(response_range)
        self.rms_response = self.measure_rms(response_range)
        self.rms_baseline = self.measure_rms(baseline)
        self.specpower()
 
    def peaktopeak(self, tr):
        tx = self.gettimeindices(tr)
        pp = np.zeros(self.waves.shape[0])
        for i in range(self.waves.shape[0]):
            pp[i] = np.max(self.waves[i][tx]) - np.min(self.waves[i][tx])
        return pp
    
    def P1N1(self):
        """
        Use Buran's peakdetect routine to find the peaks and returan a list
        of peaks. Works twice - first run finds all the positive peaks, and
        the second run finds the negative peaks that follow the positive peaks
        Note that the peaks may not be "aligned" in the sense that it is possible
        to find two positive peaks in succession without a negative peak. 
        """
        r = {}
        n = {}
        for j in range(self.waves.shape[0]):
            r[j] = peakdetect.find_np(self.sample_freq, self.waves[j,:],
                nzc_algorithm_kw={'dev': 1},
                guess_algorithm_kw={'min_latency': 2.5})
            if len(r[j]) > 0:
                n[j] = peakdetect.find_np(self.sample_freq, -self.waves[j,:],
                    nzc_algorithm_kw={'dev': 2.5},
                    guess_algorithm_kw={'min_latency': self.timebase[r[j][0]]})  # find negative peaks after positive peaks        
        self.p1n1=(r, n)

    def measure_rms(self, tr):
        tx = self.gettimeindices(tr)
        rms = np.zeros(self.waves.shape[0])
        for i in range(self.waves.shape[0]):
            rms[i] = np.std(self.waves[i][tx])
        return rms

    def gettimeindices(self, tr):
        x, = np.where((tr[0] <= self.timebase) & (self.timebase < tr[1]))
        return x
        
    def specpower(self, fr=[500., 1500.], win=[0, -1]):
        fs = 1./self.sample_rate
        # fig, ax = mpl.subplots(1, 1)
        psd = [None]*self.waves.shape[0]
        psdwindow = np.zeros(self.waves.shape[0])
        for i in range(self.waves.shape[0]):
            f, psd[i] = scipy.signal.welch(1e6*self.waves[i][win[0]:win[1]],
                fs, nperseg=256, nfft=8192, scaling='density')
            frx, = np.where((f >= fr[0]) & (f <= fr[1]))
            psdwindow[i] = np.nanmean(psd[i][frx[0]:frx[-1]])
#            ax.semilogy(f, psd[i])
        # ax.set_ylim([0.1e-4, 0.1])
        # ax.set_xlim([10., 2000.])
        # ax.set_xlabel('F (Hz)')
        # ax.set_ylabel('PSD (uV^2/Hz)')
        # mpl.show()
        self.fr = f
        self.psd = psd
        self.psdwindow = psdwindow
        return psdwindow

    def thresholds(self, waves, spls, tr=[0., 8.], SD=4.0):
        """
        Auto threshold detection:
        BMC Neuroscience200910:104  DOI: 10.1186/1471-2202-10-104
        Use last 10 msec of 15 msec window for SD estimates
        Computes SNR (max(abs(signal))/reference SD) for a group of traces
        The reference SD is the MEDIAN SD across the intensity run.
        Th
        
        """
        refwin = self.gettimeindices([15., 25.])
        sds = np.std(waves[:,refwin[0]:refwin[-1]], axis=1)
        self.median_sd = np.nanmedian(sds)
        tx = self.gettimeindices(tr)
        self.max_wave = np.max(np.fabs(waves[:, tx[0]:tx[-1]]), axis=1)
        thr, = np.where(self.max_wave >= self.median_sd*SD)  # find criteria threshold
        thrx, = np.where(np.diff(thr) == 1)  # find first contiguous point (remove low threshold non-contiguous)
        if len(thrx) > 0:
            return spls[thr[thrx[0]]]
        else:
            return np.nan

    def threshold_spec(self, waves, spls, tr=[0., 8.], SD=4.0):
        """
        Auto threshold detection:
        BMC Neuroscience200910:104  DOI: 10.1186/1471-2202-10-104
        Use last 10 msec of 15 msec window for SD estimates
        Computes SNR (max(abs(signal))/reference SD) for a group of traces
        The reference SD is the MEDIAN SD across the intensity run.
        
        MODIFIED version: criteria based on power spec
        
        """
        refwin = self.gettimeindices([15., 25.])
        sds = self.specpower(fr=[800., 1250.], win=[refwin[0], refwin[-1]])
        self.median_sd = np.nanmedian(sds)
        tx = self.gettimeindices(tr)
        self.max_wave = self.specpower(fr=[800., 1250.], win=[tx[0], tx[-1]])
        thr, = np.where(self.max_wave >= self.median_sd*SD)  # find criteria threshold
        if len(thr) > 0:
            return spls[thr[0]]
        else:
            return np.nan


class ABRWriter():
    def __init__(self, dataset, waves, filename, freqs, spls, fileindex=None):
        """
        Write one ABR file for Brad Buran's analysis program, in "EPL"
        format.
        
        TODO : fix this so that the ABR-files have names and check that it works.
        
        ABR-click-time  indicate click intenstiry series run in this
            directory, use time to distinguish file
            Example: ABR-click-0901.abr
        ABR-tone-time-freq indicate tone run number in this directory
            with frequency in Hz (6 places).
            Example: ABR-tone-0945-005460.dat
        
        """
        self.freqs = freqs
        self.spls = spls
        twaves = waves.T
        abr_filename = ('ABR-{0:s}-{1:4s}'.format(dataset['stimtype'], filename[8:12]))
        if dataset['stimtype'] == 'tonepip':
             abr_filename = abr_filename + ('{0:06d}'.format(dataset['Freqs'][fileindex]))
        abr_filename = abr_filename + '.dat' # give name an extension
        np.savetxt(abr_filename, twaves, delimiter='\t ', newline='\r\n')
        header = self.make_header(filename, fileindex)
        self.rewriteheader(abr_filename, filename, header)
        return abr_filename

    def make_header (self, filename, fileindex):
        """
        Create a (partially faked) header string for Buran's ABR analysis program, using data from our
        files.
        This mimics the Eaton Peabody data header style.
        
        Parameters
        ----------
        filename : str
            Name of the file whose information from the freqs and spl dictionaries
            will be stored here
        fileindex : int
            Index into the dict for the specific frequency.
        
        Returns
        -------
        header : str
            the whole header, ready for insertion into the file.
        
        """
        
        header1 = ":RUN-3	LEVEL SWEEP	TEMP:207.44 20120427 8:03 AM	-3852.21HR: \n"
        header2 = ":SW EAR: R	SW FREQ: " + "%.2f" % self.freqs[filename][fileindex]
        header2 += "	# AVERAGES: 512	REP RATE (/sec): 40	DRIVER: Starship	SAMPLE (usec): 10 \n"
        header3 = ":NOTES- \n"
        header4 = ":CHAMBER-412 \n"
        levs = self.spls[filename] # spls[ABRfiles[i]]
        levels = ['%2d;'.join(int(l)) for l in levs]
        header5 = ":LEVELS:" + levels + ' \n'
#            header5 = ":LEVELS:20;25;30;35;40;45;50;55;60;65;70;75;80;85;90; \n"
        header6 = ":DATA \n"
        header = header1 + header2 + header3 + header4 + header5 + header6
        return header

    def rewriteheader(self, filename, header):
        """
        Write the file header before the rest of the data, replacing any
        previous header.

        Parameters
        ----------
        filename : str
            name of the file to which the header will be prepended

        header : str
            The header string to prepend

        """
        ABRfile = open(filename, 'r')

        # remove any existing 'header' from the file, in case there are duplicate header rows in the wrong places
        datalines = []
        for line in ABRfile:
            if line[0] == ':':
                continue
            datalines.append(line)
        ABRfile.close()

        # now rewrite the file, prepending the header above the data
        ABRfile = open(filename, 'w')
        ABRfile.write(''.join(header))
        ABRfile.write('\n' .join(datalines))
        ABRfile.write(''.join('\r\r'))
        ABRfile.close()


if __name__ == '__main__':

    if len(sys.argv) > 1:
        dsname = sys.argv[1]
        mode = sys.argv[2]
    else:
        print ('Missing command arguments; call: plotABRs.py datasetname [click, tone]')
        exit(1)
    if dsname not in ABR_Datasets.keys():
        print( ABR_Datasets.keys())
        raise ValueError('Data set %s not found in our list of known datasets')
    if mode not in ['tones', 'clicks']:
        raise ValueError('Second argument must be tones or clicks')

    top_directory = os.path.join(basedir, ABR_Datasets[dsname]['dir'])
    
    dirs = [tdir for tdir in os.listdir(top_directory) if os.path.isdir(os.path.join(top_directory, tdir))]
    print( 'found dirs: ', dirs)
    if mode == 'clicks':
        if 'clickselect' in ABR_Datasets[dsname].keys():
            clicksel = ABR_Datasets[dsname]['clickselect']
        else:
            clicksel = [None]*len(dirs)
        rowlen = 8.
        m = int(np.ceil(len(clicksel)/rowlen))
        if m == 1:
            n = len(clicksel)
        else:
            n = int(rowlen)
        if m > 1:
            h = 4*m
        else:
            h = 5
        n = 2
        f, axarr = mpl.subplots(m, n, figsize=(12, h), num='Click Traces')
        f2, axarr2 = mpl.subplots(m, n, figsize=(12, h), num='Click IO Summary')
#        f3, axarr3 = mpl.subplots(m, n, figsize=(12, h))
        f4, IOax = mpl.subplots(1, 1, figsize=(6,6), num='Click IO Overlay')
        if axarr.ndim > 1:
            axarr = axarr.ravel()
        if axarr2.ndim > 1:
            axarr2 = axarr2.ravel()
        fofilename = os.path.join(top_directory, 'ClickSummary.pdf')
        nsel = len(clicksel)
        print ('Nsel: ', nsel)
        for icol, k in enumerate(range(nsel)):
            P = ABR(os.path.join(top_directory, dirs[k]), mode, info=ABR_Datasets[dsname])
            if icol == 0:
                P.summaryClick_color_map = P.makeColorMap(nsel, range(nsel))
            print ('icol: ', icol)
            P.getClickData(select=clicksel[k]) 
            P.plotClicks(select=clicksel[k], plottarget=axarr[k], superIOPlot=IOax,
                IOplot=axarr2[k], colorindex=icol)
        mpl.figure('Click Traces')
        mpl.savefig(fofilename)
        mpl.figure('Click IO Summary')
        fo2filename = os.path.join(top_directory, 'ClickIOSummary.pdf')
        mpl.savefig(fo2filename)
        mpl.figure('Click IO Overlay')
        fo4filename = os.path.join(top_directory, 'ClickIOOverlay.pdf')
        mpl.savefig(fo4filename)
        
    if mode == 'tones':
        # first one has issues: 1346 stops at 40dbspl for all traces
        # 1301+1327 seem to have the higher levels, but small responses; remainder have fragementary frequencies
        # 0901 has whole series (but, is it hte same mouse?)
        #tonesel = [['0901'], ['1446', '1505'], None, None, None, None, None, None]
        tonesel = [None]*len(dirs)
        fofilename = os.path.join(top_directory, 'ToneSummary.pdf')
        allthrs = {}
        with PdfPages(fofilename) as pdf:
            for k in range(len(tonesel)):
                P = ABR(os.path.join(top_directory, dirs[k]), mode)
                P.getToneData(select=tonesel[k]) 
                P.plotTones(select=tonesel[k], pdf=pdf)
                allthrs[dirs[k]] = P.thrs
        P.plotToneThresholds(allthrs, num='Tone Thresholds')
        tthr_filename = os.path.join(top_directory, 'ToneThresholds.pdf')
        mpl.savefig(tthr_filename)
    mpl.show()
