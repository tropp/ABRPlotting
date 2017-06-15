import sys
import os
import re
import os.path
import commands
import glob
from collections import OrderedDict
# from optparse import OptionParser
import pandas as pd
import numpy as np
import scipy.signal
import matplotlib.pyplot as mpl
#import matplotlib.colors as colors
import itertools
import seaborn as sns  # makes plot background light grey with grid, no splines. Remove for publication plots
import peakdetect  # from Brad Buran's project

from matplotlib.backends.backend_pdf import PdfPages

computer_name = commands.getoutput('scutil --get ComputerName')
if computer_name == 'Lytle':
    basedir = '/Volumes/Pegasus/ManisLab_Data3/abr_data'
elif computer_name == 'Tamalpais':
    basedir = '/Users/pbmanis/Desktop/data/ABR_DATA'
else:
    raise ValueError("Need valid computer name to set base path to data")
    
#current 10 June 2017
# This table should be converted to an excel spreadsheet organized
# with a worksheet for each data set, and within each worksheet, each subject's information
# placed in a single row. Make the table pandas-readable.

ABR_Datasets = {'NrCAM': {'dir': "Tessa/Coate-NrCAM-ABRs/", 'invert': True, 
                    'clickselect: ': [['0849'], None, None, None, None, None, None, None],
                    'toneselect': [['0901'], ['1446', '1505'], None, None, None, None, None, None],
                    'term': '\r', 'minlat': 2.2},
                 'TessaCBA': {'dir': 'Tessa/CBA', 'invert': True},
                 'TessaCBANE': {'dir': 'Tessa/CBA_NoiseExposed', 'invert': True},
                 'CNTNAP2X': {'dir': 'Tessa/CNTNAP2', 'invert': False},
                 'CNTNAP2Het': {'dir': 'Tessa/CNTNAP2_Het', 'invert': False},
                 'CNTNAP2HG': {'dir': 'Tessa/CNTNAP2_Het_GP4.3', 'invert': False},
                 'CNTNAP2KO': {'dir': 'Tessa/CNTNAP2_KO', 'invert': False},
                 'CNTNAP2WT': {'dir': 'Tessa/CNTNAP2_WT', 'invert': False},
                 'GP43Norm': {'dir': 'Tessa/GP4.3-Thy1-Normal', 'invert': False},
                 'Yong': {'dir': "Yong\'s ABRs", 'invert': True,
                     'term':'\n', 'minlat': 0.6},
                 'Jun': {'dir': "JUN\'s ABRs", 'invert': False},
                 'Ruili': {'dir': "Ruilis ABRs", 'invert': True, 'nameselect': 'CBA',
                     'term':'\n', 'minlat': 0.6},
                 'Eveleen': {'dir': "Eveleen\'s ABRs", 'invert': False},
                 } 


class ABR():
    """
    Read an ABR data set from the matlab program
    Provides routines for reading data, parsing data type (click tones),
    filtering and plotting the traces.
    
    Parameters
    ----------
    datapath : str
        Path to the datasets. The data sets are expected to be collected under this
        path into individual directories, each with the results of the ABR runs for
        one subject.
    """
    
    def __init__(self, datapath, mode='clicks', info={'invert': False, 'minlat': 0.75, 'term': '\r'}):
        
        
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
        print 'Frequency maps: '
        self.tonemaps = {}
        print self.freqs
        for f in self.freqs.keys():
            self.tonemaps[f] = {'stimtype': 'tonepip', 'Freqs': self.freqs[f], 'SPLs': self.spls[f[:13]]}
            if mode == 'tones':
                print f, self.tonemaps[f]

        print '\n    Click Intensity Runs: '
        self.clickmaps = {}
        for s in self.clicks.keys():
            self.clickmaps[s] = {'stimtype': 'click', 'SPLs': self.spls[s]}
            if mode == 'clicks':
                print '    Run: %s' % s
                print '        ', self.clickmaps[s]

    def getClickData(self, select):
        """
        Gets the click data for the current selection
        The resulting data is held in a dictionary structured as
        {mapidentity: dict of {waves, time, spls and optional marker}
        """
        select, freqs = self.adjustSelection(select)
        # get data for clicks and plot all on one plot
        self.clickdata = {}
        for i, s in enumerate(self.clickmaps.keys()):
            if select is not None:
                if s[9:] not in select:
                    continue
            print "s: ", s
            if s[0:8] == '20170419':  # these should not be here... should be in excel table 
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
        
    def plotClicks(self, plottarget=None, IOplot=None, PSDplot=None, superIOPlot=None):
        """
        Plot the click ABR intensity series, one column per subject, for one subject
        requires self.clickdata be populated
        Parameters
        ----------
        select : list of str (default : None)
            A list of the times for the datasets to plot for each subject.
            If None, then all detected click ABR runs for that subject are superimposed.
            The list looks like [['0115'], None, None]  - for each directory in order. 
        
        plottarget : matplotlib axis object
            The axis to plot the data into.

        """
        if len(self.clickdata) == 0:
            raise ValueError('plotClicks requires that self.clickdata be populated first')
        
        A = Analyzer()  # create an instance of analyzer to use here
        
        for s in self.clickdata.keys():
            waves = self.clickdata[s]['waves']
            t = self.clickdata[s]['timebase']
            spls = self.clickdata[s]['spls']
            thrs = {}
            r = {}
            n = {}
            # set minimum latency from maximal stimulus level
            # to figure out the order, since it is not stored in the data:
            first = np.std(waves[0,:])
            last = np.std(waves[-1,:])
            if first > last:
                index = -1
            else:
                index = 0
            
            # rmax = peakdetect.find_np(self.sample_freq, waves[-1,:],
            #     nzc_algorithm_kw={'dev': 1},
            #     guess_algorithm_kw={'min_latency': self.minlat})
            # min_latency = t[rmax[0]]  # latency of first + peak
            min_latency = self.minlat
            ml0 = min_latency
            for j in range(waves.shape[0])[::-1]:
                r[j] = peakdetect.find_np(self.sample_freq, waves[j,:],
                    nzc_algorithm_kw={'dev': 1},
                    guess_algorithm_kw={'min_latency': min_latency})
                min_latency = t[r[j][0]]-0.2  # new min latency is latency of peak at next higher intensity 
                if len(r[j]) > 0:
                    n[j] = peakdetect.find_np(self.sample_freq, -waves[j,:],
                        nzc_algorithm_kw={'dev': 1.5},
                        guess_algorithm_kw={'min_latency': t[r[j][0]]})  # find negative peaks after positive peaks
            A.analyze(t, waves)
            thr_spl = A.threshold_spec(waves, spls, SD=3.5)
            thrs[s] = thr_spl

            latmap = []
            spllat = []
            print '    thr_spl: %.1f  n>thr: %d' % (thr_spl, len(thrs))
            # fit a line to the latencies of the first (P1) peak for all traces
            # above threshold. The line is used to estimate the position of data
            # for measurements below threshold.
            # 
            for j in range(len(spls)):
                if spls[j] > thr_spl:
                    latmap.append(t[r[j][0]]) # get latency for first value
                    spllat.append(spls[j])
            if len(spllat) > 2:
                latp = np.polyfit(spllat, latmap, 1)
                spls2 = spls
                spls2[0] = spls2[0] - 5.
                spls2[-1] = spls2[-1] + 5.
                fitline = np.polyval(latp, spls2)
            else:
                fitline = []

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
                if len(fitline) > 2:
                    plottarget.plot(fitline, spls2, 'g-', linewidth=0.7)
                plottarget.plot([ml0, ml0], [np.min(spls), np.max(spls)], 'r--', linewidth=0.5)
                plottarget.plot([self.minlat, self.minlat], [np.min(spls), np.max(spls)], 'c--', linewidth=0.5)
                if spls[j] >= thr_spl:
                    IO[j] = 1.e6*(waves[j][r[j][0]] - waves[j][n[j][0]])
                else:
                    if len(fitline) > 2:
                        ti = int(fitline[j]/(self.sample_rate*1000.))
                        try:
                            IO[j] = 1.e6*waves[j][ti]
                        except:
                            print 'j, ti: ', j, ti, waves.shape
            
            if superIOPlot is not None:
                superIOPlot.plot(spls, IO, self.clickdata[s]['marker'])
                # print out the data for import into another plotting program
                print '*'*20
                print 'dataset: ', s
                print 't\tV'
                for i in range(len(spls)):
                    print '%.1f\t%.3f' % (spls[i], IO[i])
                print '*'*20
            
            if IOplot is not None:
                IOplot.plot(spls, 1e6*A.ppio, A.ppioMarker)
                IOplot.plot(spls, 1e6*A.rms_response, A.rmsMarker)
                IOplot.plot(spls, 1e6*A.rms_baseline, A.baselineMarker)
                ax2 = IOplot.twinx()
                ax2.plot(spls, A.psdwindow, A.psdMarker)
                ax2.tick_params('y', colors='r')
            if PSDplot is not None:
                for j in range(len(spls)):
                    PSDplot.semilogy(np.array(A.fr), np.array(A.psd[j])) # color=self.color_map[spls])
                PSDplot.set_ylim(1e-6, 0.01)
                PSDplot.set_xlim(100., 2000.)
        plottarget.set_xlim(0, 8.)
        plottarget.set_ylim(10., 115.)
        datatitle = os.path.basename(os.path.normpath(self.datapath))
        datatitle = datatitle.replace('_', '\_')
        plottarget.title.set_text(datatitle)
        plottarget.title.set_size(9)
        print""
        for s in thrs.keys():
            print "dataset: %s  thr=%.0f" % (s, thrs[s])
        print ""

    def plotTones(self, pdf=None):
        """
        Plot the tone ABR intensity series, one column per frequency, for one subject
        
        Parameters
        ----------
        select : list of str or int (default : None)
            A list of the run times for the datasets to plot for each subject.
            If None, then all detected tone ABR runs for that subject/frequency
            comtination are superimposed.
            The list looks like [['0115'], None, None]  - for each directory in order.
            Multiple runs can be indictated by [[115, 145], None, None] to include the 
            0115 and 0145 runs together (overplotted) for the first file.
        
        pdf : pdfPages object
            The pdfPages object that the plot will be appended to. Results in a multipage
            pdf, one page per subject. 

        """
        if len(self.tonemapdata) == 0:
            raise ValueError('plotTones requires that self.tonemapdata be populated first')
        
        A = Analyzer()  # create an instance of analyzer to use here
        freqs = []
        for run in self.tonemapdata.keys():  # just passing one dict... 
            freqs.extend(self.tonemapdata[run].keys())
        freqs.sort()
        f,  axarr = mpl.subplots(1, len(freqs), figsize=(12,6))
        datatitle = os.path.basename(os.path.normpath(self.datapath))[15:]
        datatitle = datatitle.replace('_', '\_')
        f.suptitle(datatitle)
        A = Analyzer()
        thrs = {}
        for i, s in enumerate(self.tonemapdata.keys()):
            for fr in self.tonemapdata[s]:  # next key is frequency
                waves = self.tonemapdata[s][fr]['waves']
                t = self.tonemapdata[s][fr]['timebase']
                spls = self.tonemapdata[s][fr]['spls']
                A.analyze(t, waves)
                thr_spl = A.threshold_spec(waves, spls, SD=3.5)
                thrs[s] = thr_spl
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
        self.cleanAxes(axarr)
        if pdf is not None:
            pdf.savefig()
            mpl.close()

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
        Remove axiessplines on top and right
        
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
             #print 'SPLfile: ', SPLfile
             rundict[run] = [float(spl.strip()) for spl in SPLfile if len(spl.strip()) > 0]  # clean up whitespace etc.
        return rundict
 
    def getFreqs(self, dir):
        """
        Return all the tonepip files in the directory. There is one tonepip file
        per frequency for its intensity run.
        
        """
        # return all the tone response runs in the directory. These are marked by having a 'khz'
        # file - but the case may change with program version, so ignore the case
        freqs = [f[:-8] for f in os.listdir(self.datapath) if re.search('kHz', f, re.IGNORECASE) is not None]
        rundict = {}
        for run in freqs:
          try:
              Freqfile = open(os.path.join(self.datapath, run+'-kHz.txt'), 'r')
          except:
              Freqfile = open(os.path.join(self.datapath, run+'-khz.txt'), 'r')
          
          rundict[run] = [float(khz.strip()) for khz in Freqfile if len(khz.strip()) > 0] # handle old data with blank line at end
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
        print dataset
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
        
        posf = pd.io.parsers.read_csv(os.path.join(self.datapath, fnamepos), delim_whitespace=True, lineterminator=self.term,
            skip_blank_lines=True, header=0)
        negf = pd.io.parsers.read_csv(os.path.join(self.datapath, fnameneg), delim_whitespace=True, lineterminator=self.term,
            skip_blank_lines=True, header=0)
        if posf.shape[0] == 0 or negf.shape[0] == 0: # malformed files
            return None
        negseries=[]
        posseries=[]
        for col in negf.columns:
            for i in range(len(negf)):
                negseries.append(negf[col][i])

        for col in posf.columns:
            for i in range(len(posf)):
                posseries.append(posf[col][i])
        
        wvfmdata = [(x + y)/2 for x, y in zip(negseries, posseries)]
        #wvfmdata = negseries  # just get one polarity
        # print fnamepos
        # print posf.shape
        # print 'col: ', col
        # print posf[col]
        # print 'len pos col: ', len(posf[col])
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
        print 'nans: ', np.argwhere(np.isnan(signal))
        if debugFlag:
            print "sfreq: %f LPF: %f HPF: %f" % (samplefreq, LPF, HPF)
        flpf = float(LPF)
        fhpf = float(HPF)
        sf = float(samplefreq)
        sf2 = sf/2.
        wp = [fhpf/sf2, flpf/sf2]
        ws = [0.5*fhpf/sf2, 2*flpf/sf2]
        if debugFlag:
            print "signalfilter: samplef: %f  wp: %f, %f  ws: %f, %f lpf: %f  hpf: %f" % (
               sf, wp[0], wp[1], ws[0], ws[1], flpf, fhpf)
        filter_b, filter_a = scipy.signal.iirdesign(wp, ws,
                gpass=1.0,
                gstop=60.0,
                ftype="ellip")
        msig = np.nanmean(signal)
        signal = signal - msig
        w = scipy.signal.lfilter(filter_b, filter_a, signal) # filter the incoming signal
        signal = signal + msig
        if debugFlag:
            print "sig: %f-%f w: %f-%f" % (np.amin(signal), np.amax(signal), np.amin(w), np.amax(w))
        return(w)


class Analyzer(object):
    """
    Provide general analysis functions for ABR waveforms
    
    """
    def __init__(self):
        self.ppioMarker = 'gs-'
        self.rmsMarker = 'bo-'
        self.psdMarker = 'r*-'
        self.baselineMarker = 'k-'
    
    def analyze(self, timebase, waves):
        """
        Main descriptive analysis. 
        Prodides the peak-to-peak, rms, baseline and mean spectral power
        for a set of ABR data
        
        Parameters
        ----------
        timebase : numpy array
            Time base for an individual trace
        waves : numpy array
            2-d array (spls x time) of ABR waveforms
        """
        
        self.sample_rate = 0.001*(timebase[1]-timebase[0])
        self.waves = waves
        self.timebase = timebase
        
        response_range = [2.0, 7] # window, in msec... 
        baseline = [0., 2.0]
        
        self.ppio = self.peaktopeak(response_range)
        self.rms_response = self.measure_rms(response_range)
        self.rms_baseline = self.measure_rms(baseline)
        self.specpower()
 
    def peaktopeak(self, tr):
        """
        Measure the peak-to-peak signal within the time window
        
        Parameters
        ----------
        tr : a 2 element list [start time, stop time]
            The time window in which to take the measure
        
        Returns
        -------
        pp : a list of the rms values across the spls in a run
        """
        tx = self.gettimeindices(tr)
        pp = np.zeros(self.waves.shape[0])
        for i in range(self.waves.shape[0]):
            pp[i] = np.max(self.waves[i][tx]) - np.min(self.waves[i][tx])
        return pp

    def measure_rms(self, tr):
        """
        Measure the rms signal within the time window
        
        Parameters
        ----------
        tr : a 2 element list [start time, stop time]
            The time window in which to take the measure
        
        Returns
        -------
        rms : a list of the rms values across the spls in a run
        """
        tx = self.gettimeindices(tr)
        rms = np.zeros(self.waves.shape[0])
        for i in range(self.waves.shape[0]):
            rms[i] = np.std(self.waves[i][tx])
        return rms

    def gettimeindices(self, tr):
        """
        Parameters
        ----------
        tr : list of 2 floats (no default)
            the two times bracketing a window
        
        Returns
        -------
        x : the indices into the array for all points within the window. The
            indices are bounded [t0...t1)
        """
        x, = np.where((tr[0] <= self.timebase) & (self.timebase < tr[1]))
        return x
        
    def specpower(self, fr=[500., 1500.], win=[0, -1]):
        """
        Compute the average power spectrum from the current waveform set
        over a window
        
        Parameters
        ----------
        fr : list of floats (default: [500., 1500.])
            frequency window for computation, in Hz
        win : list of ints (default: [0, -1])
            indices for the temporal window. Defaults are the 
            entire trace. Use result gettimeindices to get
            the index list
        
        Returns
        -------
        psdwindow : the power spectrum measured in the window across
            all spls in a run.
        """
        fs = 1./self.sample_rate
        psd = [None]*self.waves.shape[0]
        psdwindow = np.zeros(self.waves.shape[0])
        for i in range(self.waves.shape[0]):
            f, psd[i] = scipy.signal.welch(1e6*self.waves[i][win[0]:win[1]],
                fs, nperseg=256, nfft=8192, scaling='density')
            frx, = np.where((f >= fr[0]) & (f <= fr[1]))
            psdwindow[i] = np.nanmean(psd[i][frx[0]:frx[-1]])
        self.fr = f
        self.psd = psd  # save the psds to look at later if desired
        self.psdwindow = psdwindow
        return psdwindow

    def thresholds(self, waves, spls, tr=[0., 8.], SD=4.0):
        """
        Auto threshold detection:
        BMC Neuroscience200910:104  DOI: 10.1186/1471-2202-10-104
        Use last 10 msec of 25 msec window for SD estimates
        Computes SNR (max(abs(signal))/reference SD) for a group of traces
        The reference SD is the MEDIAN SD across the intensity run.
        
        Parameters
        ----------
        waves : numpy array
            2-d array, spl x time for waveforms collected from one run.
        spls : list 
            List of SPLS (floats) in the run
        tr : list of floats (default : [0., 8.])
            beginning and end time for the signal analysis window
        SD : float (default : 4.0)
            criterion to set the threshold level, in standard deviations 
        
        Returns
        -------
        spls : spls where the signal is above threshold. If no threshold was
            found, returns a np.nan
        """
        refwin = self.gettimeindices([15., 25.])
        sds = np.std(waves[:,refwin[0]:refwin[-1]], axis=1)
        self.median_sd = np.nanmedian(sds)
        tx = self.gettimeindices(tr)
        self.max_wave = np.max(np.fabs(waves[:, tx[0]:tx[-1]]), axis=1)
        self.snr = self.max_wave / self.median_sd
        thr, = np.where(self.max_wave >= self.median_sd*SD)  # find criteria threshold
        thrx, = np.where(np.diff(thr) == 1)  # find first contiguous point (remove low threshold non-contiguous)
        if len(thrx) > 0:
            return spls[thr[thrx[0]]]
        else:
            return np.nan

    def threshold_spec(self, waves, spls, tr=[0., 8.], SD=4.0):
        """
        Auto threshold detection using power spectrum:
        Use last 10 msec of 25 msec window for noise spectral estimates

        Parameters
        ----------
        waves : numpy array
            2-d array, spl x time for waveforms collected from one run.
        spls : list 
            List of SPLS (floats) in the run
        tr : list of floats (default : [0., 8.])
            beginning and end time for the signal analysis window
        SD : float (default : 4.0)
            criterion to set the threshold level, in standard deviations
            of signal power in the noise window.
        
        Returns
        -------
        spls : spls where the signal is above threshold. If no threshold was
            found, returns a np.nan
        
        """
        refwin = self.gettimeindices([15., 25.])
        sds = self.specpower(fr=[800., 1250.], win=[refwin[0], refwin[-1]])
        self.median_sd = np.nanmedian(sds)
        tx = self.gettimeindices(tr)
        self.max_wave = self.specpower(fr=[800., 1250.], win=[tx[0], tx[-1]])
        #np.max(np.fabs(waves[:, tx[0]:tx[-1]]), axis=1)
        thr, = np.where(self.max_wave >= self.median_sd*SD)  # find criteria threshold
        # thrx, = np.where(np.diff(thr) == 1)  # find first contiguous point (remove low threshold non-contiguous)
        if len(thr) > 0:
            return spls[thr[0]]
        else:
            return np.nan


class ABRWriter():
    def __init__(self, dataset, waves, filename, freqs, spls, fileindex=None):
        """
        Write one ABR file for Brad Buran's analysis program, in "EPL"
        format.
        
        TODO : fix this so that the ABR-files have proper and useful names, and check that it works.
        Note : This may not be necessary to have. The findpeaks routine in peakdetector.py from
                Buran's code is directly called above. This class is only necessary if you want to 
                use his code to get times etc.
        Filename structures:
            ABR-click-time  indicate click intenstiry series run in this
                directory, use time to distinguish file
                Example: ABR-click-0901.abr
            ABR-tone-time-freq indicate tone run number in this directory
                with frequency in Hz (6 places).
                Example: ABR-tone-0945-005460.dat

        Parameters
        ----------
        dataset : dict
            information about the dataset that will be analyzed
        
        waves : numpy array
            2-d array, spl x time for waveforms collected from one run.
        
        filename : str
            name of the file
        
        freqs : list
            List of frequencies (floats) in the filename if the stimtype is
            a tonepip
        
        spls : list 
            List of SPLS (floats) in the run
        
        fileindex : int (default: None)
            Index to the frequency to be written.
        
        Returns
        -------
        str : the filename that was written
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
        #    header5 = ":LEVELS:20;25;30;35;40;45;50;55;60;65;70;75;80;85;90; \n"
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
        print ABR_Datasets.keys()
        raise ValueError('Data set %s not found in our list of known datasets')
    if mode not in ['tones', 'clicks']:
        raise ValueError('Second argument must be tones or clicks')
    top_directory = os.path.join(basedir, ABR_Datasets[dsname]['dir'])
    
    dirs = [tdir for tdir in os.listdir(top_directory) if os.path.isdir(os.path.join(top_directory, tdir))]

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
        f, axarr = mpl.subplots(m, n, figsize=(12, h))
#        f2, axarr2 = mpl.subplots(m, n, figsize=(12, h))
#        f3, axarr3 = mpl.subplots(m, n, figsize=(12, h))
        f4, IOax = mpl.subplots(1, 1, figsize=(6,6))
        if axarr.ndim > 1:
            axarr = axarr.ravel()
        fofilename = os.path.join(top_directory, 'ClickSummary.pdf')
        for k in range(len(clicksel)):
            if 'nameselect' in ABR_Datasets[dsname].keys():  # limit directories to those with some string in the filename
                if dirs[k].find(ABR_Datasets[dsname]['nameselect']) == -1:
                    continue
            P = ABR(os.path.join(top_directory, dirs[k]), mode, info=ABR_Datasets[dsname])
            P.getClickData(select=clicksel[k]) 
            P.plotClicks(plottarget=axarr[k], superIOPlot=IOax) # , IOplot=axarr2[k])
        mpl.savefig(fofilename)
        
    if mode == 'tones':
        # first one has issues: 1346 stops at 40dbspl for all traces
        # 1301+1327 seem to have the higher levels, but small responses; remainder have fragementary frequencies
        # 0901 has whole series (but, is it hte same mouse?)
        if 'toneselect' in ABR_Datasets[dsname].keys():
            tonesel = ABR_Datasets[dsname]['toneselect']
        else:
            tonesel = [None]*len(dirs)
        fofilename = os.path.join(top_directory, 'ToneSummary.pdf')
        with PdfPages(fofilename) as pdf:
            for k in range(len(tonesel)):
                P = ABR(os.path.join(top_directory, dirs[k]), mode, info=ABR_Datasets[dsname])
                P.getToneData(select=tonesel[k]) 
                P.plotTones( pdf=pdf)

    mpl.show()
