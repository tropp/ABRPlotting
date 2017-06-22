import sys
import os
import os.path
import commands
from collections import OrderedDict
# from optparse import OptionParser
import pandas as pd
import numpy as np
import scipy.signal
import matplotlib.pyplot as mpl
#import matplotlib.colors as colors
import itertools
import seaborn as sns  # makes plot background light grey with grid, no splines. Remove for publication plots
from ABR_Analysis import peakdetect

from matplotlib.backends.backend_pdf import PdfPages

computer_name = commands.getoutput('scutil --get ComputerName')
if computer_name == 'Lytle':
    basedir = '/Volumes/Pegasus/ManisLab_Data3/abr_data'
elif computer_name == 'Tamalpais':
    #basedir = '/Users/pbmanis/Desktop/data/ABR_DATA'
    basedir = '/Volumes/Backup2B/ABR_DATA'
else:
    raise ValueError("Need valid computer name to set base path to data")
    
#current 10 June 2017
ABR_Datasets = {'NrCAM': "Tessa/Coate-NrCAM-ABRs/",
                 'TessaCBA': 'Tessa/CBA',
                 'TessaCBANE': 'Tessa/CBA_NoiseExposed',
                 'CNTNAP2X': 'Tessa/CNTNAP2',
                 'CNTNAP2Het': 'Tessa/CNTNAP2_Het',
                 'CNTNAP2HG': 'Tessa/CNTNAP2_Het_GP4.3',
                 'CNTNAP2KO': 'Tessa/CNTNAP2_KO',
                 'CNTNAP2WT': 'Tessa/CNTNAP2_WT',
                 'GP43Norm': 'Tessa/GP4.3-Thy1-Normal',
                 'Yong': "Yong\'s ABRs",
                 'Jun': "JUN\'s ABRs",
                 'RuiliCBAP40': "RuiliABRdata_2010-2015/CBA-P21-P40/",
                 'RuiliCBAP20': "RuiliABRdata_2010-2015/CBA-P10-P20/",
                 'Eveleen': "Eveleen\'s ABRs"
                 } 
                 
                 
    

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
    """
    
    def __init__(self, datapath, mode='clicks', invert=True):

        self.datapath = datapath
        self.invert = invert  # data flip... depends on "active" lead connection to vertex (false) or ear (true).
        files = [f for f in os.listdir(datapath) if os.path.isfile(os.path.join(datapath, f))]
        print '\n', datapath

        striptxt = []
        self.sample_rate = 1e-5 # standard interpolated sample rate for matlab abr program: 10 microsecnds
        self.sample_freq = 1./self.sample_rate
        self.hpf = 300.
        self.lpf = 1500.  # filter frequencies, Hz
        # self.datadate = filebase[:8]

        for f in files:
            striptxt.append(f[16:-4])

        self.spls = self.getSPLs(files)
        self.freqs = self.getFreqs(files)
        # The keys in the freqs files represent the different runs for tone maps
        # The keys in the spl files that are not in the freqs files represent click runs
        self.clicks = {}
        for s in self.spls.keys():
            if s[:13] not in [f[:13] for f in self.freqs.keys()]:
                self.clicks[s] = self.spls[s]

        # inspect the directory and get a listing of all the tone and click maps
        # that can be found.
        print 'Frequency maps: '
        self.tonemaps = {}
        for f in self.freqs.keys():
            self.tonemaps[f] = {'stimtype': 'tonepip', 'Freqs': self.freqs[f], 'SPLs': self.spls[f[:13]]}
            if mode == 'tones':
                print f, self.tonemaps[f]

        print '\nClick Intensity series: '
        self.clickmaps = {}
        for s in self.clicks.keys():
            self.clickmaps[s] = {'stimtype': 'click', 'SPLs': self.spls[s]}
            if mode == 'clicks':
                print s, self.clickmaps[s]

        # build color map where each SPL is a color (cycles over 12 levels)
        bounds = np.linspace(0, 120, 25)  # 0 to 120 db inclusive, 5 db steps
        color_labels = np.unique(bounds)
        self.color_map = self.makeColorMap(25, color_labels)
        color_labels2 = range(25)
        self.summaryClick_color_map = self.makeColorMap(25, range(25))

    def makeColorMap(self, N, labels):
        """
        create a color map of N levels using the specified labels
        
        """
        rgb_values = sns.color_palette("Set2", N)
        return dict(zip(labels, rgb_values))

        
    def plotClicks(self, select=None, plottarget=None, IOplot=None, PSDplot=None, superIOPlot=None):
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

        """
        # get data for clicks and plot all on one plot
        A = Analyzer()
        thrs = {}
        for i, s in enumerate(self.clickmaps.keys()):
            print 'select: ', select
            print 's: ', s
            if select is not None:
                if s[9:] not in select:
                    continue
            print "s: ", s
            if s[0:8] == '20170419':
                smarker = 'go-'
            elif s[0:8] in ['20170608', '20170609']:
                smarker = 'bs-'
            else:
                smarker = 'kx-'
            
            waves = self.get_combineddata(s, self.clickmaps[s])[::-1]
            t = np.linspace(0, waves.shape[1]*self.sample_rate*1000., waves.shape[1])
            spls = self.clickmaps[s]['SPLs']  # get spls
            r = {}
            n = {}
            for j in range(waves.shape[0]):
                r[j] = peakdetect.find_np(self.sample_freq, waves[j,:],
                    nzc_algorithm_kw={'dev': 1},
                    guess_algorithm_kw={'min_latency': 2.5})
                if len(r[j]) > 0:
                    n[j] = peakdetect.find_np(self.sample_freq, -waves[j,:],
                        nzc_algorithm_kw={'dev': 2.5},
                        guess_algorithm_kw={'min_latency': t[r[j][0]]})  # find negative peaks after positive peaks
            halfspl = np.max(spls)/2.
            latmap = []
            spllat = []
            for j in range(len(spls)):
                if spls[j] > halfspl:
                    latmap.append(t[r[j][0]]) # get latency for first value
                    spllat.append(spls[j])
            print spllat, latmap
            latp = np.polyfit(spllat, latmap, 1)
            fitline = np.polyval(latp, spls)
            # mpl.figure()
            # mpl.plot(spllat, latmap, 'o')
            # mpl.plot(spls, fitline, 'g-')
            #print r
            A.analyze(t, waves)
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
                    print 'ti: ', ti
                    IO[j] = 1.e6*waves[j][ti]
            
            if superIOPlot is not None:
                superIOPlot.plot(spls, IO, smarker, color=self.summaryClick_color_map[i])

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
        datatitle = os.path.basename(os.path.normpath(self.datapath))[15:]
        datatitle = datatitle.replace('_', '\_')
        plottarget.title.set_text(datatitle)
        plottarget.title.set_size(9)
        print""
        for s in thrs.keys():
            print "dataset: %s  thr=%.0f" % (s, thrs[s])
        print ""
        
    
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
        for i, s in enumerate(self.tonemaps.keys()):
            print i, s, s[9:]
            if select is not None:
                if s[9:] not in select:
                    continue
            for f in self.tonemaps[s]['Freqs']:
                if f not in freqs:
                    freqs.append(f)
        freqs.sort()
        if len(freqs) == 0:  # check of no tone pip ABR data in this directory
            return

        f, axarr = mpl.subplots(1, len(freqs), figsize=(12,6), num="Tones")
        datatitle = os.path.basename(os.path.normpath(self.datapath))[15:]
        datatitle = datatitle.replace('_', '\_')
        f.suptitle(datatitle)
        A = Analyzer()
        for i, s in enumerate(self.tonemaps.keys()):
            if select is not None:
                if s[9:] not in select:
                    continue
            thr_spls = np.zeros(len(self.tonemaps[s]['Freqs']))
            for k, fr in enumerate(self.tonemaps[s]['Freqs']):
                waves = self.get_combineddata(s, self.tonemaps[s], freq=fr)
                t = np.linspace(0, waves.shape[1]*self.sample_rate*1000., waves.shape[1])
                spls = self.tonemaps[s]['SPLs']
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
        ax = mpl.subplot(1, 1, num=num)
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
        print 'threshold lists: ', thrfrs
        # sort the threshold list
        thrs_sorted = OrderedDict(sorted(thrfrs.items(), key=lambda t: t[0]))
        frmean = np.zeros(len(thrs_sorted.keys()))
        frstd = np.zeros(len(thrs_sorted.keys()))
        print 'threshold list sorted: ', thrs_sorted
        print ('len mean: ', len(frmean))
        for i, f in enumerate(thrs_sorted):
            print 'i, f: ', i, f
            print thrs_sorted[f]
            frmean[i] = np.nanmean(thrs_sorted[f])
            frstd[i] = np.nanstd(thrs_sorted[f])
        ax.errorbar(np.array(thrs_sorted.keys())/1000., frmean, yerr=frstd, fmt='o')
        ax.set_xlim(1.8, 65.)
        xt = [2., 4., 8., 16., 32., 64.]
        mpl.xticks(xt, [str(x) for x in xt])
#        mpl.show()
        
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
    
    
    """
    def __init__(self):
        self.ppioMarker = 'gs-'
        self.rmsMarker = 'bo-'
        self.psdMarker = 'r*-'
        self.baselineMarker = 'k-'
    
    def analyze(self, timebase, waves):
        
        self.sample_rate = 0.001*(timebase[1]-timebase[0])
        self.waves = waves
        self.timebase = timebase
        
        response_range = [2.5, 7] # window, in msec... 
        baseline = [0., 2.0]
        
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
        print 'median sds: ', self.median_sd
        tx = self.gettimeindices(tr)
        self.max_wave = self.specpower(fr=[800., 1250.], win=[tx[0], tx[-1]])
        #np.max(np.fabs(waves[:, tx[0]:tx[-1]]), axis=1)
        print self.max_wave
        print self.median_sd
        thr, = np.where(self.max_wave >= self.median_sd*SD)  # find criteria threshold
        print 'thr: ', thr
#        thrx, = np.where(np.diff(thr) == 1)  # find first contiguous point (remove low threshold non-contiguous)
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
        print ABR_Datasets.keys()
        raise ValueError('Data set %s not found in our list of known datasets')
    if mode not in ['tones', 'clicks']:
        raise ValueError('Second argument must be tones or clicks')

    top_directory = os.path.join(basedir, ABR_Datasets[dsname])
    
    dirs = [tdir for tdir in os.listdir(top_directory) if os.path.isdir(os.path.join(top_directory, tdir))]

    if mode == 'clicks':
#        clicksel = [['0849'], None, None, None, None, None, None, None]
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
        f, axarr = mpl.subplots(m, n, figsize=(12, h), num='Click Traces')
        f2, axarr2 = mpl.subplots(m, n, figsize=(12, h), num='Click IO Summary')
#        f3, axarr3 = mpl.subplots(m, n, figsize=(12, h))
        f4, IOax = mpl.subplots(1, 1, figsize=(6,6), num='Click IO Overlay')
        if axarr.ndim > 1:
            axarr = axarr.ravel()
        if axarr2.ndim > 1:
            axarr2 = axarr2.ravel()
        fofilename = os.path.join(top_directory, 'ClickSummary.pdf')
        for k in range(len(clicksel)):
            P = ABR(os.path.join(top_directory, dirs[k]), mode)
            P.plotClicks(select=clicksel[k], plottarget=axarr[k], superIOPlot=IOax,
                IOplot=axarr2[k])
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
                P.plotTones(select=tonesel[k], pdf=pdf)
                allthrs[dirs[k]] = P.thrs
        P.plotToneThresholds(allthrs, num='Tone Thresholds')
        tthr_filename = os.path.join(top_directory, 'ToneThresholds.pdf')
        mpl.savefig(tthr_filename)
    mpl.show()
