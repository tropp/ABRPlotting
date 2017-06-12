import os
import os.path
# from optparse import OptionParser
import pandas as pd
import numpy as np
import scipy.signal
import matplotlib.pyplot as mpl
#import matplotlib.colors as colors
import itertools
import seaborn as sns  # makes plot background light grey with grid, no splines. Remove for publication plots

from matplotlib.backends.backend_pdf import PdfPages

#current 10 June 2017

class Pullfiles():
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
    
    def __init__(self, datapath):

        self.datapath = datapath
        files = [f for f in os.listdir(datapath) if os.path.isfile(os.path.join(datapath, f))]
        print '\n', datapath

        striptxt = []
        self.sample_rate = 1e-5 # standard interpolated sample rate for matlab abr program: 10 microsecnds
        self.sample_freq = 1./self.sample_rate
        self.hpf = 200.
        self.lpf = 3000.  # filter frequencies, Hz
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
            print f, self.tonemaps[f]

        print '\nClick Intensity series: '
        self.clickmaps = {}
        for s in self.clicks.keys():
            self.clickmaps[s] = {'stimtype': 'click', 'SPLs': self.spls[s]}
            print s, self.clickmaps[s]

        # build color map where each SPL is a color (cycles over 12 levels)
        bounds = np.linspace(0, 120, 25)  # 0 to 120 db inclusive, 5 db steps
        color_labels = np.unique(bounds)
        rgb_values = sns.color_palette("Set2", 25)
        # Map label to RGB
        self.color_map = dict(zip(color_labels, rgb_values))
        
    def plotClicks(self, select=None, plottarget=None):
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
        for i, s in enumerate(self.clickmaps.keys()):
            if select is not None:
                if s[9:] not in select:
                    continue
            waves = self.get_combineddata(s, self.clickmaps[s])
            t = np.linspace(0, waves.shape[1]*self.sample_rate*1000., waves.shape[1])
            spls = self.clickmaps[s]['SPLs']
            for j in range(len(spls))[::-1]:
                plottarget.plot(t, 2*waves[j]*1e6+spls[len(spls)-j-1], color=self.color_map[spls[len(spls)-j-1]])
        plottarget.set_xlim(0, 8.)
        plottarget.set_ylim(10., 110.)
        datatitle = os.path.basename(os.path.normpath(self.datapath))[15:]
        datatitle = datatitle.replace('_', '\_')
        plottarget.title.set_text(datatitle)
        plottarget.title.set_size(9)
    
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

        f, axarr = mpl.subplots(1, len(freqs), figsize=(12,6))
        datatitle = os.path.basename(os.path.normpath(self.datapath))[15:]
        datatitle = datatitle.replace('_', '\_')
        f.suptitle(datatitle)
        for i, s in enumerate(self.tonemaps.keys()):
            if select is not None:
                if s[9:] not in select:
                    continue
            for fr in self.tonemaps[s]['Freqs']:
                waves = self.get_combineddata(s, self.tonemaps[s], freq=fr)
                t = np.linspace(0, waves.shape[1]*self.sample_rate*1000., waves.shape[1])
                spls = self.tonemaps[s]['SPLs']
                plottarget = axarr[freqs.index(fr)]
                for j in range(len(spls))[::-1]:
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
          rundict[run] = [float(khz) for khz in Freqfile]
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

        d1 = len(wvfmdata)/len(posf[col])
        waves = np.reshape(wvfmdata,(d1,len(posf[col])))
        for i, w in enumerate(waves):
            waves[i, -1] = waves[i, -2]  # remove nan from end of waveform...
            waves[i,:] = self.filter(waves[i,:], 4, self.lpf, self.hpf, samplefreq=self.sample_freq)
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


    # def combine_polarity(self, datapath, dataset):
    #      """
    #      Read the files associated with the dataset
    #      add the + and - files together
    #      return a dict with keys by intensity of the combined waveforms.
    #
    #
    #      """
    #
    #
    #     for ABRfreq in fl1:
    #
    #         if any(('-'+ABRfreq+'.' in nff for nff in nf) and (ABRfreq in pff for pff in pf)):
    #         # nindex = nf.find(ABRfreq)
    #         # pindex = pf.find(ABRfreq)
    #         # if nindex != -1:
    #             nff = [genfile for genfile in nf if ABRfreq in genfile]
    #             print 'negative file:', nff
    #             pff = [genfile for genfile in pf if ABRfreq in genfile]
    #             print 'positive file:', pff
    #             os.chdir(datapath)
    #             negf = pd.io.parsers.read_csv(nff[0],delim_whitespace=True, lineterminator='\r',
    #                 skip_blank_lines=True, header=0)
    #             posf = pd.io.parsers.read_csv(pff[0],delim_whitespace=True, lineterminator='\r',
    #                 skip_blank_lines=True, header=0)
    #             negseries=[]
    #             posseries=[]
    #             for col in negf.columns:
    #                 print col
    #                 for i in range(len(negf)):
    #                     negseries.append(negf[col][i])
    #
    #             print 'length of negf:', len(negf[col])
    #                 #negseries.append(negf[col])
    #             for col in posf.columns:
    #                 for i in range(len(posf)):
    #                     posseries.append(posf[col][i])
    #
    #                 #posseries.append(posf[col])
    #
    #             wvfmdata = [(x + y)/2 for x, y in zip(negseries, posseries)]
    #
    #             print 'length of wvfmdata:', len(wvfmdata)
    #             print 'length of posf:', len(posf[col])
    #             print 'length of posseries:', len(posseries)
    #             print 'length of negseries:', len(negseries)
    #             d1 = len(wvfmdata)/len(posf[col])
    #             waves = np.reshape(wvfmdata,(d1,len(posf[col])))
    #             twaves = waves.T
    #
    #             denote = 'ABR-' + ABRfreq[:-4] + '-3'
    #
    #             np.savetxt(denote, twaves, delimiter='\t ', newline='\r\n')
    #         else:
    #             print 'No file with this frequency', ABRfreq
    #         ABRfiles.append(denote)
    #     return ABRfiles
    

    def writeABRfile(self, dataset, waves, filename, fileindex=None):
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
        twaves = waves.T
        abr_filename = ('ABR-{0:s}-{1:4s}'.format(dataset['stimtype'], filename[8:12])
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
    top_directory = "/Users/pbmanis/Desktop/data/ABR_DATA/Coats-NrCAM-ABRs/"
    top_directory = "/Users/pbmanis/Desktop/data/ABR_DATA/CNTNAP2/"
    dirs = [tdir for tdir in os.listdir(top_directory) if os.path.isdir(os.path.join(top_directory, tdir))]
    mode = 'tones'
#    for dirdata in dirs:
    if mode == 'clicks':
        #clicksel = [['1228'], None, None, None, None, None, None, None]
        clicksel = [None]*len(dirs)
        f, axarr = mpl.subplots(1, len(clicksel), figsize=(12, 5))
        fofilename = os.path.join(top_directory, 'ClickSummary.pdf')
        for k in range(len(clicksel)):
            P = Pullfiles(os.path.join(top_directory, dirs[k]))
            P.plotClicks(select=clicksel[k], plottarget=axarr[k])
        mpl.savefig(fofilename)
        
    if mode == 'tones':
        # first one has issues: 1346 stops at 40dbspl for all traces
        # 1301+1327 seem to have the higher levels, but small responses; remainder have fragementary frequencies
        # 0901 has whole series (but, is it hte same mouse?)
#        tonesel = [['0901'], ['1446', '1505'], None, None, None, None, None, None]
        tonesel = [None]*len(dirs)
        fofilename = os.path.join(top_directory, 'ToneSummary.pdf')
        with PdfPages(fofilename) as pdf:
            for k in range(len(tonesel)):
                P = Pullfiles(os.path.join(top_directory, dirs[k]))
                P.plotTones(select=tonesel[k], pdf=pdf)

    mpl.show()
