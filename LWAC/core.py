from obspy import Stream
import numpy as np
from scipy.signal import correlate
from scipy.signal import gaussian
from obspy.signal import detrend
from obspy.signal.invsim import cosine_taper
import cmath
import scipy.signal as ss
import numpy as np
import copy
from obspy import Trace

def auto_corr(stream,corrsize=60*60,comp1=['*HZ'],comp2=['*HZ'],allcomp=True):
    """ Auto-correlate an obspy stream object
    ...

    Parameters
    ----------
    stream : stream stream object
    corrsize : int
        size of the window to be correlated in seconds
    comp1/comp2 : list
        two equal length lists identifying the components to correlate
    allcomp : boolean
        Auto-correlate *HZ, *HE, *HN, and cross-correlate *HE and *HN for use in rotations (true/false)

    return: stream with auto-correlations
    """
    corr_stream=Stream()
    stream.merge(fill_value=0)
    if allcomp:
        comp1=['*HZ','*HN','*HE','*HE']
        comp2=['*HZ','*HN','*HE','*HN']
    for windowed_st in stream.slide(window_length=corrsize, step=corrsize,include_partial_windows=False):
        if not np.isnan(np.sum(windowed_st[0].data)):
            windowed_st.detrend(type="linear")
            windowed_st.taper(0.05)
            windowed_st=signbitnorm(windowed_st)
            for i in range(len(comp1)):
                st1=windowed_st.select(channel=comp1[i])
                st2=windowed_st.select(channel=comp2[i])
                for j in range(len(st1)):
                    tr=st1[j]
                    tr2=tr.copy()
                    tr2.stats['channel']=tr2.stats['channel']+"-"+st2[j].stats['channel']
                    cr=correlate(st1[j].data,st2[j].data,mode='full',method="fft")
                    cr=(cr/st1[j].data.std()**2)/st1[j].data.size
                    tr2.data=cr
                    corr_stream+=tr2
        return(corr_stream)

def signbitnorm(stream):
    """Sign-bit normalise an obspy stream in-place
    """
    for tr in stream:
        v=tr.data
        v[v < 0]=-1
        v[v > 0]=1
        tr.data=v
    return(stream)

setattr(Stream,'signbitnorm',signbitnorm)


def gauss_whiten(stream,gwidth=6,eps=0.001):
    """ Deconvolve an auto-correlation with a windowed version of itself (in-place)
    ...

    Parameters
    ----------
    stream : obspy stream object
    gwidth : standard deviation of the gaussian window in seconds
    eps : fraction of spectral mean used as a water level to avoid spectral holes in the deconvolution.

    return: Whitened stream
    """
    stream2=stream.copy()
    for i in range(len(stream)):
        s=stream[i].stats['sampling_rate']
        window=gaussian(stream[i].stats['npts'], std=gwidth*s)
        stream2[i].data=stream2[i].data*window
        decont=deconvolve_traces_autocorr(Stream(stream[i]),stream2[i],eps=eps)
        stream[i].data=decont[0].data
    return(stream)

setattr(Stream,'gauss_whiten',gauss_whiten)

def deconvolve_traces_autocorr(signal, divisor, eps, freq=[], residual=False):
    """ Deconvolve a time series from a set of time series.

    This function is imported from miic https://github.com/miic-sw/miic

    The function is a wrapper for the :class:`scipy.signal.convolve`
    function.
    
    :type signal: :class:`~obspy.core.stream.Stream`
    :param signal: signal from which the divisor is to be deconvolved
    :type divisor: :class:`~obspy.core.trace.Trace`
    :param divisor: time series that is to be deconvolved from signal
    :type eps: float
    :param eps: fraction of spectral mean used as a water level to 
        avoid spectral holes in the deconvolution.
    :type freq: two element array-like
    :param freq: frequency range for the estimation of the mean power
        that is scaled with ``eps`` to obtian the water level  
    :type residual: bool
    :param residual: return residual if True, defaults to False

    :rtype: obspy.core.stream
    :return: **(dcst, rst)**: decorrelated stream and residual (only
        if ``residual=True``

    """

    
    #zerotime = UTCDateTime(1971,1,1)
    # trace length is taken from signal (must be even to use real fft)
    if signal[0].stats['npts'] % 2:
        trlen = signal[0].stats['npts']+1
        delsamp=True
    else:
        trlen = signal[0].stats['npts']
        delsamp=False
    delta = divisor.stats['delta']
    
    
    # prepare divisor
    divisor.detrend(type='constant')

    divisor.trim(starttime=divisor.stats['starttime'],endtime=divisor.stats['starttime']+
                 (trlen-1)*delta,pad=True,fill_value=0,nearest_sample=False)
    # FFT divisor
    fd = np.fft.rfft(divisor.data)
    # estimate the waterlevel to stabilize deconvolution
    if freq:
        f = np.linspace(-signal[0].stats['sampling_rate']/2., signal[0].stats['sampling_rate']/2.,len(fd))
        ind = np.nonzero(np.all([f>freq[0],f<freq[1]],axis=0))
        wl = eps * np.mean((fd*fd.conj())[ind])
    else:
        wl = eps * np.mean((fd*fd.conj()))
    
    # create the output stream
    dcst = Stream()
    rst = Stream()
    for tr in signal:
        if tr.stats['sampling_rate'] != divisor.stats['sampling_rate']:
            print("Sampling rates don't match for \n %s") % tr
            continue
        
        # prepare nuerator
        tr.detrend('constant')
        tr.trim(starttime=tr.stats['starttime'], endtime=tr.stats['starttime']+
                (trlen-1)*delta,pad=True,fill_value=0,nearest_sample=False)
        taper = cosine_taper(tr.stats['npts'])
        tr.data *= taper
        # fft numerator
        sf = np.fft.rfft(tr.data)
        
        # calculate deconvolution
        fdc = sf*fd/(fd**2+wl)
        dc = np.fft.irfft(fdc)

        # template to hold results
        dctr = tr.copy()
        # propagate metadata
        dctr.stats = tr.stats#combine_stats(tr,divisor)
        if delsamp:
            dc=np.delete(dc,len(dc)-1)
        dctr.data = dc
        dctr.stats['npts'] = len(dc)
        #dctr.stats['starttime'] = tr.stats['starttime']#zerotime - (divisor.stats['starttime']-tr.stats['starttime'])
        #dctr.stats_tr1 = tr.stats
        #dctr.stats_tr2 = divisor.stats
        
        # append to output stream
        dcst.append(dctr)
        
        if residual:
            # residual
            rtr = dctr.copy()
            rtr.data = tr.data - np.fft.irfft(fdc * fd)
            # append to output stream
            rst.append(rtr)
            return (dcst, rst)
        
        return dcst


def phase_shift(stream,shift=np.pi/2):
    """ Apply a phase shift to each trace in stream (in-place)
    ...

    Parameters
    ----------
    stream : obspy stream
    shift : phase shift magnitude (radians)
    
    return : Phase shifted stream
    """
    for tr in stream:
        t=tr.times()
        signal=tr.data
        signalfft = np.fft.rfft(signal)
        newsignalfft = signalfft * cmath.rect( 1., -shift )
        newsignal = np.fft.irfft(newsignalfft)
        tr.data=newsignal
    return(stream)

setattr(Stream,'phase_shift',phase_shift)

def PWStack(rcorr,v=2,rw=False,prestack=True,psl=24):
    """Phase weighted stack
    ...

    Parameters
    ----------
    rcorr : Matrix containing traces as rows (produced from obspy stream with extract_ndarray)
    v : Order of the phase weights (0 is equal to a linear stack)
    rw : Return the phase weights (not the phase weighted stack)? (True/False)
    prestack : Apply a linear pre-stack of length psl? (True/False)
    psl : Prestack length. Number of traces to stack linearly before phase weightes stack.

    return: Single stacked array
    
    Schimmel, M., & Paulssen, H. (1997). Noise reduction and detection of weak, coherent signals through phase-weighted stacks. Geophysical Journal International, 130(2), 497-505.
    """
    if prestack:
        bl=[]
        ln=len(rcorr)
        dz=int(ln/psl)
        if dz == 0:
            bl=np.nansum(rcorr,axis=0)
        rem=ln-dz*psl
        for i in np.arange(dz):
            fin=i*psl+psl
            if fin == dz*psl:
                fin=fin+rem
            mstack=np.nansum(rcorr[i*psl:fin],axis=0)
            bl.append(mstack)
        rcorr=np.asarray(bl)
        if dz == 0:
            return(rcorr)
    analytic_signal=ss.hilbert(rcorr)
    instantaneous_phase = np.unwrap(np.angle(analytic_signal))
    ip=np.exp(instantaneous_phase*1j)
    stack=abs(np.nansum(ip,axis=0)/len(rcorr))
    linstack=np.nansum(rcorr,axis=0)/len(rcorr)   
    pws=(stack**v)*linstack
    pw=stack**v
    if rw:
        return(pw)
    else:
        return(pws)



def stack_autocorr(st,stacklength='all',pw=True,v=2,prestack=True,psl=24,rw=False,half=True):
    """Stack auto-correlation (stream)
    ...

    Parameters
    ----------
    st : obspy stream
    stacklength : integer (number of traces), or 'all' (to stack all traces), or 'month' (for monthly stacks).
    pw : Apply a phase weighted stack? (True/False)
    v : Order of the phase weights (0 is equal to a linear stack)
    prestack : Apply a linear pre-stack of length psl? (True/False)
    psl : Prestack length. Number of traces to stack linearly before phase weightes stack.
    rw : Return the phase weights (not the phase weighted stack)? (True/False)
    half : Half the auto-correlation so its not mirrored? (True/False)

    return: obspy stream containing stacks equal to value defined in 'stacklength'
    """
    stream=st.copy()
    stream.sort(keys=["starttime"])
    stream=half_autocorr(stream,flip=False,half=half)
    stream2=Stream()
    start=copy.deepcopy(stream[0].stats['starttime'])
    end=copy.deepcopy(stream[len(stream)-1].stats['endtime'])
    curtime=copy.deepcopy(start)
    while curtime < end:
        swin=copy.deepcopy(curtime)
        if stacklength == 'month':
            ewin=copy.deepcopy(swin)
            if ewin.month < 12:
                try:
                    ewin.month=ewin.month+1
                except:
                    ewin.day=ewin.day-3
                    ewin.month=ewin.month+1
            else:
                 ewin.year=ewin.year+1
                 ewin.month=1
        else:
            if stacklength == 'all':
                ewin=end
            else:
                ewin=copy.deepcopy(curtime+stacklength*60*60*24)
        sl=stream.slice(swin,ewin)
        curtime=copy.deepcopy(ewin)
        if len(sl) > 0:
            slr= extract_ndarray(sl,smartpad=False)
            if len(slr) > 0:
                if len(slr) == 1:
                    stt=slr[0]
                else:
                    if pw:
                        stt=PWStack(slr,v=v,psl=psl,rw=rw,prestack=prestack)
                    else:
                        stt=np.nansum(slr,axis=0)
                sltr=Trace(stt)
                sltr.stats=sl[0].stats
                stream2+=sltr
    return(stream2)


def half_autocorr(stream,flip=False,half=True):
    nds=extract_ndarray(stream)
    if flip:
        ndf=np.flip(nds,axis=1)
    else:
        ndf=nds
    if half:
        top=np.ceil(len(ndf[0])/2.0)
        ndf=ndf[:,:int(top)]
    for i in range(len(stream)):
        stream[i].data=ndf[i]
    return(stream)

def extract_ndarray(stream,time=False,tadj=[],smartpad=True):
    """Take obspy steam and turn into matrix of traces (for later stacking or plotting)
    ...

    Parameters
    ----------
    stream : obspy stream
    time : Return matrix of times? (True/False)
    tadj : If time=True, adjust the time of each trace by some constant (array, list)
    smartpad : Pad short auto-correlations so each trace has the same length. The auto-correlation peak remains at the centre. (True/False)
    
    return : Numpy array
    """
    if smartpad:
        stream=smart_pad(stream)
    lls=max([len(x.data) for x in stream])
    arry=np.ndarray(shape=(len(stream),lls), dtype=float, order='F')*0
    for i in range(len(stream)):
        if time:
            if len(tadj) == 0:
                tadd=0
            else:
                tadd=tadj[i]
            arry[i]=stream[i].times() +tadd
        else:
            if len(arry[i]) == len(stream[i].data):
                arry[i]=stream[i].data
    return(arry)



def smart_pad(stream):
    lls=np.asarray([len(x.data) for x in stream])
    if max(lls) != min(lls):
        ind=np.arange(len(lls))
        gd=ind[lls == max(lls)]
        bd=ind[lls != max(lls)]
        cmparry=stream[gd[0]].data
        maxind=np.arange(len(cmparry))[cmparry == max(cmparry)][0]
        for i in bd:
            ex=stream[i].data
            fdif=max(lls)-lls[i]
            dif=int(fdif/2)
            l1=np.zeros(dif)
            ex2=np.concatenate((l1,ex,l1))
            rem=fdif-dif*2
            if rem == 1:
                maxindex=np.arange(len(ex2))[ex2 == max(ex2)][0]
                if maxindex == maxind:
                    ex2=np.concatenate((ex2,[0]))
                else:
                    ex2=np.concatenate(([0],ex2))
            stream[i].data=ex2
    return(stream)