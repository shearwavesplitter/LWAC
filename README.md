# LWAC 0.1.1
Light Weight Auto-Correlator. For calculating auto-correlations of seismic ambient noise or coda. 

## Publication / Reference

Stefan Mroczek, Frederik Tilmann, Joint ambient noise autocorrelation and receiver function analysis of the Moho, *Geophysical Journal International*, Volume 225, Issue 3, June 2021, Pages 1920–1934, [https://doi.org/10.1093/gji/ggab065](https://doi.org/10.1093/gji/ggab065)

```bibtex
@article{mroczek2021acor,
  title={Joint ambient noise autocorrelation and receiver function analysis of the Moho},
  author={Mroczek, Stefan and Tilmann, Frederik},
  journal={Geophysical Journal International},
  volume={225},
  number={3},
  pages={1920--1934},
  year={2021},
  publisher={Oxford University Press}
}
```

## Installation

```bash
pip install -e /home/user/path/to/LWAC
```

Requirements: [Obspy](https://docs.obspy.org/), its dependencies, and python 3.6 or greater. 

## Example

```python
from obspy.core import read
from obspy.signal import filter
from LWAC.core import *
#read with no arguments returns an example earthquake (for the example we ignore the actual content). We assume that the response has been removed
st=read()
#Sign bit normalisation
st.signbitnorm()
#Perform the auto-correlation. Because this example is just one event we adjust the window size to the duration of the trace
corr_st=auto_corr(st,30)
#The EHE-EHN correlation is only used for rotations so we can remove it
corr_st.remove(corr_st.select(channel='EHE-EHN')[0])
#Whiten with a gaussian window
corr_st.gauss_whiten(gwidth=3,eps=.01)
#Taper and filter
corr_st.taper(0.03)
corr_st.filter("bandpass",freqmin=0.3,freqmax=1,zerophase=True,corners=4)
#Apply pi/2 phase-shift
corr_st=phase_shift(corr_st,np.pi/2)
#Stack the auto-correlations
#Normally this is done on a stream containing only one component and more than three traces
#Because there are less than 24 traces only a linear stack is applied and the mirrored part of the trace is removed
stack=stack_autocorr(corr_st)
#Plot it
stack.plot()
```