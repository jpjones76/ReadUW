# ReadUW
Matlab/Octave reader for the old University of Washington (UW) seismic data format.

Please place the file in your Matlab path and type "help ReadUW" at a command prompt. :)

## Notes
* Unprocessed UW format waveform data won't work unless the data file contains a valid channel structure.
* The UW data format had two versions, UW-1 and UW-2; this only works with the latter. However, converters from UW-1 to UW-2 exist; contact me if you need help locating one.
* Many thanks to prof. Steve Malone (Univ. Washington) for providing the original UW data format source code.
