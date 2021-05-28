"Extract_Regions_with_BDP2" Fiji macro:

Extracts regions from images that are too big to be opened in non virtual
mode in Fiji:

1. Asks for input image path and output folder;
2. Opens the image using Bio-formats Dialog box: 
   If too big to be opened in real mode, the image should be opened in virtual mode;
3. In an infinite loop:
   For each region to be extracted, asks to:
   - draw a rectangle and to select first slice, channel and time-point if any;
   - select last slice, channel and timepoint;
   Validating last slice, channel and timepoint with the SHIFT key down terminates
   the loop and closes the image.
4. Saves ROIs from "ROI Manager" as a zip to output folder;
5. Opens the image using BigDataProcessor2;
6. For each region (XYZCT domain) defined previously:
   - Crops the image to a new viewer;
   - Saves the crop to output folder.

Tested with following input types:
 - Metamorph Multi Dimensional series (input file: .nd);
 - Leica Image Files (.lif);
 - Zeiss (.czi) pyramidal images need the crop domains to be defined in the
   highest resolution image (#1). Each crop results in pyramidal image file
   containing a number of resolutions depending on the XY size of the crop- 
   regions. The resolution sequence of output images is 1, 1/2, 1/4, ...;
 - Single plane (Z*C*T = 1) TIFFs and compressed ImageJ TIFFs (.zip).

Notes:
1. In case of slow reactivity of Fiji: to terminate the crop-domains definition
   the SHIFT key must be bept down until the image closes.
   If not, a new domain definition may be initiated.
2. For one dimensional TIFFs, only one "Action Required" window opens per crop-region
3. The macro launches the "Memory" window if not open and "MemoryMonitor_Launcher.jar" 
   is found in the "plugins" folder (does nothing if not found).
   The "Memory" window stays open after execution if it was open before, closes otherwise.
   "MemoryMonitor_Launcher.jar" can be downloaded here:
   https://github.com/MicPhotonIGBMC/ImageJ-Macros/blob/master/Metamorph%20Utilities/MemoryMonitor_Launcher.jar

Known problems:
- Repeated use of the macro results in an increasing memory occupation because it cannot
  close BDP2 viewers (no viewer-closing macro command found in BDP2 0.5.7).
  Closing the viewers manually saves some memory, but becovery of its basal level can be
  obtained only by restarting Fiji.
