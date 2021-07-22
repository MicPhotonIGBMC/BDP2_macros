"Extract_Regions_with_BDP2_20":

Extracts regions from images that are too big to be opened in non virtual
mode in Fiji:

1. Choose input format and working-mode:
   Check "¤ Crop whole CZT domain" to crop all Channels, Z and T (Simple mode);
   Otherwise (Full mode), you can define different C, Z and T begins and ends for each crop.
2. Choose input image and output folder;
3. In the Bio-formats Dialog box choose virtual mode;
4. In the infinite crop-regions loop:
   For each region:
   ¤ In Full mode:
     - draw a rectangle, select first slice, channel and time-point;
     - validate by OK;
     - select last slice, channel and timepoint;
       validate by OK to define another region;
                or SHIFT-OK to terminate;
   ¤ In Simple mode :
     - draw a rectangle;
     - validate by OK to define another region;
                or SHIFT-OK to terminate;
5. Crop-ROIs are saved as a zip to output folder;
6. The image is opened using BigDataProcessor2;
7. Each region is cropped in a new viewer and saved to output folder;

Tested with following input types:
 - Metamorph Multi Dimensional series (input file: .nd);
 - Leica Image Files (.lif);
 - Zeiss (.czi) pyramidal images need the crop domains to be defined in the
   highest resolution image (#1). Each crop results in an pyramidal image file
   containing a number of resolutions depending on the XY size of the crop- 
   regions. The resolution sequence of output images is 1, 1/2, 1/4, ...;
 - Single plane (Z*C*T = 1) TIFFs, PNG.

Notes:
1. In case Fiji is slow, the SHIFT key must be kept down until the image closes
   to terminate the crop-regions definition to avoid the creation of a new one.
2. The macro launches the "Memory" window if not open and "MemoryMonitor_Launcher.jar" 
   is found in the "plugins" folder, works without otherwise.
   "MemoryMonitor_Launcher.jar" can be downloaded here:
   https://github.com/MicPhotonIGBMC/ImageJ-Macros/blob/master/Metamorph%20Utilities/MemoryMonitor_Launcher.jar

Known problems:
- In "Single Plane TIFF" image-type, the output image names are those of the folder
  containing the images. Use "TIFF stack" or "METAMORPH, CZI, LIF, HDF5 ..." instead.
- Repeated use of the macro results in an increasing memory occupation because it cannot
  close the viewers (no viewer-closing macro command found in BDP2 0.5.7).
  Closing the viewers manually saves some memory, but recovery of its basal level requires
  to restart Fiji.
