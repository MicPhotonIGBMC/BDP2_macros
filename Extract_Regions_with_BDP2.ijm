/** Extract_Regions_with_BDP2_
 * Marcel Boeglin 
 * boeglin@igbmc.fr
 * May - July 2021
 */

// 2021-07-16 : BigDataProcessor2 no longer working with ND files (.nd)

/**
Extracts regions from images that are too big to be opened in non virtual
mode in Fiji:
1. Choose input format and working-mode:
   Check "¤ Crop whole CZT domain" to crop all Channels, Z and T (Simple mode);
   Otherwise (Full mode), you can define different C, Z and T begins and ends for
   each crop.
2. Choose input image and output folder;
3. If Bio-formats Dialog box opens, check "¤ Use virtual stack";
4. Regions creation loop:
   ¤ Full mode:
     - draw a rectangle, select first slice, channel and time-point;
     - validate by OK;
     - select last slice, channel and timepoint;
       validate by OK to define another region;
                or SHIFT-OK to terminate;
   ¤ Simple mode :
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
   highest resolution image (#1). Each crop results in pyramidal image file
   containing a number of resolutions depending on the XY size of the crop- 
   regions. The resolution sequence of output images is 1, 1/2, 1/4, ...;
 - Single plane (Z*C*T = 1) TIFFs, PNG.

Notes:
1. Termination the crop-regions creation:
   Keep SHIFT key down until the image closes to avoid creation of a new region.
2. The macro launches the "Memory" window if not open and "MemoryMonitor_Launcher.jar" 
   is found in the "plugins" folder, works without otherwise.
   "MemoryMonitor_Launcher.jar" can be downloaded here:
   https://github.com/MicPhotonIGBMC/ImageJ-Macros/blob/master/Metamorph%20Utilities/MemoryMonitor_Launcher.jar
 */

/* Known problems:
 * 1.
 * In case of multiple series in the input file, the macro may process another series
 * than the one selected in the Bio-Formats window.
 * This is due to the fact that the series-names in Bio-Formats may differ strongly 
 * according to the experiment type, leading to a false series-index determination
 * 2.
 * BDP2 does not open properly time-series from LIF files:
 * generally opens only the first time-point sometimes repeated in a few frames
 *
*/
var macroname = "Extract_Regions_with_BDP2_";
var version = 25;
var author = "Marcel Boeglin";
var email = "boeglin@igbmc.fr";

var width, height, channels, slices, frames;
var stackSize;
//keep memory of previous slider positions to reuse in next crop-region
var previousMinC, previousMinZ, previousMinT;
var previousMaxC, previousMaxZ, previousMaxT;

var crop_whole_CZT_domain = true;

var minx, miny, minz, minc, mint;
var maxx, maxy, maxz, maxc, maxt;

var memoryWasOpen = isOpen("Memory");

var inputFormats = newArray("nd, CZI, LIF, HDF5 ...",
							"TIFF stack",
							"Single plane TIFF");
var inputFormat = "nd, CZI, LIF, HDF5 ...";
var inputDir, inputFileName, extension, inputPath;
var inputImage, inputImageName;
var isMetamorph;
var inputImageTitle;
var outputDir;
var seriesindex;
var seriesName;

var BigDataViewerAlreadyOpen = true;

//Macro BEGIN

IJ.redirectErrorMessages();
run("Bio-Formats Macro Extensions");

roiManager("Associate", "true");
roiManager("Centered", "false");
roiManager("UseNames", "false");

execute();
logPath = outputDir+inputImageTitle+"_LOG.txt";
selectWindow("Log");
saveAs("text", logPath);
if (memoryWasOpen) exit;
closeMemoryMonitor();

//Macro END


function getInputFormat_And_Working_Mode() {
	Dialog.create("Input format and working mode");
//	Dialog.addCheckbox("BigDataViewer is already open", BigDataViewerAlreadyOpen);
//	BDP2 opens a new viewer even if "viewingmodality=[Show in current viewer] "
	Dialog.addChoice(inputFormat, inputFormats, inputFormat);
	Dialog.addCheckbox("Crop_whole_CZT_domain", crop_whole_CZT_domain);
	Dialog.show();
//	BigDataViewerAlreadyOpen = Dialog.getCheckbox();
	inputFormat = Dialog.getChoice();
	crop_whole_CZT_domain = Dialog.getCheckbox();
	print("\ninputFormat = "+inputFormat);
	print("crop_whole_CZT_domain = "+crop_whole_CZT_domain);
}

function getInputPath()  {
	inputPath = File.openDialog("Open Image as a Virtual stack");
	inputDir = File.getDirectory(inputPath);
	inputFileName = File.getNameWithoutExtension(inputPath);
	extension = substring(inputPath, lastIndexOf(inputPath, "."));
}

function execute() {
	print("\\Clear");
	getInputFormat_And_Working_Mode();
	closeExceptionWindows();
	close("*");
	print(macroname+version+".ijm");
	print(macroDescription());
	print(author);
	print(email);
	getInputPath();
	print("\ninputDir = "+inputDir);
	print("inputFileName = "+inputFileName);
	print("extension = "+extension);
	isMetamorph = (extension==".nd");
	print("isMetamorph = "+isMetamorph);
	lowercaseExtension = toLowerCase(extension);
	outputDir = getDirectory("Destination Directory ");
	print("outputDir = "+outputDir);

	msg = "\nPyramidal images:\n"+
		"Open highest resolution:"+
		"\nFor single series or multi-series 1st series:"+
		"\nZeiss CZI :  #1"+
		"\nHDF5 :         _0\n ";
	print(msg);

	launchMemoryMonitor();

	//Open image
	if (inputFormat=="nd, CZI, LIF, HDF5 ...") {//MAY FAIL
		seriesIndex = 0;
		open(inputPath);//opens Bio-Formats window
		title = getTitle();
		print("title = "+title);
		imageID = getImageID();
		if (indexOf(title, " - ") >= 0)  {
			//seriesName = substring(title, lastIndexOf(title, " - ")+3);
			seriesName = substring(title, indexOf(title, " - ")+3);
			print("seriesName = "+seriesName);
			if (isMetamorph) // multiple stage positions
				seriesindex = getSeriesIndex(imageID);
			else {
				seriesindex = findSeriesIndex(inputDir,
					inputFileName+extension, seriesName);
			}
		}
		//exit;
	}
	else if	(inputFormat=="TIFF stack") {
		//print(extension);
		if (toLowerCase(extension)==".zip") exit("ZIP files not supported");
		run("TIFF Virtual Stack...", "open=["+inputPath+"]");
	}
	else if	(inputFormat=="Single plane TIFF") {
/*
		print("\ninputDir = "+inputDir);
		print("inputFileName+extension = "+inputFileName+extension+"\n ");
		run("Image Sequence...", "select=["+inputDir+"] "+
			" filter=["+inputFileName+extension+"] count=1 sort use");
		getDimensions(w, h, c, d, t);
		if (c * d * t > 1) exit("1D image required");
*/
		//image title may be image path instead of image name
		run("Bio-Formats Importer", "open=["+inputDir+inputFileName+extension+"] "+
			"color_mode=Composite rois_import=[ROI manager] view=Hyperstack "+
			"stack_order=XYCZT use_virtual_stack series_1");
	
		rename(inputFileName+extension);
	}
	if (nImages<1) exit("No image: Macro aborted");
	imageID = getImageID();
/*
	if (isMetamorph)
		seriesindex = getSeriesIndex(imageID);
	else
		seriesindex = findSeriesIndex(inputDir, inputFileName+extension, inputImageTitle);
*/
//	seriesindex = findSeriesIndex(inputDir, inputFileName+extension, seriesName);

	Stack.getDimensions(width, height, channels, slices, frames);
	print(width);
	stackSize = channels*slices*frames;
	inputImageTitle = getTitle();
	print("seriesindex : "+seriesindex);
	print("inputImageTitle : "+inputImageTitle);
	print("inputFileName+extension : "+inputFileName+extension);
	inputImageTitle = replace(inputImageTitle, "\"", "");
	inputImageTitle = replace(inputImageTitle, "/", "-");

	if (inputImageTitle==inputFileName+extension) {
		//inputImageTitle = seriesName;
		inputImageTitle += "_"+seriesName;
	}

	imageID = getImageID();
	if (channels>1)
		Stack.setDisplayMode("composite");
	s = slices/2; if (s<1) s=1; if (slices>1) Stack.setSlice(s);
	t = frames/2; if (t<1) t=1; if (frames>1) Stack.setFrame(t);
	if (channels>1) {
		for (c=1; c<=channels; c++) {
			Stack.setChannel(c); run("Enhance Contrast", "saturated=0.25");
		}
	}
	roiManager("Associate", "false");
	roiManager("Centered", "false");
	roiManager("UseNames", "false");
	roiManager("reset");
	setTool("rectangle");
	nmax = 40;
	previousMinC=1; previousMinZ=1; previousMinT=1;
	previousMaxC=channels; previousMaxZ=slices; previousMaxT=frames;
	if (stackSize>1)
		Stack.setPosition(previousMinC, previousMinZ, previousMinT);
	makeRectangle(3*width/8, 3*height/8, width/4, height/4);
	Roi.setStrokeColor(0, 255, 0);
	i=0;
	while (true) {
		if (++i > nmax) break;
		if (isKeyDown("shift\")")) break;
		addRegionToManager();
	}
//	close();
	roiManager("deselect");
	roiManager("save", outputDir+inputImageTitle+"_Crop-Rois.zip");
	//Process Crop-Rois from Roi Manager
	nregions = roiManager("count")/2;//each region is defined by 2 Rois


//PROBLEMES POTENTIELS : pas les bonnes limites des domaines de Crop

	//ESSAI : ne pas fermer l'image a traiter pour ne pas avoir a 
	//en creer une autre forcement non virtuelle (PROBLEMEs aleatoires de bornes de crop
	close();//close image used for Crop-Rois drawing


	//IMAGE AUXILIAIRE pour pouvoir afficher les rois (ACTIVER ROI DANS RoiManager NECESSITE UNE IMAGE)
	//newImage("HyperStack", "8-bit", width, height, 1, 1, 1);//w,h,c,z,t
	//if false z or t range in output, use command below:
//	newImage("HyperStack", "8-bit", width, height, channels, slices, frames);
	//newImage("HyperStack", "8-bit", 1, 1, channels, slices, frames);
	//newImage("HyperStack", "8-bit", width, height, 1, 1, 1);
	newImage("Untitled", "8-bit", width, height, 1);//ESSAI
	tmpid = getImageID();
	size = nregions;
	minx = newArray(size); miny = newArray(size); minz = newArray(size);
	minc = newArray(size); mint = newArray(size);
	maxx = newArray(size); maxy = newArray(size); maxz = newArray(size);
	maxc = newArray(size); maxt = newArray(size);
	for (r=0; r<nregions; r++) {
		print("");
		//i = 2*r;
		roiManager("select", 2*r);
		Roi.getBounds(minX, minY, w, h);
		//print("minX = "+minX);
		//print("minY = "+minY);
		//Roi.getBounds(minx[r], miny[r], w, h);
		minx[r] = minX; miny[r] = minY;
		//Roi.getPosition(channel, slice, frame);
		Roi.getPosition(mincPlus1, minzPlus1, mintPlus1);
		minc[r] = mincPlus1-1;
		minz[r] = minzPlus1-1;
		mint[r] = mintPlus1-1;
		print("minx="+minx[r]+"  miny="+miny[r]+"  minz="+minz[r]+
			"  minc="+minc[r]+"  mint="+mint[r]);
		roiManager("select", 2*r+1);
		Roi.getPosition(maxcPlus1, maxzPlus1, maxtPlus1);
		maxc[r] = maxcPlus1-1;
		maxz[r] = maxzPlus1-1;
		maxt[r] = maxtPlus1-1;
		maxx[r] = minx[r]+w;
		maxy[r] = miny[r]+h;
		print("maxx="+maxx[r]+"  maxy="+maxy[r]+"  maxz="+maxz[r]+
			"  maxc="+maxc[r]+"  maxt="+maxt[r]);
	}
	selectImage(tmpid);
	close();
	print("inputPath:");
	print(inputPath);
/*
	print("\nBigDataViewerAlreadyOpen = "+BigDataViewerAlreadyOpen);
	if (BigDataViewerAlreadyOpen) {//opens a new viewer anyway
		print("viewingmodality=[Show in current viewer] ");
		run("BDP2 Open Bio-Formats...", "viewingmodality=[Show in current viewer] "+
			"enablearbitraryplaneslicing=true file=["+inputPath+"] seriesindex="+
			seriesindex);
	}
	else
*/
	{
		print("\n \nviewingmodality=[Show in new viewer] ");
		print("seriesindex = "+seriesindex);
		run("BDP2 Open Bio-Formats...", "viewingmodality=[Show in new viewer] "+
			"enablearbitraryplaneslicing=true file=["+inputPath+"] seriesindex="+
			seriesindex);
	}

	/* Resulting BigDataViewer can't be closed using an ImageJ command because
	 * it's not added to the ImageJ windows using
	 * ij.WindowManager.addWindow(java.awt.Frame win) 
	 * or
	 * ij.WindowManager.addWindow(java.awt.Window win)
	 */
	inputimage = inputFileName;
	print("\ninputimage = "+inputimage);

	for (r=0; r<nregions; r++) {
		print("");
		print("minx="+minx[r]+"  miny="+miny[r]+"  minz="+minz[r]+
			"  minc="+minc[r]+"  mint="+mint[r]);
		print("maxx="+maxx[r]+"  maxy="+maxy[r]+"  maxz="+maxz[r]+
			"  maxc="+maxc[r]+"  maxt="+maxt[r]);
		outputimagename = inputImageTitle+"_Crop_"+r;

		// BDP2 Crop throws an exception if input image is a TIFF stack
		// but output looks OK

		run("BDP2 Crop...", "inputimage=["+inputimage+
			"] outputimagename=["+outputimagename+
			"] viewingmodality=[Show in current viewer]"+
			" minx="+minx[r]+" miny="+miny[r]+" minz="+minz[r]+
				" minc="+minc[r]+" mint="+mint[r]+""+
			" maxx="+maxx[r]+" maxy="+maxy[r]+" maxz="+maxz[r]+
				" maxc="+maxc[r]+" maxt="+maxt[r]+" ");

/*
//ORIGINAL
		run("BDP2 Crop...", "inputimage=["+inputimage+
			"] outputimagename=["+outputimagename+
			"] viewingmodality=[Show in new viewer]"+
			" minx="+minx[r]+" miny="+miny[r]+" minz="+minz[r]+
				" minc="+minc[r]+" mint="+mint[r]+""+
			" maxx="+maxx[r]+" maxy="+maxy[r]+" maxz="+maxz[r]+
				" maxc="+maxc[r]+" maxt="+maxt[r]+" ");
*/

// Commentarize for Debug
		run("BDP2 Save As...", "inputimage=["+outputimagename+
			"] directory=["+outputDir+"] numiothreads=1 numprocessingthreads=4"+
			" filetype=[BigDataViewerXMLHDF5] saveprojections=false"+
			" savevolumes=true tiffcompression=[None]"+
			" tstart="+mint[r]+" tend="+mint[r]);
// END Commentarize for Debug

		closeBDPViewer(outputimagename);//Doesn't work
	}
	closeBDPViewer(inputimage);//Doesn't work
}

function addRegionToManager() {
	if (crop_whole_CZT_domain) {
	//region begin :
		Stack.setPosition(1, slices/2, frames/2);
		msg = "Draw a rectangle"+
			"\nPress OK to validate and add another Region"+
			"\nPress Shift-OK to finish";
		Roi.setStrokeColor("green");
		waitForUser(msg);
		Roi.setPosition(1, 1, 1);
		roiManager("add");
		roiManager("select", roiManager("count")-1);
		roiNum = (roiManager("count")-1)/2 + 1;
		roiNumStr = String.pad(roiNum, 2);
		roiManager("Rename", roiNumStr+"_BEGIN");
		RoiManager.setGroup(0);
		roiManager("Set Color", "green");
		roiManager("Set Line Width", 0);
		roiManager("deselect");
		print("\nframes = "+frames+"\n ");
	//region end :
		Roi.setPosition(channels, slices, frames);
		roiManager("add");
		roiManager("select", roiManager("count")-1);
		roiManager("Rename", roiNumStr+"_END");
		RoiManager.setGroup(0);
		roiManager("Set Color", "red");
		roiManager("Set Line Width", 0);
	}
	else {//set CZT domain for each Roi
	//region begin :
		Stack.setPosition(previousMinC, previousMinZ, previousMinT);
		msg = "Draw a rectangle, Select 1st Slice, Channel and Frame\n"+
			"Press OK to validate";
		if (stackSize<2)
			msg = "Draw a rectangle\n"+
				"Press OK to add a new Region\nPress Shift-OK to finish";
		Roi.setStrokeColor("green");
		waitForUser(msg);
		Stack.getPosition(channel, slice, frame);
		print("BEGIN : channel="+channel+"  "+"slice="+slice+"  "+"frame="+frame);
		Roi.setPosition(channel, slice, frame);
		previousMinC=channel; previousMinZ=slice; previousMinT=frame;
		roiManager("add");
		roiManager("select", roiManager("count")-1);
		roiNum = (roiManager("count")-1)/2 + 1;
		roiNumStr = String.pad(roiNum, 2);
		roiManager("Rename", roiNumStr+"_BEGIN");
		RoiManager.setGroup(0);
		roiManager("Set Color", "green");
		roiManager("Set Line Width", 0);
	//region end :
		Roi.setStrokeColor("red");
		roiManager("deselect");
		Stack.setPosition(previousMaxC, previousMaxZ, previousMaxT);
		if (stackSize>1)
			waitForUser("Select last Slice and Frame of region\n"+
				"DO NOT REMOVE ROI\n"+
				"Press OK to add a new Region\nPress Shift-OK to finish");
		Stack.getPosition(channel, slice, frame);
		print("END : channel="+channel+"  "+"slice="+slice+"  "+"frame="+frame);
		Roi.setPosition(channel, slice, frame);
		previousMaxC=channel; previousMaxZ=slice; previousMaxT=frame;
		roiManager("add");
		roiManager("select", roiManager("count")-1);
		roiManager("Rename", roiNumStr+"_END");
		RoiManager.setGroup(0);
		roiManager("Set Color", "red");
		roiManager("Set Line Width", 0);
	}
}


/**
 * Doesn't work because BDP2 viewers are not ImageJ windows
 */
function closeBDPViewer(imagename) {//Doesn't work
	return;
	nonimageWindows = getList("window.titles");
	if (nonimageWindows.length<1) return;
	print("\nCurrent Non Image Windows:");
	for (i=0; i<nonimageWindows.length; i++) {
		wintitle = nonimageWindows[i];
		print(wintitle);
		//if (wintitle != imagename) continue;
		if (wintitle == imagename) {
			selectWindow(imagename);
			run("Close");
		}
	}
}

function closeExceptionWindows() {
	nonimageWindows = getList("window.titles");
	if (nonimageWindows.length<1) return;
	print("\nCurrent Non Image Windows:");
	//Array.print(nonimageWindows);
	for (i=0; i<nonimageWindows.length; i++) {
		wintitle = nonimageWindows[i];
		print(wintitle);
		//if (wintitle != imagename) continue;
		if (wintitle == "Exception") {
			selectWindow(wintitle);
			run("Close");
		}
	}
}

/**
 * Current version works only with Metamorph Multi Dimensional Series
 * Returns 0 for other input formats.
 */
function getSeriesIndex(imageID) {
	position = "0";
	selectImage(imageID);
	title = getTitle();
	if (isMetamorph) {
		if (indexOf(title, "Stage")<0) return 0;
		index = -1;
		if (indexOf(title, "Stage") >=0)
			index = indexOf(title, "Stage");
		str = substring(title, index);
		print("str = "+str);
		index2 = indexOf(str, " \"");
		print("index2 = "+index2);
		position = substring(str, 5, index2);
		print("position = "+position);
		return parseInt(position) - 1;
	}
	return parseInt(position);
}

/** Requires run("Bio-Formats Macro Extensions");
 * 'file' filename with extension
 * Adapted from 'extract_series_from_leica_lif.ijm' macro
 * from Erwan Grandgirard and Bertrand Vernay 
 * Does not work with HDF5 and Metamorph ND files
 **/
function findSeriesIndex(dir, file, seriesname) {
	print("findSeriesIndex(dir, file),\nfile = "+file);
	failed = "findSeriesIndex(dir, file) not adapted for HDF5 or ND files";
	if (endsWith(file, ".h5")) {print(failed); return;}
	if (endsWith(file, ".nd")) {print(failed); return;}
	path=dir+file;
	print("file = "+file);
	print("seriesname = "+seriesname);
	Ext.setId(path);
	Ext.getCurrentFile(file);
	Ext.getSeriesCount(seriesCount);//gets the number of series in input file
	//print("Processing the file = " + file);
// See:
//http://imagej.1557.x6.nabble.com/multiple-series-with-bioformats-importer-td5003491.html
	//while next size is a fraction of current (1/3 for CZI, 1/2 for BDP2 HDF5
	//the images belon to the same pyramidal series with different resolutions
	index = 0;
	Ext.getUsedFileCount(count);
	print("usedFileCount = "+count);
	for (j=0; j<seriesCount; j++) {
		//print("Extracting Series "+j+1+" / "+seriesCount);
		Ext.setSeries(j);
		Ext.getSeriesName(seriesName);
		print("seriesName = "+seriesName);
		print("index = "+index);
//		Ext.getUsedFileCount(count);
		Ext.getCurrentFile(currentFile);
		//print("file = "+currentFile);
		Ext.getSizeX(sizeX);
        Ext.getSizeY(sizeY);
		//print("sizeX = "+sizeX+"    sizeY = "+sizeY);
		//seriesName = replace(seriesName, "\"", "");
		//if (seriesName == imageTitle) break;
		//if (seriesName == seriesname) break;
		//if (seriesName == seriesname) return j;
		if (endsWith(seriesname, seriesName)) return j;
		index ++;
	}
	print("findSeriesIndex(dir, file, imageTitle) END");
	//inputImageTitle = seriesName;//case only 1 series in input file
	return index;
}

/*
 * from extract_series_from_leica_lif_split.ijm,
 * Erwan Grandgirard and Bertrand Vernay
 */
function chooseImageToProcess(path) {
	//from extract_series_from_leica_lif_split.ijm, 
	Ext.setId(path);
	Ext.getCurrentFile(path);
	Ext.getSeriesCount(seriesCount); // this gets the number of series
	if (seriesCount==1) 
	print("Processing the file = " + fileToProcess);
	for (j=0; j<seriesCount; j++) {
        Ext.setSeries(j);
        Ext.getSeriesName(seriesName);
		run("Bio-Formats Importer", "open=&path color_mode=Default "+
			"view=Hyperstack stack_order=XYCZT series_"+j+1); 
		fileNameWithoutExtension = File.nameWithoutExtension;
		//print(fileNameWithoutExtension);
		saveAs("tiff", dir2+fileNameWithoutExtension+"_"+seriesName+".tif");
		run("Close");
	}
}

function macroDescription() {
	s = "This macro uses BigDataProcessor2 to extract regions \n"+
	"from big images to make them openable in ImageJ.";
	return s;
}

function launchMemoryMonitor() {
	if (isOpen("Memory")) return;
	pluginsDir = getDirectory("plugins");
	if (findFile(pluginsDir, "MemoryMonitor_Launcher.jar"))
		run("MemoryMonitor Launcher");
	else
		print("\nLaunch \"Monitor Memory...\" in a macro"+
			"\nrequires \"MemoryMonitor_Launcher.jar\""+
			"\nto be installed in plugins folder\n ");
}

function closeMemoryMonitor() {
	if (!isOpen("Memory")) return;
	selectWindow("Memory");
	run("Close");
}

function getFiles(dir) {
	list = getFileList(dir);
	if (list.length==0) {
		showMessage(macroName,"Input folder:\n"+dir+"\nseems empty");
		exit();
	}
	j=0;
	list2 = newArray(list.length);
	for (i=0; i<list.length; i++) {
		s = list[i];
		if (File.isDirectory(dir+s)) continue;
		skip = false;
		for (k=0; k<projPrefixes.length; k++)
			if (startsWith(s, projPrefixes[k])) {
				skip = true;
				break;
			}
		if (skip) continue;
		list2[j++] = s;
	}
	if (j<1) {
		showMessage(macroName,"Input folder:\n"+dir+
			"\nseems not to contain Metamorph images");
		exit();
	}
	for (i=0; i<list2.length; i++) {
		list2[i] = toString(list2[i]);
	}
	list2 = Array.trim(list2, j);
	return Array.sort(list2);
}

/** Returns TIFF and STK files contained in 'list' matching fileFilter and 
	not matching excludingFilter */
function filterList(list, fileFilter, excludingFilter) {
	list2 = newArray(list.length);
	j=0;
	for (i=0; i<list.length; i++) {
		s = list[i];
		if (fileFilter!="" && indexOf(s, fileFilter)<0) continue;
		if (excludingFilter!="" && indexOf(s, excludingFilter)>=0) continue;
		s2 = toLowerCase(s);
		ext = getExtension(s);
		if (!endsWith(s2, ".tif") && !endsWith(s2, ".stk")) continue;
		list2[j++] = s;
	}
	if (j<1) {
		showMessage(macroName,
				"Input folder seems not to contain TIFF or STK files "+
				"matching "+fileFilter);
		exit();
	}
	for (i=0; i<list2.length; i++) {
		list2[i] = toString(list2[i]);
	}
	list2 = Array.trim(list2, j);
	list2 = Array.sort(list2);
	return list2;
}

function ndFileNames(filenames) {
	dbg = true;
	nFiles = filenames.length;
	ndNames = newArray(nFiles);
	j=0;
	for (i=0; i<nFiles; i++) {
		fname = filenames[i];
		if (endsWith(fname, ".nd"))
			ndNames[j++] = substring(fname,0,lastIndexOf(fname,".nd"));
	}
	if (dbg) for (i=0; i<j; i++) print("ndNames["+i+"] = "+ndNames[i]);
	return Array.trim(ndNames, j);
}

function findFile(dir, filename) {
	lst = getFileList(dir);
	for (i=0; i<lst.length; i++) {
		if (File.isDirectory(""+dir+lst[i]))
			findFile(""+dir+lst[i], filename);
		else {
			if (lst[i]==filename) return true;
		}
	}
	return false;
}

//80 chars:
//23456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789