/** Extract_Regions_with_BDP2_
 * Authors: 
 * Marcel Boeglin: boeglin@igbmc.fr
 * Bertand Vernay: vernayb@igbmc.fr
 * May 2021
 */

/**
 * Procedure:
 * 1. Open image as a virtual stack;
 * 2. Multi-Channel images: select best channel to choose regions to be extracted 
 *    or activate Composite mode;
 * 3. Infinite loop of:
 *    Draw a rectangle surounding the region to be extracted;
 *    WaitForUser 1:
 * 	  Select first time-point and z-position to be extracted;
 *    The Begin Roi is automatically added to Roi Manager
 *    WaitForUser 2:
 *    Select last time-point and z-position to be extracted;
 *    The End Roi is automatically added to Roi Manager
 * 4. Break loop by pressing "OK" with "Shift" key down
 * 5. Close virtual image;
 * 6. Save RoiSet to output folder
 * 7. Open image with BigDataProcessor2
 * 8. For each Roi couple in Roi Manager:
 *    - get xmin, ymin, zmin, cmin, tmin, xmax, ymax, zmax, cmax tmax
 *    - run("BDP2 Crop...", "inputimage=..., xmin, ymin, zmin, cmin tmin, 
 *      xmax, ymax, zmax, cmax, tmax); in a new window
 *    - save crop to output folder
 *    - close crop
 */


var macroname = "Extract_Regions_with_BDP2_";
var version = 08;
var copyRight = "Authors: Marcel Boeglin - Bertrand Vernay,  May 2021";
var email = "e-mail: boeglin@igbmc.fr";

var width, height, channels, slices, frames;

//keep memory of previous slider positions
var previousMinC, previousMinZ, previousMinT;
var previousMaxC, previousMaxZ, previousMaxT;

var minx, miny, minz, minc, mint;
var maxx, maxy, maxz, maxc, maxt;

/*
//not used
var doRangeAroundMedianSlice = false;
var rangeAroundMedianSlice = 50;// % of stack-size
var timepoints;
var firstTimePoint = 0;//in BDP, time range is [0, timepoints-1]
var lastTimePoint;
var medianTimepoint;
var rangeAroundMedianTimePoint = 50;// % of timepoints
var doRangeAroundMedianTimePoint = false;// if false do all time points
var resize = false;
var resizeFactor = 0.5;
*/
var memoryWasOpen = isOpen("Memory");
var inputPath, inputDir, inputFile, inputFileName;
var inputImage, inputImageName;
var isMultiSeries;
var isMetamorph;
var inputImageTitle;
var outputDir;
var seriesindex;


//Macro BEGIN

run("Bio-Formats Macro Extensions");
execute();
if (memoryWasOpen) exit;
closeMemoryMonitor();

//Macro END


function getInputPath() {
	inputPath = File.openDialog("Open Image as a Virtual stack");
	inputDir = File.getDirectory(inputPath);
	inputFile = File.getName(inputPath);
	inputFileName = File.getNameWithoutExtension(inputPath);
}

function execute() {
	print("\\Clear");
	close("*");
	print(macroname+version+".ijm");
	print(macroDescription());
	print(copyRight);
	print(email);
	getInputPath();
	print("\ninputDir = "+inputDir);
	print("inputFile = "+inputFile);
	print("inputFileName = "+inputFileName);
	isMetamorph = endsWith(inputFile, ".nd");
	print("isMetamorph = "+isMetamorph);
	
	outputDir = getDirectory("Destination Directory ");
	print("outputDir = "+outputDir);
	launchMemoryMonitor();
	open(inputPath);
	imageID = getImageID();
	//if (isMetamorph)
		seriesindex = getSeriesIndex(imageID);
	if (nImages<1) {showMessage("No image: Macro aborted"); return;}
	inputImageTitle = getTitle();
	print("inputImageTitle : "+inputImageTitle);
	inputImageTitle = replace(inputImageTitle, "\"", "");
	imageID = getImageID();
	Stack.setDisplayMode("composite");
	Stack.getDimensions(width, height, channels, slices, frames);
	s = slices/2; if (s<1) s=1; Stack.setSlice(s);
	t = frames/2; if (t<1) t=1; Stack.setFrame(t);
	for (c=1; c<=channels; c++) {
		Stack.setChannel(c); run("Enhance Contrast", "saturated=0.25");
	}
	roiManager("Associate", "false");
	roiManager("Centered", "false");
	roiManager("UseNames", "false");
	roiManager("reset");
	setTool("rectangle");
	nmax = 20;
	previousMinC=1; previousMinZ=1; previousMinT=1;
	previousMaxC=channels; previousMaxZ=slices; previousMaxT=frames;
	Stack.setPosition(previousMinC, previousMinZ, previousMinT);
	makeRectangle(3*width/8, 3*height/8, width/4, height/4);
	i=0;
	while (true) {
		if (++i > nmax) break;
		if (isKeyDown("shift")) break;
		wait(1000);
		addRegionToManager();
		wait(1000);
	}
	//close();
	roiManager("deselect");
	roiManager("save", outputDir+inputImageTitle+"_Crop-Rois.zip");
	
	//Process Crop-Rois from Roi Manager
	nregions = roiManager("count")/2;//each region is defined by 2 Rois
	close();//close image used for Crop-Rois drawing
	//newImage("Tmp", "8-bit", 1, 1, 1);// -> false RoiX and RoiY
	newImage("Tmp", "8-bit", width, height, 1);
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
		//print("minx="+minx[r]+"  miny="+miny[r]+"  minz="+minz[r]+
		//	"  minc="+minc[r]+"  mint="+mint[r]);

		roiManager("select", 2*r+1);
		Roi.getPosition(maxcPlus1, maxzPlus1, maxtPlus1);
		maxc[r] = maxcPlus1-1;
		maxz[r] = maxzPlus1-1;
		maxt[r] = maxtPlus1-1;
		maxx[r] = minx[r]+w;
		maxy[r] = miny[r]+h;
		//print("maxx="+maxx[r]+"  maxy="+maxy[r]+"  maxz="+maxz[r]+
		//	"  maxc="+maxc[r]+"  maxt="+maxt[r]);
	}
	selectImage(tmpid);
	close();

	print("inputPath");
	print(inputPath);


//	Can't be closed using an ImageJ command
	run("BDP2 Open Bio-Formats...", "viewingmodality=[Show in new viewer] enablearbitraryplaneslicing=true file="
	+inputPath+" seriesindex="+seriesindex);

	/*
	run("BDP2 Open Bio-Formats...", "viewingmodality=[Do not show] enablearbitraryplaneslicing=true file="
	+inputPath+" seriesindex="+seriesindex);
	//image can be closed only by closing Fiji
	*/
	//Attempts to get image or viewer title to be able to close it
	//inputimage = getTitle();//not an ImageJ image
	inputimageWindowTitle = getInfo("window.title");
	print("inputimageWindowTitle = "+inputimageWindowTitle);
	/*
	 * Conclusion: When a new viewer is created in BDP2, it's not added to the ImageJ windows using
	 * ij.WindowManager.addWindow(java.awt.Frame win) 
	 * or
	 * ij.WindowManager.ddWindow(java.awt.Window win)
	 * and can't therefore not be closed using any ImageJ command
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
		run("BDP2 Crop...", "inputimage=["+inputimage+
			"] outputimagename=["+outputimagename+
			"] viewingmodality=[Show in new viewer]"+
			" minx="+minx[r]+" miny="+miny[r]+" minz="+minz[r]+
				" minc="+minc[r]+" mint="+mint[r]+""+
			" maxx="+maxx[r]+" maxy="+maxy[r]+" maxz="+maxz[r]+
				" maxc="+maxc[r]+" maxt="+maxt[r]+" ");
		run("BDP2 Save As...", "inputimage=["+outputimagename+
			"] directory=["+outputDir+"] numiothreads=1 numprocessingthreads=4"+
			" filetype=[BigDataViewerXMLHDF5] saveprojections=false"+
			" savevolumes=true tiffcompression=[None]"+
			" tstart="+mint[r]+" tend="+mint[r]);
		closeBDPViewer(outputimagename);
	}
	closeBDPViewer(inputimage);
}

/**
 * Doesn't work because BDP2 viewers are not ImageJ windows
 */
function closeBDPViewer(imagename) {
	nonimageWindows = getList("window.titles");
	if (nonimageWindows.length<1) return;
	print("\nCurrent Non Image Windows:");
	//Array.print(nonimageWindows);
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

function getSeriesIndex(imageID) {
	position = "0";
	selectImage(imageID);
	title = getTitle();
	if (isMetamorph) {
		index = -1;
		if (indexOf(title, "Stage") >=0)
			index = indexOf(title, "Stage");
		str = substring(title, index);
		print("str = "+str);
		index2 = indexOf(str, " \"");
		print("index2 = "+index2);
		//position = substring(str, 0+8, index2);
		position = substring(str, 5, index2);
		print("position = "+position);
		return parseInt(position) - 1;
	}
	return parseInt(position);
}

function chooseImageToProcess(path) {
	Ext.setId(path);
	Ext.getCurrentFile(fileToProcess);
	Ext.getSeriesCount(seriesCount); // this gets the number of series
	if (seriesCount==1) 
	print("Processing the file = " + fileToProcess);
	for (j=0; j<seriesCount; j++) {
        Ext.setSeries(j);
        Ext.getSeriesName(seriesName);
		run("Bio-Formats Importer", "open=&path color_mode=Default view=Hyperstack stack_order=XYCZT series_"+j+1); 
		fileNameWithoutExtension = File.nameWithoutExtension;
		//print(fileNameWithoutExtension);
		saveAs("tiff", dir2+fileNameWithoutExtension+"_"+seriesName+".tif");
		run("Close");
	}
}

function addRegionToManager() {
//region begin
	Stack.setPosition(previousMinC, previousMinZ, previousMinT);
		wait(1000);
	waitForUser("Draw a rectangle, Select 1st Slice, Channel and Frame\n"+
		"Press OK to validate");
	Stack.getPosition(channel, slice, frame);
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
//region end
	roiManager("deselect");
	Stack.setPosition(previousMaxC, previousMaxZ, previousMaxT);
		wait(1000);
	waitForUser("Select last Slice and Frame of region\n"+
		"DO NOT REMOVE ROI\n"+
		"Press OK to validate\nPress Shift-OK to finish");
	Stack.getPosition(channel, slice, frame);
	Roi.setPosition(channel, slice, frame);
	previousMaxC=channel; previousMaxZ=slice; previousMaxT=frame;
	roiManager("add");
	roiManager("select", roiManager("count")-1);
	roiManager("Rename", roiNumStr+"_END");
	RoiManager.setGroup(0);
	roiManager("Set Color", "red");
	roiManager("Set Line Width", 0);
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
		showMessage("Launch \"Monitor Memory...\" in a macro"+
			"\nrequires \"MemoryMonitor_Launcher.jar\""+
			"\nto be installed in plugins folder");
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
			"\nseems not to contain Metamorph images to merge");
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
//23456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890