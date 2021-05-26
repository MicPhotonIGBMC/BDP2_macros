/*
//from BDP recorder:
run("BDP2 Open Bio-Formats...", "viewingmodality=[Show in new viewer] enablearbitraryplaneslicing=true file=D:\\Extract_Regions_with_BigDataProcessor2\\Celine_20180627_partie\\test-4.nd seriesindex=0");
run("BDP2 Crop...", "inputimage=[test-4] outputimagename=[test-4-crop] viewingmodality=[Show in current viewer] minx=24 miny=90 minz=0 minc=0 mint=0 maxx=325 maxy=507 maxz=10 maxc=1 maxt=7 ");
run("BDP2 Save As...", "inputimage=[test-4-crop] directory=[C:\\Users_Data\\Celine\\20180627_partie_XMLHDF5_crop\\] numiothreads=1 numprocessingthreads=4 filetype=[BigDataViewerXMLHDF5] saveprojections=false savevolumes=true tiffcompression=[None] tstart=0 tend=7 ");
//record other viewing modalities;
//loop over seriesindex from 0 to npositions-1
*/

//run("BDP2 Open Bio-Formats...", "viewingmodality=[Show in new viewer] enablearbitraryplaneslicing=true file=D:\\Extract_Regions_with_BigDataProcessor2\\Celine_20180627_partie\\test-4.nd seriesindex=0");

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
var version = 03;
var copyRight = "Author: Marcel Boeglin - May 2021";
var email = "e-mail: boeglin@igbmc.fr";

var doCrop = true;
var cropRoisImagesDir;
var dir1;//input folder
var dir2;//output folder
var dir3;//folder containing images with crop-rois in their overlay

var doRangeAroundMedianSlice = false;
var rangeAroundMedianSlice = 50;// % of stack-size
/*
var firstSlice = 0;//in BDP, slice range is [0, slices-1]
var lastSlice = slices-1;
*/
var width, height, channels, slices, frames;

//to keep memory of previous slider positions
var previousMinC, previousMinZ, previousMinT;
var previousMaxC, previousMaxZ, previousMaxT;

var minx=0, miny=0, minz=0, minc=0, mint=0;
var maxx, maxy, maxz, maxc, maxt;

//inutile
var timepoints;
var firstTimePoint = 0;//in BDP, time range is [0, timepoints-1]
var lastTimePoint = timepoints-1;

//Non utilise, a revoir : 
//var doRangeFrom_t0 = false;
//var rangeFrom_t0 = 50;// % of timepoints
//remplacement par:

var medianTimepoint = timepoints/2;
var rangeAroundMedianTimePoint = 50;// % of timepoints
var doRangeAroundMedianTimePoint = false;// if false do all time points

var resize = false;
var resizeFactor = 0.5;

var inputPath, inputDir, inputFile, inputFileName;
var inputImage, inputImageName;
var inputImageTitle;
var outputDir;

//Macro BEGIN

run("Bio-Formats Macro Extensions");
execute();

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
/*
	dir3 = cropRoisImagesDir;
	if (!File.exists(dir3)) dir3 = "";
*/
	print(macroname+version+".ijm");
	print(macroDescription());
	print(copyRight);
	print(email);
/*
	print("dir1 = "+dir1);
	print("dir2 = "+dir2);
	print("dir3 = "+dir3);
	print("cropRoisImagesDir = "+cropRoisImagesDir);
*/
	getInputPath();
	print("\ninputDir = "+inputDir);
	print("inputFile = "+inputFile);
	print("inputFileName = "+inputFileName);
	
	outputDir = getDirectory("Destination Directory ");
	print("outputDir = "+outputDir);

	open(inputPath);
	if (nImages<1) {showMessage("No image: Macro aborted"); return;}
	inputImageTitle = getTitle();
	print("inputImageTitle : "+inputImageTitle);
	inputImageTitle = replace(inputImageTitle, "\"", "");
	imageID = getImageID();
	Stack.setDisplayMode("composite");
	Stack.getDimensions(width, height, channels, slices, frames);
	s = slices/2;
	if (s<1) s=1;
	Stack.setSlice(s);
	t = frames/2;
	if (t<1) t=1;
	Stack.setFrame(t);
	for (c=1; c<=channels; c++) {
		Stack.setChannel(c);
		run("Enhance Contrast", "saturated=0.35");
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
		addRegionToManager();
	}
	//close();
	roiManager("deselect");
	roiManager("save", outputDir+inputImageTitle+"_Crop-Rois.zip");
	
	//Process Crop-Rois from Roi Manager
	nregions = roiManager("count")/2;//each region is defined by 2 Rois
	//utiliser le nom du fichier ouvert avec Bioformat
	//Probleme: comment recuperer l'indice de la position ouverte en mode virtuel

/*
	//Open image to be cropped
	//setBatchMode(false);
	run("BDP2 Open Bio-Formats...",
		"viewingmodality=[Show in new viewer] enablearbitraryplaneslicing=true "+
		"file=D:\\Extract_Regions_with_BigDataProcessor2\\Celine_20180627_partie\\test-4.nd"+
		" seriesindex=0");
*/
/*
	run("BDP2 Open Bio-Formats...",
		"viewingmodality=[Show in new viewer] enablearbitraryplaneslicing=true "+
		"["+inputPath+"] seriesindex=0");
*/
	//Close virtual input image and Replace with a real image having 1 slice,
	//width=1, height=1 to get Roi coordinates
	close();//close image used for Crop-Rois drawing
	//newImage("Tmp", "8-bit", width, height, 1);
	newImage("Tmp", "8-bit", 1, 1, 1);//works also
	tmpid = getImageID();
	for (r=0; r<nregions; r++) {
		print("");
		//i = 2*r;
		roiManager("select", 2*r);
		Roi.getBounds(minx, miny, w, h);
		//Roi.getPosition(channel, slice, frame);
		Roi.getPosition(mincPlus1, minzPlus1, mintPlus1);
		minc = mincPlus1-1;
		minz = minzPlus1-1;
		mint = mintPlus1-1;
		print("minx="+minx+"  miny="+miny+"  minz="+minz+"  minc="+minc+"  mint="+mint);

		roiManager("select", 2*r+1);
		Roi.getPosition(maxcPlus1, maxzPlus1, maxtPlus1);
		maxc = maxcPlus1-1;
		maxz = maxzPlus1-1;
		maxt = maxtPlus1-1;
		maxx = minx+w;
		maxy = miny+h;
		print("maxx="+maxx+"  maxy="+maxy+"  maxz="+maxz+"  maxc="+maxc+"  maxt="+maxt);
	}
	selectImage(tmpid);
	close();

print("inputPath");
print(inputPath);

//OK, mais comment ouvrir une serie autre que 0
//on peut le deduire du titre de l'image et du fichier .nd
run("BDP2 Open Bio-Formats...", "viewingmodality=[Show in new viewer] enablearbitraryplaneslicing=true file="
	+inputPath+" seriesindex=0");//OK, mais comment ouvrir une serie autre que 0
}

//run("Bio-Formats Importer", "open=&path color_mode=Default view=Hyperstack stack_order=XYCZT series_"+j+1);
//run("Bio-Formats", "[open=D:/Extract_Regions_with_BigDataProcessor2/Celine_20180627_partie/test-4.nd] color_mode=Default view=Hyperstack stack_order=XYCZT series_1");
//run("Bio-Formats", "open=D:/Extract_Regions_with_BigDataProcessor2/Celine_20180627_partie/test-4.nd color_mode=Default view=Hyperstack stack_order=XYCZT use_virtual_stack series_0");//Problem
//run("Bio-Formats", "open=D:/Extract_Regions_with_BigDataProcessor2/Celine_20180627_partie/test-4.nd color_mode=Custom rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT use_virtual_stack series_1 series_0_channel_0_red=0 series_0_channel_0_green=255 series_0_channel_0_blue=0 series_0_channel_1_red=255 series_0_channel_1_green=0 series_0_channel_1_blue=0");
//probleme SCI java
//var ; initializeSciJavaParameters ( ) ; run ( "Bio-Formats" , "open=D:/Extract_Regions_with_BigDataProcessor2/Celine_20...
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

function getDirs() {
	dir1 = getDirectory("Source Directory ");
	if (!File.exists(dir1)) return false;
	if (doCrop) {
		cropRoisImagesDir = getDirectory("Directory containing Crop-Rois");
		if (!File.exists(cropRoisImagesDir)) return false;
	}
	dir2 = getDirectory("Destination Directory ");
	while (dir2==dir1) {
		showMessage("Destination must be different from source");
		dir2 = getDirectory("Choose Destination Directory ");
	}
	if (!File.exists(dir2)) return false;
	return true;
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

function prefixLess(roiImages) {
	n = roiImages.length;
	list2 = newArray(n);
	for (i=0; i<n; i++) {
		s = roiImages[i];
		for (j=0; j<projPrefixes.length; j++) {
			prefix  = projPrefixes[j];
			if (startsWith(s, prefix)) {
				list2[i] = substring(s, lengthOf(prefix), lengthOf(s));
				break;
			}
		}
	}
	return list2;
}

/** Adds ROIs from roiImage corresponding to current series and position
	to ROI Manager.
	roiImages: list of filnames in folder dir1+"ImagesWithRois"
	seriesName: name of currently processed series
	channels: array of current series channel-strings. channel-string may be ""
	position: name of currently processed position: "_s1", "_s2", etc
	In case of several images with rois matching a given series or position
	the used one is first in the 'prefixLessRoiImages' list */
function getRoisFromImage(prefixLessRoiImages, seriesName, channels, position) {
	print("getRoisFromImage(roiImages, "+seriesName+
			", channels, position = "+position+")");
	print("seriesName: "+seriesName);
	print("position = "+position);
	roiManager("reset");
	for (i=0; i<prefixLessRoiImages.length; i++) {
		prefixLessRoiImage = prefixLessRoiImages[i];
		if (!startsWith(prefixLessRoiImage, seriesName)) continue;
		if (indexOf(prefixLessRoiImage, position)<0) continue;
		pos = getPositionString(prefixLessRoiImage);
		if (pos!=position) continue;
		seriesLessName = substring(prefixLessRoiImage, lengthOf(seriesName), 
				lengthOf(prefixLessRoiImage));
		for (j=0; j<channels.length; j++) {
			print("channels["+j+"] = "+channels[j]);
			channel = channels[j];
			if (channel!="" && !startsWith(seriesLessName, channel))
				continue;
			channelLessName = substring(seriesLessName, lengthOf(channel),
					lengthOf(seriesLessName));
			if (startsWith(channelLessName, position)) {
				print("getting rois from\n"+roiImages[i]+"\n ");
				open(dir1+"ImagesWithRois"+File.separator+roiImages[i]);
				run("To ROI Manager");
				close();
				return true;
			}
		}
	}
	return false;
}

function getPositionString(fileName) {
	if (matches(fileName, ".*_s\\d{1,3}.*")) {
		str = substring(fileName, lastIndexOf(fileName, "_s"));
		//print("str = "+str);
		if (indexOf(str, ".")>=0)
			str = substring(str, 0, lastIndexOf(str, "."));
		//print("str = "+str);
		if (matches(str, ".*_t\\d{1,5}.*"))
			return substring(str, 0, lastIndexOf(str, "_t"));
	}
	return "";
}

function reduceSeriesToRoiImages() {
	print("\n\nreduceSeriesToRoiImages()");
	nSeries = seriesNames.length;
	seriesNames2 = newArray(nSeries);
	n2 = 0;
	for (j=0; j<nSeries; j++) {
		seriesName = seriesNames[j];
		for (i=0; i<roiImages.length; i++) {
			prefixLessName = prefixLessRoiImages[i];
			print("prefixLessName = "+prefixLessName);
			if (!startsWith(prefixLessName, seriesName)) continue;
			if (belongsToSeries(prefixLessName, seriesName)) {
				seriesNames2[n2++] = seriesName;
				break;
			}
		}
	}
	//print("nSeries2 = "+n2);
	return Array.trim(seriesNames2, n2);
}


function belongsToSeries(prefixLessRoiImageName, seriesName) {
	seriesLessName = substring(prefixLessRoiImageName, lengthOf(seriesName), 
			lengthOf(prefixLessRoiImageName));
	print("seriesLessName = "+seriesLessName);
	if (startsWith(seriesLessName, "_w")) return true;
	if (startsWith(seriesLessName, "_s")) return true;
	if (startsWith(seriesLessName, "_t")) return true;
	if (startsWith(toLowerCase(seriesLessName), ".tif")) return true;
	if (startsWith(toLowerCase(seriesLessName), ".stk")) return true;
	return false;
}

/* Returns array of TIF and STK imagenames in imagesWithRoisDir */
function roiImageList(imagesWithRoisDir) {
	l = getFileList(imagesWithRoisDir);
	nf = l.length;
	l2 = newArray(nf);
	n = 0;
	for (i=0; i<nf; i++) {
		f = l[i];
		if (File.isDirectory(f)) continue;
		if (!endsWith(toUpperCase(f), ".TIF") &&
			!endsWith(toUpperCase(f), ".STK")) continue;
		open(imagesWithRoisDir+File.separator+f);
		nr = Overlay.size;
		close();
		if (nr<1) continue;
		l2[n++] = f;
	}
	return Array.trim(l2, n);
}

//80 chars:
//23456789 123456789 123456789 123456789 123456789 123456789 123456789 1234567890