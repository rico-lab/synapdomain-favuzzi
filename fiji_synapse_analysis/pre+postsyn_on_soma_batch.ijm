/*
 * Macro template to process multiple images in a folder
 */

input = getDirectory("Input directory");
output = getDirectory("Output directory");

Dialog.create("File type");
Dialog.addString("File suffix: ", ".tif", 5);
Dialog.show();
suffix = Dialog.getString();

processFolder(input);
setBatchMode(true);
function processFolder(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + list[i]))
			processFolder("" + input + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	};
};
function processFile(input, output, file) {
	open(input + file);
//Start script
saveAs("Tiff", "/Users/Emilia/Desktop/output/picture-1.tif");
// NeuN/soma
run("Duplicate...", "title=picture1-3.tif duplicate channels=3");
selectWindow("picture1-3.tif");
// define PV threshold
title = "PV Threshold";
  width=1000; height=1000;
  Dialog.create("Please select PV threshold");
  Dialog.addNumber("PV Min Threshold:", 25);
    Dialog.show();
  BC_TR = Dialog.getNumber();
run("Subtract Background...", "rolling=50");
run("Despeckle");
run("Gaussian Blur...", "sigma=0.5");
run("Smooth");
run("Enhance Contrast...", "saturated=0.5");
run("RGB Color");
run("Color Threshold...");
// Color Thresholder 1.47v
// Autogenerated macro, single images only!
min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
min[0]=0;
max[0]=255;
filter[0]="pass";
min[1]=0;
max[1]=255;
filter[1]="pass";
min[2]=25;
max[2]=255;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop")  run("Invert");
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++){
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename(a);
// Colour Thresholding-------------

//setTool("polygon") find ROI;
setOption("BlackBackground", true);
run("Make Binary");
run("Analyze Particles...", "size=100-Infinity show=Masks display include summarize");
run("Create Selection");
run("Enlarge...", "enlarge=-2.5");
run("Enlarge...", "enlarge=2.5");
waitForUser("CHECK ROI");
roiManager("Add");
roiManager("Measure");
run("Fill");
run("Clear Outside");
saveAs("Tiff", "/Users/emilia/Desktop/output/mask1.tif");
close("\\Others");
// NeuN end
// Syt2
open("/Users/emilia/Desktop/output/picture-1.tif");
selectWindow("picture-1.tif");
run("Duplicate...", "title=picture1-1.tif duplicate channels=1");
selectWindow("picture1-1.tif");
run("Subtract Background...", "rolling=50");
run("Despeckle");
run("Gaussian Blur...", "sigma=0.5");
run("Smooth");
run("Enhance Contrast...", "saturated=0.5");
run("RGB Color");
run("Color Threshold...");
// Color Thresholder 1.47v
// Autogenerated macro, single images only!
min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
min[0]=0;
max[0]=255;
filter[0]="pass";
min[1]=0;
max[1]=255;
filter[1]="pass";
min[2]=110;
max[2]=255;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop")  run("Invert");
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++){
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename(a);
// Colour Thresholding-------------
run("Analyze Particles...", "size=0.10-infinity circularity=0.00-1.00 show=Masks");
run("Invert");
run("Watershed");
run("Invert");
saveAs("Tiff", "/Users/emilia/Desktop/output/mask2.tif");
close("\\Others");
// Syt2 end
// Gephyrin
open("/Users/emilia/Desktop/output/picture-1.tif");
selectWindow("picture-1.tif");
run("Duplicate...", "title=picture1-2.tif duplicate channels=2");
selectWindow("picture1-2.tif");
run("Subtract Background...", "rolling=50");
run("Despeckle");
run("Gaussian Blur...", "sigma=0.5");
run("Smooth");
run("Enhance Contrast...", "saturated=0.5");
run("RGB Color");
run("Color Threshold...");
// Color Thresholder 1.47v
// Autogenerated macro, single images only!
min=newArray(3);
max=newArray(3);
filter=newArray(3);
a=getTitle();
run("HSB Stack");
run("Convert Stack to Images");
selectWindow("Hue");
rename("0");
selectWindow("Saturation");
rename("1");
selectWindow("Brightness");
rename("2");
min[0]=0;
max[0]=255;
filter[0]="pass";
min[1]=0;
max[1]=255;
filter[1]="pass";
min[2]=110;
max[2]=255;
filter[2]="pass";
for (i=0;i<3;i++){
  selectWindow(""+i);
  setThreshold(min[i], max[i]);
  run("Convert to Mask");
  if (filter[i]=="stop")  run("Invert");
}
imageCalculator("AND create", "0","1");
imageCalculator("AND create", "Result of 0","2");
for (i=0;i<3;i++){
  selectWindow(""+i);
  close();
}
selectWindow("Result of 0");
close();
selectWindow("Result of Result of 0");
rename(a);
// Colour Thresholding-------------
run("Analyze Particles...", "size=0.10-infinity circularity=0.00-1.00 show=Masks");
run("Invert");
run("Watershed");
run("Invert");
saveAs("Tiff", "/Users/emilia/Desktop/output/mask3.tif");
close("\\Others");
// Geph end
open("/Users/emilia/Desktop/output/mask2.tif");
open("/Users/emilia/Desktop/output/mask3.tif");
run("Merge Channels...", "c1=mask2.tif c2=mask3.tif create");
// select colocalizing pre+post
run("Invert");
run("Stack to RGB");
run("8-bit");
//run("Threshold...");
setThreshold(10, 255);
run("Convert to Mask");
// Add band around cell & counts synapses
selectWindow("Composite (RGB)");
roiManager("Select", 0);
run("Enlarge...", "enlarge=-0.5");
run("Make Band...", "band=0.7");
run("Analyze Particles...", "size=0.03-Infinity show=Masks");
// analysis of particles - need to invert image for analysis
run("Analyze Particles...", "size=0.02-Infinity display summarize");
selectWindow("Composite (RGB)");
saveAs("Tiff", output + file);
// clear screen
close("\\Others"); 
run("Close All");
// counts Syt2 puncta
open("/Users/emilia/Desktop/output/mask2.tif");
run("8-bit");
run("Invert");
//run("Threshold...");
setThreshold(10, 255);
run("Convert to Mask");
// Add band around cell
roiManager("Select", 0);
run("Enlarge...", "enlarge=-0.1");
run("Make Band...", "band=0.7");
run("Analyze Particles...", "size=0.03-Infinity show=Masks");
// analysis of particles - need to invert image for analysis
run("Analyze Particles...", "size=0.02-Infinity display summarize");
run("Close All");
// counts Geph clusters
open("/Users/emilia/Desktop/output/mask3.tif");
run("8-bit");
run("Invert");
//run("Threshold...");
setThreshold(10, 255);
run("Convert to Mask");
// Add band around cell
roiManager("Select", 0);
run("Enlarge...", "enlarge=-0.6");
run("Make Band...", "band=0.7");
run("Analyze Particles...", "size=0.03-Infinity show=Masks");
// analysis of particles - need to invert image for analysis
run("Analyze Particles...", "size=0.02-Infinity display summarize");
waitForUser("click 'OK'");
roiManager("Delete");
	print("Processing: " + input + file);
	print("Saving to: " + output);
}
setBatchMode(false);
print("THIS THIS THE END");