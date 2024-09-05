import qupath.lib.gui.dialogs.Dialogs
//import javax.swing.*
import qupath.lib.regions.ImagePlane
import qupath.lib.objects.PathDetectionObject
import qupath.lib.objects.PathCellObject
import qupath.lib.objects.classes.PathClass
import qupath.lib.objects.classes.PathClassFactory

def plane = ImagePlane.getPlane(0, 0)
def myProject = getProject()
def imageData = getCurrentImageData()
hierarchy = imageData.getHierarchy()
pixelSize = imageData.getServer().getPixelCalibration().getAveragedPixelSizeMicrons()

// Specify the delimiter used in your file (e.g., ',' for CSV, '\t' for tab-separated)
def columnX = 'x'
def columnY = 'y'
int t = 0
int z = 0
def delimiter = ','
File baseDir = new File(buildFilePath(PROJECT_BASE_DIR))

println("Please choose the file that contains multiple classifications in columns.")
File mxClFile = Dialogs.promptForFile("Choose Multiple Clusters File", baseDir, "csv")
print(mxClFile)

// Function to read and list the first row of a file
def listFirstRow(file, delimiter) {
    try {
        def firstRow = file.text.split('\n')[0]
        //println("First Row: $firstRow")
        return firstRow.split(delimiter)
    } catch (Exception e) {
        println("Error reading the file: ${e.message}")
        return null
    }
}
// Example usage
def header = listFirstRow(mxClFile, delimiter)
//println(header)
def clColumn = Dialogs.showChoiceDialog("Choose Cluster", "Choose Cluster", header, header[0])
print("Selected: "+clColumn)


resetDetectionClassifications()


originalList = []
cIndx = header.findIndexOf { it == clColumn }
xIndx = header.findIndexOf { it == columnX }
yIndx = header.findIndexOf { it == columnY }

mxClFile.splitEachLine(delimiter) {fields ->
    if(fields[xIndx] != columnX) {  // need to skip header
        cX = Float.valueOf(fields[xIndx].trim())
        cY = Double.parseDouble(fields[yIndx])
        
        phenoClass = fields[cIndx]
        cellHere = PathObjectTools.getObjectsForLocation(hierarchy, cX/pixelSize, cY/pixelSize,  t, z, -1).find{it.isDetection()}
        if( cellHere != null) {
            //print(phenoClass)
            cellHere.setPathClass(getPathClass(phenoClass))
            originalList << phenoClass
        } else {
            print(" >> NULL LOCATION  x="+cX+" y="+cY+"   toName="+phenoClass)
            
        }
        
    }
}
   
fireHierarchyUpdate()

def newUniqueList = originalList.unique()
listOfClasses = [getPathClass(null)]
newUniqueList.each {
   listOfClasses <<  getPathClass(it)
    }
// Delete the default classes
def pathClasses = getQuPath().getAvailablePathClasses()
pathClasses.setAll(listOfClasses)


project.syncChanges()