// This script allows for loading cluster labels for the entire project.
// It requires some configuration for file paths/names and column indices (see comments)

import qupath.lib.regions.ImagePlane
import qupath.lib.regions.ImagePlane
import qupath.lib.roi.EllipseROI;
import qupath.lib.objects.PathDetectionObject
import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathCellObject
import qupath.lib.objects.PathObjects

import qupath.lib.analysis.DistanceTools
import static qupath.lib.gui.scripting.QPEx.*

def plane = ImagePlane.getPlane(0, 0) 
def myProject = getProject()

// Modify with the path to your cluster directory
predDir = /path\to\cluster\directory/

myProject.getImageList().each{
    String imgName = it.getImageName().split("\\.")[0]
    def imageData = it.readImageData()
    def hierarchy = imageData.getHierarchy()
    def server = imageData.getServer()
    
    def cal = server.getPixelCalibration()
    def calibx = cal.pixelWidth
    def caliby = cal.pixelHeight
    
    println(imgName+" calX="+calibx+" CalY="+caliby)
    
    // Modify according to file naming pattern for clusters
    def prefFH = new File(predDir + "/seurat_clusters_" + imgName +".csv" )
    if( prefFH.exists() ){
            println("Have Predictions for : "+imgName)
    } else {
        println("  >> SKIP: "+imgName)
        return
    }
    
    trainingObjects = []
    predLines = prefFH.readLines()
    predLines.remove(0)
    predLines.each { String line ->
        lineFromFile = line.split(",")
        
        // Modify with indices of X, Y, and cluster label columns
        annoX = Float.parseFloat(lineFromFile[0]) * (1/calibx)
        annoY = Float.parseFloat(lineFromFile[1]) * (1/caliby)
        phenoClass = lineFromFile[2]
        
        size = 3.0
        def roi = new EllipseROI(annoX-size/2,annoY-size/2,size,size, ImagePlane.getDefaultPlane())
        trainingObjects << PathObjects.createDetectionObject(roi, PathClass.fromString(phenoClass))
    }
    hierarchy.addPathObjects(trainingObjects)
    trainingObjects.each{
        it.getMeasurementList().putMeasurement("GroundTruth",1)
    }

    fireHierarchyUpdate()
    it.saveImageData(imageData)
    
    
    ///Loop all cells
    hierarchy.getDetectionObjects().eachWithIndex{ cell, idx ->
        roi = cell.getROI()
        if(cell.getMeasurementList().containsNamedMeasurement("GroundTruth")){return}
        ptCell = hierarchy.getObjectsForROI(qupath.lib.objects.PathDetectionObject,roi).find{it.getMeasurementList().containsNamedMeasurement("GroundTruth")}
         if(ptCell){
            cell.setPathClass(ptCell.getPathClass())
        } else {
           cell.setPathClass(null)
        }
    }

    // Wipeout
    hierarchy.removeObjects(hierarchy.getDetectionObjects().findAll{it.getMeasurementList().containsNamedMeasurement("GroundTruth")}, true)
    fireHierarchyUpdate()
    
    it.saveImageData(imageData)
    imageData.getServer().close() 
    println("")
}

// Changes should now be reflected in the project directory
myProject.syncChanges()





