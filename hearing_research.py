# History:
#   2018.10.03  farhad.qureshi@ucalgary.ca    Created
#
# Description:
#   Auditory Neuroscience - Hearing Research
#
#
# Usage:
#   python hearing_research.py


import Tkinter
from tkFileDialog import askopenfilename
import scipy.io as sio
import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy, numpy_to_vtk
import os

"""The software is built for visualizing auditory neural data comprising of
    spikes or action potentials recorded in response to complex acoustic
    stimulation. The visualization maps the data on 3 domains: frequency,
    intensity and time, the 3 fundamentals of acoustic signals."""

class Values(Tkinter.Tk):
    """docstring for Values"""
    def __init__(self, parent):
        Tkinter.Tk.__init__(self,parent)
        self.parent = parent
        self.initialize()
    def initialize(self):

        self.wm_title("Auditory Neuroscience") #Makes the window title
        # Tkinter.config(background = "#FFFFFF")
        self.config(background = "#FFFFFF")
        #Draw frames
        leftFrame = Tkinter.Frame(self, width=200, height = 600)
        leftFrame.grid(row=0, column=0, padx=10, pady=2)
        rightFrame = Tkinter.Frame(self, width=200, height = 600)
        rightFrame.grid(row=0, column=1, padx=10, pady=2)
        Tkinter.Label(leftFrame, text="Options:").grid(row=0, column=0, padx=10, pady=2)
    #
        btnFrame = Tkinter.Frame(leftFrame, width=200, height = 200)
        btnFrame.grid(row=1, column=0, padx=10, pady=2)
        #Putting a heading "Auditory neuroscience"
        lbl = Tkinter.Label(self, text="Hearing Research", font=("Arial Bold", 20))
        lbl.grid(column=0, row=200)

    # Create labels
        self.lbl_filename = Tkinter.Label(rightFrame, text='Filename: ')
        self.lbl_filename.grid(row=1, column=0, padx=10, pady=2)

        self.lbl_thresh = Tkinter.Label(rightFrame, text='Set threshold')
        self.lbl_thresh.grid(row=4, column=0, padx=10, pady=2)
        self.set_thresh = Tkinter.Entry(rightFrame,width=10)
        self.set_thresh.grid(column=1, row=4)

        # lbl_BF = Tkinter.Label(rightFrame, text='BF: ')
        # lbl_BF.grid(row=5, column=0, padx=10, pady=2)
        #
        # lbl_CF = Tkinter.Label(rightFrame, text='CF: ')
        # lbl_CF.grid(row=6, column=0, padx=10, pady=2)
        #
        # lbl_MT = Tkinter.Label(rightFrame, text='MT: ')
        # lbl_MT.grid(row=7, column=0, padx=10, pady=2)

        btn = Tkinter.Button(btnFrame, text="Load Data", command=self.read_data)
        btn.grid(row=0, column=0, padx=10, pady=2)

        btn_1 = Tkinter.Button(btnFrame, text="Volume Rendering",
        command = self.create_volren)
        btn_1.grid(row=0, column=1, padx=10, pady=2)
        #
        btn_2 = Tkinter.Button(btnFrame, text="Surface Rendering",
        command = self.create_surfren)
        btn_2.grid(row=0, column=3, padx=10, pady=2)

        btn_3 = Tkinter.Button(btnFrame, text="Slicing",
        command = self.create_slicing)
        btn_3.grid(column=4, row=0,padx=10, pady=2)

    def create_volren(self):
        """Creates the volume rendering visualization"""
        Matrix_3D = self.matrix
        Matrix_3D = self.set_scale(Matrix_3D,100)
        matrix_vtk,mIn,mAx,scalar_range = self.preprocess_visualization(Matrix_3D)
        matrix_vtk = self.castfunction(matrix_vtk)
        volume = self.volumerendering(matrix_vtk,mIn,mAx,scalar_range)
        self.create_renwin(volume,700,700)
    def create_surfren(self):
        """Creates the surface rendering visualization"""
        Matrix_3D = self.matrix
        Matrix_3D = self.set_scale(Matrix_3D,100)
        matrix_vtk,mIn,mAx,scalar_range = self.preprocess_visualization(Matrix_3D)
        matrix_vtk = self.castfunction(matrix_vtk)
        iso_val = self.set_thresh.get()  # Get the text value entered for threshold
        surface = self.surfacerendering(matrix_vtk,float(iso_val),mIn,mAx,scalar_range)
        self.create_renwin(surface,700,700)
    def create_slicing(self):
        """Creates the slicing visualization"""
        Matrix_3D = self.matrix
        Matrix_3D = self.set_scale(Matrix_3D,100)
        matrix_vtk,mIn,mAx,scalar_range = self.preprocess_visualization(Matrix_3D)
        matrix_vtk = self.castfunction(matrix_vtk)
        slices = self.slicing(matrix_vtk,mIn,mAx,scalar_range)
        self.create_renwin(slices,1600,700)
    def read_data(self):
        """This function loads a dialog box and asks user to select the 3D Matrix"""
        #Read 3D Matrices from Matlab file
        filenamex = askopenfilename() #show an "Open" dialog box and return the
                                    #path to the selected file
        base = os.path.basename(filenamex)#Seperating out the directory from
                                            #the filename
        File = os.path.splitext(base)[0]#Seperating out the
                                            #extension from the filename
    #The name of the matrix is always the number of repetitions used if 2 rep then
    #it would be strf2rep so extracting out the name of the matrix from the filename
        repLoc = File.find("rep") #Finds where in the name is rep
        matrixName = File[repLoc-1:repLoc+3] ##Extracts the repetition number and
                                            #'rep' from the name of the file
        matrixName = 'strf'+matrixName #Adds strf at the start which makes the name
                                    #of the file in matlab
        filename = File+"ftc" #for FTC
        filename1 = File+"iSTRF" # for iSTRF
        mAtrix = File+".mat"
        Matlab_Data = sio.loadmat(mAtrix)
        Matrix_3D=Matlab_Data[matrixName]
        self.lbl_filename.configure(text="Filename: "+mAtrix)
        self.lbl_thresh.configure(text="Set threshold: ",font=("Arial Bold", 10))
        #input for threshold on gui
        #global matrix
        self.matrix = Matrix_3D
        return Matrix_3D
    def set_scale(self,Matrix_3D,i):
        """
        This function scales up the data for increasing resolution of
        visualizations
        Matrix_3D = 3D matrix
        i = number by which to scale up
        """
        return Matrix_3D*i
    def preprocess_visualization(self,Matrix_3D):

        """The function takes in the 3D matrix filters and smooths it and
        returns the filtered matrix and the scalar range for setting up
        the color map"""

        vtk_data = numpy_to_vtk(num_array = Matrix_3D.transpose(2,1,0).ravel(),
     deep = True, array_type = vtk.VTK_FLOAT)
#---------------------------------------------------------
#-------------------------------
        img_vtk = vtk.vtkImageData() # Converting vtk data to vtk image data format
        img_vtk.SetDimensions(Matrix_3D.shape)
        #xd,yd,zd =  Matrix_3D.shape
        #img_vtk.SetDimensions(xd,yd-100,zd)
        img_vtk.GetPointData().SetScalars(vtk_data)
    #-------------------------------
        flipXFilter =vtk.vtkImageFlip()
        flipXFilter.SetFilteredAxis(0)
        flipXFilter.SetInputData(img_vtk)
        flipXFilter.Update()
    #---------------------------------------------------------
        ##Filtering data
        filter_image_m = vtk.vtkImageMedian3D()
        filter_image_m.SetInputData(flipXFilter.GetOutput())
        filter_image_m.SetKernelSize(2,2,2)
        filter_image_m.Update()
        medicalImage_median = filter_image_m.GetOutput()

        filter_image = vtk.vtkImageGaussianSmooth()
        filter_image.SetInputData(medicalImage_median)
        filter_image.SetStandardDeviation(2)
        filter_image.SetRadiusFactors(1,1,1)
        filter_image.SetDimensionality(3)
        filter_image.Update()
        medicalImage = filter_image.GetOutput()
        # Finding min max range
        Range = medicalImage.GetScalarRange()
        mIn = Range[0]
        mAx = Range[1]
        scalar_range = mAx-mIn
        scalar_range = scalar_range/4
        return medicalImage,mIn,mAx,scalar_range
    def castfunction(self,medicalImage):
        """The function converts the filtered image to unsigned short format for
        volume and surface rendering visualizations only"""
    # Converting image into unsigned short for the raycast function
        castFilter = vtk.vtkImageCast()
        castFilter.SetInputData(medicalImage)
        castFilter.SetOutputScalarTypeToUnsignedShort()
        castFilter.Update()
        medicalImageUnsignedShort = castFilter.GetOutput()
        return medicalImageUnsignedShort
    def volumerendering(self,medicalImageUnsignedShort,mIn,mAx,scalar_range):
        """
        Takes in the unsignedshort image and puts it in the volumerendering
        pipeline. min, max and scalar range will be used to define the
            color map

        mIn = float (minimum value in scalar range)
        mAx = float (max value in the scalar range)
        scalar_range = the scalar range of the matrix

        """
        # Using MIPF ray cast mapper
        rayCastFunction = vtk.vtkVolumeRayCastMIPFunction()
        rayCastFunction.SetMaximizeMethodToOpacity()
        #---------------------------------------------------------
        #Mapper
        volumeMapper = vtk.vtkVolumeRayCastMapper()
        volumeMapper.SetInputData(medicalImageUnsignedShort)
        volumeMapper.SetVolumeRayCastFunction(rayCastFunction)
        #---------------------------------------------------------
        # Defining color map
        volumeColor = vtk.vtkColorTransferFunction()
        volumeColor.AddRGBSegment(mIn, 0.0, 0.0, 0.,scalar_range+mIn+0.1,  0.0, 0.0, 1)
        volumeColor.AddRGBSegment(2*scalar_range+mIn,  0, 1, 0.7,3*scalar_range +
         mIn, 1.0, 1.0, 0.0)
        volumeColor.AddRGBPoint(mAx, 1, 0, 0)
        #---------------------------------------------------------
        #Defining opacity
        volumeScalarOpacity = vtk.vtkPiecewiseFunction()
        volumeScalarOpacity.AddPoint(mIn,0.5)
        volumeScalarOpacity.AddPoint(mAx,1)
        volumeGradientOpacity = vtk.vtkPiecewiseFunction()
        volumeGradientOpacity.AddPoint(mIn,0.1)
        volumeGradientOpacity.AddPoint(mAx,1)
        #---------------------------------------------------------
        #Putting everying together
        volumeProperty = vtk.vtkVolumeProperty()
        volumeProperty.SetColor(volumeColor)
        volumeProperty.SetScalarOpacity(volumeScalarOpacity)
        volumeProperty.SetGradientOpacity(volumeGradientOpacity)
        volumeProperty.SetInterpolationTypeToLinear()
        volumeProperty.ShadeOn()
        volumeProperty.SetAmbient(0.3)
        volumeProperty.SetDiffuse(0.1)
        volumeProperty.SetSpecular(0.2)
        #---------------------------------------------------------
        volume = vtk.vtkVolume()
        volume.SetMapper(volumeMapper)
        volume.SetProperty(volumeProperty)
        #---------------------------------------------------------
        return volume
    def surfacerendering(self,medicalImageUnsignedShort,isovalue,mIn,mAx,scalar_range):
        """
        The function takes in the unsignedshort matrix and applies the surface-
        rendering pipeline. Isovalue defines the contour level (0 - 1) 0.5 is
        usually selected for auditory neural data as threshold.min, max and scalar
        range will be used to define the color map

        isovalue = float (0.0 - 1.0)
        mIn = float (minimum value in scalar range)
        mAx = float (max value in the scalar range)
        scalar_range = the scalar range of the matrix

        """
        isovalu = self.set_thresh.get()
        rayCastFunctionIsoSurf = vtk.vtkVolumeRayCastIsosurfaceFunction()
        rayCastFunctionIsoSurf.SetIsoValue(mAx*isovalue)
        volumeMapperIso = vtk.vtkVolumeRayCastMapper()
        volumeMapperIso.SetInputData(medicalImageUnsignedShort)
        volumeMapperIso.SetVolumeRayCastFunction(rayCastFunctionIsoSurf)
        #---------------------------------------------------------
        # Defining color map
        volumeColor = vtk.vtkColorTransferFunction()
        volumeColor.AddRGBSegment(mIn, 0.0, 0.0, 0.,scalar_range+mIn+0.1,  0.0, 0.0, 1)
        volumeColor.AddRGBSegment(2*scalar_range+mIn,  0, 1, 0.7,3*scalar_range +
         mIn, 1.0, 1.0, 0.0)
        volumeColor.AddRGBPoint(mAx, 1, 0, 0)
        #---------------------------------------------------------
        #Putting everying together
        volumeProperty = vtk.vtkVolumeProperty()
        volumeProperty.SetColor(volumeColor)
        volumeProperty.SetInterpolationTypeToLinear()
        volumeProperty.ShadeOn()
        volumeProperty.SetSpecular(0.6)
        #---------------------------------------------------------
        volumeIso = vtk.vtkVolume()
        volumeIso.SetMapper(volumeMapperIso)
        volumeIso.SetProperty(volumeProperty)
        return volumeIso
    def slicing(self,preprocessed_vis,mIn,mAx,scalar_range):
        """
        The function takes in the vtk data converted matrix and performs slicing.
        min, max and scalar range will be used to define the color map
        mIn = float (minimum value in scalar range)
        mAx = float (max value in the scalar range)
        scalar_range = the scalar range of the matrix

        """
        volumeColor = vtk.vtkColorTransferFunction()
        volumeColor.AddRGBSegment(mIn, 0.0, 0.0, 0.,scalar_range+mIn+0.1,  0.0, 0.0, 1)
        volumeColor.AddRGBSegment(2*scalar_range+mIn,  0, 1, 0.7,3*scalar_range +
         mIn, 1.0, 1.0, 0.0)
        volumeColor.AddRGBPoint(mAx, 1, 0, 0)
        SliceMapper = vtk.vtkImageResliceMapper()
        SliceMapper.SetInputData(preprocessed_vis)
        #IF all 4 below are commented slice begins at 0 position
        SliceMapper.SliceAtFocalPointOn() ##Slice begins from the center
        #SliceMapper.SetJumpToNearestSlice() # Select a particular slice (0-Z)
        SliceMapper.SliceFacesCameraOn() ## Change slice position on xyz keys
        SliceMapper.BorderOff()
        #---------------------------------------------------------
        #Image property for source bianry image and registered image
        SliceProperty = vtk.vtkImageProperty()
        SliceProperty.SetColorWindow(500)
        SliceProperty.SetColorLevel(1000)
        SliceProperty.SetLookupTable(volumeColor)
        # SliceProperty.UseLookupTableScalarRangeOn()
        SliceProperty.SetAmbient(0.0)
        SliceProperty.SetDiffuse(1.0)
        SliceProperty.SetOpacity(1.0)
        SliceProperty.SetInterpolationTypeToLinear()
        #---------------------------------------------------------
        #Actors
        SliceActor = vtk.vtkImageSlice()
        SliceActor.SetMapper(SliceMapper)
        SliceActor.SetProperty(SliceProperty)
        return SliceActor
    def draw_axes(self):
        """The function adds 3D axes for visualizations"""
        matrix_vtk = self.preprocess_visualization(self.matrix)[0]

        IntensityLabels=[]
        Intensity=[]
        maxInten=76
        for n in range(0,matrix_vtk.GetDimensions()[2]):

            Intensity.append(maxInten)
            maxInten=maxInten-2

        Intensity=np.sort(Intensity)
        mi=Intensity[0]
        mX=max(Intensity)
        o=[]
        iRange= ((mX-mi)/2) + 1
        for m in range(0,len(Intensity),5):

            o=Intensity[m]
            p=str(o)
            IntensityLabels.append(p)
    #--------------------------------------------------------- #
    #This version of vtk does not support log axis so making a pseudo log axis for
    # Frequency axis

        Freqs=(" "," 40"," "," ","35 "," "," "," ","30 "," "," "," ","25 ",
        " "," "," "," ","20 "," "," "," "," "," "," "," 15"," "," "," "," "," "," ",
        " "," "," ","10 "," "," "," "," "," "," "," "," "," "," "," "," "," "," 5")

        #Freqs=("5","10","15","20","25","30","35","40")
        FreqArray=vtk.vtkStringArray()
        for i in range(0,len(Freqs)):
            FreqArray.InsertNextValue(Freqs[i])


        FrequencyAxisActor =vtk.vtkAxisActor()
        FrequencyAxisActor.SetPoint1(0,0,0)
        FrequencyAxisActor.SetPoint2(matrix_vtk.GetDimensions()[0],0,0)
        FrequencyAxisActor.SetBounds(matrix_vtk.GetBounds())
        FrequencyAxisActor.SetAxisTypeToX()
        FrequencyAxisActor.SetTitle("Frequency (kHz) ")
        FrequencyAxisActor.SetCalculateTitleOffset(0)
        FrequencyAxisActor.SetTitleScale(3)
        FrequencyAxisActor.SetLabelScale(5)
        FrequencyAxisActor.SetAxisType(0)
        FrequencyAxisActor.SetTickVisibility(0)
        FrequencyAxisActor.SetSaveTitlePosition(1)
        FrequencyAxisActor.GetTitleTextProperty().SetColor(0, 0, 0 )
        FrequencyAxisActor.GetLabelTextProperty().SetColor(0, 0, 0 )
        FrequencyAxisActor.GetAxisLinesProperty().SetColor(0,0,0)
        FrequencyAxisActor.SetTickLocationToOutside()
        FrequencyAxisActor.SetLabels(FreqArray)
        FrequencyAxisActor.MinorTicksVisibleOff()
        FrequencyAxisActor.SetMajorStart(0,0)
        FrequencyAxisActor.SetHorizontalOffsetYTitle2D(-10)
        FrequencyAxisActor.SetDeltaMajor(0,1)
        FrequencyAxisActor.SetCamera(self.ren.GetActiveCamera())

        IntensityArray=vtk.vtkStringArray()
        for i in range(0,len(IntensityLabels)):
            IntensityArray.InsertNextValue(IntensityLabels[i])

        IntensityAxisActor =vtk.vtkAxisActor()
        IntensityAxisActor.SetPoint1(0,0,0)
        IntensityAxisActor.SetPoint2(0,0,matrix_vtk.GetDimensions()[2])
        IntensityAxisActor.SetBounds(matrix_vtk.GetBounds())
        IntensityAxisActor.SetTickLocationToBoth()
        IntensityAxisActor.SetAxisTypeToZ()
        IntensityAxisActor.SetLabels(IntensityArray)
        IntensityAxisActor.SetTitle("Intensity (dB SPL)")
        IntensityAxisActor.SetTitleScale(2)
        IntensityAxisActor.SetLabelScale(1.5)
        IntensityAxisActor.GetTitleTextProperty().SetColor(0, 0, 0 )
        IntensityAxisActor.GetLabelTextProperty().SetColor(0, 0, 0 )
        IntensityAxisActor.GetAxisLinesProperty().SetColor(0,0,0)
        IntensityAxisActor.SetTickVisibility(2)
        IntensityAxisActor.SetTitleVisibility(1)
        IntensityAxisActor.SetMajorTickSize(1)
        #IntensityAxisActor.SetHorizontalOffsetYTitle2D(-10)
        #IntensityAxisActor.SetTickLocationToInside()
        IntensityAxisActor.SetDeltaMajor(2,5)
        IntensityAxisActor.MinorTicksVisibleOff()
        IntensityAxisActor.SetCamera(self.ren.GetActiveCamera())

        Time=("0","20 ","40 ","60","80","100","120","140","160","180","200")
        TimeArray=vtk.vtkStringArray()
        for i in range(0,len(Time)):
            TimeArray.InsertNextValue(Time[i])

        TimeAxisActor =vtk.vtkAxisActor()
        TimeAxisActor.SetPoint1(0,0,0)
        TimeAxisActor.SetPoint2(0,matrix_vtk.GetDimensions()[1],0)
        TimeAxisActor.SetBounds(matrix_vtk.GetBounds())
        TimeAxisActor.SetTickLocationToBoth()
        TimeAxisActor.SetAxisTypeToY()
        TimeAxisActor.SetLabels(TimeArray)
        TimeAxisActor.SetTitle("Time (ms)")
        TimeAxisActor.GetAxisLinesProperty().SetColor(0,0,0)
        TimeAxisActor.SetCalculateTitleOffset(0)
        TimeAxisActor.SetTitleScale(2)
        TimeAxisActor.SetLabelScale(3)
        TimeAxisActor.GetTitleTextProperty().SetColor(0., 0., 0. )
        TimeAxisActor.GetLabelTextProperty().SetColor(0., 0., 0.)
        TimeAxisActor.SetTickVisibility(1)
        TimeAxisActor.SetTitleVisibility(1)
        TimeAxisActor.SetMajorTickSize(1)
        TimeAxisActor.SetRange(0,200)

        TimeAxisActor.SetTickLocationToInside()
        TimeAxisActor.SetDeltaMajor(1,20)
        #TimeAxisActor.SetAxisPosition(0)
        TimeAxisActor.MinorTicksVisibleOff()
        TimeAxisActor.SetCamera(self.ren.GetActiveCamera())

        return TimeAxisActor, FrequencyAxisActor, IntensityAxisActor
    def create_renwin(self,data,x,y):
        """The function takes in the volume of data and the window size (x,y)
            as int to display the visualization """

        self.ren = vtk.vtkRenderer()
        renWin = vtk.vtkRenderWindow()
        renWin.AddRenderer(self.ren)
        iren = vtk.vtkRenderWindowInteractor()
        iren.SetRenderWindow(renWin)
        TimeAxisActor, FrequencyAxisActor, IntensityAxisActor = self.draw_axes()
        self.ren.AddActor(TimeAxisActor)
        self.ren.AddActor(FrequencyAxisActor)
        self.ren.AddActor(IntensityAxisActor)
        self.ren.AddVolume(data)
        camera =  self.ren.GetActiveCamera()
        c = data.GetCenter()
        camera.SetFocalPoint(c[0], c[1], c[2])
        camera.SetPosition(c[0] + 400, c[1], c[2])
        camera.SetViewUp(0, 0, -1)
        style = vtk.vtkInteractorStyleImage()
        style.SetInteractionModeToImage3D()
        iren.SetInteractorStyle(style)
        cam1 = self.ren.GetActiveCamera()
        renWin.Render()
        renWin.SetSize(x,y)
        self.ren.SetBackground(0.6,0.6,0.6)
        iren.Start()


#---------------------------------------------------------

if __name__ == '__main__':
    app = Values(None)
    app.mainloop() #this will run until it closes
