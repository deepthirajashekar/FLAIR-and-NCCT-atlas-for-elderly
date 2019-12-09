import vtk 

# prepare a output image:
mask = vtk.vtkImageData()
bounds = [0]*6  
mask.SetSpacing(1, 1, 1)
mask.SetDimensions(193, 229, 193)
mask.SetExtent(0,192,0,228,0,192) #[96 132 -78], [289 361 115]
mask.SetOrigin(0,0,0)
mask.AllocateScalars(vtk.VTK_UNSIGNED_CHAR, 1)

# fill the imahe with foreground voxels:
inval = 255 
outval = 0
count = mask.GetNumberOfPoints()
for i in range(count):
	mask.GetPointData().GetScalars().SetTuple1(i,inval)


# read mesh file 
reader = vtk.vtkSTLReader()
reader.SetFileName('origin-vtk-friendly-smoothed.stl')
reader.Update()

# sweep polydata: 
extruder = vtk.vtkLinearExtrusionFilter() 
extruder.SetInputData(reader.GetOutput())
extruder.SetScaleFactor(1.)
extruder.SetExtrusionTypeToNormalExtrusion()
extruder.SetVector(0, 0, 1)
extruder.Update()


# poly data to image stencil: 
pol2stenc = vtk.vtkPolyDataToImageStencil()
pol2stenc.SetTolerance(0)
pol2stenc.SetInputConnection(reader.GetOutputPort())
pol2stenc.SetOutputOrigin(0,0,0)
pol2stenc.SetOutputSpacing(1, 1, 1)
pol2stenc.SetOutputWholeExtent(mask.GetExtent())
pol2stenc.Update()

# cut the corresponding white image and set the background:
imgstenc = vtk.vtkImageStencil()
imgstenc.SetInputData(mask)
imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
imgstenc.ReverseStencilOff()
imgstenc.SetBackgroundValue(outval)
imgstenc.Update()

# write out the image
imageWriter = vtk.vtkNIFTIImageWriter()
imageWriter.SetFileName("origin-vtk-friendly-smoothed.nii.gz")
imageWriter.SetInputConnection(imgstenc.GetOutputPort())
imageWriter.Write()
