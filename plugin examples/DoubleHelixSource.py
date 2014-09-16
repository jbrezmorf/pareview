# An example Programmable Filter from http://www.paraview.org/Wiki/Python_Programmable_Filter#Programmable_Source:_Double_Helix

Name = 'DoubleHelixSource'
Label = 'Double Helix Source'
Help = 'An example that draws a double helix with connecting lines (like DNA).'

NumberOfInputs = 0
OutputDataType = 'vtkPolyData'

Properties = dict(
  Number_of_points = 80,
  Length = 8.0,
  Rounds = 3.0,
  Phase_shift = 1.5
  )


def RequestData():

    #This script generates a helix double.
    #This is intended as the script of a 'Programmable Source'
    import math

    numPts = Number_of_points # Points along each Helix
    length = float(Length) # Length of each Helix
    rounds = float(Rounds) # Number of times around
    phase_shift = math.pi/Phase_shift # Phase shift between Helixes

    #Get a vtk.PolyData object for the output
    pdo = self.GetPolyDataOutput()

    #This will store the points for the Helix
    newPts = vtk.vtkPoints()
    for i in range(0, numPts):
       #Generate Points for first Helix
       x = i*length/numPts
       y = math.sin(i*rounds*2*math.pi/numPts)
       z = math.cos(i*rounds*2*math.pi/numPts)
       newPts.InsertPoint(i, x,y,z)

       #Generate Points for second Helix. Add a phase offset to y and z.
       y = math.sin(i*rounds*2*math.pi/numPts+phase_shift)
       z = math.cos(i*rounds*2*math.pi/numPts+phase_shift)
       #Offset Helix 2 pts by 'numPts' to keep separate from Helix 1 Pts
       newPts.InsertPoint(i+numPts, x,y,z)

    #Add the points to the vtkPolyData object
    pdo.SetPoints(newPts)

    #Make two vtkPolyLine objects to hold curve construction data
    aPolyLine1 = vtk.vtkPolyLine()
    aPolyLine2 = vtk.vtkPolyLine()

    #Indicate the number of points along the line
    aPolyLine1.GetPointIds().SetNumberOfIds(numPts)
    aPolyLine2.GetPointIds().SetNumberOfIds(numPts)
    for i in range(0,numPts):
       #First Helix - use the first set of points
       aPolyLine1.GetPointIds().SetId(i, i)
       #Second Helix - use the second set of points
       #(Offset the point reference by 'numPts').
       aPolyLine2.GetPointIds().SetId(i,i+numPts)

    #Allocate the number of 'cells' that will be added.
    #Two 'cells' for the Helix curves, and one 'cell'
    #for every 3rd point along the Helixes.
    links = range(0,numPts,3)
    pdo.Allocate(2+len(links), 1)

    #Add the poly line 'cell' to the vtkPolyData object.
    pdo.InsertNextCell(aPolyLine1.GetCellType(), aPolyLine1.GetPointIds())
    pdo.InsertNextCell(aPolyLine2.GetCellType(), aPolyLine2.GetPointIds())

    for i in links:
       #Add a line connecting the two Helixes.
       aLine = vtk.vtkLine()
       aLine.GetPointIds().SetId(0, i)
       aLine.GetPointIds().SetId(1, i+numPts)
       pdo.InsertNextCell(aLine.GetCellType(), aLine.GetPointIds())
