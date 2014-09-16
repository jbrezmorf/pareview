

'''

'''


Name = 'EclipseReader'
Label = 'Eclipse EGRID and UNRST reader'
Help = 'Read mesh in EGRID format and time dependent data from UNRST output files.'

NumberOfInputs = 0
OutputDataType = 'vtkPolyData'

Properties = dict(
  FileNameBase = ''
  )


'''
Search for given 'keyword' in the given binary file 'f'.
Seeks 'f' to the first byte after the keyword. 
Returns -1 on fail.
'''
def skip_after_keyword(f, keyword):
    data_str=""
    while True:
        new_data=f.read(1024*4)
        if (not new_data):
            return -1
        data_str += new_data 
        pos=data_str.find(keyword)
        if (pos >= 0):
            f.seek(-(len(data_str)-pos),os.SEEK_CUR)
            return 1
        else:
            data_str=data.str[len(data_str) - len(keyword):]
                

'''
Return dictionary with data from the header. Possible keys:

'''
def read_egrid_header(f):
    return


'''
Reads EGRID file with name given by parameter 'filename' and
store it into (empty) output object given by 'pdo' parameter.
Returns updated pdo object.
'''
def read_egrid_mesh(pdo, filename):
    f=file.open(filename, 'rb')
    header_data=read_egrid_header(f)
    



def RequestData():

    #This script generates a helix double.
    #This is intended as the script of a 'Programmable Source'
    import math

    #Get a vtk.PolyData object for the output
    pdo = self.GetPolyDataOutput()
    read_egrid_mesh(pdo,FileNameBase+'.EGRID')
    
    numPts = Number_of_points # Points along each Helix
    length = float(Length) # Length of each Helix
    rounds = float(Rounds) # Number of times around
    phase_shift = math.pi/Phase_shift # Phase shift between Helixes


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









if __name__ == "__main__":
    import sys
    fib(int(sys.argv[1]))