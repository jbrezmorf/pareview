import vtk
import numpy

def to_str(string):
      string=string[0:8]  
      i_str=numpy.fromstring(string, dtype="uint64")
      return i_str[0]


def make_rect(coords):
    xmin,ymin,xmax,ymax=coords

    pts = vtk.vtkPoints()
    pts.InsertPoint(0, xmin, ymin, 0)
    pts.InsertPoint(1, xmax, ymin, 0)
    pts.InsertPoint(2, xmax, ymax, 0)
    pts.InsertPoint(3, xmin, ymax, 0)
    rect = vtk.vtkCellArray()
    rect.InsertNextCell(5)
    rect.InsertCellPoint(0)
    rect.InsertCellPoint(1)
    rect.InsertCellPoint(2)
    rect.InsertCellPoint(3)
    rect.InsertCellPoint(0)
    
    output = vtk.vtkPolyData()
    output.SetPoints(pts)
    output.SetLines(rect)

    labels=vtk.vtkStringArray()
    labels.InsertNextValue("one")
    labels.InsertNextValue("two")
    labels.InsertNextValue("three")
    labels.InsertNextValue("four")
    labels.SetName("labels")
    output.GetPointData().AddArray(labels)
    
    tenths=vtk.vtkTypeInt64Array()
    tenths.InsertNextValue(to_str("10abcdefgh"))
    tenths.InsertNextValue(to_str("20abcdefgh"))
    tenths.InsertNextValue(to_str("30abcdefgh"))
    tenths.InsertNextValue(to_str("40abcdefgh"))
    tenths.SetName("tenths")
    output.GetPointData().AddArray(tenths)
    
    return output
    
output.SetBlock(0, make_rect( (0,0,1,1) ) )    
output.SetBlock(1, make_rect( (1,1,2,2) ) )

'''
block=output.GetBlock(0)
pd=block.GetPointData()
index=vtk.mutable(-1)
array=pd.GetAbstractArray("labels", index)
print "index: ", index
print array
string_array=vtk.vtkStringArray.SafeDownCast(array)
print string_array
'''