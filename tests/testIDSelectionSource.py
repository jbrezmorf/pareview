#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()
import vtk


# create a new 'Table To Points'
'''
datacsv = CSVReader(FileName=['/home/jb/workspace/pareview/test/data.csv'])
source = TableToPoints(Input=datacsv)
source.XColumn = 'x'
source.YColumn = 'y'
source.ZColumn = 'z'
'''

with open("./labeled_poly_source.py", 'r') as content_file:
    script = content_file.read()

source=ProgrammableSource()
source.OutputDataSetType=0 # PolyData
source.Script=script

rep = Show()




n_points=source.GetDataInformation().DataInformation.GetNumberOfPoints()
selection=IDSelectionSource()
IDs = []
for i in range(n_points):
    IDs.append(0L)
    IDs.append(long(i))
selection.IDs = IDs

source.SetSelectionInput(0,selection,0)

# set active source
SetActiveSource(source)

rep.SelectionPointLabelVisibility = 1
rep.SelectionPointFieldDataArrayName = 'label'
rep.SelectionPointLabelFormat = '%s'
#rep.SelectionPointSize = 0
#rep.SelectionPointLabelColor = [0,0,0]

ResetCamera()
