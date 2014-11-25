import paraview.simple as ps
import numpy as np

merge_groups_script="""
import vtk
self_input=self.GetInput()
self_output=self.GetOutput()
assert(self_input is not None)
assert(self_input.IsA("vtkMultiBlockDataSet"))
assert(self_output.IsA("vtkPolyData"))
if self_input.GetNumberOfBlocks() >0:
    groups=self_input.GetBlock(1)
    assert(groups.IsA("vtkMultiBlockDataSet"))
    append=vtk.vtkAppendPolyData()
    iterator = groups.NewIterator()
    iterator.InitTraversal()
    while not iterator.IsDoneWithTraversal():
        block=iterator.GetCurrentDataObject()
        assert(block.IsA("vtkPolyData"))
        append.AddInputData(block)
        iterator.GoToNextItem()
    append.Update()
    output.ShallowCopy(append.GetOutput())
"""

class WrongInput(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# get and check active source
source=ps.GetActiveSource()
if source is None:
    raise WrongInput('No source selected.')
info=source.GetDataInformation().DataInformation
if info is None:
    raise WrongInput('Source has no information, uncomplete.')
composite_class=info.GetCompositeDataClassName()
if composite_class is None or not composite_class=='vtkMultiBlockDataSet':   
    raise WrongInput('Source produce wrong type of data. MultiBlockDataSet required.')

# make a filter to collect all well lines in order
# to workaround a bug that do not porduce selsection labels 
# on composite datasets
merged=ps.ProgrammableFilter()
merged.Script=merge_groups_script
merged.OutputDataSetType=0 # PolyData
merged.Input=source

# update pipeline in order to get total number of points
ps.SetActiveSource(merged)
ps.UpdatePipeline()

# make selection where to show labels
n_points=merged.GetDataInformation().DataInformation.GetNumberOfPoints()
selection=ps.IDSelectionSource()
IDs = []
for i in range(n_points):
     # select only non-head points
    IDs.append(0L)
    IDs.append(long(i))
selection.IDs = IDs
selection.FieldType=1 # select point
merged.SetSelectionInput(0,selection,0)

# get representation and select its selection properties
ps.SetActiveSource(merged)
rep=ps.Show()
rep.SelectionPointFieldDataArrayName = 'label'
rep.SelectionPointLabelColor = [0,1,0]
rep.SelectionPointLabelFormat = "%s"
# try to do not displsy selection since 
# they are displaied with an offset 
rep.SelectionPointSize=0.001
rep.SelectionOpacity=0
rep.SelectionPointLabelVisibility = 1

# try to get correct orientation (Z pointing down)
camera=GetActiveCamera()
position=np.array(camera.GetPosition())
norm=np.linalg.norm(position)
camera.SetPosition(position[0], -norm, 0)

view=GetActiveView()
view.CameraViewUp=[0,0,-1]
Render()