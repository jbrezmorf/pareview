import vtk

assert(inputs[0].IsA("vtkMultiBlockDataSet"))
assert(output.IsA("vtkPolyData"))

groups=inputs[0].GetBlock(1)
assert(groups.IsA("vtkMultiBlockDataSet"))

#out_point_array=output.GetPoints().GetData()

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
    
    
    
