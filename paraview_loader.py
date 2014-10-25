import paraview.simple as ps
import vtk

with open("eclipse_reader_data.py","rt") as f:
    request_data_script=f.read()
    
with open("eclipse_reader_info.py","rt") as f:
    request_information_script=f.read()

program=vtk.vtkPythonProgrammableFilter()
program.SetOutputDataSetType(13) # vtkMultiblockDataSet
program.SetScript(request_data_script)
program.SetInformationScript(request_information_script)
#help(program)
program.SetParameter("FileName","test.unrst")
program.Update()

output.ShallowCopy(program.GetOutput())
