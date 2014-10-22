import paraview.simple as ps

with open("eclipse_reader.py","rt") as f:
    request_data_script=f.read()
    
with open("eclipse_reader_info.py","rt") as f:
    request_information_script=f.read()

program=ps.ProgrammableSource()
program.OutputDataSetType=13 # vtkMultiblockDataSet
program.Script=request_data_script
program.ScriptRequestInformation=request_information_script
program.UpdatePipeline()

output.ShallowCopy(program.GetOutput())
