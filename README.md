pareview
========

Python Paraview plugin implementing reader for results produced by Eclipse simulator.

Files:
  eclipse_reader_template.xml - XML template without the script body
  eclipse_reader_main.py - main body of the reader, this can be copy-pasted into a programmable source (for testing)
  make_plugin_xml.py - script that creates final XML of the plugin (adding script into the template)
  makefile
  test.egrid - test mesh
  test_eclipse_reader.py - unit tests


Eclipse file extensions:
DBG - debug output
EGRID - egrid mesh file
INIT - initial file
INSPEC
MSG - message output
PRT - TXT, main printer output
RFT - vector file (graphing data)
RSM - TXT, Run Summary
RSSPEC - unformatted index for restart files
SMSPEC - summary specification file
         (root file name, number of vector data - wells, positions of wells (cell or region number))
UNRST - unified restart file
UNSMRY - unified summary file (vector data at ministeps)
