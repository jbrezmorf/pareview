pareview
========

Python Paraview plugin implementing reader for results produced by Eclipse simulator.

Files:
  README.md - readme and basic documentation
  eclipse_reader.xml - encapsuled plugin file
  show_wells.py - macro to show labels of wells
  
  eclipse_reader_template.xml - XML template without the script body
  eclipse_reader_data.py - script implementing RequestData, just call particular method in Eclipse class
  eclipse_reader_info.py - RequestInformation script, Eclipse class
  make_plugin_xml.py - script that creates final XML of the plugin (adding script into the template)
  makefile
  test.egrid - test mesh
  test.unrst - test data
  test_eclipse_reader.py - unit tests

Installation
============

The plugin was tested with ParaView version 4.2.
1) Put eclipse_reader.xml into lib/paraview-4.2/plugins of your 
   paraview installation.
2) Start paraview.
3) In the menu choose:  Tools -> Manage Plugins ... 
   You should see "eclipse_reader" plugin in the list. 
   If not continue with following:
   Choose "Load New...".
   In "Files of type:", choose "XML".
   Select your copy of eclipse_reader.xml.
   Expand eclipse_reader item in the list of plugins. And
   Check "Auto Load"
4) Use Macros -> Add new macro ... to select "show_wells.py" file.
   This step also copy the file into ParaView directory.

Usage
=====
To read Eclipse result data use File -> Open ...
and select either "*.egrid" or "*.unrst" file. 

Wait a moment until Python libraries are loaded. 
Use green "Apply" button. To actually read the data.

In order to show well's labels use Macro -> show_wells. 
This creates two additional filters in the pipeline and change some 
view settings in order to show the labels. These can be deleted at any time
without influence on the main data.

Set Decimation
==============
In order to have fast response, Paraview use decimation. The limit on data size 
for using the decimation is quite low and can be enlarged in 
Edit -> Settings... -> "Render View" tab.
Set the first slider name "Interactive Rendering Options" to size of the data in MB.

Elevate data
============
To elevate the data you have two options:
1) Use "Transform" filter: Filters -> Aplhabetical -> Transform, and set last component of Scale.
   This creates elevated copy of your data. With Z-coordinates scaled.
2) In Properties panel, Properties tab use wheel icon to show advanced properties.
   Scroll down to "Transforming" section ans set last component of Scale.
   This just change the view of your data, so Z-coordinates are preserved.


=====================================================================
Development notes:

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
