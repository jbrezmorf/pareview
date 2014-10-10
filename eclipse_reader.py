import os
import numpy as np

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

class EclipseIO :

    '''
    Set definitions of format headers with proper endian (default is big-endian).
    Call with endian='<i4' for choosing little-endian.
    
    '''
    def __init__(self, endian='>i4'):
        '''
        dtype specification for every header is stored under 
        its keyword in the egrid_header record.
        '''
        self.egrid_header={} 
        i4=self.i4=endian
        real=self.real='>f'
               
        self.egrid_header['FILEHEAD']=np.dtype([ 
          ('version', i4),
          ('release_year',i4),
          ('reserved_00',i4),
          ('backward_compatibility_version',i4),  # smallest version that can read the file
          ('grid_type',i4), # 0 - corner point; 1 - unstructured; 2 - hybrid
          ('dual_porosity_model', i4), # 0 - single porosity; 1 - dual porosity; 2 - dual permeability
          ('original_grid_format', i4),
          ('not_used',i4,(93,))])
        
        self.egrid_header['MAPUNITS']=np.dtype([('units', 'S1')])
        
        self.egrid_header['MAPAXES']=np.dtype([
          ('y_axis_end',real, (2,)), # x,y coordinate tuple
          ('origin',real, (2,)),
          ('x_axis_end',real, (2,))
          ])
        
        self.egrid_header['GRIDUNIT']=np.dtype([('units', 'S2')])
        
        # grid header for corner point grids
        self.egrid_header['GRIDHEAD']=np.dtype([
          ('grid_type',i4),
          ('dimensions',i4, (3,)), # (nx,ny,nz)
          ('LGR_idx',i4),
          ('not_used_1',i4,(19,)),
          ('numres',i4),
          ('nseg',i4),
          ('ntheta',i4),
          ('host_box',i4,(6,)), # (lower i, j, k, upper i, j, k)
          ('not_used_2',i4,(67,))
          ])

        self.egrid_header['BOXORIG ']=np.dtype([('origin_xyz',i4, (3,))])
        
        self.egrid_header['COORDSYS']=np.dtype([
          ('size', i4), 
          ('type','S4'),
          ]) 
        
        self.egrid_header['reservoir']=np.dtype([
          ('lower_k_bound', i4), 
          ('upper_k_bound',i4),
          ('incomplete_circle',i4),
          ('isolate_reservoir',i4),
          ('lower_lateral_bound',i4),
          ('upper_lateral_bound',i4)
          ]) 
        
        self.egrid_header['ENDGRID ']=np.dtype([])
        
        
    '''
    Search given 'keyword' in the given binary file 'f'.
    Seeks file 'f' to the first byte of  the keyword. 
    Returns -1 on fail.
    Note that keywords for eclipse files are always upper case 8 character long.
    '''
    def skip_after_keyword(self,f, keyword):
        data_str=""
        while True:
            new_data=f.read(1024*4)
            if (not new_data):
                return -1
            data_str += new_data 
            pos=data_str.find(keyword)
            if (pos >= 0):
                f.seek(-(len(data_str)-pos)+len(keyword),os.SEEK_CUR)
                return 1
            else:
                data_str=data_str[len(data_str) - len(keyword):]
           
    '''
    Return dictionary with data from the header. Possible keys:

    '''
    def read_header(self,keyword):
        in_keyword=self.grid_file.read(8)
        if (keyword==in_keyword):
          self.grid_file.read(8) # skip size and type

          dtype=self.egrid_header[keyword]
          if ( dtype==np.dtype([]) ):
              return None

          self.grid_file.read(8) # 0x10 block_code
          header=np.fromfile(self.grid_file, dtype=dtype, count=1)
          self.grid_file.read(8) # block_code 0x10
          return header[0]
        else: 
          self.grid_file.seek(-8, os.SEEK_CUR)
          return None
        
    def read_array(self,keyword,elem_type,count):    
        in_keyword=self.grid_file.read(8)
        if (keyword==in_keyword):
          self.grid_file.read(8) # skip size and type
          self.grid_file.read(8) # 0x10 block_code
          array=np.fromfile(self.grid_file, dtype=elem_type, count=count)
          self.grid_file.read(8) # block_code 0x10
          return array
        else: 
          self.grid_file.seek(-8, os.SEEK_CUR)
          return None

    '''
    Reads EGRID file with name given by parameter 'filename' and
    store it into (empty) output object given by 'pdo' parameter.
    Returns updated pdo object.
    '''
    def read_egrid_mesh(self,pdo, filename):
        self.grid_file=open(filename, 'rb')
        # skip one int
        self.grid_file.read(4)
        
        filehead=self.read_header( "FILEHEAD")
        mapunits=self.read_header( "MAPUNITS")        
        mapaxes=self.read_header("MAPAXES")
        
        gridhead=self.read_header("GRIDHEAD")
        (nx,ny,nz) = gridhead['dimensions']
        nlines=6*(nx+1)*(ny+1)*gridhead['numres']

        boxorig=self.read_header("BOXORIG ")       
        lines=self.read_array("COORD   ", np.dtype(self.real), nlines)
        reservoirs=self.read_array("COORDSYS", self.egrid_header['reservoir'], count=gridhead['numres'])
        corners=self.read_array("ZCORN   ", np.dtype(self.real), count=8*nx*ny*nz)
        activecells=self.read_array("ACTNUM  ", np.dtype(self.i4), count=nx*ny*nz) 
        # 0-inactive, 1-active, 2-active fracture, 3-active matrix and fracture 
        
        coarsening=self.read_array("CORSNUM ", np.dtype(self.i4), count=nx*ny*nz) 
        hostcells=self.read_array("HOSTNUM ", np.dtype(self.i4), count=nx*ny*nz) 
        self.read_header("ENDGRID ")
        
        assert(lines != None)
        assert(corners != None)
        
        print "LINES\n",lines
        print "CORNERS\n",corners
        
        self.grid_file.close()


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