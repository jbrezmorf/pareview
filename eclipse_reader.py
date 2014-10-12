
 
 import os
 import numpy as np
 from vtk.util import numpy_support
 
 
 
 VTK_HEXAHEDRON=12
 
 
 class EclipseIO :
     running_reader=False
     FileName=None
 
     '''
     Set definitions of format headers with proper endian (default is big-endian).
     Call with endian='<i4' for choosing little-endian.
     
     '''
     def __init__(self, endian='>'):
         '''
         dtype specification for every header is stored under 
         its keyword in the egrid_header record.
         '''
         self.egrid_header={} 
         i4=self.i4=endian+'i4'
         real=self.real=endian+'f'
               
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
     def read_egrid_mesh(self, output, filename):
         self.grid_file=open(filename, 'rb')
         # skip one int
         self.grid_file.read(4)
         
         filehead=self.read_header( "FILEHEAD")
         mapunits=self.read_header( "MAPUNITS")        
         mapaxes=self.read_header("MAPAXES")
         
         gridhead=self.read_header("GRIDHEAD")
         (nx,ny,nz) = gridhead['dimensions']
         nlines=(nx+1)*(ny+1)*gridhead['numres']
 
         boxorig=self.read_header("BOXORIG ")       
         lines=self.read_array("COORD   ", np.dtype(self.real), 6*nlines)
         reservoirs=self.read_array("COORDSYS", self.egrid_header['reservoir'], count=gridhead['numres'])
         z_corners=self.read_array("ZCORN   ", np.dtype(self.real), count=8*nx*ny*nz)
         activecells=self.read_array("ACTNUM  ", np.dtype(self.i4), count=nx*ny*nz) 
         # 0-inactive, 1-active, 2-active fracture, 3-active matrix and fracture 
         
         coarsening=self.read_array("CORSNUM ", np.dtype(self.i4), count=nx*ny*nz) 
         hostcells=self.read_array("HOSTNUM ", np.dtype(self.i4), count=nx*ny*nz) 
         self.read_header("ENDGRID ")
         
         assert(lines != None)
         assert(z_corners != None)
         lines.shape=(nlines,6)
         
         # Create corresponding VTK mesh in pdo object
         #print "LINES\n",lines
         #print "CORNERS\n",z_corners
         
         self.grid_file.close()
         
         output.corners=np.empty(3*8*nx*ny*nz, dtype=float)
         output.cells=np.empty(9*nx*ny*nz, dtype=int)
         
         i_point=0
         i_corner=0
         i_cell=0
         
         # eclipse coordinate system:  
         #           / y 
         #  x  <---|/
         #         |
         #         v z
         
         
         # local coordinates (x,y,z)
         hexahedron_local_vtx=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
         for ix in xrange(nx) :
             for iy in xrange(ny) :
                 for iz in xrange(nz) :
                     # print "CELL = ", i_cell
                     output.cells[i_cell]=8 # number of vertices
                     i_cell+=1
                     # set corners of one cell and cell indices to points
                     
                     # cell vertical lines are at (ix,ix+1) x ( iy, iy+1) in coords
                     # cell z coords are at 
                     
                     for i_vtx in xrange(8):
                         output.cells[i_cell]=i_point
                         i_cell+=1
                         i_point+=1
                         
                         loc_x,loc_y,loc_z=hexahedron_local_vtx[i_vtx]
                         i_zcoord=2*ix+loc_x + 2*nx*(2*iy+loc_y)+ 2*nx*2*ny*(2*iz+loc_z)
                         z_coord= z_corners[i_zcoord]
                         line=lines[(ix+loc_x) + (nx+1)*(iy+loc_y)]
                         top=(line[0], line[1])
                         z_top=line[2]
                         bot=(line[3], line[4])
                         z_bot=line[5]
                         t=(z_coord-z_bot)/(z_top-z_bot)
                         (x_coord,y_coord)= top*t + bot*(1-t)
                         output.corners[i_corner] = x_coord
                         i_corner+=1
                         output.corners[i_corner] = y_coord
                         i_corner+=1
                         output.corners[i_corner] = z_coord
                         i_corner+=1
                         
                         # print "    vtx: ", i_vtx, x_coord, y_coord, z_coord
         
         output.corners.shape=(8*nx*ny*nz, 3)                  
         output.points=vtk.vtkPoints()
         output.points.SetData(numpy_support.numpy_to_vtk(output.corners)) # 8*nx*ny*nz (x,y,z)
         output.SetPoints(output.points)
         
         output.cell_array = vtk.vtkCellArray()
         output.cell_array.SetCells(nx*ny*nz, numpy_support.numpy_to_vtkIdTypeArray(output.cells)) # nx*ny*nz (n,8*i_point)
         output.SetCells(VTK_HEXAHEDRON, output.cell_array) 
         
 
 
 def main():
     if (not EclipseIO.running_reader):
         EclipseIO.running_reader=True
         
         if (EclipseIO.FileName==None):
             EclipseIO.FileName=self.GetProgressText()
             self.SetProgressText(None)
             if (EclipseIO.FileName==None):
                 raise IOError("No input filename.")
 
         #print "f:", FileName
         # for debugging
         #if (FileName==None or FileName==""):
         #    FileName='/home/jb/workspace/pareview/test.egrid'
         #print "f:", FileName
 
         output = self.GetOutput()
         io=EclipseIO()
         io.read_egrid_mesh(output,EclipseIO.FileName)
         
         #from paraview import servermanager
         #help(servermanager)
         #help(servermanager.GetRenderView() )
 
         #print servermanager.GetRenderView()
         #view=servermanager.CreateRenderView()
         #view.ResetCamera() 
         import paraview.simple 
         paraview.simple.ResetCamera() 
 
 if __name__ == '__main__':
     main()  
 