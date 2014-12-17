import threading
import paraview.simple 
import os
import sys
import numpy as np
from vtk.util import numpy_support
import vtk
from contextlib import contextmanager
from collections import namedtuple
import traceback

'''
Type of numpy arrays passed to VTK through numpy_support.
'''
vtk_int=np.int64
vtk_double=np.float64
'''
Simple timer.
'''
import time

class timewith():
    def __init__(self, name=''):
        self.name = name
        self.start = time.time()

    @property
    def elapsed(self):
        return time.time() - self.start

    def checkpoint(self, name=''):
        print '{timer} {checkpoint} took {elapsed} seconds'.format(
            timer=self.name,
            checkpoint=name,
            elapsed=self.elapsed,
        ).strip()

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.checkpoint('finished')
        pass

'''
Class for reading Eclipse files and producing Paraview datasets.
Attributes:
  self.egrid_header - directory with type definitions for various eclipse headers
  self.out - output object of the Programmable Filter, to avoid deep copies we store all
      auxiliary data directly here
  
  self.times - times of the restart file; set by RequestInformation
  self.file_names - file names named tuple with fields egrid and unrst. Set by SetFileName 
  self.grid - egrid data; set by read_egrid
  self.restart - array of data for time steps in UNRST file; set by read_restart
  self.output - VTK output object, set by ExtractFilterOutput
'''
class EclipseIO :
  
    '''
    Running guard
    '''
    @contextmanager
    def running_guard(self, persistent_object):
        #print persistent_object, hasattr(persistent_object,'Running')
        if not hasattr(persistent_object,'Running'):
            #print "init false"
            persistent_object.Running=False
         
        try:
            running=persistent_object.Running
            persistent_object.Running=True
            yield not running
        finally:
            #print "set false"
            persistent_object.Running=False
  
    VTK_HEXAHEDRON=12
    VTK_POLY_DATA=0

    '''
    Set definitions of format headers with proper endian (default is big-endian).
    Call with endian='<' for choosing little-endian.    
    '''
    def __init__(self, endian='>'):
       
        '''
        dtype specification for every header is stored under 
        its keyword in the egrid_header record.
        '''        
        self.types={'INTE':endian+'i4',
                    'DOUB':endian+'f8',
                    'LOGI':endian+'i4',
                    'CHAR':'a8',
                    'REAL':endian+'f',
                    'MESS':'a1'} # no data ??     

        self.block_spec={} 
             
        self.block_spec['FILEHEAD']=[int,
          'version',
          'release_year',
          'reserved_00',
          (4,'backward_compatibility_version'), # smallest version that can read the file
          'grid_type', # 0 - corner point; 1 - unstructured; 2 - hybrid
          'dual_porosity_model', # 0 - single porosity; 1 - dual porosity; 2 - dual permeability
          'original_grid_format',
          ]
        
        self.block_spec['MAPUNITS']=[str,'units']
        
        self.block_spec['MAPAXES ']=[float,          
          (0,'y_axis_end',2), # x,y coordinate tuple
          (0,'origin', 2),
          (0,'x_axis_end',2)
          ]
        
        self.block_spec['GRIDUNIT']=[str,
           (0,'units', 2)]
        
        # grid header for corner point grids
        self.block_spec['GRIDHEAD']=[int,
          'grid_type',
          (0,'dimensions',3), # (nx,ny,nz)
          'LGR_idx',
          (25,'numres'),
          'nseg',
          'ntheta',
          (0,'host_box',6) # (lower i, j, k, upper i, j, k)          
          ]

        self.block_spec['BOXORIG ']=[int,(0,'origin_xyz',3)]
               
        self.block_spec['reservoir']=[int,
          'lower_k_bound', 
          'upper_k_bound',
          'incomplete_circle',
          'isolate_reservoir',
          'lower_lateral_bound',
          'upper_lateral_bound',
          ] 
               
        self.block_spec['SEQNUM  ']=[int,'file_sequence_number']
        
        self.block_spec['INTEHEAD']=[int,
          'creation_time',
          (3,'units_type'), # 1-metric, 2-field, 3-lab
          (9,'dimensions',3),
          (12,'n_active_cells'),
          (15,'i_phase'), # [1-oil, water, oil/water, gas, oil/gas. gas/water, oil/water/gas]
          (17,'n_wells'),                                                
          (18,'n_max_completitions_per_well'),
          (20,'n_max_wells_per_group'),
          (21,'n_max_groups'),
          (25,'n_data_per_well'), # in IWELL array
          (28,'n_words_per_well'), # n of 8-char words in ZWELL array
          (33,'n_data_per_completition'), # in ICON array
          (37,'n_data_per_group'), # in IGRP array
          (65,'date',3), # date of the report time
          (95,'program_id'),
          (176,'n_max_segmented_wells'),
          (177,'n_max_segments_per_well'),
          (179,'n_data_per segment') # in ISEG array
          ]
        
        self.block_spec['LOGIHEAD']=[int,
          (4,'radial_model_flag_300'), # for eclipse 300
          (5,'radial_model_flag_100'), # for eclipse 100
          (15,'dual_porosity_flag'),
          (31,'coal_bed_methane_flag')
          ]                                     
        
        self.block_spec['DOUBHEAD']=[float,'time_in_days'] # day of the report step

        self.block_spec['well_data']=[int,
          (0,'wellhead_pos_ijk',3),  # well position; cell indices i,j,k
          'n_connections',              
          'i_group',
          'well_type', #1-producer; 2-oil injection; 3-water injection; 4-gass injection
          (11,'well_status'), # >0 open; <=0 shut
          (43,'i_LGR'),
          (49,'friction_flag'),
          (71, 'segment_well_number') # =0 for ordinary wells
          ]

        self.block_spec['completion']=[int,
          'connection_index', # -IC if no in current LGR
          (0,'coordinates',3), # ijk ???
          (6,'status'), # connection status >0 open, <=0 shut
          (14, 'penetration_direction'), # 1=x, 2=y, 3=z, 4=fractured in x, 5=fractured in y
          (15, 'segment') # segment containing connection, 0- for multi segment wells
          ]
        
        # Every field is tuple (in_label, out_label, unit)
        solution_fields=[
          ("PRESSURE", "Pressure",""),
          ("SWAT",      "WaterSat",""),
          ("SGAS",      "GasSat",""),
          ("SOIL",      "OilSat",""),
          ("RS",        "GasOilRatio",""),
          ("RV",        "OilGasRatio",""),
          ("OWC",       "OilWaterContact",""),
          ("OGC",       "OilGasContact",""),
          ("GWC",       "GasWaterContact",""),
          ("OILAPI",    "OilAPI",""),
          ("FIPOIL",    "OilFIP",""),
          ("FIPGAS",    "GasFIP",""),
          ("FIPWAT",    "WaterFIP",""),
          ("OIL-POTN",  "OilPotential",""),
          ("GAS-POTN",  "GasPotential",""),
          ("WAT-POTN",  "WaterPotential",""),
          ("POLYMER",   "PolymerConc",""),
          ("PADS",      "PolymerAdsorbedConc",""),
          ("XMF",       "LiquidMoleFrac",""),
          ("YMF",       "VaporMoleFrac",""),
          ("ZMF",       "TotalMoleFrac",""),
          ("SSOL",      "SolventSat",""),
          ("PBUB",      "BubblePressure",""),
          ("PDEW",      "DewPressure",""),
          ("SURFACT",   "SurfaceInteraction",""),
          ("SURFADS",   "AdsorbedSurfactant",""),
          ("SURFMAX",   "MaxSurfactantConc",""),
          ("SURFCNM",   "SurfactantCapilaryNumber",""),
          ("GGI",       "GI_InjectedGasRatio",""),
          ("WAT-PRES",  "WaterPressure",""),
          ("WAT_PRES",  "WaterPressure",""),
          ("GAS-PRES",  "GasPressure",""),
          ("GAS_PRES",  "GasPressure",""),
          ("OIL-VISC",  "OilViscosity",""),
          ("OIL_VISC",  "OilViscosity",""),
          ("VOIL",      "OilViscosity",""),
          ("WAT-VISC",  "WaterViscosity",""),
          ("WAT_VISC",  "WaterViscosity",""),
          ("VWAT",      "WaterViscosity",""),
          ("GAS-VISC",  "GasViscosity",""),
          ("GAS_VISC",  "GasViscosity",""),
          ("VGAS",      "GasViscosity",""),
          ("OIL-DEN",   "OilDensity",""),
          ("OIL_DEN",   "OilDensity",""),
          ("WAT-DEN",   "WaterDensity",""),
          ("WAT_DEN",   "WaterDensity",""),
          ("GAS-DEN",   "GasDensity",""),
          ("GAS_DEN",   "GasDensity",""),
          ("DRAINAGE",  "DrainageRegionNumber","")
          ]
        self.solution_fields={}
        for item in solution_fields:
            (key,label,unit)=item
            self.solution_fields[key]=(label,unit)
    
    '''
    Delete whole class as it is after __init__. Namely delete all VTK objects.
    This is called only when file name has changed.
    '''
    def __del__(self):
        #? should I delete these?
        #del self.output.corners
        #del self.output.cells
        #del self.output.points
        del self.file_names
        pass

    '''
    Get output object from the filter, check its type and store in self.
    '''
    def ExtractFilterOutput(self, program_filter):
        if not hasattr(self, "output"):
            multi_output = program_filter.GetOutput()
            if not multi_output.IsA("vtkMultiBlockDataSet"):
                print "Wrong output data type. Should be vtkMultiBlockDataSet."
                raise SystemExit
            self.output=multi_output        
        
        
    '''
    Separate filename base and creates filename for egrid file and restart
    file.
    '''
    def SetFileName(self, file_name):
        (base, ext)=os.path.splitext(file_name)
        egrid=base+".egrid"
        unrst=base+".unrst"
        if not os.path.isfile(egrid):
            egrid=None
        
        if not os.path.isfile(unrst):
            unrst=None
        FileNames=namedtuple('FileNames', 'egrid unrst')        
        self.file_names=FileNames(egrid, unrst)
    
    
    '''
    Search given 'keyword' in the given binary file 'f'.
    Seeks file 'f' to the first byte of  the keyword. 
    Returns -1 on fail.
    Note that keywords for eclipse files are always upper case 8 character long.
    '''
    @staticmethod
    def skip_to_keyword(f, keyword):
        data_str=""
        while True:
            new_data=f.read(1024*4)
            if (not new_data):
                return -1
            data_str += new_data 
            pos=data_str.find(keyword)
            if (pos >= 0):
                f.seek(-(len(data_str)-pos),os.SEEK_CUR)
                return 1
            else:
                data_str=data_str[len(data_str) - len(keyword):]
                
    
    '''
    Without parameters read the block and return tuple consisting of the 
    block name and block data.
    
    If the keyword parameter is given check that the block has correct keyword.
    If yes, read the block to the array of appropriate type and return it.
    If no, seek back and return None.
    
    Arrays are stored in named block using Fortran UNFORMATED data format (unfortunately not standardized).
    It seems that the format is as follows:
    position    bytes           meaning
    0x0000      8               blockname
    0x0008      4               number of stored items
    0x000c      4               type signature of the items
    0x0010      4               0x10 ??
    0x0014      4               N - subblock size in bytes
    ...         N 
                4               N - subblock size in bytes
    ...                         repetition of subblocks up to total number of items
                4               0x10 ??
    '''
    def read_array(self, f, keyword=None, check_type=None, check_size=None, optional=False):            
        in_keyword=f.read(8)
        if (keyword is not None and keyword!=in_keyword):
            if optional:
                f.seek(-8, os.SEEK_CUR)
                return None
            else:
                print "Missing obligatory block: ", keyword, "have block: ", in_keyword
                raise SystemExit

        size=np.fromfile(f, dtype=np.dtype(self.types['INTE']), count=1 )[0]
        elem_type=np.fromfile(f, dtype=np.dtype('S4'), count=1 )[0]
        elem_dtype=np.dtype(self.types[elem_type]) # possibly catch KeyError exception
        
        if (check_type is not None): 
            assert(self.types[check_type]==elem_dtype)
        if (check_size is not None):
            assert(check_size==size)
        f.read(4) # 0x10 block_code
        
        array=np.array([], dtype=elem_dtype)
        while (array.size < size):
            n_bytes=np.fromfile(f, dtype=np.dtype(self.types['INTE']), count=1 )[0]
            string=f.read(n_bytes)
            array=np.append(array, np.fromstring(string, dtype=elem_dtype))
            n_bytes_end=np.fromfile(f, dtype=np.dtype(self.types['INTE']), count=1 )[0]             
            assert(n_bytes==n_bytes_end)                         
        f.read(4) # block_code 0x10
        
        if (keyword is None):
            return (in_keyword, array)
        else:  
            return array
             
    
    
    
    '''
    Numpy array wrapper which check validity of indices and slices 
    and return None for an index out of the array.
    Further provides conversion for a dictionary.
    '''
    class Array:
      
        '''
        Construct from a one dim (numpy) array
        '''
        def __init__(self, np_array):
            self.array=np_array
            
        def __getitem__(self, key):
            if (isinstance(key, int)):
                if (key>len(self.array)): return None
            elif (isinstance(key,slice)):
                if (key.stop>len(self.array)): return None
            return self.array[key]  

        '''
        Read a block by read_array and convert it to dictionary
        using specification is self.header[keyword]. One block
        specification is an array, where the first is singleton of type of 
        produced dictionary items. Elements from input array are processed along with 
        items of specification array assigning the values from the input array to the keys
        given by the specification array. Specification may contain also tuple (position, key),
        to give explicit position of the key in the input array. Position must be greater then current 
        position in the input array effectively skip some of its elements.
        Position is numbered from 1. 
        
        Examples:
        'name' - set key 'name' to current value
        (20,'name') - jump forward to position 20 and set that value to the key 'name'
        (0,'name') - same as simply 'name'
        (0,'name',3) - form numpy array from 3 following values and assign to the key 'name'
        (10,'name',3) - jump forward and form a numpy array ...
        
        Vector items are numpy arrays.
        '''
        def make_dict(self, spec):
            in_array=self.array
            if (in_array is None or len(in_array) is 0):
                return None
            array_type=spec[0]  
            assert( np.can_cast(type(in_array[0]), array_type) 
                   or array_type == bool)  
            i_in_array=0
            
            result_dict={}
            for item in spec[1:]:
                vector_size=None # default are scalar values
                if type(item) is tuple:
                    position=item[0]-1
                    if (position>0): 
                        assert(position >= i_in_array)
                        i_in_array=position               
                    key=item[1]
                    if len(item)==3:
                        vector_size=item[2]
                else:
                    key=item
                              
                              
                if (i_in_array >= len(in_array)):
                    new_value=None
                else:    
                    if vector_size is None:
                        # scalar value
                        new_value=in_array[i_in_array]
                        i_in_array+=1
                    else:
                        # vector value
                        new_value=in_array[i_in_array:i_in_array+vector_size]
                        i_in_array+=vector_size                
                # print "Item: ", key, i_in_array, new_value
                result_dict[key]=array_type(new_value)
                
            return result_dict        

    
    
    '''
    Read a block by read_array and convert it to dictionary
    using specification is self.header[keyword]. ..
    '''
    def read_dict(self,f,keyword, optional=False):
        spec=self.block_spec[keyword]
        in_array=self.read_array(f,keyword, optional=optional)
        return self.Array(in_array).make_dict(spec)        

    '''
    Read a block by read_array and connects individual character octets into an
    array of strings.
    '''
    def read_string_array(self,f, keyword, n_strings, n_octets_per_string, optional=False):            
        octets=self.read_array(f, keyword, 'CHAR', n_strings*n_octets_per_string, optional)  
        octets.shape=(n_strings, n_octets_per_string)
        
        strings=n_strings*[""]
        for i in xrange(n_strings):
            tmp_str="".join(octets[i])
            strings[i]=tmp_str.strip()
            
        return strings    
    
    '''
    Reads EGRID file with name given by parameter 'filename' and
    store it into (empty) output object given by 'pdo' parameter.
    Returns updated pdo object.
    '''
    def read_egrid(self):
        if not self.file_names.egrid:
            print "No grid file. ABORT."
            raise SystemExit
          
        with open(self.file_names.egrid, 'rb') as f:
            # skip one int
            f.read(4)
            grid={}
            
            grid['filehead']=self.read_dict(f, "FILEHEAD")
            grid['mapunits']=self.read_dict(f, "MAPUNITS", optional=True)        
            grid['mapaxes']=self.read_dict(f,"MAPAXES ", optional=True)
            grid['grid_unit']=self.read_dict(f,"GRIDUNIT", optional=True)
            
            grid['gridhead']=self.read_dict(f,"GRIDHEAD")
            
            numres=grid['gridhead']['numres']
            #print "Numres: ", numres
            (nx,ny,nz) = grid['dimensions'] = grid['gridhead']['dimensions']            
            nlines=(nx+1)*(ny+1)*numres

            grid['boxorig']=self.read_dict(f,"BOXORIG ", optional=True)       
            grid['lines']=self.read_array(f,"COORD   ", 'REAL', 6*nlines)
            grid['lines'].shape=(nlines,6)
            
            res_data=self.read_array(f,"COORDSYS", 'INTE', 6*grid['gridhead']['numres'], optional=True)
            if (res_data):
                res_data.shape=(grid['gridhead']['numres'], 6)
                grid['reservoirs']=res_data
            
            grid['z_corners']=self.read_array(f,"ZCORN   ", 'REAL', 8*nx*ny*nz)
            active=self.read_array(f,"ACTNUM  ", 'INTE', nx*ny*nz, optional=True)
            if active is not None:             
                grid['active_cells']=(active==1) # create boolean array
                assert(active.size == nx*ny*nz)
            else:
                grid['active_cells']=None
            # 0-inactive, 1-active, 2-active fracture, 3-active matrix and fracture 
            
            grid['coarsening']=self.read_array(f,"CORSNUM ", 'INTE', nx*ny*nz, optional=True) 
            grid['hostcells']=self.read_array(f,"HOSTNUM ", 'INTE', nx*ny*nz, optional=True) #LGR only
            self.read_array(f,"ENDGRID ")
            
            assert(grid['lines'] != None)
            assert(grid['z_corners'] != None)
            
            #print "LINES\n",lines
            #print "CORNERS\n",z_corners
            self.grid=grid
        

    '''
    Create corresponding VTK mesh in given output object, that
    should be vtkUnstructuredGrid. Takes data from self.grid 
    and assumes that they persist. New data
    are created in provided object output.
    '''
    def create_grid(self):
        if not hasattr(self, "grid"):
            print "Grid data not loaded."
            raise SystemExit
          
        if hasattr(self, "vtk_grid"):
            return self.vtk_grid
          
        output=vtk.vtkUnstructuredGrid()

        # eclipse coordinate system:  
        #           / 
        #  x  <---|/
        #         |
        #         v z
        
        grid=self.grid
        (nx,ny,nz)=grid['dimensions']
        lines=grid['lines']
        z_corners=grid['z_corners']
        
        output.Reset()
        
        # number lines
        i_line=np.arange(lines.shape[0]).reshape((ny+1,nx+1))
        # distribute line numbers to corners in one X, Y layer
        # make them broadcastable to Z direction
        i_line=np.repeat( np.repeat(i_line, 2, axis=0), 2, axis=1)[1:-1, 1:-1].reshape(1,-1)
        assert(i_line.shape[1] == 2*nx*2*ny)
        # separate Z direction to use broadcasting of line numbers
        z_corners.shape=(-1, 2*nx*2*ny)
        assert(z_corners.shape[0]==2*nz)
        
        # interpolate X,Y of corners on lines
        z_line_top=lines[i_line,2]
        z_line_bot=lines[i_line,5]
        t= (z_corners - z_line_top)/(z_line_bot-z_line_top)
        x_corners=lines[i_line,0] - t*(lines[i_line,3] - lines[i_line,0])
        y_corners=lines[i_line,1] - t*(lines[i_line,4] - lines[i_line,1])
        
        # permute corners to individual cells
        corners=np.vstack( (x_corners.flatten(), y_corners.flatten(), z_corners.flatten())  )
        # permute corners on one cell into VTK numbering 
        local_points=np.arange(24).reshape(2,4,3)[:,[0,1,3,2],:]
        corn=np.transpose(corners.reshape(3,nz,2,ny,2,nx,2), axes=(1,3,5,2,4,6,0)).reshape(nx*ny*nz, 24)[:,local_points]
        # compute cell centers for wells
        self.cell_centers=np.mean(corn.reshape(nz,ny,nx,8,3), axis=3)
        
        # possibly remove non-active cells
        active=grid['active_cells']
        if active is not None:
            output.corners=corn[active, :]
        else:
            output.corners=corn
        n_cells=output.corners.shape[0]            
        output.corners.shape=(8*n_cells, 3)                  
        
        # create cell to corner map
        output.cells=np.empty((n_cells,9), dtype=vtk_int)
        output.cells[:,0]=8
        output.cells[:,1:9]=np.arange(8*n_cells).reshape((-1,8))               
        output.cells.shape=(-1)        
        
        output.points=vtk.vtkPoints()
        output.points.SetData(numpy_support.numpy_to_vtk(output.corners)) # 8*nx*ny*nz (x,y,z)
        output.SetPoints(output.points)
        
        output.cell_array = vtk.vtkCellArray()
        output.cell_array.SetCells(n_cells, numpy_support.numpy_to_vtkIdTypeArray(output.cells)) # nx*ny*nz (n,8*i_point)
        output.SetCells(self.VTK_HEXAHEDRON, output.cell_array) 
        
        self.vtk_grid=output
        
        return self.vtk_grid
        
       
        
    '''
    Read a restart file as an array of times
    self.restart[ step1, step2, ...]
    '''
    def read_restart(self):
        if not self.file_names.unrst:
            return
        
        with open(self.file_names.unrst, 'rb') as f:
            # skip one int
            f.read(4)
            self.restart=[]
            
            while (1):
                one_step={}
                one_step['seq_num']=self.read_dict(f,'SEQNUM  ')
                
                one_step['head']=self.read_dict(f,'INTEHEAD')
                
                one_step['head'].update(self.read_dict(f,'LOGIHEAD'))
                #print one_step['head']
                data=self.read_array(f,'DOUBHEAD') # !! much more complex then described, skip
                #print data
                #one_step['head'].update()
                
                n_groups=one_step['head']['n_max_groups']
                group_data_size=one_step['head']['n_data_per_group']
                n_wells_in_group=one_step['head']['n_max_wells_per_group']
                group_i_data=self.read_array(f,'IGRP    ', 'INTE', n_groups*group_data_size)
                group_i_data.shape=(n_groups, group_data_size)
                
                # first create array of groups with raw data
                group_data=[]
                for i_data in group_i_data:
                    one_group_data=self.Array(i_data)
                    data={}
                    childs=one_group_data[0:n_wells_in_group]
                    n_childs=one_group_data[n_wells_in_group]
                    data['childs']=childs[0:n_childs]
                    data['type']=one_group_data[n_wells_in_group+26]
                    # 0-well_group, 1-node_group (childs are groups), 2- satellite group, 3-slave group
                    data['level']=one_group_data[n_wells_in_group+27]
                    data['parent_group']=one_group_data[n_wells_in_group+28]                    
                    #print "====="
                    #print data
                    #print "-----"
                    #print one_group_data[n_wells_in_group+29:]
                    group_data.append(data)
                
                self.read_array(f,'SGRP    ')
                self.read_array(f,'XGRP    ')
                
                # undocumented block with group names
                n_words_per_group_name=5 # seems that this is fixed constant                
                name_array=self.read_string_array(f,'ZGRP    ', n_groups, n_words_per_group_name)
                for i_group in xrange(n_groups):                    
                    group_data[i_group]['name']=name_array[i_group]
                    #print "=====================", i_group
                    #print group_data[i_group]
                
                # create a tree from the raw group data
                #i_root=None
                #for i_grp in xrange(n_groups):
                #    data=group_data[i_grp]
                #    parent=data['parent_group']
                #    if parent==0 and data['group_level']==0:   
                #        i_root=parent-1
                #    group_data['subgroups']    
                
                one_step['group_data']=group_data 
                
                '''
                n_seg_wells=one_step['n_max_segmented_wells']
                n_seg_per_well=one_step['n_max_segments_per_well']
                segment_data_size=one_step['n_data_per_segment']
                segments_data=self.read_array(f,'ISEG    ', np.types['INTE'], n_seg_wells*n_seg_per_well*segment_data_size)
                segments_data.shape(n_seg_wells,n_seg_per_well,segment_data_size)
                for well in segments_data:
                    for segment in well:
                        outlet_segment_number=segment[1]
                        branch_for_segment=segment[3]
                '''
                
                n_wells=one_step['head']['n_wells']
                well_data_size=one_step['head']['n_data_per_well']                
                wells_data=self.read_array(f,'IWEL    ', 'INTE', n_wells*well_data_size)
                wells_data.shape=(n_wells, well_data_size)
                wells=[]
                for well_data in wells_data:
                    wells.append(self.Array(well_data)
                                .make_dict(self.block_spec['well_data']))
                
                self.read_array(f,'SWEL    ')
                self.read_array(f,'XWEL    ')
                
                name_array=self.read_string_array(f,'ZWEL    ', n_wells, one_step['head']['n_words_per_well'])
                for i_well in xrange(n_wells):                    
                    wells[i_well]['name']=name_array[i_well]
                
                n_completion=one_step['head']['n_max_completitions_per_well']
                n_per_completion=one_step['head']['n_data_per_completition']
                completion_data=self.read_array(f,'ICON    ', 'INTE', n_wells*n_completion*n_per_completion)
                completion_data.shape=(n_wells, n_completion, n_per_completion)
                for i_well in xrange(n_wells):
                    completions=[]
                    for compl in completion_data[i_well]:
                        completions.append(self.Array(compl)
                                          .make_dict(self.block_spec['completion']))
                    wells[i_well]['completions']=completions
                            
                one_step['wells']=wells
                
                #self.read_array(f,'SCON    ')
                #self.read_array(f,'XCON    ')
                #self.read_array(f,'DLYTIM  ')
                #self.read_array(f,'IAAQ    ')
                #self.read_array(f,'SAAQ    ')
                #self.read_array(f,'XAAQ    ')
                #self.read_array(f,'ICAQNUM ')
                #self.read_array(f,'ICAQ    ')
                #self.read_array(f,'SCAQNUM ')
                #self.read_array(f,'SCAQ    ')
                #self.read_array(f,'ACAQNUM '))
                #self.read_array(f,'ACAQ    '))                        
                
                # hidden contains names of solution arrays that 
                # are just internal and should not be post-processed
                # self.skip_to_keyword(f,'HIDDEN')
                # self.read_array(f,'HIDDEN  ') # skip
                #self.read_array(f,'ZTRACER ') # skip for now

                # read data fields  
                self.skip_to_keyword(f,'STARTSOL')
                data1=f.read(24) # 'STARTSOL',0x0,'MESS',0x10,0x10
                
                key=""
                fields=[]
                while (1):
                    (key,array)=self.read_array(f)
                    if (key == 'ENDSOL  '):
                        f.seek(-8,os.SEEK_CUR)
                        break
                    key=key.strip()
                    vtk_array=self.make_data_set(key, array)
                    if vtk_array:
                        fields.append( vtk_array )      
                one_step['fields']=fields    
                self.restart.append(one_step)
                
                #buff=str(f.read(4))
                #print "buff:", buff
                #f.seek(-4,os.SEEK_CUR)
                if (self.skip_to_keyword(f,'SEQNUM  ') == -1):
                    break 
            # end step loop
            
        
        
    '''
    Read only DOUBHEAD sections to retrieve output times.
    We also check SEQNUM section for consecutive sequence.
    '''
    def read_restart_times(self):
        if not self.file_names.unrst:
            self.times=[0]
            return
          
        with open(self.file_names.unrst, "rb") as f:
            i_time=0
            times=[]
            while (1):
                eof=self.skip_to_keyword(f, 'SEQNUM')
                if (eof==-1):
                    break
                i_time+=1  
                seq_num=self.read_dict(f,'SEQNUM  ')['file_sequence_number']
                #if (i_time != seq_num):
                #    print "Wrong sequence number", seq_num, " at position ", i_time
                #    raise AssertionError
                eof=self.skip_to_keyword(f, 'DOUBHEAD')
                if (eof==-1):
                    print "No DOUBHEAD section after SEQNUM section."
                    raise AssertionError
                times.append(self.read_dict(f,'DOUBHEAD')['time_in_days']) 
        self.times=times
                

    '''
    '''
    def add_dataset_to_multiblock(self, multiblock, dataset, name):
        n=multiblock.GetNumberOfBlocks()
        multiblock.SetBlock(n, dataset)
        multiblock.GetMetaData(n).Set(multiblock.NAME(), name)
        
        
    '''
    Add new PointData to given vkUnstructuredGrid with given name
    and data in a numpy array. Assumes a grid 
    '''
    def make_data_set(self, key, np_array):
        # check that key is known
        if self.solution_fields.has_key(key):
            name=self.solution_fields[key][0]
        else:
            return None

        in_dtype=np_array.dtype
        if issubclass(in_dtype.type, np.number):
            if issubclass(in_dtype.type, np.integer):
                fixed_np_array=np_array.astype(np.dtype(vtk_int))
                new_array=numpy_support.numpy_to_vtkIdTypeArray(fixed_np_array, deep=True)
            else:
                new_array=numpy_support.numpy_to_vtk(np_array.astype(np.dtype('d')), deep=True)
        else:
            return None  
        new_array.SetName(name)
        return new_array

    '''
    Create vtkPolyData representing the wells.
    Input - restart data for current time step.
    '''
    def make_wells(self, one_step):   
        n_total_points=0
        
        wells=one_step['wells']        
        groups=one_step['group_data']
        n_groups=len(groups)
        
        out=vtk.vtkMultiBlockDataSet()
        for i_group in xrange(n_groups):
            group=groups[i_group]
            if not group['type'] == 0: # non-well group
                continue
                          
            # check well indices
            group_well_ids=[]
            for i_well in group['childs']:
                #if not wells[i_well-1]['i_group']==i_group+1:
                #    print "Warning: well of a group do not point back", group['childs'], i_group, wells[i_well-1]['i_group']
                group_well_ids.append(i_well-1)
                    
            n_wells=len(group_well_ids)
            if n_wells==0:
                continue
            
            points=np.empty( (n_wells, 2, 3), dtype=self.cell_centers.dtype)
            lines=np.empty( (n_wells, 3), dtype=vtk_int)
            labels=vtk.vtkStringArray()
            labels.SetName("label")
            for i_well in xrange(n_wells):                
                well=wells[ group_well_ids[i_well] ]
                head_pos=well['wellhead_pos_ijk']-1
                well_type=well['well_type'] #1-producer; 2-oil injection; 3-water injection; 4-gass injection
                well_status=well['well_status'] # >0 open; <=0 shut
                name=well['name']
                #connections=well['completions']
                '''
                for conect in connections:
                    pos=conect['coordinates'] # ijk cell 
                    status=conect['status'] # connection status >0 open, <=0 shut
                '''
                # need cell centers for all cells
                #print self.cell_centers.shape
                #print head_pos 
                head=self.cell_centers[head_pos[2]-1, head_pos[1]-1, head_pos[0],0:3]
                
                #print i_well
                #print head
                points[i_well,0,0:3]=head
                labels.InsertNextValue("")
                
                head_shift=head-np.array([0,0,100])
                points[i_well,1,0:3]=head_shift
                lines[i_well,0]=2
                lines[i_well,1]=2*i_well
                lines[i_well,2]=2*i_well+1
                labels.InsertNextValue(name)
                i_well+=1
            
            points.shape=(-1,3)
            lines.shape=(-1)

            group_block=vtk.vtkPolyData()
            
            vtk_points=vtk.vtkPoints()
            vtk_points.SetData(numpy_support.numpy_to_vtk(points, deep=True))
            group_block.SetPoints(vtk_points)
            n_total_points+=points.shape[0]
            
            point_cells=vtk.vtkCellArray()   
            point_cells.SetCells(n_wells, numpy_support.numpy_to_vtkIdTypeArray(lines, deep=True))
            group_block.SetLines(point_cells)
            
            group_block.GetPointData().AddArray(labels)
            self.add_dataset_to_multiblock(out, group_block, group['name'])
            self.n_wells_poly_data_points=n_total_points

        return out   


        
    '''
    Make datasets from all fields on input.
    '''
    def set_all_data_sets(self, one_step, output):
        for vtk_field in one_step['fields']:
            name=vtk_field.GetName()
            cell_data=output.GetCellData()
            output_array=cell_data.GetArray(name)
            #help(cell_data.GetArray)
            if output_array:
                # this possibly could be done better
                cell_data.RemoveArray(name)
                cell_data.AddArray(vtk_field)
            else:
                cell_data.AddArray(vtk_field)
            
    
    
    '''
    Setting information about the filter output.
    '''
    def RequestInformation(self, program_filter):
        try:
            with self.running_guard(program_filter) as r:
              if r:
                self.read_restart_times()
                self.np_times=np.array(self.times)
                
                executive=program_filter.GetExecutive()
                #help(executive.__class__)
                out_info=executive.GetOutputInformation(0)
                #out_info=executive.GetOutputInformation().GetInformationObject(0)
                out_info.Remove(executive.TIME_STEPS())
                for time in self.times:
                    out_info.Append(executive.TIME_STEPS(), time)
                    #out_info.Append(vtkStreamingDemandDrivePipeline.TIME_STEPS(), time)
                out_info.Remove(executive.TIME_RANGE())
                out_info.Append(executive.TIME_RANGE(), self.times[0])
                out_info.Append(executive.TIME_RANGE(), self.times[-1])
                #print out_info
                #print "Times:", times
        except BaseException:
            print "== Eclipse Reader Exception (RequestInfo) =="
            (et, ex, tr)=sys.exc_info()
            print "Exception: ", et, ex
            traceback.print_tb(tr)
                
            
    '''
    Get timestep to which we should set the data on the grid.
    '''
    def GetUpdateTimeStep(self, program_filter):
        # get requested timestep or the frist one if not present
        executive = program_filter.GetExecutive()
        out_info = executive.GetOutputInformation(0)
        #print out_info
        if not out_info.Has(executive.UPDATE_TIME_STEP()):
            return self.times[0]
        else:
            return out_info.Get(executive.UPDATE_TIME_STEP())

            
    '''
    Setting data to the filter output
    '''
    def RequestData(self, program_filter):
        try:
            with self.running_guard(program_filter) as r:
              if r:
                if not hasattr(self, "grid"):
                    self.read_egrid() 
                
                self.ExtractFilterOutput(program_filter)
                timestep=self.GetUpdateTimeStep(program_filter)
                
                # optionaly create the grid
                if not self.output.GetBlock(0):                              
                    self.output.SetBlock(0, self.create_grid() )
                    self.output.GetMetaData(0).Set(self.output.NAME(), "eclipse grid");                    
                
                # optionaly read data
                if self.file_names.unrst:
                    if not hasattr(self, "restart"):
                        self.read_restart()
                
                    #find time
                    i_time=np.abs(self.np_times-timestep).argmin()
                    timestep=self.times[i_time]
                    # make datasets        
                    self.set_all_data_sets(self.restart[i_time], self.output.GetBlock(0))
                    # mark correct timestep 
                    self.output.GetInformation().Set(self.output.DATA_TIME_STEP(), timestep)

                    # create well groups  block  
                    self.add_dataset_to_multiblock(self.output, 
                                                   self.make_wells(self.restart[i_time]),
                                                   "eclipse well groups")
                    

                # possibly reset camera to get correct view
                # need working recursion prevention
                # Finally it is done automaticaly but we may want other orientation of the system.
                #print "reset camera"
                #paraview.simple.ResetCamera()
                
                # Create BlockSelection (do not work)
                '''
                group_block_ids=[]                
                iterator = self.output.NewIterator()
                iterator.InitTraversal()
                while not iterator.IsDoneWithTraversal():
                    obj_type=iterator.GetCurrentDataObject().GetDataObjectType()
                    if  obj_type== self.VTK_POLY_DATA:
                        group_block_ids.append( iterator.GetCurrentFlatIndex() )
                    iterator.GoToNextItem()
                
                print group_block_ids
                self.group_block_ids=group_block_ids

                #sel_source.FieldType="point"
                
                #help(self.programmable_filter)
                #self.programmable_filter.SetSelectionInput(0, sel_source, 0)
                
                #paraview.simple.UpdatePipeline()
                source=paraview.simple.FindSource(self.reader_name)
                paraview.simple.SetActiveSource(source)
                
                sel_source=paraview.simple.BlockSelectionSource()
                sel_source.Blocks=group_block_ids
                
                rep=paraview.simple.Show()
                rep.SelectionPointLabelVisibility = 1
                rep.SelectionPointFieldDataArrayName = 'labels'
                rep.SelectionPointLabelFormat = '%s'
                #rep.SelectionPointSize = 0
                #rep.SelectionPointLabelColor = [1,1,1]
                print "done RequestData"
                '''
                #Merge to vtkPolyData
                              
                '''
                #print "find source:", self.reader_name
                #paraview.simple.UpdatePipeline()
                source=paraview.simple.FindSource(self.reader_name)
                print source
                
                #print "set active"
                

                #print "make merged"
                merged=paraview.simple.ProgrammableFilter()
                merged.Script=self.merge_groups_script
                merged.OutputDataSetType=0 # PolyData
                merged.Input=source
                
                
                print "make selection"                
                #n_points=merged.GetDataInformation().DataInformation.GetNumberOfPoints()
                n_points=self.n_wells_poly_data_points
                selection=paraview.simple.IDSelectionSource()
                IDs = []
                for i in range(n_points):
                    IDs.append(0L)
                    IDs.append(long(i))
                selection.IDs = IDs
                selection.FieldType=1
                merged.SetSelectionInput(0,selection,0)

                

                print "Show"  
                paraview.simple.SetActiveSource(merged)
                print paraview.simple.Show()
                reps=paraview.simple.GetRepresentations()
                for (rep_name, rep_id) in reps:
                    print (rep_name, rep_id)
                    if rep_name=='GeometryRepresentation2':
                        key=(rep_name, rep_id)
                rep=reps[key]
                print rep
                
                
                print "setup selection"  
                rep.SelectionPointFieldDataArrayName
                rep.SelectionPointFieldDataArrayName = 'label'
                rep.SelectionPointLabelColor = [0,0,0]
                rep.SelectionPointLabelFormat = "%s"
                rep.SelectionPointLabelVisibility = 1
                #rep.SelectionPointSize = 0
                #

                #print "run ResetCamera"
                #paraview.simple.ResetCamera()
                print "done"
                '''
                
        except BaseException:
            print "== Eclipse Reader Exception =="
            (et, ex, tr)=sys.exc_info()
            print "Exception: ", et, ex
            traceback.print_tb(tr)
            
            
            
            # clean output
            for i_block in range(self.output.GetNumberOfBlocks()):
                self.output.RemoveBlock(i_block)
                


      
if (not hasattr(self,"info_done_event")):
    self.info_done_event=threading.Event()
      
# main RequestInformation script
# every internal data are stored in self.code
if hasattr(self, "code"):
    del self.code

self.code=EclipseIO()
self.code.programmable_filter=self
self.code.SetFileName(FileName)
self.code.RequestInformation(self)
# indicator that RequestInformation is done
self.info_done_event.set()


 