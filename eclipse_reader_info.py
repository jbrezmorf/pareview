import threading
import paraview.simple 
import os
import sys
import numpy as np
from vtk.util import numpy_support
import vtk
from contextlib import contextmanager
from collections import namedtuple




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
            #print persistent_object, persistent_object.Running
            if not persistent_object.Running:
                #print "set true"
                persistent_object.Running=True
                yield
            else:
                pass
                #print "skip"
        finally:
            #print "set false"
            persistent_object.Running=False
  
    VTK_HEXAHEDRON=12

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
        
        self.block_spec['MAPAXES']=[float,          
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
          (0,'wellhead_pos_ijk',3),
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
    Check that next block has correct keyword.
    If yes, read the block to the array of appropriate type and return it.
    If no, seek back and return None.
    '''
    def read_array(self, f, keyword=None, check_type=None, check_size=None):            
        in_keyword=f.read(8)
        if (keyword is not None and keyword!=in_keyword):
            f.seek(-8, os.SEEK_CUR)
            return None
          

        size=np.fromfile(f, dtype=np.dtype(self.types['INTE']), count=1 )[0]
        elem_type=np.fromfile(f, dtype=np.dtype('S4'), count=1 )[0]
        elem_dtype=np.dtype(self.types[elem_type]) # possibly catch KeyError exception
        
        if (check_type is not None): 
            assert(self.types[check_type]==elem_dtype)
        if (check_size is not None):
            assert(check_size==size)
        f.read(8) # 0x10 block_code
        array=np.fromfile(f, dtype=elem_dtype, count=size)
        f.read(8) # block_code 0x10
        
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
        (0,'name',3) - form array from 3 following values and assign to the key 'name'
        (10,'name',3) - jump forward and form array ...
        '''
        def make_dict(self, spec):
            in_array=self.array
            if (in_array is None or len(in_array) is 0):
                return None
            assert( np.can_cast(type(in_array[0]), spec[0]) 
                   or spec[0] == bool)  
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
                result_dict[key]=new_value    
                
            return result_dict        

    
    
    '''
    Read a block by read_array and convert it to dictionary
    using specification is self.header[keyword]. ..
    '''
    def read_dict(self,f,keyword):
        spec=self.block_spec[keyword]
        in_array=self.read_array(f,keyword)
        return self.Array(in_array).make_dict(spec)        
              
    
    
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
            grid['mapunits']=self.read_dict(f, "MAPUNITS")        
            grid['mapaxes']=self.read_dict(f,"MAPAXES")
            
            grid['gridhead']=self.read_dict(f,"GRIDHEAD")
            (nx,ny,nz) = grid['dimensions'] = grid['gridhead']['dimensions']
            
            nlines=(nx+1)*(ny+1)*grid['gridhead']['numres']

            grid['boxorig']=self.read_dict(f,"BOXORIG ")       
            grid['lines']=self.read_array(f,"COORD   ", 'REAL', 6*nlines)
            grid['lines'].shape=(nlines,6)
            
            res_data=self.read_array(f,"COORDSYS", 'INTE', 6*grid['gridhead']['numres'])
            if (res_data):
                res_data.shape=(grid['gridhead']['numres'], 6)
                grid['reservoirs']=res_data
            
            grid['z_corners']=self.read_array(f,"ZCORN   ", 'REAL', 8*nx*ny*nz)
            grid['active_cells']=self.read_array(f,"ACTNUM  ", 'INTE', nx*ny*nz) 
            # 0-inactive, 1-active, 2-active fracture, 3-active matrix and fracture 
            
            grid['coarsening']=self.read_array(f,"CORSNUM ", 'INTE', nx*ny*nz) 
            grid['hostcells']=self.read_array(f,"HOSTNUM ", 'INTE', nx*ny*nz) 
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
    def create_grid(self, output):
        if not output.IsA("vtkUnstructuredGrid"):
            print "Creating grid in wrong dataset. Should be vtkUnstructuredGrid"
            raise SystemExit
        
        grid=self.grid
        (nx,ny,nz)=grid['dimensions']
        lines=grid['lines']
        z_corners=grid['z_corners']
        
        output.Reset()
        output.corners=np.empty(3*8*nx*ny*nz, dtype=float)
        output.cells=np.empty(9*nx*ny*nz, dtype=int)
        
        i_point=0
        i_corner=0
        i_cell=0
        
        # eclipse coordinate system:  
        #           / 
        #  x  <---|/
        #         |
        #         v z
        
        
        # local coordinates (x,y,z)
        hexahedron_local_vtx=[(0,0,0),(1,0,0),(1,1,0),(0,1,0),(0,0,1),(1,0,1),(1,1,1),(0,1,1)]
        for iz in xrange(nz) :
            for iy in xrange(ny) :
                for ix in xrange(nx) :
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
        output.SetCells(self.VTK_HEXAHEDRON, output.cell_array) 
        
          
        
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
                groups_data=self.read_array(f,'IGRP    ', 'INTE', n_groups*group_data_size)
                groups_data.shape=(n_groups, group_data_size)
                groups=[]
                for group_data in groups_data:
                    one_group_data=self.Array(group_data)
                    group={}
                    group['childs']=one_group_data[0:n_wells_in_group-1]
                    group['n_childs']=one_group_data[n_wells_in_group]
                    group['group_type']=one_group_data[n_wells_in_group+26]
                    # 0-well_group, 1-node_group (childs are groups), 2- satellite group, 3-slave group
                    group['group_level']=one_group_data[n_wells_in_group+27]
                    group['parent_group']=one_group_data[n_wells_in_group+28]
                    
                    groups.append(group)
                one_step['groups']=groups    
                self.read_array(f,'SGRP    ')
                self.read_array(f,'XGRP    ')
                self.read_array(f,'ZGRP    ')
                
                
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
                
                n_words_per_well=one_step['head']['n_words_per_well']
                well_names_data=self.read_array(f,'ZWEL    ', 'CHAR', n_wells*n_words_per_well)  
                well_names_data.shape=(n_wells, n_words_per_well)
                for i_well in xrange(n_wells):
                    name="".join(well_names_data[i_well])
                    wells[i_well]['name']=name.strip()
                
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
                
                self.read_array(f,'SCON    ')
                self.read_array(f,'XCON    ')
                self.read_array(f,'DLYTIM  ')
                self.read_array(f,'IAAQ    ')
                self.read_array(f,'SAAQ    ')
                self.read_array(f,'XAAQ    ')
                self.read_array(f,'ICAQNUM ')
                self.read_array(f,'ICAQ    ')
                self.read_array(f,'SCAQNUM ')
                self.read_array(f,'SCAQ    ')
                
                self.read_array(f,'HIDDEN  ') # skip
                self.read_array(f,'ZTRACER ') # skip for now

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
                print "restart size:", len(self.restart)
                
                #buff=str(f.read(4))
                #print "buff:", buff
                #f.seek(-4,os.SEEK_CUR)
                if (self.skip_to_keyword(f,'SEQNUM  ') == -1):
                    print "EOF"
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
                if (i_time != seq_num):
                    print "Wrong sequence number", seq_num, " at position ", i_time
                    raise AssertionError
                eof=self.skip_to_keyword(f, 'DOUBHEAD')
                if (eof==-1):
                    print "No DOUBHEAD section after SEQNUM section."
                    raise AssertionError
                times.append(self.read_dict(f,'DOUBHEAD')['time_in_days']) 
        self.times=times
                
        
        
    '''
    Add new PointData to given vkUnstructuredGrid with given name
    and data in a numpy array. Assumes a grid 
    '''
    def make_data_set(self, key, np_array):
        # check tht key is known
        if self.solution_fields.has_key(key):
            name=self.solution_fields[key][0]
        else:
            return None

        in_dtype=np_array.dtype
        if issubclass(in_dtype.type, np.number):
            if issubclass(in_dtype.type, np.integer):
                fixed_np_array=np_array.astype(np.dtype('i'))
                new_array=numpy_support.numpy_to_vtkIdTypeArray(fixed_np_array, deep=True)
            else:
                new_array=numpy_support.numpy_to_vtk(np_array.astype(np.dtype('d')), deep=True)
        else:
            return None  
        new_array.SetName(name)
        return new_array


        
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
        with self.running_guard(program_filter):
            
            self.read_restart_times()
            self.np_times=np.array(self.times)
            
            executive=program_filter.GetExecutive()
            help(executive.__class__)
            out_info=executive.GetOutputInformation(0)
            #out_info=executive.GetOutputInformation().GetInformationObject(0)
            out_info.Remove(executive.TIME_STEPS())
            for time in self.times:
                out_info.Append(executive.TIME_STEPS(), time)
                #out_info.Append(vtkStreamingDemandDrivePipeline.TIME_STEPS(), time)
            out_info.Remove(executive.TIME_RANGE())
            out_info.Append(executive.TIME_RANGE(), self.times[0])
            out_info.Append(executive.TIME_RANGE(), self.times[-1])
            print out_info
            #print "Times:", times
            
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
        with self.running_guard(program_filter):
            self.ExtractFilterOutput(program_filter)
            timestep=self.GetUpdateTimeStep(program_filter)
            print "TimeStep: ", timestep
            
            # optionaly create the grid
            grid_block=self.output.GetBlock(0)
            if not grid_block:
                self.read_egrid()
                grid_block=vtk.vtkUnstructuredGrid()  
                self.create_grid(grid_block)
                self.output.SetBlock(0, grid_block)
            
            if hasattr(self, 'sphere'):
                del self.sphere
            self.sphere=vtk.vtkSphereSource()
            self.sphere.SetCenter((0,0,timestep))
            self.sphere.Update()
            self.output.SetBlock(1, self.sphere.GetOutput())                            
            
            #optionaly read data
            if not hasattr(self,"restart"):
                self.read_restart()
            
            # set data
            
            #find time
            i_time=np.abs(self.np_times-timestep).argmin()
            print "I time:", i_time
            timestep=self.times[i_time]
            print "New time:", timestep
            print "restart size:", len(self.restart)
            # make datasets        
            self.set_all_data_sets(self.restart[i_time], grid_block)
            # mark correct timestep 
            self.output.GetInformation().Set(self.output.DATA_TIME_STEP(), timestep)
            
            # possibly reset camera to get correct view
            # need working recursion prevention
            print "reset camera"
            #paraview.simple.ResetCamera()



      
if (not hasattr(self,"info_done_event")):
    self.info_done_event=threading.Event()
      
# main RequestInformation script
# every internal data are stored in self.code
if hasattr(self, "code"):
    del self.code


FileName="test.egrid"    
self.code=EclipseIO()   
self.code.SetFileName(FileName)
self.code.RequestInformation(self)
# indicator that RequestInformation is done
self.info_done_event.set()

 