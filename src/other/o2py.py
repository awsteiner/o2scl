"""
  -------------------------------------------------------------------

  Copyright (C) 2006-2016, Andrew W. Steiner
  
  This file is part of O2scl.
  
  O2scl is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  O2scl is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with O2scl. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
"""
"""
  Todos:
  - Be able to close HDF5 file and retain data set?
"""

import getopt, sys, h5py, math, os, hashlib
import matplotlib.pyplot as plot
from matplotlib.colors import LinearSegmentedColormap
import urllib.request
import numpy

class cloud_file:

    force_subdir=False
    env_var=''
    username=''
    password=''
    verbose=1
    
    def md5(fname):
        """
        Compute the md5 hash of the specified file. This function
        reads 4k bytes at a time and updates the hash for each
        read in order to prevent from having to read the entire
        file in memory at once.
        """
        hash_md5 = hashlib.md5()
        with open(fname, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def download_file(self,data_dir,fname_orig,url,mhash):
        force_subdir_val=self.force_subdir
        self.force_subdir=False
        fname=self.download_file_subdir(data_dir,'',fname_orig,url,
                                        mhash)
        self.force_subdir=force_subdir_val
        return fname

    def download_file_subdir(self,data_dir,subdir_orig,fname_orig,url,
                             mhash):
        """
        This function attempts to find a file named 'fname_orig' in
        subdirectory 'subdir_orig' of the data directory 'data_dir'. If
        'data_dir' is empty, it attempts to set it equal to the value of
        the environment variable 'env_var'. If that environment variable
        is not present, the user is prompted for the correct data
        directory. If the file is not found, then this function uses curl
        (or wget if curl was unsuccessful) to download the file from
        'url'. If this process was successful at finding or downloading
        the file, then the full filename is returned. Otherwise, an
        exception is thrown.
        
        Warning: this function has several potential security issues 
        and should not be used without due care.
    
        """
        # First obtain the data directory
        method=''
        if data_dir!='':
            method='arg'
        else:
            if self.env_var in os.environ:
                data_dir=os.environ[self.env_var]
                method='ev'
            if data_dir=='':
                data_dir=input('Data directory not set. Enter data directory: ')
                if data_dir!='':
                    method='ui'
        if data_dir=='' or method=='':
            raise ValueError('Failed to obtain data directory.')
        if method=='arg':
            if verbose>0:
                print('Data directory set (by function argument) to:',data_dir)
        elif method=='ev':
            if verbose>0:
                print('Data directory set (by environment variable) to:',
                      data_dir)
        else:
            if verbose>0:
                print('Data directory set (by user input) to:',data_dir)
    
        # Now test for the existence of the subdirectory and create it if
        # necessary
        subdir=''
        if self.force_subdir==True:
            subdir=data_dir+'/'+subdir_orig
            if os.path.isdir(subdir)==False:
                if verbose>0:
                    print('Directory not found and force_subdir is true.')
                    print('Trying to automatically create using "mkdir -p"')
                cmd='mkdir -p '+subdir
                ret=os.system(cmd)
                if ret!=0:
                    raise FileNotFoundError('Correct subdirectory does '+
                                            'not exist and failed to create.')
        else:
            subdir=data_dir

        # The local filename
        fname=subdir+'/'+fname_orig
        hash_match=False
        if os.path.isfile(fname)==True:
            mhash2=mda5(fname)
            if mhash==mhash2:
                hash_match=True
            elif verbose>0:
                print('Hash of file',fname,'did not match',mhash)
        elif verbose>0:
            print('Could not find file',fname)
            
        # Now download the file if it's not already present
        if hash_match==False or os.path.isfile(fname)==False:
            response=input('Hash did not match or data file '+fname+
                           ' not found. Download (y/Y/n/N)? ')
            ret=1
            if response=='y' or response=='Y':
                if verbose>0:
                    print('Trying two download:')
                urllib.request.urlretrieve(url,fname)
                ret=0
            if ret==0:
                mhash2=mda5(fname)
                if mhash!=mhash2:
                    raise ConnectionError('Downloaded file but '+
                                          'has does not match.')
            if ret!=0:
                raise ConnectionError('Failed to obtain data file.')
    
        # Return the final filename
        return fname

class hdf5_reader:

    list_of_dsets=[]
    search_type=''

    def hdf5_is_object_type(self,name,obj):
        """
        If object 'obj' named 'name' is of type 'search_type', then add
        that name to 'list_of_dsets'
        """
        # Convert search_type to a bytes object
        search_type_bytes=bytes(self.search_type,'utf-8')
        if isinstance(obj,h5py.Group):
            if 'o2scl_type' in obj.keys():
                o2scl_type_dset=obj['o2scl_type']
                if o2scl_type_dset.__getitem__(0) == search_type_bytes:
                    self.list_of_dsets.append(name)
        return

    def h5read_first_type(self,fname,loc_type):
        """
        Read the first object of type 'loc_type' from file named 'fname'
        """
        del self.list_of_dsets[:]
        self.search_type=loc_type
        file=h5py.File(fname,'r')
        file.visititems(self.hdf5_is_object_type)
        if len(self.list_of_dsets)==0:
            str='Could not object of type '+loc_type+' in file '+fname+'.'
            raise RuntimeError(str)
        return file[self.list_of_dsets[0]]

    def h5read_name(self,fname,name):
        """
        Read object named 'name' from file named 'fname'
        """
        file=h5py.File(fname,'r')
        obj=file[name]
        o2scl_type_dset=obj['o2scl_type']
        loc_type=o2scl_type_dset.__getitem__(0)
        return (obj,loc_type)
    
    def h5read_type_named(self,fname,loc_type,name):
        """
        Read object of type 'loc_type' named 'name' from file named 'fname'
        """
        del self.list_of_dsets[:]
        self.search_type=loc_type
        file=h5py.File(fname,'r')
        file.visititems(self.hdf5_is_object_type)
        if name in self.list_of_dsets:
            return file[name]
        str='No object of type '+loc_type+' named '+name+' in file '+fname+'.'
        raise RuntimeError(str)
        return
    

# This is probably best replaced by get_str_array() below
#
# def parse_col_names(dset):
#     nc=dset['nc'].__getitem__(0)
#     nw=dset['nw'].__getitem__(0)
#     counter=dset['counter']
#     data=dset['data']
#     clist=[]
#     k=0
#     for i in range(0,nw):
#         column=''
#         for j in range(0,counter[i]):
#             column=column+str(unichr(data[k]))
#             k=k+1
#         clist.append(column)
#     return clist

def default_plot(lmar=0.14,bmar=0.12,rmar=0.04,tmar=0.04):
    """
    Common plot defaults
    """
    plot.rc('text',usetex=True)
    plot.rc('font',family='serif')
    plot.rcParams['lines.linewidth']=0.5
    fig=plot.figure(1,figsize=(6.0,6.0))
    fig.set_facecolor('white')
    ax=plot.axes([lmar,bmar,1.0-lmar-rmar,1.0-tmar-bmar])
    ax.minorticks_on()
    ax.tick_params('both',length=12,width=1,which='major')
    ax.tick_params('both',length=5,width=1,which='minor')
    plot.grid(False)
    return (fig,ax)
    
def get_str_array(dset):
    """
    Extract a string array from O2scl HDF5 dataset dset as a list
    """
    nw=dset['nw'][0]
    nc=dset['nc'][0]
    data=dset['data']
    counter=dset['counter']
    char_counter=1
    word_counter=0
    list=[]
    col=''
    for ix in range(0,nc):
        # Skip empty strings in the array
        done=0
        while done==0:
            if word_counter==nw:
                done=1
            elif counter[word_counter]==0:
                word_counter=word_counter+1
                list.append('')
            else:
                done=1
        col=col+str(chr(data[ix]))
        if char_counter==counter[word_counter]:
            list.append(col)
            col=''
            word_counter=word_counter+1
            char_counter=1
        else:
            char_counter=char_counter+1
    # We're done with the characters, but there are some blank
    # strings left. Add the appropriate blanks at the end.
    while word_counter<nw:
        list.append('')
        word_counter=word_counter+1
    return list
    
def parse_arguments(argv,verbose=0):
    """
    Parse command
    """
    list=[]
    unproc_list=[]
    if verbose>1:
        print('Number of arguments:', len(argv), 'arguments.')
        print('Argument List:', str(argv))
    ix=1
    while ix<len(argv):
        if verbose>1:
            print('Processing index',ix,'with value',argv[ix],'.')
        # Find first option, at index ix
        initial_ix_done=0
        while initial_ix_done==0:
            if ix==len(argv):
                initial_ix_done=1
            elif argv[ix][0]=='-':
                initial_ix_done=1
            else:
                if verbose>1:
                     print('Adding',argv[ix],' to unprocessed list.')
                unproc_list.append(argv[ix])
                ix=ix+1
        # If there is an option, then ix is its index
        if ix<len(argv):
            list_one=[]
            # Strip single and double dashes
            cmd_name=argv[ix][1:]
            if cmd_name[0]=='-':
                cmd_name=cmd_name[1:]
            # Add command name to list
            list_one.append(cmd_name)
            if verbose>1:
                print('Found option',cmd_name,'at index',ix)
            # Set ix_next to the next option, or to the end if
            # there is no next option
            ix_next=ix+1
            ix_next_done=0
            while ix_next_done==0:
                if ix_next==len(argv):
                    ix_next_done=1
                elif argv[ix_next][0]=='-':
                    ix_next_done=1
                else:
                    if verbose>1:
                        print('Adding '+argv[ix_next]+' with index '+
                              str(ix_next)+' to list for '+cmd_name)
                    list_one.append(argv[ix_next])
                    ix_next=ix_next+1
            list.append(list_one)
            ix=ix_next
    return (list,unproc_list)
                    
class plotter:

    h5r=hdf5_reader()
    logx=0
    logy=0
    logz=0
    xtitle=''
    ytitle=''
    xlo=0
    xhi=0
    xset=0
    ylo=0
    yhi=0
    yset=0
    zlo=0
    zhi=0
    zset=0
    verbose=1
    dset=0
    axes=0
    fig=0
    canvas_flag=0
    dtype=''
    cmap='jet'
    colbar=0

    def myreds(self):
        cdict={'red': ((0.0,1.0,1.0),(1.0,1.0,1.0)),
               'green': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'blue': ((0.0,1.0,1.0),(1.0,0.0,0.0))}
        myreds=LinearSegmentedColormap('MyReds',cdict)
        plot.register_cmap(cmap=myreds)
        self.cmap='MyReds'

    def mygreens(self):
        cdict={'red': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'green': ((0.0,1.0,1.0),(1.0,1.0,1.0)),
               'blue': ((0.0,1.0,1.0),(1.0,0.0,0.0))}
        mygreens=LinearSegmentedColormap('MyGreens',cdict)
        plot.register_cmap(cmap=mygreens)
        self.cmap='MyGreens'

    def myblues(self):
        cdict={'red': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'green': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'blue': ((0.0,1.0,1.0),(1.0,1.0,1.0))}
        myblues=LinearSegmentedColormap('MyBlues',cdict)
        plot.register_cmap(cmap=myblues)
        self.cmap='MyBlues'
        
    def contour_plot(self,level,**kwargs):
        if self.dtype!='vector<contour_line>':
            print('Wrong type for contour_plot.')
            return
        if self.verbose>2:
            print('contour_plot',level,kwargs)
        if self.canvas_flag==0:
            self.canvas()
        n_lines=self.dset['n_lines'][0]
        for i in range(0,n_lines):
            line_level=self.dset['line_'+str(i)+'/level'][0]
            if abs(level-line_level) < 1.0e-7:
                if self.logx==1:
                    if self.logy==1:
                        plot.loglog(self.dset['line_'+str(i)+'/x'],
                                    self.dset['line_'+str(i)+'/y'],**kwargs)
                    else:
                        plot.semilogx(self.dset['line_'+str(i)+'/x'],
                                      self.dset['line_'+str(i)+'/y'],**kwargs)
                else:
                    if self.logy==1:
                        plot.semilogy(self.dset['line_'+str(i)+'/x'],
                                      self.dset['line_'+str(i)+'/y'],**kwargs)
                    else:
                        plot.plot(self.dset['line_'+str(i)+'/x'],
                                  self.dset['line_'+str(i)+'/y'],**kwargs)
        return
 
    def points(self,colx,coly,**kwargs):
        if self.dtype!='table':
            print('Wrong type for plot.')
            return
        if self.verbose>2:
            print('plot',colx,coly,kwargs)
        if self.canvas_flag==0:
            self.canvas()
        if self.logx==1:
            if self.logy==1:
                plot.loglog(self.dset['data/'+colx],
                            self.dset['data/'+coly],lw=0,marker='.',**kwargs)
            else:
                plot.semilogx(self.dset['data/'+colx],
                              self.dset['data/'+coly],lw=0,marker='.',**kwargs)
        else:
            if self.logy==1:
                plot.semilogy(self.dset['data/'+colx],
                              self.dset['data/'+coly],lw=0,marker='.',**kwargs)
            else:
                plot.plot(self.dset['data/'+colx],
                          self.dset['data/'+coly],lw=0,marker='.',**kwargs)
        if self.xset==1:
            plot.xlim([self.xlo,self.xhi])
        if self.yset==1:
            plot.ylim([self.ylo,self.yhi])
        return

    def plot(self,colx,coly,**kwargs):
        
        if self.force_bytes(self.dtype)==b'table':
            if self.verbose>2:
                print('plot',colx,coly,kwargs)
            if self.canvas_flag==0:
                self.canvas()
            if self.logx==1:
                if self.logy==1:
                    plot.loglog(self.dset['data/'+colx],
                                self.dset['data/'+coly],**kwargs)
                else:
                    plot.semilogx(self.dset['data/'+colx],
                                  self.dset['data/'+coly],**kwargs)
            else:
                if self.logy==1:
                    plot.semilogy(self.dset['data/'+colx],
                                  self.dset['data/'+coly],**kwargs)
                else:
                    plot.plot(self.dset['data/'+colx],
                              self.dset['data/'+coly],**kwargs)
            if self.xset==1:
                plot.xlim([self.xlo,self.xhi])
            if self.yset==1:
                plot.ylim([self.ylo,self.yhi])
        elif self.force_bytes(self.dtype)==b'hist':
            size=dset['size'][0]
            bins=dset['bins']
            weights=dset['weights']
            rmode=dset['rmode'][0]
            reps=bins[0:size-1]
            for i in range(0,size):
                reps[i]=(bins[i]+bins[i+1])/2
            if self.logx==1:
                if self.logy==1:
                    plot.loglog(reps,weights,**kwargs)
                else:
                    plot.semilogx(reps,weights,**kwargs)
            else:
                if self.logy==1:
                    plot.semilogy(reps,weights,**kwargs)
                else:
                    plot.plot(reps,weights,**kwargs)
            if self.xset==1:
                plot.xlim([self.xlo,self.xhi])
            if self.yset==1:
                plot.ylim([self.ylo,self.yhi])
            return
        return

    def plot1(self,col,**kwargs):
        if self.force_bytes(self.dtype)!=b'table':
            print('Wrong type for plot1.')
            return
        if self.verbose>2:
            print('plot1',col,kwargs)
        if self.canvas_flag==0:
            self.canvas()
        tlist=range(1,len(self.dset['data/'+col])+1)
        if self.logx==1:
            if self.logy==1:
                plot.loglog(tlist,self.dset['data/'+col],**kwargs)
            else:
                plot.semilogx(tlist,self.dset['data/'+col],**kwargs)
        else:
            if self.logy==1:
                plot.semilogy(tlist,self.dset['data/'+col],**kwargs)
            else:
                plot.plot(tlist,self.dset['data/'+col],**kwargs)
        if self.xset==1:
            plot.xlim([self.xlo,self.xhi])
        if self.yset==1:
            plot.ylim([self.ylo,self.yhi])
        return

    def plot1m(self,col,files,**kwargs):
        if self.verbose>2:
            print('plot1m',col,kwargs)
        if self.canvas_flag==0:
            self.canvas()
        for i in range(0,len(files)):
            self.dset=self.h5r.h5read_first_type(files[i],'table')
            self.dtype='table'
            tlist=range(1,len(self.dset['data/'+col])+1)
            if self.logx==1:
                if self.logy==1:
                    plot.loglog(tlist,self.dset['data/'+col],**kwargs)
                else:
                    plot.semilogx(tlist,self.dset['data/'+col],**kwargs)
            else:
                if self.logy==1:
                    plot.semilogy(tlist,self.dset['data/'+col],**kwargs)
                else:
                    plot.plot(tlist,self.dset['data/'+col],**kwargs)
        if self.xset==1:
            plot.xlim([self.xlo,self.xhi])
        if self.yset==1:
            plot.ylim([self.ylo,self.yhi])
        return
    
    def plotm(self,colx,coly,files,**kwargs):
        if self.verbose>2:
            print('plotm',colx,coly,kwargs)
        if self.canvas_flag==0:
            self.canvas()
        for i in range(0,len(files)):
            self.dset=self.h5r.h5read_first_type(files[i],'table')
            self.dtype='table'
            if self.logx==1:
                if self.logy==1:
                    plot.loglog(self.dset['data/'+colx],
                                self.dset['data/'+coly],**kwargs)
                else:
                    plot.semilogx(self.dset['data/'+colx],
                                  self.dset['data/'+coly],**kwargs)
            else:
                if self.logy==1:
                    plot.semilogy(self.dset['data/'+colx],
                                  self.dset['data/'+coly],**kwargs)
                else:
                    plot.plot(self.dset['data/'+colx],
                              self.dset['data/'+coly],**kwargs)
        if self.xset==1:
            plot.xlim([self.xlo,self.xhi])
        if self.yset==1:
            plot.ylim([self.ylo,self.yhi])
        return
    
    def hist(self,col,**kwargs):
        if self.verbose>2:
            print('hist',kwargs)
        if self.canvas_flag==0:
            self.canvas()
        if self.force_bytes(self.dtype)==b'table':
            plot.hist(self.dset['data/'+col],**kwargs)
        else:
            print('Wrong type',self.dtype,'for hist()')
        return

    def hist2d(self,colx,coly,**kwargs):
        if self.verbose>2:
            print('hist',colx,coly,kwargs)
        if self.canvas_flag==0:
            self.canvas()
        plot.hist2d(self.dset['data/'+colx],self.dset['data/'+coly],**kwargs)
        return

    def line(self,x1,y1,x2,y2,**kwargs):
        if self.verbose>2:
            print('Line',x1,y1,x2,y1)
        if self.canvas_flag==0:
            self.canvas()
        plot.plot([x1,x2],[y1,y2],**kwargs)
        return

    def reset_xlimits(self):
        self.xset=0
        return

    def xlimits(self,xlo,xhi):
        self.xlo=xlo
        self.xhi=xhi
        self.xset=1
        if self.canvas_flag==1:
            plot.xlim([xlo,xhi])
        return

    def reset_ylimits(self):
        self.yset=0
        return

    def ylimits(self,ylo,yhi):
        self.ylo=ylo
        self.yhi=yhi
        self.yset=1
        if self.canvas_flag==1:
            plot.ylim([ylo,yhi])
        return

    def canvas(self):
        if self.verbose>2:
            print('Canvas')
        # Default o2mpl plot
        (self.fig,self.axes)=default_plot()
        # Plot limits
        if self.xset==1:
            plot.xlim([self.xlo,self.xhi])
        if self.yset==1:
            plot.ylim([self.ylo,self.yhi])
        # Titles
        if self.xtitle!='':
            plot.xlabel(self.xtitle,fontsize=16)
        if self.ytitle!='':
            plot.ylabel(self.ytitle,fontsize=16)
        self.canvas_flag=1
        return

    def move_labels(self):
        # Move tick labels
        for label in self.axes.get_xticklabels():
            t=label.get_position()
            t2=t[0],t[1]-0.01
            label.set_position(t2)
            label.set_fontsize(16)
        for label in self.axes.get_yticklabels():
            t=label.get_position()
            t2=t[0]-0.01,t[1]
            label.set_position(t2)
            label.set_fontsize(16)
        return

    def show(self):
        plot.show()
        return

    def save(self,filename):
        if self.verbose>0:
            print('Saving as',filename,'.')
        plot.savefig(filename)
        return

    def read(self,filename):
        if self.verbose>0:
            print('Reading file',filename,'.')
        self.dset=self.h5r.h5read_first_type(filename,'table')
        self.dtype='table'
        return

    def read_type(self,filename,loc_type):
        if self.verbose>0:
            print('Reading object of type',loc_type,
                  'in file',filename,'.')
        self.dset=self.h5r.h5read_first_type(filename,loc_type)
        self.dtype=loc_type
        return

    def read_name(self,filename,name):
        if self.verbose>0:
            print('Reading object named',name,'in file',filename,'.')
        atuple=self.h5r.h5read_name(filename,name)
        self.dset=atuple[0]
        self.dtype=atuple[1]
        return

    def force_bytes(self,obj):
        """
        In cases where we're unsure whether or not obj is
        a string or bytes object, we ensure it's a bytes
        object by converting if necessary.
        """
        if isinstance(obj,numpy.bytes_)==False:
            return bytes(obj,'utf-8')
        return obj
    
    def list(self):
        if self.force_bytes(self.dtype)==b'table':
            col_list=get_str_array(self.dset['col_names'])
            if self.verbose>2:
                print('-----------------------')
            unit_list=[]
            unit_flag=self.dset['unit_flag'][0]
            print('unit_flag',unit_flag)
            if self.verbose>2:
                print('unit_flag:',unit_flag)
            if unit_flag==1:
                unit_list=get_str_array(self.dset['units'])
                if self.verbose>2:
                    print('-----------------------')
                    print('unit_list:',unit_list)
            print(len(col_list),'columns.')
            for ix in range(0,len(col_list)):
                if unit_flag:
                    print(str(ix)+'. '+col_list[ix]+' ['+unit_list[ix]+']')
                else:
                    print(str(ix)+'. '+col_list[ix])
            print(self.dset['nlines'][0],'lines.')
            if self.verbose>2:
                print('Done in list')
        elif self.force_bytes(self.dtype)==b'table3d':
            sl_list=get_str_array(self.dset['slice_names'])
            print(len(sl_list),'slices.')
            for ix in range(0,len(sl_list)):
                print(str(ix)+'. '+sl_list[ix])
            xgrid=self.dset['xval'].value
            ygrid=self.dset['yval'].value
            lxgrid=len(xgrid)
            lygrid=len(ygrid)
            print('X-grid start: '+str(xgrid[0])+' end: '+
                  str(xgrid[lxgrid-1])+' size '+str(lxgrid))
            print('Y-grid start: '+str(ygrid[0])+' end: '+
                  str(ygrid[lygrid-1])+' size '+str(lygrid))
        else:
            print('Cannot list type',self.dtype)
        return

    def den_plot(self,slice_name,**kwargs):
        if self.force_bytes(self.dtype)==b'table3d':
            name='data/'+slice_name
            sl=self.dset[name].value
            sl=sl.transpose()
            xgrid=self.dset['xval'].value
            ygrid=self.dset['yval'].value
            if self.canvas_flag==0:
                self.canvas()
            if self.logx==1:
                for i in range(0,len(xgrid)):
                    xgrid[i]=math.log(xgrid[i],10)
            if self.logy==1:
                for i in range(0,len(ygrid)):
                    ygrid[i]=math.log(ygrid[i],10)
            lx=len(xgrid)
            ly=len(ygrid)
            plot.imshow(sl,cmap=self.cmap,interpolation='nearest',
                        origin='lower',
                        extent=[xgrid[0]-(xgrid[1]-xgrid[0])/2,
                                xgrid[lx-1]+(xgrid[lx-1]-xgrid[lx-2])/2,
                                ygrid[0]-(ygrid[1]-ygrid[0])/2,
                                ygrid[ly-1]+(ygrid[ly-1]-ygrid[ly-2])/2],
                        aspect='auto',**kwargs)
            if self.colbar>0:
                plot.colorbar()
                
        else:
            print('Cannot density plot object of type',self.dtype)
        return

    def set(self,name,value):
        if self.verbose>0:
            print('Set',name,'to',value)
        if name=='logx':
            self.logx=int(value)
        elif name=='logy':
            self.logy=int(value)
        elif name=='xtitle':
            self.xtitle=value
        elif name=='ytitle':
            self.ytitle=value
        elif name=='xlo':
            self.xlo=float(value)
            self.xset=1
        elif name=='xhi':
            self.xhi=float(value)
            self.xset=1
        elif name=='xset':
            self.xset=int(value)
        elif name=='ylo':
            self.ylo=float(value)
            self.yset=1
        elif name=='yhi':
            self.yhi=float(value)
            self.yset=1
        elif name=='yset':
            self.yset=int(value)
        elif name=='zlo':
            self.zlo=float(value)
            self.zset=1
        elif name=='zhi':
            self.zhi=float(value)
            self.zset=1
        elif name=='zset':
            self.zset=int(value)
        elif name=='verbose':
            self.verbose=int(value)
        elif name=='cmap':
            self.cmap=value
        elif name=='colbar':
            self.colbar=int(value)
        else:
            print('No variable named',name)
        return

    def help(self,arg=''):
        if arg=='canvas':
            print('canvas()')
        elif arg=='plot':
            print('---------------------------------------------------------')
            print('plot(x,y,**kwargs)')
            print('---------------------------------------------------------')
            print('Useful kwargs:')
            print('color or c                   matplotlib color')
            temp=str('dashes                       '+
                     'sequence of on/off ink in points')
            print(temp)
            temp=str("fillstyle                    "+
                     "['full' | 'left' | 'right' | 'bottom' | 'top'")
            print(temp)
            print("                              | 'none']")
            temp=str("label                        "+
                     "string or anything printable with '%s' conversion.")
            print(temp)
            temp=str("linestyle or ls              "+
                     "['-' | '--' | '-.' | ':' | 'None' | ' ' | '']")
            print(temp)
            temp=str("                              "+
                     "and any drawstyle in combination with a")
            print(temp)
            temp=str("                              "+
                     "linestyle, e.g., 'steps--'.")
            print(temp)
            print("linewidth or lw              float value in points")
            print('marker                       marker type')
            print('markeredgecolor or mec       matplotlib color')
            print('markeredgewidth or mew       float value in points')
            print('markerfacecolor or mfc       matplotlib color')
            print('markerfacecoloralt or mfcalt matplotlib color')
            print('markersize or ms             float')
            print('---------------------------------------------------------')
        elif arg=='den_plot':
            print('den_plot(slice_name)')
        elif arg=='get':
            print('get(name)')
            print('')
            print('logx,logy,xtitle,ytitle,xlo,xhi,xset,ylo,yhi,yset')
            print('zlo,zhi,zset,verbose,cmap,colbar')
        elif arg=='line':
            print('line(x1,y1,x2,y2,color,lstyle)')
        elif arg=='list':
            print('list()')
        elif arg=='move_labels':
            print('move_labels()')
        elif arg=='read':
            print('read(filename)')
        elif arg=='read_name':
            print('read_name(filename,name)')
        elif arg=='read_type':
            print('read_type(filename,type)')
        elif arg=='reset_xlimits':
            print('reset_xlimits()')
        elif arg=='reset_ylimits':
            print('reset_ylimits()')
        elif arg=='set':
            print('set(name,value)')
            print('')
            print('logx,logy,xtitle,ytitle,xlo,xhi,xset,ylo,yhi,yset')
            print('zlo,zhi,zset,verbose,cmap')
        elif arg=='show':
            print('show()')
        elif arg=='text':
            print('text(x,y,string,color)')
        elif arg=='xlimits':
            print('xlimits(xlow,xhigh)')
        elif arg=='ylimits':
            print('ylimits(ylow,yhigh)')
        else:
            print('Methods:')
            print('canvas()')
            print('contour_plot(level,**kwargs)')
            print('den_plot(slice_name,**kwargs)')
            print('get(name)')
            print('hist([col],**kwargs)')
            print('hist2d(colx,coly,**kwargs)')
            print('line(x1,y1,x2,y2,**kwargs)')
            print('list()')
            print('move_labels()')
            print('parse_argv(argv)')
            print('plot(colx,coly,**kwargs)')
            print('read(filename)')
            print('read_name(filename,name)')
            print('read_type(filename,type)')
            print('reds()')
            print('reset_xlimits()')
            print('reset_ylimits()')
            print('set(name,value)')
            print('show()')
            print('text(x,y,string,**kwargs)')
            print('xlimits(xlow,xhigh)')
            print('ylimits(ylow,yhigh)')
        return

    def text(self,tx,ty,str,**kwargs):
        if self.canvas_flag==0:
            self.canvas()
        self.axes.text(tx,ty,str,transform=self.axes.transAxes,
                       fontsize=16,va='center',ha='center',**kwargs)
        return

    def get(self,name):
        if name=='logx':
            print('The value of logx is',self.logx)
        elif name=='logy':
            print('The value of logy is',self.logy)
        elif name=='xtitle':
            print('The value of xtitle is',self.xtitle)
        elif name=='ytitle':
            print('The value of ytitle is',self.ytitle)
        elif name=='xlo':
            print('The value of xlo is',self.xlo)
        elif name=='xhi':
            print('The value of xhi is',self.xhi)
        elif name=='xset':
            print('The value of xset is',self.xset)
        elif name=='ylo':
            print('The value of ylo is',self.ylo)
        elif name=='yhi':
            print('The value of yhi is',self.yhi)
        elif name=='yset':
            print('The value of yset is',self.yset)
        elif name=='zlo':
            print('The value of zlo is',self.zlo)
        elif name=='zhi':
            print('The value of zhi is',self.zhi)
        elif name=='zset':
            print('The value of zset is',self.zset)
        elif name=='verbose':
            print('The value of verbose is',self.verbose)
        elif name=='cmap':
            print('The value of cmap is',self.cmap)
        else:
            print('No variable named',name)
        return

    def parse_argv(self,argv):
        if self.verbose>2:
            print('Number of arguments:', len(argv), 'arguments.')
            print('Argument List:', str(argv))
        ix=0
        while ix<len(argv):
            if self.verbose>2:
                print('Processing index',ix,'with value',argv[ix],'.')
            # Find first option, at index ix
            initial_ix_done=0
            while initial_ix_done==0:
                if ix==len(argv):
                    initial_ix_done=1
                elif argv[ix][0]=='-':
                    initial_ix_done=1
                else:
                    if self.verbose>2:
                         print('Incrementing ix')
                    ix=ix+1
            # If there is an option, then ix is its index
            if ix<len(argv):
                cmd_name=argv[ix][1:]
                if self.verbose>2:
                    print('Found option',cmd_name,'at index',ix)
                # Set ix_next to the next option, or to the end if
                # there is no next option
                ix_next=ix+1
                ix_next_done=0
                while ix_next_done==0:
                    if ix_next==len(argv):
                        ix_next_done=1
                    elif argv[ix_next][0]=='-':
                        ix_next_done=1
                    else:
                        if self.verbose>2:
                            print('Incrementing ix_next')
                        ix_next=ix_next+1
                # Now process the option
                if cmd_name=='set':
                    if self.verbose>2:
                        print('Process set.')
                    if ix_next-ix<3:
                        print('Not enough parameters for set option.')
                    else:
                        self.set(argv[ix+1],argv[ix+2])
                elif cmd_name=='get':
                    if self.verbose>2:
                        print('Process get.')
                    if ix_next-ix<2:
                        print('Not enough parameters for get option.')
                    else:
                        self.get(argv[ix+1])
                elif cmd_name=='text':
                    if self.verbose>2:
                        print('Process text.')
                    if ix_next-ix<4:
                        print('Not enough parameters for text option.')
                    else:
                        self.text(argv[ix+1],argv[ix+2],argv[ix+3])
                elif cmd_name=='read':
                    if self.verbose>2:
                        print('Process read.')
                    if ix_next-ix<2:
                        print('Not enough parameters for read option.')
                    else:
                        self.read(argv[ix+1])
                elif cmd_name=='read-name':
                    if self.verbose>2:
                        print('Process read-name.')
                    if ix_next-ix<3:
                        print('Not enough parameters for read-name option.')
                    else:
                        self.read_name(argv[ix+1],argv[ix+2])
                elif cmd_name=='read-type':
                    if self.verbose>2:
                        print('Process read-type.')
                    if ix_next-ix<3:
                        print('Not enough parameters for read-type option.')
                    else:
                        self.read_type(argv[ix+1],argv[ix+2])
                elif cmd_name=='den-plot':
                    if self.verbose>2:
                        print('Process den-plot.')
                    if ix_next-ix<2:
                        print('Not enough parameters for den-plot option.')
                    else:
                        self.den_plot(argv[ix+1])
                elif cmd_name=='read-name':
                    if self.verbose>2:
                        print('Process read.')
                    if ix_next-ix<3:
                        print('Not enough parameters for read option.')
                    else:
                        self.read_name(argv[ix+1],argv[ix+2])
                elif cmd_name=='read-type':
                    if self.verbose>2:
                        print('Process read.')
                    if ix_next-ix<3:
                        print('Not enough parameters for read option.')
                    else:
                        self.read_type(argv[ix+1],argv[ix+2])
                elif cmd_name=='xlimits':
                    if self.verbose>2:
                        print('Process xlimits.')
                    if ix_next-ix<3:
                        print('Not enough parameters for xlimits option.')
                    else:
                        self.xlimits(float(argv[ix+1]),float(argv[ix+2]))
                elif cmd_name=='ylimits':
                    if self.verbose>2:
                        print('Process ylimits.')
                    if ix_next-ix<3:
                        print('Not enough parameters for ylimits option.')
                    else:
                        self.ylimits(float(argv[ix+1]),float(argv[ix+2]))
                elif cmd_name=='plot':
                    if self.verbose>2:
                        print('Process plot.')
                    if ix_next-ix<3:
                        print('Not enough parameters for plot option.')
                    elif ix_next-ix<4:
                        self.plot(argv[ix+1],argv[ix+2])
                    else:
                        print('plot parse_argv',argv[ix+1],argv[ix+2],
                              argv[ix+3])
                        self.plot(argv[ix+1],argv[ix+2],eval(argv[ix+3]))
                elif cmd_name=='points':
                    if self.verbose>2:
                        print('Process points.')
                        print(ix,ix_next)
                    if ix_next-ix<3:
                        print('Not enough parameters for points option.')
                    elif ix_next-ix<4:
                        self.points(argv[ix+1],argv[ix+2])
                    else:
                        print('points parse_argv',argv[ix+1],argv[ix+2],
                              argv[ix+3])
                        self.points(argv[ix+1],argv[ix+2],eval(argv[ix+3]))
                elif cmd_name=='plot1':
                    if self.verbose>2:
                        print('Process plot1.')
                    if ix_next-ix<2:
                        print('Not enough parameters for plot1 option.')
                    else:
                        self.plot1(argv[ix+1])
                elif cmd_name=='plot1m':
                    if self.verbose>2:
                        print('Process plot1m.')
                    if ix_next-ix<2:
                        print('Not enough parameters for plot1m option.')
                    else:
                        files=[]
                        for i in range(ix+2,len(argv)):
                            if argv[i][0]!='-':
                                files.append(argv[i])
                            else:
                                i=len(argv)
                        self.plot1m(argv[ix+1],files)
                elif cmd_name=='plotm':
                    if self.verbose>2:
                        print('Process plotm.')
                    if ix_next-ix<3:
                        print('Not enough parameters for plot1 option.')
                    else:
                        files=[]
                        for i in range(ix+3,len(argv)):
                            if argv[i][0]!='-':
                                files.append(argv[i])
                            else:
                                i=len(argv)
                        self.plotm(argv[ix+1],argv[ix+2],files)
                elif cmd_name=='hist':
                    if self.verbose>2:
                        print('Process hist.')
                    if ix_next-ix<2:
                        print('Not enough parameters for hist option.')
                    else:
                        self.hist(argv[ix+1])
                elif cmd_name=='save':
                    if self.verbose>2:
                        print('Process save.')
                    if ix_next-ix<2:
                        print('Not enough parameters for save option.')
                    else:
                        plot.savefig(argv[ix+1])
                elif cmd_name=='line':
                    if self.verbose>2:
                        print('Process line.')
                    if ix_next-ix<5:
                        print('Not enough parameters for line option.')
                    else:
                        # Attempt to include keyword arguments.
                        # Doesn't work yet.
                        self.line(argv[ix+1],argv[ix+2],argv[ix+3],argv[ix+4])
                elif cmd_name=='list':
                    if self.verbose>2:
                        print('Process list.')
                    self.list()
                elif cmd_name=='type':
                    if self.verbose>2:
                        print('Process type.')
                    print(self.dtype)
                elif cmd_name=='move-labels':
                    if self.verbose>2:
                        print('Process move-labels.')
                    self.move_labels()
                elif cmd_name=='help':
                    if self.verbose>2:
                        print('Process help.')
                    if ix_next-ix<2:
                        self.help()
                    else:
                        self.help(argv[ix+1])
                elif cmd_name=='show':
                    if self.verbose>2:
                        print('Process show.')
                    self.show()
                elif cmd_name=='canvas':
                    if self.verbose>2:
                        print('Process canvas.')
                    self.canvas()
                else:
                    print('No option named',cmd_name)
                # Increment to the next option
                ix=ix_next
            if self.verbose>2:
                print('Going to next.')
        return
