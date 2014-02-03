"""
  -------------------------------------------------------------------

  Copyright (C) 2006-2013, Andrew W. Steiner
  
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

import matplotlib.pyplot as plot
import h5py 
import math
from matplotlib.colors import LinearSegmentedColormap

list_of_dsets=[]
search_type=''

def hdf5_is_object_type(name,obj):
    if isinstance(obj,h5py.Group):
         if 'o2scl_type' in obj.keys():
            o2scl_type_dset=obj['o2scl_type']
            if o2scl_type_dset.__getitem__(0) == search_type:
                list_of_dsets.append(name)

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

def h5read_first_type(fname,loc_type):
    del list_of_dsets[:]
    global search_type
    search_type=loc_type
    file=h5py.File(fname,'r')
    file.visititems(hdf5_is_object_type)
    if len(list_of_dsets)==0:
        str='Could not object of type '+loc_type+' in file '+fname+'.'
        raise RuntimeError(str)
    return file[list_of_dsets[0]]

def h5read_name(fname,name):
    file=h5py.File(fname,'r')
    obj=file[name]
    o2scl_type_dset=obj['o2scl_type']
    loc_type=o2scl_type_dset.__getitem__(0)
    return (obj,loc_type)

def h5read_type_named(fname,loc_type,name):
    del list_of_dsets[:]
    global search_type
    search_type=loc_type
    file=h5py.File(fname,'r')
    file.visititems(hdf5_is_object_type)
    if name in list_of_dsets:
        return file[name]
    str='No object of type '+loc_type+' named '+name+' in file '+fname+'.'
    raise RuntimeError(str)
    return

def default_plot(lmar=0.14,bmar=0.12,rmar=0.05,tmar=0.05):
    """
    Works on mac, but not on riddler...
    plot.rc('text',usetex=True)
    """
    plot.rc('font',family='serif')
    plot.rcParams['lines.linewidth']=0.5
    fig=plot.figure(1,figsize=(6.4,6.4))
    fig.set_facecolor('white')
    ax=plot.axes([lmar,bmar,1.0-lmar-rmar,1.0-tmar-bmar])
    ax.minorticks_on()
    ax.tick_params('both',length=12,width=1,which='major')
    ax.tick_params('both',length=5,width=1,which='minor')
    plot.grid(False)
    return ax
    
def get_str_array(dset):
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
        col=col+str(unichr(data[ix]))
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
    
class plotter:
    
    logx=0
    logy=0
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
    canvas_flag=0
    dtype=''
    cmap='jet'

    def myreds(self):
        cdict={'red': ((0.0,1.0,1.0),(1.0,1.0,1.0)),
               'green': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'blue': ((0.0,1.0,1.0),(1.0,0.0,0.0))}
        myreds=LinearSegmentedColormap('MyReds',cdict)
        plot.register_cmap(cmap=myreds)
        cmap='MyReds'

    def mygreens(self):
        cdict={'red': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'green': ((0.0,1.0,1.0),(1.0,1.0,1.0)),
               'blue': ((0.0,1.0,1.0),(1.0,0.0,0.0))}
        mygreens=LinearSegmentedColormap('MyGreens',cdict)
        plot.register_cmap(cmap=mygreens)
        cmap='MyGreens'

    def myblues(self):
        cdict={'red': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'green': ((0.0,1.0,1.0),(1.0,0.0,0.0)),
               'blue': ((0.0,1.0,1.0),(1.0,1.0,1.0))}
        myblues=LinearSegmentedColormap('MyBlues',cdict)
        plot.register_cmap(cmap=myblues)
        cmap='MyBlues'
        
    def contour_plot(self,level,**kwargs):
        if self.dtype!='vector<contour_line>':
            print 'Wrong type for contour_plot.'
            return
        if self.verbose>2:
            print 'contour_plot',level,kwargs
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
        n_lines=self.dset['n_lines'][0]
        print n_lines,'lines.'
        for i in range(0,n_lines):
            line_level=self.dset['line_'+str(i)+'/level'][0]
            print 'level for line',i,' is ',line_level
            if abs(level-line_level) < 1.0e-10:
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
 
    def hist_plot(self,**kwargs):
        if self.dtype!='hist':
            print 'Wrong type for hist_plot.'
            return
        if self.verbose>2:
            print 'hist_plot',kwargs
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
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
        return
        

    def plot(self,colx,coly,**kwargs):
        if self.verbose>2:
            print 'plot',colx,coly,kwargs
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
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
        return

    def plot1(self,col,**kwargs):
        if self.verbose>2:
            print 'plot1',col,kwargs
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
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
        return

    def hist(self,col,**kwargs):
        if self.verbose>2:
            print 'hist',col,kwargs
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
        plot.hist(self.dset['data/'+col],**kwargs)
        return

    def hist2d(self,colx,coly,**kwargs):
        if self.verbose>2:
            print 'hist',colx,coly,kwargs
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
        plot.hist2d(self.dset['data/'+colx],self.dset['data/'+coly],**kwargs)
        return

    def line(self,x1,y1,x2,y2,**kwargs):
        if self.verbose>2:
            print 'Line',x1,y1,x2,y1
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
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
            print 'Canvas'
        # Default o2mpl plot
        self.axes=default_plot()
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

    def read(self,filename):
        self.dset=h5read_first_type(filename,'table')
        self.dtype='table'
        return

    def read_type(self,filename,loc_type):
        self.dset=h5read_first_type(filename,loc_type)
        self.dtype=loc_type
        return

    def read_name(self,filename,name):
        atuple=h5read_name(filename,name)
        self.dset=atuple[0]
        self.dtype=atuple[1]
        return

    def list(self):
        if self.dtype=='table':
            col_list=get_str_array(self.dset['col_names'])
            if self.verbose>2:
                print '-----------------------'
            unit_list=[]
            unit_flag=self.dset['unit_flag'][0]
            print 'unit_flag ',unit_flag
            if self.verbose>2:
                print 'unit_flag:',unit_flag
            if unit_flag==1:
                unit_list=get_str_array(self.dset['units'])
                if self.verbose>2:
                    print '-----------------------'
                    print 'unit_list:',unit_list
            print len(col_list),'columns.'
            for ix in range(0,len(col_list)):
                if unit_flag:
                    print str(ix)+'. '+col_list[ix]+' ['+unit_list[ix]+']'
                else:
                    print str(ix)+'. '+col_list[ix]
            print self.dset['nlines'][0],'lines.'
            if self.verbose>2:
                print 'Done in list'
        elif self.dtype=='table3d':
            sl_list=get_str_array(self.dset['slice_names'])
            print len(sl_list),'slices.'
            for ix in range(0,len(sl_list)):
                print str(ix)+'. '+sl_list[ix]
            xgrid=self.dset['xval'].value
            ygrid=self.dset['yval'].value
            lxgrid=len(xgrid)
            lygrid=len(ygrid)
            print('X-grid start: '+str(xgrid[0])+' end: '+
                  str(xgrid[lxgrid-1])+' size '+str(lxgrid))
            print('Y-grid start: '+str(ygrid[0])+' end: '+
                  str(ygrid[lygrid-1])+' size '+str(lygrid))
        else:
            print 'Cannot list type',self.dtype
        return

    def den_plot(self,slice_name,**kwargs):
        if self.dtype=='table3d':
            name='data/'+slice_name
            sl=self.dset[name].value
            sl=sl.transpose()
            xgrid=self.dset['xval'].value
            ygrid=self.dset['yval'].value
            if self.canvas_flag==0:
                self.canvas()
                self.canvas_flag=1
            if self.logx==1:
                for i in range(0,len(xgrid)):
                    xgrid[i]=math.log(xgrid[i],10)
            if self.logy==1:
                for i in range(0,len(ygrid)):
                    ygrid[i]=math.log(ygrid[i],10)
            plot.imshow(sl,cmap=self.cmap,interpolation='nearest',
                        origin='lower',
                        extent=[xgrid[0],xgrid[len(xgrid)-1],ygrid[0],
                                ygrid[len(ygrid)-1]],aspect='auto',**kwargs)
        else:
            print 'Cannot density plot object of type',self.dtype
        return

    def set(self,name,value):
        if self.verbose>1:
            print 'Set',name,'to',value
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
        elif name=='xhi':
            self.xhi=float(value)
        elif name=='xset':
            self.xset=int(value)
        elif name=='ylo':
            self.ylo=float(value)
        elif name=='yhi':
            self.yhi=float(value)
        elif name=='yset':
            self.yset=int(value)
        elif name=='zlo':
            self.zlo=float(value)
        elif name=='zhi':
            self.zhi=float(value)
        elif name=='zset':
            self.zset=int(value)
        elif name=='verbose':
            self.verbose=int(value)
        else:
            print 'No variable named',name
        return

    def help(self,arg=''):
        if arg=='canvas':
            print 'canvas()'
        elif arg=='plot':
            print '---------------------------------------------------------'
            print 'plot(x,y,**kwargs)'
            print '---------------------------------------------------------'
            print 'Useful kwargs:'
            print 'color or c                   matplotlib color'
            temp=str('dashes                       '+
                     'sequence of on/off ink in points')
            print temp
            temp=str("fillstyle                    "+
                     "['full' | 'left' | 'right' | 'bottom' | 'top'")
            print temp
            print "                              | 'none']"
            temp=str("label                        "+
                     "string or anything printable with '%s' conversion.")
            print temp
            temp=str("linestyle or ls              "+
                     "['-' | '--' | '-.' | ':' | 'None' | ' ' | '']")
            print temp
            temp=str("                              "+
                     "and any drawstyle in combination with a")
            print temp
            temp=str("                              "+
                     "linestyle, e.g., 'steps--'.")
            print temp
            print "linewidth or lw              float value in points"
            print 'marker                       marker type'
            print 'markeredgecolor or mec       matplotlib color'
            print 'markeredgewidth or mew       float value in points'
            print 'markerfacecolor or mfc       matplotlib color'
            print 'markerfacecoloralt or mfcalt matplotlib color'
            print 'markersize or ms             float'
            print '---------------------------------------------------------'
        elif arg=='den_plot':
            print 'den_plot(slice_name)'
        elif arg=='get':
            print 'get(name)'
        elif arg=='line':
            print 'line(x1,y1,x2,y2,color,lstyle)'
        elif arg=='list':
            print 'list()'
        elif arg=='move_labels':
            print 'move_labels()'
        elif arg=='read':
            print 'read(filename)'
        elif arg=='read_name':
            print 'read_name(filename,name)'
        elif arg=='read_type':
            print 'read_type(filename,type)'
        elif arg=='reset_xlimits':
            print 'reset_xlimits()'
        elif arg=='reset_ylimits':
            print 'reset_ylimits()'
        elif arg=='set':
            print 'set(name,value)'
        elif arg=='show':
            print 'show()'
        elif arg=='text':
            print 'text(string,x,y,color)'
        elif arg=='xlimits':
            print 'xlimits(xlow,xhigh)'
        elif arg=='ylimits':
            print 'ylimits(ylow,yhigh)'
        else:
            print 'Methods:'
            print 'canvas()'
            print 'contour_plot(level,**kwargs)'
            print 'den_plot(slice_name,**kwargs)'
            print 'get(name)'
            print 'hist(col,tbins,tcolor)'
            print 'hist_plot(**kwargs)'
            print 'hist2d(colx,coly,**kwargs)'
            print 'line(x1,y1,x2,y2,**kwargs)'
            print 'list()'
            print 'move_labels()'
            print 'parse_argv(argv)'
            print 'plot(colx,coly,**kwargs)'
            print 'read(filename)'
            print 'read_name(filename,name)'
            print 'read_type(filename,type)'
            print 'reds()'
            print 'reset_xlimits()'
            print 'reset_ylimits()'
            print 'set(name,value)'
            print 'show()'
            print 'text(string,x,y,**kwargs)'
            print 'xlimits(xlow,xhigh)'
            print 'ylimits(ylow,yhigh)'
        return

    def text(self,str,tx,ty,**kwargs):
        if self.canvas_flag==0:
            self.canvas()
            self.canvas_flag=1
        self.axes.text(tx,ty,str,transform=self.axes.transAxes,
                       fontsize=16,verticalalignment='top',**kwargs)
        return

    def get(self,name):
        if name=='logx':
            print 'logx is',self.logx
        elif name=='logy':
            print 'logy is',self.logy
        elif name=='xtitle':
            print 'xtitle is',self.xtitle
        elif name=='ytitle':
            print 'ytitle is',self.ytitle
        elif name=='xlo':
            print 'xlo is',self.xlo
        elif name=='xhi':
            print 'xhi is',self.xhi
        elif name=='xset':
            print 'xset is',self.xset
        elif name=='ylo':
            print 'ylo is',self.ylo
        elif name=='yhi':
            print 'yhi is',self.yhi
        elif name=='yset':
            print 'yset is',self.yset
        elif name=='zlo':
            print 'zlo is',self.zlo
        elif name=='zhi':
            print 'zhi is',self.zhi
        elif name=='zset':
            print 'zset is',self.zset
        elif name=='verbose':
            print 'verbose is',self.verbose
        else:
            print 'No variable named',name
        return

    def parse_argv(self,argv):
        if self.verbose>2:
            print 'Number of arguments:', len(argv), 'arguments.'
            print 'Argument List:', str(argv)
        ix=0
        while ix<len(argv):
            if self.verbose>2:
                print 'Processing index',ix,'with value',argv[ix],'.'
            # Find first option, at index ix
            initial_ix_done=0
            while initial_ix_done==0:
                if ix==len(argv):
                    initial_ix_done=1
                elif argv[ix][0]=='-':
                    initial_ix_done=1
                else:
                    if self.verbose>2:
                         print 'Incrementing ix'
                    ix=ix+1
            # If there is an option, then ix is its index
            if ix<len(argv):
                cmd_name=argv[ix][1:]
                if self.verbose>2:
                    print 'Found option',cmd_name,'at index',ix
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
                            print 'Incrementing ix_next'
                        ix_next=ix_next+1
                # Now process the option
                if cmd_name=='set':
                    if self.verbose>2:
                        print 'Process set.'
                    if ix_next-ix<3:
                        print 'Not enough parameters for set option.'
                    else:
                        set(argv[ix+1],argv[ix+2])
                elif cmd_name=='get':
                    if self.verbose>2:
                        print 'Process get.'
                    if ix_next-ix<2:
                        print 'Not enough parameters for get option.'
                    else:
                        get(argv[ix+1])
                elif cmd_name=='text':
                    if self.verbose>2:
                        print 'Process text.'
                    if ix_next-ix<4:
                        print 'Not enough parameters for text option.'
                    else:
                        text(argv[ix+1],argv[ix+2],argv[ix+3])
                elif cmd_name=='read':
                    if self.verbose>2:
                        print 'Process read.'
                    if ix_next-ix<2:
                        print 'Not enough parameters for read option.'
                    else:
                        self.read(argv[ix+1])
                elif cmd_name=='den-plot':
                    if self.verbose>2:
                        print 'Process den-plot.'
                    if ix_next-ix<2:
                        print 'Not enough parameters for den-plot option.'
                    else:
                        self.den_plot(argv[ix+1])
                elif cmd_name=='read-name':
                    if self.verbose>2:
                        print 'Process read.'
                    if ix_next-ix<3:
                        print 'Not enough parameters for read option.'
                    else:
                        self.read_name(argv[ix+1],argv[ix+2])
                elif cmd_name=='read-type':
                    if self.verbose>2:
                        print 'Process read.'
                    if ix_next-ix<3:
                        print 'Not enough parameters for read option.'
                    else:
                        self.read_type(argv[ix+1],argv[ix+2])
                elif cmd_name=='plot':
                    if self.verbose>2:
                        print 'Process plot.'
                    if ix_next-ix<3:
                        print 'Not enough parameters for plot option.'
                    else:
                        self.line(argv[ix+1],argv[ix+2])
                elif cmd_name=='line':
                    if self.verbose>2:
                        print 'Process line.'
                    if ix_next-ix<5:
                        print 'Not enough parameters for line option.'
                    else:
                        self.line(argv[ix+1],argv[ix+2],argv[ix+3],argv[ix+4])
                elif cmd_name=='list':
                    if self.verbose>2:
                        print 'Process list.'
                    self.list()
                elif cmd_name=='move-labels':
                    if self.verbose>2:
                        print 'Process move-labels.'
                    self.move_labels()
                elif cmd_name=='canvas':
                    if self.verbose>2:
                        print 'Process canvas.'
                    self.canvas()
                else:
                    print 'No option named',cmd_name
                # Increment to the next option
                ix=ix_next
            if self.verbose>2:
                print 'Going to next.'
        return
