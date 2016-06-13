#Plot data from MST processed excel spreadsheet
#VJS 1/2016
    
######
#Plot from processed data
#VJS 1/2016
def rdMSTpcsv(procfile):
    '''
    Read in data from an MST processed csv file, exported from excel 
    (first sheet only)
    Input:
        procfile:     String of path to csv Proc file (exported from excel)
    Output:  ** These will contain nans**
        sbdepth:      Array with the core sub bottom depths 
        pvel:         Array with the p-wave velocity
        rho:          Array with the density
        msu:          Array with the magnetic susceptibility
        units:        Array with units for the respective output variables
        fileinfo:     String with information about datafile
    '''
    
    import numpy as np
    
    #First, read in units:
    mstp=open(procfile,'r')
    #Read in information about file and date of creation
    mstpdat=mstp.readline()
    #It's all in one line if it's a csv, so split it:
    fileinfo=mstpdat.split('\r')[0].split(',')[0]
    udat=mstpdat.split('\r')[4].split(',')[11:23]    
    #Close file
    mstp.close()
    
    #Make units array:
    units=[udat[0],udat[5],udat[6],udat[8]]
    
    #Read in data:
    dat=np.genfromtxt(procfile,delimiter=',',skip_header=6)
    #Assign to variables:
    sbdepth=dat[:,11]
    pvel=dat[:,16]
    rho=dat[:,17]
    msu=dat[:,19]


    return sbdepth,pvel,rho,msu,units,fileinfo
    
##############################
def rdax(axfile):
    #VJS 1/2016
    '''
    Read a text file with axis limits
    Input:
        axfile:     String to path of the input axfile.  
                    Format:
                        vmin=#
                        vmax=#
                        rmin=#
                        rmax=#
                        msmin=#
                        msmax=#
    Output:
        Returns pax: pax=np.array([[vmin,vmax],[rmin,rmax],[msmin,msmax]])
    '''
    
    import numpy as np
    
    #Open ax file
    ax=open(axfile,'r')
    
    #Read in vmin:
    vminl=ax.readline()
    vmaxl=ax.readline()
    rminl=ax.readline()
    rmaxl=ax.readline()
    msminl=ax.readline()
    msmaxl=ax.readline()
    #Close
    ax.close()
    
    #Convert:
    vmin=np.float(vminl.split('=')[1].split('\n')[0])
    vmax=np.float(vmaxl.split('=')[1].split('\n')[0])
    rmin=np.float(rminl.split('=')[1].split('\n')[0])
    rmax=np.float(rmaxl.split('=')[1].split('\n')[0])
    msmin=np.float(msminl.split('=')[1].split('\n')[0])
    msmax=np.float(msmaxl.split('=')[1].split('\n')[0])
    
    pax=np.array([[vmin,vmax],[rmin,rmax],[msmin,msmax]])
    
    return pax
    
    
##############################
def pltMST(mstfile,pdffile):
    '''
    Plot data from an output, processed, MST file 
    Data to plot: P-wave speed, magnetic susceptibility, and density
    Input:
        pdffile:    String with output pdf file
        pltparams:  
    Output: 
        Prints a pdf to pdffile
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    
    #Read units from MST file:
    mstread=open(mstfile,'r')
    #Read in information about file and date of creation
    fileinfo=mstread.readline()
    
    #Skip blank and processing parameters lines:
    for i in range(14):
        mstread.readline()
    
    #Read in the units, now on row 16:
    units_dos=mstread.readline()
    #remove the tabs, and the newline entry at the end of the line:
    units=units_dos.split('\t')
    units[-1]=units[-1].split('\r')[0]
    
    #Close file
    mstread.close()
    ######
    
    #Get title/core name:
    corename=fileinfo.split(' ')[4].split('_')[1]

    
    #Read in data - it's actually processes, but here call it raw
    rdat=np.genfromtxt(mstfile,delimiter='\t',skip_header=17)
    
    #Remove calibration pieces:
    cal_inds=np.where(np.isnan(rdat[:,1]))
    cal_end=np.max(np.max(cal_inds))+1
    dat=rdat[cal_end:-1,:]
    
    cdep=dat[:,0]
    pvel=dat[:,5]
    rho=dat[:,6]
    ms=dat[:,7]
    
    #Find indices that are NOT nan...can plot these:
    pind=np.where(np.isfinite(pvel))
    rind=np.where(np.isfinite(rho))
    mind=np.where(np.isfinite(ms))
    
    #Get axes:
    #For y, always plot 0 to ymax:
    ymin=np.max(cdep)
    ymax=0
    ytick=np.arange(ymin,ymax,0.2)
    
    #For pwave velocity:
    vstd=np.std(pvel)
    vm=np.mean(pvel)
    vxmin=vm-1*vstd
    vxmax=vm+1*vstd
    #Get ticks:
    vtmin=vxmin-np.remainder(vxmin,100)
    vtmax=vxmax-np.remainder(vxmax,100)
    
    vax=[vtmin,vtmax,ymin,ymax]
    vtick=np.arange(vtmin,vtmax+100,300)
    
    #Density:
    #Std dev:
    rstd=np.std(rho)
    rm=np.mean(rho)
    rxmin=np.round(rm-3.5*rstd,1)
    rxmax=np.round(rm+3.5*rstd,1)
    rax=[rxmin,rxmax,ymin,ymax]
    #Get ticks
    rtick=np.arange(rxmin,rxmax,0.3)

    
    #Magnetic susceptibility:
    msstd=np.std(ms)
    msm=np.mean(ms)
    msxmin=np.round(msm-3.5*msstd,1)
    msxmax=np.round(msm+3.5*msstd,1)
    msax=[msxmin,msxmax,ymin,ymax]
    #Get ticks
    mtick=np.arange(msxmin,msxmax,1)
    
    ####
    #INitiate Plotting!!
    f=plt.figure(num=1,figsize=(7,10))
    
    #Plot P-wave velocity
    plt.subplot(131)
    plt.scatter(pvel[pind],cdep[pind],verts='inverse',color='blue')
    plt.plot(pvel[pind],cdep[pind],linewidth=1,color='blue')
    plt.axis(vax)
    plt.xticks(vtick)
    #plt.yticks(ytick)
    plt.grid(which='major')
    plt.ylabel('Core depth ('+units[0]+')')
    plt.xlabel('('+units[5]+')')
    plt.title('Pwave vel')
    
    #Density
    plt.subplot(132)
    plt.scatter(rho[rind],cdep[rind],color='red')
    plt.plot(rho[rind],cdep[rind],linewidth=1,color='red')
    plt.axis(rax)
    plt.xticks(rtick)
    plt.grid(which='major')
    plt.xlabel('('+units[6]+')')
    plt.title('Density')
    
    #Magnetic susceptibility
    plt.subplot(133)
    plt.scatter(ms[mind],cdep[mind],color='green')
    plt.plot(ms[mind],cdep[mind],linewidth=1,color='green')
    plt.axis(msax)
    plt.xticks(mtick)
    plt.grid(which='major')
    plt.xlabel('('+units[7]+')')
    plt.title('Magnetic susc.')
    
    #Label
    plt.suptitle('Core '+corename,fontsize=14)
    
    #Save figure
    plt.savefig(pdffile)
    plt.clf()    
    
    
#########################
#VJs 1/21/16
def pltProc(procfile,pdffile,pax):
    '''
    Plot data from a Processed csv file (from excel)
    Input:
        procfile:   String with path to the procfile
        pdffile:    STring with path to the pdffile
        pax:         Array with axes:
                    [[vmin,vmax],[rmin,rmax],[msumin,msumax]]  
    Output:
        Prints pdf to pdffile
    '''
    
    import numpy as np
    import matplotlib.pyplot as plt
    import MSTtools as mst
    
    #Get data:
    cdep,pvel,rho,msu,units,fileinfo=mst.rdMSTpcsv(procfile)
    
    #Get cruise name and corename:
    cruisename=fileinfo.split(' ')[4].split('_')[0]
    corename=fileinfo.split(' ')[4].split('_')[1].split('.')[0]
    
    #Find indices that are NOT nan...can plot these:
    pind=np.where(np.isfinite(pvel))
    rind=np.where(np.isfinite(rho))
    mind=np.where(np.isfinite(msu))
    
    #Get axes for plotting:
    ymax=0
    cind=np.where(np.isfinite(cdep))
    ymin=np.max(cdep[cind])
    
    #Velocity
    vmin=pax[0,0]
    vmax=pax[0,1]
    vax=[vmin,vmax,ymin,ymax]
    step=(vmax-vmin)/4
    vtick=np.around(np.linspace(vmin,vmax,4),0)
        
    #Rho
    rmin=pax[1,0]
    rmax=pax[1,1]
    rax=[rmin,rmax,ymin,ymax]
    rtick=np.around(np.linspace(rmin,rmax,4),1)
    
    ##mass mag sus
    msmin=pax[2,0]
    msmax=pax[2,1]
    msax=[msmin,msmax,ymin,ymax]
    mstick=np.around(np.linspace(msmin,msmax,4),1)   
    
        
    
    #Plot:
    ####
    #INitiate Plotting!!
    plt.clf()
    f=plt.figure(num=1,figsize=(7,10))
    
    #Plot P-wave velocity
    plt.subplot(131)
    plt.scatter(pvel[pind],cdep[pind],verts='inverse',color='blue')
    plt.plot(pvel[pind],cdep[pind],linewidth=1,color='blue')
    plt.axis(vax)
    plt.xticks(vtick)
    #plt.yticks(ytick)
    plt.grid(which='major')
    plt.ylabel('Core depth ('+units[0]+')')
    plt.xlabel('('+units[1]+')')
    plt.title('Pwave vel')
    
    #Density
    plt.subplot(132)
    plt.scatter(rho[rind],cdep[rind],color='red')
    plt.plot(rho[rind],cdep[rind],linewidth=1,color='red')
    plt.axis(rax)
    plt.xticks(rtick)
    plt.grid(which='major')
    plt.xlabel('('+units[2]+')')
    plt.title('Density')
    
    #Magnetic susceptibility
    plt.subplot(133)
    plt.scatter(msu[mind],cdep[mind],color='green')
    plt.plot(msu[mind],cdep[mind],linewidth=1,color='green')
    plt.axis(msax)
    plt.xticks(mstick)
    plt.grid(which='major')
    plt.xlabel('('+units[3]+')')
    plt.title('Magnetic susc.')
    
    #Label
    plt.suptitle('Core '+corename,fontsize=14)
    
    #Save figure
    plt.savefig(pdffile)
    plt.clf()
    
   
        
#####################
#Plot trigger core and piston core together
#Vjs 1/2016
def pltJpcTc(jpcfile,tcfile,pdffile,pax):
    '''
    Plot a trigger core and piston core on the same plot
    Input:
        jpcfile:     String for path of jpc csv file
        tcfile:      String for path of tc csv file
        pdffile:     STring for path of output pdffile
        pax:         Array with axes:
                    [[vmin,vmax],[rmin,rmax],[msumin,msumax]]  
    Output:
        Prints to pdffile with tri plot: P-wave velocity, rho, mag sus.
    '''
            
    import numpy as np
    import matplotlib.pyplot as plt
    import MSTtools as mst
    
    #Get data:
    cdepJ,pvelJ,rhoJ,msuJ,unitsJ,fileinfoJ=mst.rdMSTpcsv(jpcfile)
    cdepT,pvelT,rhoT,msuT,unitsT,fileinfoT=mst.rdMSTpcsv(tcfile)
    
    #Get cruise name and corename:
    cruisename=fileinfoJ.split(' ')[4].split('_')[0]
    corename=fileinfoJ.split(' ')[4].split('_')[1].split('.')[0]
    
    #Find indices that are NOT nan...can plot these:
    pindJ=np.where(np.isfinite(pvelJ))
    rindJ=np.where(np.isfinite(rhoJ))
    mindJ=np.where(np.isfinite(msuJ))
    pindT=np.where(np.isfinite(pvelT))
    rindT=np.where(np.isfinite(rhoT))
    mindT=np.where(np.isfinite(msuT))        
        
        
    #Get axes for plotting:
    ymax=0
    ymin=max(np.max(cdepJ),np.max(cdepT))
    
    #Velocity
    vmin=pax[0,0]
    vmax=pax[0,1]
    vax=[vmin,vmax,ymin,ymax]
    step=(vmax-vmin)/4
    vtick=np.around(np.linspace(vmin,vmax,4),0)
        
    #Rho
    rmin=pax[1,0]
    rmax=pax[1,1]
    rax=[rmin,rmax,ymin,ymax]
    rtick=np.around(np.linspace(rmin,rmax,4),1)
    
    ##mass mag sus
    msmin=pax[2,0]
    msmax=pax[2,1]
    msax=[msmin,msmax,ymin,ymax]
    mstick=np.around(np.linspace(msmin,msmax,4),1)        
    
    #Plot:
    ####
    #INitiate Plotting!!
    plt.figure(num=1,figsize=(10,7))
    
    #Plot P-wave velocity
    plt.subplot(131)
    #JPC
    plt.scatter(pvelJ[pindJ],cdepJ[pindJ],verts='inverse',color='blue')
    plt.plot(pvelJ[pindJ],cdepJ[pindJ],linewidth=1,color='blue')
    #TC
    plt.scatter(pvelT[pindT],cdepT[pindT],verts='inverse',color='red')
    plt.plot(pvelT[pindT],cdepT[pindT],linewidth=1,color='red')
    #Axes
    plt.axis(vax)
    plt.xticks(vtick)
    plt.grid(which='major')
    plt.ylabel('Core depth ('+unitsJ[0]+')')
    plt.xlabel('('+unitsJ[1]+')')
    plt.title('Pwave vel')
    
    #Density
    plt.subplot(132)
    #JPC
    plt.scatter(rhoJ[rindJ],cdepJ[rindJ],verts='inverse',color='blue')
    plt.plot(rhoJ[pindJ],cdepJ[pindJ],linewidth=1,color='blue')
    #TC
    plt.scatter(rhoT[rindT],cdepT[rindT],verts='inverse',color='red')
    plt.plot(rhoT[rindT],cdepT[rindT],linewidth=1,color='red')
    #Axes
    plt.axis(rax)
    plt.xticks(rtick)
    plt.grid(which='major')
    plt.xlabel('('+unitsJ[2]+')')
    plt.title('Density')
    
    #Magnetic susceptibility
    plt.subplot(133)
    #JPC
    plt.scatter(msuJ[mindJ],cdepJ[mindJ],verts='inverse',color='blue',label='JPC')
    plt.plot(msuJ[mindJ],cdepJ[mindJ],linewidth=1,color='blue')
    #TC
    plt.scatter(msuT[mindT],cdepT[mindT],verts='inverse',color='red',label='TC')
    plt.plot(msuT[mindT],cdepT[mindT],linewidth=1,color='red')
    #Axes
    plt.axis(msax)
    plt.xticks(mstick)
    plt.grid(which='major')
    plt.xlabel('('+unitsJ[3]+')')
    plt.title('Magnetic susc.')
    
    #Label
    plt.suptitle('Core TC/'+corename,fontsize=14)
    plt.legend(bbox_to_anchor=(0.99,0.2),loc='upper right')
    
    #Save figure
    plt.savefig(pdffile)
    plt.clf()
    
#Get axes automatically:
#VJS 1/2016
      

