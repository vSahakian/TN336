#Batchfile to print MST data for descriptions
#VJS 1/2016
#In terminal, copy all massMS.out files into a massMS directory, like:
#cp TN*/*massMS.out massMSout/.

import numpy as np
import MSTtools as mt
from glob import glob
from os import path


##*****CHANGE ME*****##
#Make sure these directories exist first!  
#Make sure there is a slash at the end of the directory string!

##Directory from which you are taking the massMS.out files:
gdir='/Users/vjsahakian/Documents/Cruise/TN336/MST/massMSout/'
##Directory where you want the pdfs:
pdir='/Users/vjsahakian/Documents/Cruise/TN336/MST/pdfs/'

##*****CONTINUE...*****##



#Make the globfile
gfile=glob(gdir+'*.out')

#Loop over files in here to plot:
for i in range(len(gfile)):
    #Get info for file names:
    cruiseName=path.split(gfile[i])[1].split('_')[0]
    coreName=path.split(gfile[i])[1].split('_')[1]
    
    #INput and output files:
    mstFile=gfile[i]
    pdfFile=pdir+cruiseName+coreName+'.pdf'
    
    #Plot:
    mt.pltMST(mstFile,pdfFile)
    print cruiseName+coreName
    
    
    