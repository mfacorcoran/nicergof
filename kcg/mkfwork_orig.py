#
import sys
import numpy as np
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table
from astropy.io import fits

#mkffile=sys.argv[1]
#print("MKF file = ", mkffile)
#fits.info(mkffile)
mkffile = '/Volumes/SXDC/Data/NICER/kp_model/h1743322kp.mkf'
mkf = Table.read(mkffile, hdu=1)

#print(mkf.columns)

#print(mkf['COR_SAX'].unit)
#print(mkf['KP'].unit)



###################################################################################
#
# Read Library
#
###################################################################################
#librarydirectory="/export/home/kgendrea/work/SPACE_WEATHER/BKGROUND_DATA/gen10/"
librarydirectory = '/Volumes/SXDC/Data/NICER/kp_model/gen10/'
kpmax=7
nn=0
lib=[]
flist=[]
for c in range(15):
    c1=c
    c2=c1+1
    for k in range(7):
        k1=k-0.5
        k2=k1+1
        file=librarydirectory+"cor"+str(c1)+"_"+str(c2)+"_kp"+str(k1)+"_"+str(k2)+".normalized"
        flist.append(file)
#        print file
        data=[]
        f=open(file,'r')
        for line in f:
            line = line.strip()
            columns = line.split()
            data.append(float(columns[1]))
            libe=[c1,c2,k1,k2,data]
#            lib.append(libe)
        lib.append(data)
#        print("A", nn, c,k1, data[90])
        nn=nn+1
#        print("A", data[77])


print("HEY HEY", lib[9][77], lib[0][77])
#print(len(lib))


###################################################################################
#
# Build background spectra
#
###################################################################################
bkgd=np.zeros(1501)
BKGD_EST_PI_0000_0025=[]
BKGD_EST_PI_0035_0200=[]
BKGD_EST_PI_0200_0800=[]
BKGD_EST_PI_0800_1200=[]
BKGD_EST_PI_1200_1500=[]

i=0
length=len(mkf['COR_SAX'])
#print("Length = ", length)
while i < length:
#    print(mkf['COR_SAX'][i],i,length)
# note kludge to handle unlikely cases of too high a COR or KP.. we would normally throw out too high a KP,
    cor=min(mkf['COR_SAX'][i],13.99)
    kp=min(mkf['KP'][i],5.9)
    kpi=round(kp)
    corin=round(cor-0.5)
    indexa=int(corin*kpmax+kpi)
    x=cor-corin-0.5
    indexb=int((corin+np.sign(x))*kpmax+kpi)
    # lib is an array of length 105; each element is a list containing the normalized spectra (1501 channels)
    # example:
    # lib[0] is the spectrum cor0_1_kp-0.5_0.5.normalized
    # lib[1] is the spectrum cor0_1_kp0.5_1.5.normalized etc
    vala=lib[indexa]
    valb=lib[indexb]
    valc=np.array(vala)+(np.array(valb)-np.array(vala))*x*np.sign(x)
    bkgd=np.array(bkgd)+np.array(valc)
# nomalize to number of FPMs used in background field library
    B_PI_0000_0025=sum(valc[0:25])/49
    B_PI_0035_0200=sum(valc[35:200])/49
    B_PI_0200_0800=sum(valc[200:800])/49
    B_PI_0800_1200=sum(valc[800:1200])/49
    B_PI_1200_1500=sum(valc[1200:1500])/49
    BKGD_EST_PI_0000_0025.append(B_PI_0000_0025)
    BKGD_EST_PI_0035_0200.append(B_PI_0035_0200)
    BKGD_EST_PI_0200_0800.append(B_PI_0200_0800)
    BKGD_EST_PI_0800_1200.append(B_PI_0800_1200)
    BKGD_EST_PI_1200_1500.append(B_PI_1200_1500)

#    print("B", indexa, cor,kpi, vala[90])
#    print(cor,corin,kp,indexa,indexb)
    print mkf['TIME'][i], mkf['COR_SAX'][i], mkf['KP'][i], \
        mkf['FPM_XRAY_PI_0035_0200'][i]+mkf['FPM_XRAY_PI_0200_0800'][i],  \
        B_PI_0035_0200+B_PI_0200_0800
    i += 1
#######################################################
#
# now write the estimated background rates to the mkf file
#
#########################################################3





#########################################################3
#
#Write ascii file of built up estimated background spectrum... note- may need to deadtime correct
#
#################################################3
file="predicted_background.qdp"
f=open(file,'wb')
for i in range(1501):
    output=str(i)+" "+str(bkgd[i])+"\n"
    f.write(output)
f.close()


#print("Length = ", length)
#print("XX ", bkgd[30:50])
