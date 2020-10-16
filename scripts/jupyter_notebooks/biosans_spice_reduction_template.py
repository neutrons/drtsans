
# coding: utf-8

# # Quick q reduction

# Import Libraries

# In[2]:


from mantid.simpleapi import *
import numpy as np
from drtsans.mono import gpsans as sans 
from drtsans.mono.gpsans import attenuation_factor
from drtsans.mono.absolute_units import empty_beam_scaling
import drtsans
from drtsans.mask_utils import apply_mask, load_mask
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D  # noqa E402
from drtsans.mono.gpsans import convert_to_q
from drtsans.iq import BinningMethod, BinningParams, determine_1d_linear_bins  # noqa E402
from drtsans.api import subtract_background
from drtsans.iq import bin_all
from drtsans.plots import plot_IQmod, plot_IQazimuthal
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages


# Set up bkgd, center, flood data to be used.
# Note:  This notebook only works on one configuration at a time

# In[3]:


# Empty cell changer 
ipts=26004
Bkgd=[]#low_q 12102 #midq 12099
centers=[13078] #low q 12102 #midq 12100
sample=[13078]# low q 12077 #midq 12011
thickness=0.1
nbins=100
linear_binning=False
output_dir = f'/HFIR/CG2/shared/UserAcceptance/LDS_center/straight_through/masking'
flood_file='/HFIR/CG2/shared/drt_sensitivity/sens_c489_bar.nxs'
sample_temp='p1'
for subfolder in ['1D','2D']:
    output_folder=os.path.join(output_dir,subfolder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
image_args={'vmin':0}


# Next allows users to see a quick 2D images for debugging or quality check. Users might want to play with title (last line)

# In[62]:


def show_instr(ws, detector='main'):
    pixel_size_x=0.55
    pixel_size_y=0.45
    if detector == 'main':
        npix = 192
    elif detector == 'wing':
        npix = 160
    else:
        npix = 192 + 160
    data=ws.extractY().reshape(-1,8,256).T
    data2=data[:,[0,4,1,5,2,6,3,7],:]
    data2=data2.transpose().reshape(-1,256)
    Z=np.ma.masked_where(data2<1,data2)
 #   Z=data2

    x=np.arange(npix)
    y=np.arange(256)
    Z = np.log(np.transpose(Z))
#    Z=np.transpose(Z)
    scale=20
    plt.figure(figsize=(pixel_size_x*npix/scale, pixel_size_y*256/scale))
#    plt.pcolor(Z[:, :npix],vmin=-1e3,vmax=3e4)
    plt.pcolor(Z[:, :npix])
    plt.colorbar()
 #   plt.title("Field {} T".format(meta_data))
  #  filename=os.path.join(output_dir,f'figure_{runs}_2D.png')
 #   fig.savefig(filename)


# In[34]:


#defining 2D Gaussian fititng functions

def Gaussian2D(x1, x2, amp, wid1, wid2, cen1, cen2):
    if wid1<=0:
        print('wid1 =',wid1)
    if wid2<=0:
        print('wid2 =',wid2)
    val = (amp/(np.sqrt(2*np.pi)*np.sqrt(wid1*wid2))) * np.exp(-((x1-cen1)**2/(2*wid1**2)+(x2-cen2)**2/(2*wid2**2)))
    return val

def _buildModel(self, **kwargs):
    model = lmfit.Model(self._function, independent_vars=["x1", "x2"],
                        param_names=["amp", "wid", "cen1", "cen2"])
    return model

def fit(self, data, freeX, **kwargs):
    freeX = np.asarray(freeX, float)
    model = self._buildModel(**kwargs)
    params = self._generateModelParams(model, **kwargs)

    model.fit(data, x1=freeX[0], x2=freeX[1], params=params)


# Find beam center for main detector

# In[86]:


import lmfit

#loading beam center data
center_filename = f"/HFIR/CG2/IPTS-{ipts}/nexus/CG2_{centers[0]}.nxs.h5"
ws=sans.load_events(center_filename,output_workspace='ws_center', pixel_calibration=True)
#LoadInstrument(ws, InstrumentName='CG2', RewriteSpectraMap=False)
ws=sans.transform_to_wavelength(ws)
ws=drtsans.process_uncertainties.set_init_uncertainties(ws) 
#MaskAngle(ws, MinAngle=0.2)
sans.solid_angle_correction(ws)
drtsans.apply_sensitivity_correction(ws,flood_file,min_threshold=0.5,max_threshold=1.5)
MaskBTP(ws,Pixel="1-70,186-256")
xc, yc = sans.find_beam_center(ws) #I use COM to help give good starting position for 2D Gaussian but it works with giving it 0,0

print(xc,yc)
show_instr(ws,'main')
#fitting 2D gaussian to center data
x=[]
y=[]
intes=[]
intes_err=[]
keep=[]
si=ws.spectrumInfo()
for spectrum_number in range(ws.getNumberHistograms()):
    pos=si.position(spectrum_number)
    is_masked=si.isMasked(spectrum_number)
    x.append(pos.X())
    y.append(pos.Y())
    keep.append(not is_masked and  np.isfinite(ws.readY(spectrum_number)[0]))
    intes.append(ws.readY(spectrum_number)[0])
    intes_err.append(ws.readE(spectrum_number)[0])
    
x=np.array(x)
y=np.array(y)
intes=np.array(intes)
intes_err=np.array(intes_err)

x=x[keep]
y=y[keep]
intes=intes[keep]
intes_err=intes_err[keep]



ws2=ws.clone() # used for seeing if 2D Gaussian model gives good data before fitting (not needed for implementation) but good for testing
idx=np.arange(ws.getNumberHistograms())
idx=idx[keep]

#for i,j in enumerate(idx):
 #   ws2.setY(int(j),np.array([model[i]]))

show_instr(ws2)

#absolute scaling
sans.center_detector(ws,xc,yc)
scale_factor, scale_factor_error = attenuation_factor(ws)



# In[89]:



model = lmfit.Model(Gaussian2D, independent_vars=["x1", "x2"],
                        param_names=["amp", "wid1","wid2", "cen1", "cen2"])

params=lmfit.Parameters()
params.add("amp",value=ws.extractY().max())
params.add('wid1',value=0.01,min=0.005)#width in x
params.add('wid2',value=0.01,min=1e-10)#width in y
params.add('cen1',value=0)
params.add('cen2',value=0)
results=model.fit(intes, x1=x, x2=y,weights=1/intes_err, params=params)

print(results.fit_report())

fitted=model.eval(results.params,x1=x,x2=y)


for i,j in enumerate(idx):
    ws2.setY(int(j),np.array([fitted[i]]))
ws3=ws.clone()

for i,j in enumerate(idx):
    ws3.setY(int(j),np.array([intes[i]]))
show_instr(ws2) #fitted results in 2D form
show_instr(ws3) #data array

