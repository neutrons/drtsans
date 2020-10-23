import pytest
from mantid.simpleapi import MaskBTP
import os
import numpy as np
from drtsans.mono import gpsans as sans 
from drtsans.mono.gpsans import attenuation_factor
import drtsans
import lmfit
import pytest

#defining 2D Gaussian fititng functions

def Gaussian2D(x1, y1, amp, sigma_x, sigma_y, theta, x0, y0):
    a = np.cos(theta)**2/(2.*sigma_x**2) + np.sin(theta)**2/(2.*sigma_y**2)
    b = -np.sin(2.*theta)/(4.*sigma_x**2) + np.sin(2.*theta)/(4.*sigma_y**2)
    c = np.sin(theta)**2/(2.*sigma_x**2) + np.cos(theta)**2/(2.*sigma_y**2)
    amplitude = amp /(np.sqrt(2.*np.pi)*np.sqrt(sigma_x*sigma_y))
    val = amplitude * np.exp(-(a*(x1-x0)**2 + 2.*b*(x1-x0)*(y1-y0) + c*(y1-y0)**2))
    return val

def test_gaussian_fit():
    flood_file='/HFIR/CG2/shared/drt_sensitivity/sens_c489_bar.nxs'
    #Find beam center for main detector
    #loading beam center data
    center_filename = "/HFIR/CG2/IPTS-26004/nexus/CG2_13078.nxs.h5"
    ws=sans.load_events(center_filename,output_workspace='ws_center', pixel_calibration=True)
    ws=sans.transform_to_wavelength(ws)
    ws=drtsans.process_uncertainties.set_init_uncertainties(ws) 
    sans.solid_angle_correction(ws)
    drtsans.apply_sensitivity_correction(ws,flood_file,min_threshold=0.5,max_threshold=1.5)
    MaskBTP(ws,Pixel="1-70,186-256")
    xc, yc = sans.find_beam_center(ws) #I use COM to help give good starting position for 2D Gaussian but it works with giving it 0,0

    print(xc,yc)

    #fitting 2D gaussian to center data
    x=[]
    y=[]
    intes=[]
    intes_err=[]
    keep=[]
    si=ws.spectrumInfo()
    for i,it in enumerate(si):
        pos=si.position(i)
        is_masked=si.isMasked(i)
        x.append(pos.X())
        y.append(pos.Y())
        keep.append(not is_masked and np.isfinite(ws.readY(i)[0]))
        intes.append(ws.readY(i)[0])
        intes_err.append(ws.readE(i)[0])
    
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

    #absolute scaling
    sans.center_detector(ws,xc,yc)
    scale_factor, scale_factor_error = attenuation_factor(ws)

    model = lmfit.Model(Gaussian2D, independent_vars=["x1", "y1"],
                        param_names=["amp", "sigma_x","sigma_y", "theta", "x0", "y0"])

    params=lmfit.Parameters()
    params.add("amp",value=ws.extractY().max())
    params.add('sigma_x',value=0.01,min=1.e-10)#width in x
    params.add('sigma_y',value=0.01,min=1.e-10)#width in y
    params.add('theta',value=0.1, min=0., max=np.pi)
    params.add('x0',value=0.)
    params.add('y0',value=0.)
    results=model.fit(intes, x1=x, y1=y,weights=1/intes_err, params=params)

    print(results.fit_report())

    fitted=model.eval(results.params,x1=x,y1=y)

    for i,j in enumerate(idx):
        ws2.setY(int(j),np.array([fitted[i]]))
    ws3=ws.clone()

    for i,j in enumerate(idx):
        ws3.setY(int(j),np.array([intes[i]]))

if __name__ == '__main__':
    pytest.main([__file__])
