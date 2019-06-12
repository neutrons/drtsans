from __future__ import print_function

from mantid.plots.helperfunctions import get_spectrum
from mantid.simpleapi import SaveCanSAS1D


def save_ascii_1D(wksp, logName, fileName):
    q, intensity, sigIntensity, dq = get_spectrum(wksp, 0, True,
                                                  True, True)
    f = open(fileName, "w+")
    f.write('# reduction log: ' + logName + '\n')
    f.write('#Q (1/A)        I (1/cm)        dI (1/cm)       dQ (1/A)\n')
    for i in range(len(intensity)):
        f.write('{:.6f}\t'.format(q[i]))
        f.write('{:.6f}\t'.format(intensity[i]))
        f.write('{:.6f}\t'.format(sigIntensity[i]))
        f.write('{:.6f}\n'.format(dq[i]))
    f.close()


def save_xml_1D(wksp, logName, fileName):
    SaveCanSAS1D(InputWorkspace=wksp, Process=logName, Filename=fileName)


def save_ascii_2D(q2, q2x, q2y, logName, fileName):
    f = open(fileName, "w+")
    f.write('# reduction log: ' + logName + '\n')
    f.write('#Qx (1/A)       Qy (1/A)        I (1/cm)        dI (1/cm)'
            + '       dQx (1/A)       dQy (1/A)\n')
    f.write('ASCII data\n\n')
    for i in range(len(q2.readY(0))):
        for j in range(q2.getNumberHistograms()):
            qy = float(q2.getAxis(1).label(j))
            x = 0.5*(q2.readX(j)[i]+q2.readX(j)[i+1])
            f.write('{:.6f}\t'.format(x))
            f.write('{:.6f}\t'.format(qy))
            f.write('{:.6f}\t'.format(q2.readY(j)[i]))
            f.write('{:.6f}\t'.format(q2.readE(j)[i]))
            f.write('{:.6f}\t'.format(q2x.readY(j)[i]))
            f.write('{:.6f}\n'.format(q2y.readY(j)[i]))
    f.close()
