import tempfile
from os.path import join
from ornl.sans.save_2d import save_nist_dat, save_nexus
from mantid.simpleapi import LoadNexus
from ornl.settings import unique_workspace_name

def test_save_nist_dat(refd):
    filename = join(refd.new.eqsans,'test_save_output/EQSANS_68200_iq.nxs')
    ws = LoadNexus(filename)
    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_nist_dat(ws, tmp.name)
        save_nist_dat(ws,'/tmp/save_nist_dat.DAT')
        #assert tmp.name == ''

def test_save_nexus(refd):
    filename = join(refd.new.eqsans,'test_save_output/EQSANS_68200_iq.nxs')
    ws = LoadNexus(filename)
    with tempfile.NamedTemporaryFile('r+') as tmp:
        save_nexus(ws, 'EQSANS 68200', tmp.name)
        save_nexus(ws, 'EQSANS 68200', '/tmp/save_nexus.nxs')
        #assert tmp.name == ''

