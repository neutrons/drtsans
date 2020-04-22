# Test gpsans.api.load_all_files()
# ... ...
# ... ...
from drtsans.mono.gpsans import load_all_files
from drtsans.mono.gpsans import load_and_split


def test_load_all_files():

    load_all_files(reduction_input=balbla,
                   prefix='TestCase1',
                   load_params=None,
                   path=Whatever)

    for ws_name in mtd.getObjectNames():
        print(ws_name)