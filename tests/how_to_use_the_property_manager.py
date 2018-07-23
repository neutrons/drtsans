from __future__ import print_function

from mantid.kernel import PropertyManagerDataService, PropertyManager

pm = PropertyManager()
pmds = PropertyManagerDataService.add("pm_name", pm)
pm.declareProperty("p1","v1")
pm.declareProperty("p2","v2")

for k,v in zip(pm.keys(), pm.values()):
    print("{} -> ".format(k), end="")
    try:
        print(v.value)
    except:
        print(v)