#python thermodata.py */*.json
import json
import sys
from dfttk.analysis.ywplot import Myjsonout

record = []
for jj in range (1,len(sys.argv)):
    print (sys.argv[jj])
    with open (sys.argv[jj], "r") as fp:
        orec = json.load (fp)
        record.extend(orec)

with open ("ExptData.json", "w") as out:
    Myjsonout(record,out)

print ("\n", len(record), "records handled.\n")

with open ("ExptData.json", "r") as fp:
    orec = json.load (fp)

