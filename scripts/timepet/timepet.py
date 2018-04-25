#!/usr/bin/env python
from datetime import datetime
import sys

# run as ./timepet.py [path to IPE PET???.ESMF_Logfile]
# will print out the length of time between each `sub-update_IPE finished` call

d1 = ""
flag = True

for line in open(sys.argv[1]):
	if "sub-update_IPE finished" in line or flag and ">>>IPM" in line:
		d2 = datetime.strptime(line[:19], "%Y%m%d %H%M%S.%f")
		if d1 is not "":
			print round(float(abs((d2 - d1).total_seconds()))/60,4)
		d1 = d2
		flag = False
