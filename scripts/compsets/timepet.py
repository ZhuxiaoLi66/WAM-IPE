#!/usr/bin/env python
from datetime import datetime
import sys

if len(sys.argv) == 3:
	d1 = datetime.strptime(sys.argv[1], "%Y%m%d %H%M%S")
	d2 = datetime.strptime(sys.argv[2], "%Y%m%d %H%M%S")
	print round(float(abs((d2 - d1).total_seconds()))/60,2)
