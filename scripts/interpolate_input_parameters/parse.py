#!/usr/bin/env python
import xml.etree.ElementTree as ET
import sys

input = 'wam_input2.xsd'
output = 'wam_input.asc'

tree = ET.parse(input)
root = tree.getroot()

f = open(output,'w')

f.write('Issue Date          '+tree.find('issue-date').text+"\n")
f.write('F10 81 Day Avg      '+tree.find('f10-81-avg-currentday').text+" \n")
f.write('Flags:  0=Forecast, 1=Estimated, 2=Observed \n\n')

f.write(' Date_Time                   F10           Kp       F10Flag      KpFlag   \n')
f.write('-----------------------------------------------------------------------   \n')

for child in root.findall('data-item'):
	f.write('{0}{1:>12}{2:>12}{3:>12}{4:>12}'.format(child.get('time-tag'),child.find('f10').text,child.find('kp').text,child.find('f10-flag').text,child.find('kp-flag').text)+' \n')
