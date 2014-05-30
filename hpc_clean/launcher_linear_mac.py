#!/usr/bin/env python
# encoding: utf-8
"""
launcher.py

Created by dzt on 2014-05-20.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import subprocess

def main():
	makeall()	
	loops()
	
#make all versions
def makeall(): 
	print '*** Make process starting for all versions ***'
	dirs = ['hpc_linear_mm']
	options = ['','opti']
	for path in dirs:
		os.chdir(path) #move into the desired path
		for opt in options:
			shellcom = "make " + opt
			shellcom = shellcom.split()
			common(shellcom)
		os.chdir("..") #move back to the previous path
	print '*** Make process finished ***'
	
def common(args):
	popen = subprocess.Popen(args, stdout=subprocess.PIPE)
	popen.wait()
	output = popen.stdout.read()
	print output
	
def loops():
	print '*** Starting loops ***'
	execs = ['hpc_linear_mm/bin/hpc_linear_mm']
	optimization = [' ','-op ']
	#dims = ['192','960']
	#configs = [' 0 0', ' 0 1', ' 0 2', ' 1 0', ' 1 1', ' 1 2',] #first digist is for i/o on/off, second digit is for load_func low/medium/high
	dimconfigs = ['192 0 0','960 1 2']
	
	for exe in execs:
		for op in optimization:
			for dimconfig in dimconfigs:
				temp_args = " ./" + exe + op + dimconfig
				args = temp_args.split()
				common(args)
	
	print '*** Ending loops ***'
	
if __name__ == '__main__':
	main()

