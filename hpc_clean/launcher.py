#!/usr/bin/env python
# encoding: utf-8
"""
launcher.py

Created by Jian Banzai on 2014-05-06.
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
	dirs = ['hpc_mpi_mm','hpc_mpi_mm-fast','hpc_mpi_cannon','hpc_mpi_farm']
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
	nprocs = ['5','17']
	execs = ['hpc_mpi_mm/bin/hpc_mpi_mm', 'hpc_mpi_mm-fast/bin/hpc_mpi_mm-fast', 'hpc_mpi_cannon/bin/hpc_mpi_cannon', 'hpc_mpi_farm/bin/hpc_mpi_farm']
	optimization = [' ','-op ']
	dims = ['16','64','256']
	configs = [' 0 0', ' 0 1', ' 0 2', ' 1 0', ' 1 1', ' 1 2',] #first digist is for i/o on/off, second digit is for load_func low/medium/high
	for nproc in nprocs:
		for exe in execs:
			for op in optimization:
				for dim in dims:
					for config in configs:
						temp_args = "mpiexec -n " + nproc + " ./" + exe + op + dim + config
						args = temp_args.split()
						common(args)
	print '*** Ending loops ***'
	
if __name__ == '__main__':
	main()

