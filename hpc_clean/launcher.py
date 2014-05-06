#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Jian Banzai on 2014-05-06.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import subprocess

def main():
	#args = ("bin/bar", "-c", "somefile.xml", "-d", "text.txt", "-r", "aString", "-f", "anotherString")
	#args = ("mpiexec", "-n 10", "./hpc_mpi_mm-fast/bin/hpc_mpi_mm-fast", "64", "0", "0")
	#Or just:
	#args = "bin/bar -c somefile.xml -d text.txt -r aString -f anotherString".split()
	args = "mpiexec -n 10 ./hpc_mpi_mm-fast/bin/hpc_mpi_mm-fast 1920 0 0".split()
	popen = subprocess.Popen(args, stdout=subprocess.PIPE)
	popen.wait()
	output = popen.stdout.read()
	print output

if __name__ == '__main__':
	main()

