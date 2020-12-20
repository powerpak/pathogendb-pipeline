#!/usr/bin/env python

# -*- coding: utf-8 -*-

"""
@Author: Ajay Obla
@Created: 07-24-2020
@Description: This script updates assembly table in PDB with mec cassette annotation and spa type
"""

import mysql.connector
import os
import sys
import argparse

parser = argparse.ArgumentParser(description='add spa and mec annotation in PDB')
parser.add_argument('--spa', "-s", help="spa type output file")
parser.add_argument('--mec', "-m", help="mec output file")
parser.add_argument('--assembly_id', "-a", help="Assembly ID")
args = parser.parse_args()

path = os.path.expanduser('~') + '/.my.cnf'
with open(path) as cnf_file:
	for line in cnf_file:
		if line.startswith('user='):
			user = line.rstrip()[5:]
		if line.startswith('password='):
			pw = line.rstrip()[9:]
		if line.startswith('host='):
			host = line.rstrip()[5:]
		if line.startswith('database='):
			database = line.rstrip()[9:]

db = mysql.connector.connect(host=host,user=user,passwd=pw,db=database)

spa_result_file=args.spa
mec_result_file=args.mec
assembly_id=args.assembly_id

try:
	with open(spa_result_file,'r') as spa_result:
		for i, line in enumerate(spa_result):
			try:
				if i==1:					
					spa_type=line.rsplit('\t')[2]
					#if not spa_type.startswith('t'):
					#	spa_type='novel'
			except IndexError:
				spa_type='none'
	spa_type=spa_type.strip()
	cur2=db.cursor()
	cur2.execute("UPDATE `tAssemblies` SET spa_type ='%s' where assembly_ID='%s'" % (spa_type,assembly_id))

except IOError:
	print "Spa Type result file not found"

try:
	with open(mec_result_file,'r') as mec_result:
		for i, line in enumerate(mec_result):
			if line.startswith('The input organism'):
					if 'MSSA' in line:
						mec_cassette_type='none'
					if 'MRSA' in line:
						mec_cassette_type=1	
			elif line.startswith('Predicted SCCmec element:'):
					mec_cassette_type=line.split(':')[1]
			elif line.startswith('Alert!'):
					mec_cassette_type='multiple'
	mec_cassette_type=mec_cassette_type.strip()
	cur2=db.cursor()
	cur2.execute("UPDATE `tAssemblies` SET mec_cassette ='%s' where assembly_ID='%s'" % (mec_cassette_type,assembly_id))
except IOError:
	print "MEC Cassette result file not found"

db.commit()
	

