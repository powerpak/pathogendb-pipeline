#!/usr/bin/env python

'''
Author: Ajay Obla, Luis Cunha
Date: 03/21/2018
Description: Script to scrap mslt information from pubmlst.org
Usage:  fetch_mslt.py <file.fasta>
'''


from bs4 import BeautifulSoup
import csv
import sys
import re
import os
import argparse
import json

def fetch_results(inputfile,mlstdb,output):
    #Fetch the result from the server
	curl_command = """(echo -n '{"base64":true,"sequence": "'; base64 """ + inputfile + """ ; echo '"}') | curl --proxy http://mgmt01.chimera.hpc.mssm.edu:3128 -s -H "Content-Type: application/json" -X POST "http://rest.pubmlst.org/db/pubmlst_""" + mlstdb + """_seqdef/schemes/1/sequence" -d @-"""
	
	
	stream = os.popen(curl_command)
	json_string = stream.read()
	json_data = json.loads(json_string)

	
	
	try:
		ST=json_data['fields']['ST']
	except KeyError:
		ST='No Match'
	try:
		mlst_clade=json_data['fields']['mlst_clade']
	except KeyError:
		mlst_clade='No Match'
	

	with open(output, "w") as f:
		
		f.write('\t'.join(['Allele','Length','Startposition','Endposition'])+'\n'+'\n')
		for gene in sorted(json_data['exact_matches'].keys()):
			allele_id=json_data['exact_matches'][gene][0]['allele_id']
			f.write(gene+':'+allele_id+'\n')

		f.write('\n'+'\n'+'MLST'+'\n'+'\n')
		f.write('ST'+'\t'+ST+'\n')
		f.write('mlstclade'+'\t'+mlst_clade)
			

def usage():
	print "\n\tScript to get mslt information from pubmlst.org.\n"
	print "\tRequirement: fasta input file and MLST database (Examples: 'cdifficile' for C.diff or 'mlst' for MRSA).\n" 
	print "\tUsage: fetch_mslt.py [-h] [--fasta FASTA] [--output OUTPUT] [--mlst DATABASE]\n\n"
	sys.exit(0)


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='scrap pubmlst.org for mslt')
	parser.add_argument('--fasta', "-f", help="Input fasta file ")
	parser.add_argument('--output', "-o", help="Output file with MLST info ")
	parser.add_argument('--mlst', "-m", help="MLST database to use. Examples: 'cdifficile' for C.diff or 'mlst' for MRSA.")
	args = parser.parse_args()
	if args.fasta==None or args.output==None or args.mlst==None:
		usage()
	if args.fasta!=None and args.output!=None and args.mlst!=None:  
		os.environ["CURL_CA_BUNDLE"] = "/etc/ssl/certs/ca-bundle.crt"
		fastaInput = args.fasta
		output = args.output
		MLSTdb = args.mlst

	fetch_results(fastaInput,MLSTdb,output)
    	

#Code to process legacy pubMLST curl command - Not in use
"""
def fetch_results(inputfile, mlstdb):
    #Fetch the result from the server

    curl_command = 'curl --proxy proxy.mgmt.hpc.mssm.edu:8123 --form "fasta_upload=@' + inputfile + '" --form "db=pubmlst_' + mlstdb + '_seqdef" --form "page=sequenceQuery" --form "locus=SCHEME_1" --form "order=locus" -F "submit=submit" -F "no_ajax=1" "https://pubmlst.org/perl/bigsdb/bigsdb.pl"'

	

    stream = os.popen(curl_command)
    html = stream.read()

 
    #load the result into a dom representation and parse the tables
    #there's 3 tables in the result. we are interested in the last 2
    dom = BeautifulSoup(html, "lxml")
    
    tables = dom.find_all("table")
    if len(tables)==3:
      alleles = tables[1]
      MSLT = tables[2]
    elif len(tables)==2:
      alleles = tables[0]
      MSLT=tables[1]
    elif len(tables)==1:
      alleles = tables[0]
      MSLT=""
    else:
      print "No match found\n"
      sys.exit(0)
    return (alleles, MSLT)

def  convertTabl2tsv(allelesTable, mlstTable, output):
    #convert from table to csv
    #handle new lines and whitespaces
    t1 = allelesTable
    t1 = re.sub('\n', '\t', t1.prettify())
    t1_ = BeautifulSoup(t1, "lxml")
   
    t2_=False
    if mlstTable!="": 
    	t2 = mlstTable
    	t2 = re.sub('\n', '\t', t2.prettify())
    	t2_ = BeautifulSoup(t2, "lxml")
    
    with open("out.tmp", "wb") as f:
        writer = csv.writer(f)
        h = t1_.find_all("th")
        row = [elem.text.encode('utf-8') for elem in h[:4]]
        writer.writerow(row)
    
        for tr in t1_.find_all("tr"):
           tds = tr.find_all('td')
           row = [elem.text.encode('utf-8') for elem in tds[:4]]
           writer.writerow(row)
    
        f.write("\n\n")
        f.write("MLST\n\n")
    
        if t2_:
        	for tr in t2_.find_all("tr"):
            		tds = tr.find_all('th')
            		c1 = [elem.text.encode('utf-8') for elem in tds[:4]]
            		tds = tr.find_all('td')
            		c2 = [elem.text.encode('utf-8') for elem in tds[:4]]
            		row = c1 + c2
            		writer.writerow(row)
        else:
		row = ["no match found"]
		#writer.writerow(row)
 		f.write("no_match_found")
    f=open(output, "w")
    f.close()
    with open(output, "wb") as f:
         g= open("out.tmp", "r")
         lines=g.readlines()
         g.close()
         for line in lines:
                 l=re.sub(" ", "", line)
                 l=re.sub("\t,\t", "\t", l)
                 l=re.sub("^\t", "", l)
                 l=l.strip()
                 f.write(l)
                 f.write("\n")
    
         os.remove("out.tmp")

"""
