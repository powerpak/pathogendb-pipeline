#!/usr/bin/env python

'''
Author: Luis Cunha
Date: 02/03/2015
Description: Script to scrap mslt information from pubmlst.org
Usage:  fect_mslt.py <file.fasta>
'''


from bs4 import BeautifulSoup
import csv
from poster.encode import multipart_encode
from poster.streaminghttp import register_openers
import urllib2
import sys
import re
import os
import argparse

register_openers()


def fetch_results(file):
    #Fetch the result from the server
    params = {'fasta_upload': open(file, "rb"), 'db': 'pubmlst_cdifficile_seqdef', "page":"sequenceQuery", "locus":"SCHEME_1", "submit":"submit" }
    datagen, headers =multipart_encode(params)
    upload_url = "http://pubmlst.org/perl/bigsdb/bigsdb.pl"
   
    #os.environ['http_proxy'] = 'proxy.mgmt.hpc.mssm.edu:8123' 
    #proxy = urllib2.ProxyHandler({'http': 'proxy.mgmt.hpc.mssm.edu:8123'})    
    #opener = urllib2.build_opener(proxy)
    #urllib2.install_opener(opener)
    
    request = urllib2.Request(upload_url, datagen, headers)
    #result = urllib2.urlopen(upload_url)
    result = urllib2.urlopen(request)
    html=result.read()
 
    #load the result into a dom representation and parse the tables
    #there's 3 tables in the result. we are interested in the last 2
    dom = BeautifulSoup(html)
    
    tables = dom.find_all("table")
    alleles = tables[1]
    MLST = tables[2]
    return (alleles, MLST)




def  convertTabl2tsv(allelesTable, mlstTable, output):
    #convert from table to csv
    #handle new lines and whitespaces
    t1 = allelesTable
    t1 = re.sub('\n', '\t', t1.prettify())
    t1_ = BeautifulSoup(t1)
    
    t2 = mlstTable
    t2 = re.sub('\n', '\t', t2.prettify())
    t2_ = BeautifulSoup(t2)
    
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
    
    
        for tr in t2_.find_all("tr"):
            tds = tr.find_all('th')
            c1 = [elem.text.encode('utf-8') for elem in tds[:4]]
            tds = tr.find_all('td')
            c2 = [elem.text.encode('utf-8') for elem in tds[:4]]
            row = c1 + c2
            writer.writerow(row)
    
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


def usage():
    print "\n\tScript to get mslt information from pubmlst.org.\n"
    print "\tRequirement: fasta input file\n" 
    print "\tUsage: fetch_mslt.py [-h] [--fasta FASTA] [--output OUTPUT]\n\n"
    sys.exit(0)


if __name__ == "__main__":

    #if len(sys.argv)!=2:
    #    usage()
    parser = argparse.ArgumentParser(description='scrap pubmlst.org for mslt')
    parser.add_argument('--fasta', "-f", help="input fasta file ")
    parser.add_argument('--output', "-o", help="input fasta file ")
    args = parser.parse_args()
    if args.fasta==None or args.output==None:
    	usage()
    if args.fasta!=None and args.output!=None:  
    	fastaInput = args.fasta
    	output = args.output

    	allelesTable, mlstTable = fetch_results(fastaInput)
    	convertTabl2tsv(allelesTable, mlstTable, output)


