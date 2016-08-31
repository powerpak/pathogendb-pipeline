import sys
import MySQLdb
import re
import argparse
from ftplib import FTP
import collections

parser=argparse.ArgumentParser(description='FTP script to send PacBio reads by ftp to SRA')
parser.add_argument('-p', '--password', help='database password')
parser.add_argument('-u', '--user', help='database user')
parser.add_argument('-t', '--host', help='database host')
parser.add_argument('-d', '--database', help='database name')
parser.add_argument('-f', '--ftp_password', help='ftp password')
parser.add_argument('-s', '--service', help='ftp service name')
parser.add_argument('-n', '--ncbi', help='ncbi ftp address')
parser.add_argument('-r', '--directory', help='ftp directory')
parser.add_argument('-o', '--organism', help='organism name')
parser.add_argument('-i', '--isolate_file', help='isolate file')
args=parser.parse_args()

def queryPathogenDB(organism, isolateID_list):
    results=[]
    db=MySQLdb.connect(host=args.host, db=args.database, user=args.user, passwd=args.password)
    cur=db.cursor()
    for extractID in isolateID_list:
        cur.execute("select tIsolates.isolate_ID,tSequencing_runs.sequence_run_ID,tSequencing_runs.run_data_link from tSequencing_runs join tExtracts on tSequencing_runs.extract_ID=tExtracts.extract_ID join tStocks on tExtracts.stock_ID=tStocks.stock_ID join tIsolates on tStocks.isolate_ID=tIsolates.isolate_ID join tOrganisms on tIsolates.organism_ID=tOrganisms.organism_ID where tOrganisms.full_name=\'"+organism+"\' and tExtracts.extract_ID=\'"+extractID+"\'")
        for row in cur.fetchall():
            results.append(row)
        db.close()
        return results

def fill_data(results):
    run_data={}
    run_data=collections.defaultdict()
#    print results
    for i in (0, len(results)-1):
        run_data[results[i][0]]={}
    for i in (0, len(results)-1):
        run_data[results[i][0]][results[i][1]]=results[i][2]
    return run_data

isolate_ID_list=[]
fh=open(args.isolate_file, 'r')
for line in fh:
    isolate_ID_list.append(line.rstrip())
fh.close()
results=queryPathogenDB(args.organism, isolate_ID_list)

run_data=fill_data(results)





for isolate in run_data.keys():
    for exp in run_data[isolate].keys():
        url=re.split("/", run_data[isolate][exp])
        fh=open("/sc/orga/projects/pacbio/userdata_permanent/jobs/"+str(url[-1])[:3]+"/"+url[-1]+"/input.fofn", 'r')
        for line in fh.readlines():
            ftp=FTP(args.ncbi)
            ftp.login(args.service, args.ftp_password)
            ftp.cwd(args.directory)
            reads=line.rstrip().split("/")
            fh=open(line.rstrip(),'r')
            ftp.storbinary("STOR "+reads[-1], fh)
            ftp.close()
