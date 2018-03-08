from ftplib import FTP
import sys
import argparse

parser=argparse.ArgumentParser(description='FTP script to send XML to SRA by ftp')
parser.add_argument('-p', '--password', help='ftp password')
parser.add_argument('-s', '--service', help='ftp service')
parser.add_argument('-n', '--ncbi', help='ncbi ftp address')
parser.add_argument('-d', '--directory', help='ftp directory')
args=parser.parse_args()

ftp=FTP(args.ncbi)
ftp.login(args.service, args.password)
ftp.cwd(args.directory)
fh=open("/sc/orga/work/attieo02/pathogendb-pipeline/scripts/submission.xml", 'r')
ftp.storlines("STOR submission.xml", fh)
#ftp.storlines("STOR submit.ready", "/sc/orga/work/attieo02/pathogendb-pipeline/scripts/submit.ready")

ftp.close()
