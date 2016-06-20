'''This code generates an XML file for submission to NCBI from a flat file 
with BioProject information and the PathogenDB database.'''
import sys
import collections
import MySQLdb
import re
import argparse
parser=argparse.ArgumentParser(description='XML generator to generate XML from PathogenDB database for GenBank submission')
parser.add_argument('-r','--release_date', help='Release date for GenBank submission')
parser.add_argument('-n','--bioproject_name', help='Name of BioProject')
parser.add_argument('-b','--bioproject_title', help='Title of BioProject')
parser.add_argument('-m','--medical_relevance', help='Medical Relevance')
parser.add_argument('-t','--type', help='Mono or Multiisolate')
parser.add_argument('-o','--organism_name', help='Organism Name')
parser.add_argument('-i','--isolates', help='Isolates File')
args=parser.parse_args()

'''Queries pathogendb to get the necessary data for the XML'''
def queryPathogenDB(organism,isolateID_list):
    results=[]
    db=MySQLdb.connect(host="db.hpc.mssm.edu", db="vanbah01_pathogens", user="pathogendb_ro", passwd="avbikCog3")
    cur=db.cursor()
    for isolateID in isolateID_list:
        cur.execute("select tIsolates.isolate_ID, tIsolates.collection_sourceA, tIsolates.collection_sourceB, tOrganisms.full_name, tIsolates.collection_date, tSequencing_runs.sequence_run_ID, tSequencing_runs.sequencing_platform, tSequencing_runs.read_length, tSequencing_runs.paired_end, tSequencing_runs.run_data_link from tSequencing_runs join tExtracts on tSequencing_runs.extract_ID=tExtracts.extract_ID join tStocks on tExtracts.stock_ID=tStocks.stock_ID join tIsolates on tStocks.isolate_ID=tIsolates.isolate_ID join tOrganisms on tIsolates.organism_ID=tOrganisms.organism_ID where tOrganisms.full_name=\'"+organism+"\' and tIsolates.isolate_ID=\'"+isolateID+"\'")
        for row in cur.fetchall():
            results.append(row)
    db.close()
    return results

'''Fills the dictionaries with the data from PathogenDB'''
def fill_data(results):
    bs_organism={}
    bs_collection_sourceA={}
    bs_collection_sourceB={}
    bs_collection_date={}
    exp_platform=collections.defaultdict()
    exp_readlength=collections.defaultdict()
    exp_paired_end=collections.defaultdict()
    exp_run_data=collections.defaultdict()
    for i in range(0,len(results)):
        bs_organism[results[i][0]]=results[i][3]
        bs_collection_sourceA[results[i][0]]=results[i][1]
        bs_collection_sourceB[results[i][0]]=results[i][2]
        bs_collection_date[results[i][0]]=results[i][4]
        exp_platform[results[i][0]]={}
        exp_readlength[results[i][0]]={}
        exp_paired_end[results[i][0]]={}
        exp_run_data[results[i][0]]={}
    for i in range(0, len(results)):
        exp_platform[results[i][0]][results[i][5]]=results[i][6]
        exp_readlength[results[i][0]][results[i][5]]=results[i][7]
        exp_paired_end[results[i][0]][results[i][5]]=results[i][8]
        exp_run_data[results[i][0]][results[i][5]]=results[i][9]
    return (bs_organism, bs_collection_sourceA, bs_collection_sourceB, bs_collection_date, exp_platform, exp_readlength, exp_paired_end, exp_run_data)


bp=[]
isolate_ID_list=[]
#datafile=sys.argv[1]
#isolates_list=sys.argv[2]
#fh=open(datafile, 'r')
#for line in fh:
#    if(line.rstrip().split("\t")[0]=='<BioProject>'):
#        BioProject_flag=1
#    if(line.rstrip().split("\t")[0]=='</BioProject>'):
#        BioProject_flag=0
#    data_list=line.rstrip().split("\t")
#    if(BioProject_flag==1):
#        bp.append(data_list)
#fh.close()
fh1=open(args.isolates, 'r')
for line in fh1:
    isolate_ID_list.append(line.rstrip())
fh1.close()

results=queryPathogenDB(args.organism_name,isolate_ID_list)
(bs_organism, bs_collection_sourceA, bs_collection_sourceB, bs_collection_date,exp_platform, exp_readlength, exp_paired_end,exp_run_data)=fill_data(results)
print "<Submission xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"../../xml/portal/submission.xsd\">"
print "<Description>"
print "<Comment> BP(1.0)+BS(1.0)+SRA</Comment>"
print "<Organization role=\"owner\" type=\"institute\">"
print "<Name>Icahn School of Medicine at Mount Sinai</Name>"
print "<Contact email=\"oliver.attie@mssm.edu\">"
print "<Name>"
print "<First>Oliver</First>"
print "<Last>Attie</Last>"
print "</Name>"
print "</Contact>"
print "</Organization>"
print "<Hold release_date=\'"+args.release_date+"\'></Hold>"
print "</Description>"
print "<Action>"
print "<AddData target_db=\"BioProject\">"
print "<Data content_type=\"xml\">"
print "<XmlContent>"
print "<Project schema_version=\"2.0\">"
print "<ProjectID>"
print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+args.bioproject_name+"</SPUID>"
print "</ProjectID>"
print "<Descriptor>"
print "<Title>"+args.bioproject_title+"</Title>"
print "<Description><p>"+args.bioproject_name+"</p></Description>"
print "<Relevance>"
print "<Medical>Yes</Medical>"
print "</Relevance>"
print "</Descriptor>"
print "<ProjectType>"
print "<ProjectTypeSubmission sample_scope=\"eMonoisolate\">"
print "<Organism>"
print "<OrganismName>"+args.organism_name+"</OrganismName>"
print "</Organism>"
print "<IntendedDataTypeSet>"
print "<DataType>genome sequencing</DataType>"
print "</IntendedDataTypeSet>"
print "</ProjectTypeSubmission>"
print "</ProjectType>"
print "</Project>"
print "</XmlContent>"
print "</Data>"
print "<Identifier>"
print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+args.bioproject_name+"</SPUID>"
print "</Identifier>"
print "</AddData>"
print "</Action>"
print "<Action>"
print "<AddData target_db=\"BioSample\">"
print "<Data content_type=\"xml\">"
print "<XmlContent>"
print "<BioSample schema_version=\"2.0\">"
print "<SampleId>"
for isolate in bs_organism.keys():
    print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+isolate+"</SPUID>"
    print "</SampleId>"
    print "<Descriptor>"
    if(len(bs_collection_sourceA[isolate])>0):
        print "<Title>"+bs_organism[isolate]+" sample from "+bs_collection_sourceA[isolate]+"</Title>"
    else:
        print "<Title>"+bs_organism[isolate]+" sample from "+bs_collection_sourceB[isolate]+"</Title>"
    print "</Descriptor>"
    print "<Organism>"
    print "<OrganismName>"+bs_organism[isolate]+"</OrganismName>"
    print "</Organism>"
    print "<Package>Pathogen.env.1.0</Package>"
    print "<Attributes>"
    print "<Attribute attribute_name=\"strain\">"+isolate+"</Attribute>"
    print "<Attribute attribute_name=\"collected_by\">ISMMS</Attribute>"
    print "<Attribute attribute_name=\"collection_date\">"+str(bs_collection_date[isolate])+"</Attribute>"
    if(len(bs_collection_sourceA[isolate])>0):
        print "<Attribute attribute_name=\"isolation_source\">"+bs_collection_sourceA[isolate]+"</Attribute>"
    else:
        print "<Attribute attribute_name=\"isolation_source\">"+bs_collection_sourceB[isolate]+"</Attribute>"
    print "<Attribute attribute_name=\"geo_loc_name\">USA:New York City</Attribute>"
    print "<Attribute attribute_name=\"lat_lon\">40N 73W</Attribute>"
    print "</Attributes>"
    print "</BioSample>"
    print "</XmlContent>"
    print "</Data>"
    print "<Identifier>"
    print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+isolate+"</SPUID>"
    print "</Identifier>"
    print "</AddData>"
    print "</Action>"
    for exp in exp_platform[isolate].keys():
        print "<Action>"
        print "<AddFiles target_db=\"SRA\">"
        url=re.split("/", exp_run_data[isolate][exp])
        fh=open("/sc/orga/projects/pacbio/userdata_permanent/jobs/"+str(url[-1])[:3]+"/"+url[-1]+"/input.fofn", 'r')

        for line in fh.readlines():
            reads=line.split("/")
            print "<File file_path=\""+reads[-1]+"\">"
#        print runs
#        print range(0, len(runs[i-1][j-1].keys())+1)
#        print runs[0][0][1][1][0]
#        for k in range(1, len(runs[i-1][j-1].keys())+1):
#            print "\t\t<File file_path=\""+str(runs[i-1][j-1][i][1][0])+"\">"
            print "<DataType>generic-data</DataType>"
            print "</File>"
        if(exp_platform[isolate][exp]=='Pacbio'):
            print "<Attribute name=\"instrument_model\">PacBio RSII</Attribute>"
        else:
            print "<Attribute name=\"instrument_model\">Illumina HiSeq 2500</Attribute>"
        print "<Attribute name=\"library_name\">"+bs_organism[isolate]+" "+isolate+"</Attribute>"
        print "<Attribute name=\"library_strategy\">WGS</Attribute>"
        print "<Attribute name=\"library_source\">GENOMIC</Attribute>"
        print "<Attribute name=\"library_selection\">RANDOM</Attribute>"
        if(exp_paired_end[isolate][exp]=="No"):
            print "<Attribute name=\"library_layout\">FRAGMENTED</Attribute>"
        else:
            print "<Attribute name=\"library_layout\">PAIRED_END</Attribute>"
        print "<AttributeRefId name=\"BioProject\">"
        print "<RefId>"
        print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+args.bioproject_name+"</SPUID>"
        print "</RefId>"
        print "</AttributeRefId>"
        print "<AttributeRefId name=\"BioSample\">"
        print "<RefId>"
        print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+isolate+"</SPUID>"
        print "</RefId>"
        print "</AttributeRefId>"
        print "<Identifier>"
        print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+exp+"</SPUID>"
        print "</Identifier>"
        print "</AddFiles>"
        print "</Action>"
print "</Submission>"
