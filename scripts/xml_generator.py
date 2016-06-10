'''This code generates an XML file for submission to NCBI from a flat file 
with BioProject information and the PathogenDB database.'''
import sys
import collections
import MySQLdb
import re

'''Queries pathogendb to get the necessary data for the XML'''
def queryPathogenDB():
    results=[]
    db=MySQLdb.connect(host="db.hpc.mssm.edu", db="vanbah01_pathogens", user="pathogendb_ro", passwd="avbikCog3")
    cur=db.cursor()
    cur.execute("select tIsolates.isolate_ID, tIsolates.collection_sourceA, tIsolates.collection_sourceB, tOrganisms.full_name, tIsolates.collection_date, tSequencing_runs.sequence_run_ID, tSequencing_runs.sequencing_platform, tSequencing_runs.read_length, tSequencing_runs.paired_end, tSequencing_runs.run_data_link from tSequencing_runs join tExtracts on tSequencing_runs.extract_ID=tExtracts.extract_ID join tStocks join tIsolates on tStocks.isolate_ID=tIsolates.isolate_ID join tOrganisms on tIsolates.organism_ID=tOrganisms.organism_ID where tOrganisms.organism_ID=88 limit 10")
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
datafile=sys.argv[1]

fh=open(datafile, 'r')
for line in fh:
    if(line.rstrip().split("\t")[0]=='<BioProject>'):
        BioProject_flag=1
    if(line.rstrip().split("\t")[0]=='</BioProject>'):
        BioProject_flag=0
    data_list=line.rstrip().split("\t")
    if(BioProject_flag==1):
        bp.append(data_list)

results=queryPathogenDB()
(bs_organism, bs_collection_sourceA, bs_collection_sourceB, bs_collection_date,exp_platform, exp_readlength, exp_paired_end,exp_run_data)=fill_data(results)
print "<Submission xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"../../xml/portal/submission.xsd\">"
print "<Description>"
print "<Comment> BP(1.0)+BS(1.0)+SRA</Comment>"
print "<Name>Icahn School of Medicine at Mount Sinai</Name>"
print "<Contact email=\"oliver.attie@mssm.edu\">"
print "<Name>"
print "<First>Oliver</First>"
print "<Last>Attie</Last>"
print "</Name>"
print "</Organization>"
print "<Hold release_date=\'"+bp[1][0]+"\'></Hold>"
print "</Description>"
print "<Action>"
print "<AddData target_db=\"BioProject\">"
print "<Data content_type=\"xml\">"
print "<XmlContent>"
print "<Project schema_version=\"2.0\">"
print "<ProjectID>"
print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+bp[2][0]+"</SPUID>"
print "</ProjectID>"
print "<Descriptor>"
print "<Title>"+bp[3][0]+"</Title>"
print "<Description><p>"+bp[2][0]+"</p></Description>"
print "<Relevance>"
print "<Medical>Yes</Medical>"
print "</Relevance>"
print "</Descriptor>"
print "<ProjectType>"
print "<ProjectTypeSubmission sample_scope=\"eMonoisolate\">"
print "<Organism>"
print "<OrganismName>"+bp[6][0]+"</OrganismName>"
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
print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+bp[2][0]+"</SPUID>"
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
    print "<Title>"+bs_organism[isolate]+" sample from "+bs_collection_sourceA[isolate]+"/"+bs_collection_sourceB[isolate]+"</Title>"
    print "</Descriptor>"
    print "<Organism>"
    print "<OrganismName>"+bs_organism[isolate]+"</OrganismName>"
    print "</Organism>"
    print "<Package>Pathogen.env.1.0</Package>"
    print "<Attributes>"
    print "<Attribute attribute_name=\"strain\">"+isolate+"</Attribute>"
    print "<Attribute attribute_name=\"collected_by\">ISMMS</Attribute>"
    print "<Attribute attribute_name=\"collection_date\">"+str(bs_collection_date[isolate])+"</Attribute>"
    print "<Attribute attribute_name=\"isolation_source\">"+bs_collection_sourceA[isolate]+"</Attribute>"
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
#    print reads
#    print range(1,len(reads[i].keys())+1)
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
        print "<Attribute name=\"instrument_model\">"+exp_platform[isolate][exp]+"</Attribute>"
        print "<Attribute name=\"library_name\">"+bs_organism[isolate]+" "+isolate+"</Attribute>"
        print "<Attribute name=\"library_strategy\">WGS</Attribute>"
        print "<Attribute name=\"library_source\">GENOMIC</Attribute>"
        print "<Attribute name=\"library_selection\">RANDOM</Attribute>"
        print "<Attribute name=\"library_layout\">FRAGMENTED</Attribute>"
        print "<AttributeRefId name=\"BioProject\">"
        print "<RefId>"
        print "<SPUID spuid_namespace=\"ISMMS_PSP\">"+isolate+"</SPUID>"
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
