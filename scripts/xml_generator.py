import sys
bp=[]
bs=[]
reads=[]
BioProject_flag=0
BioSample_flag=0
Reads_flag=0
datafile=sys.argv[1]

fh=open(datafile, 'r')
for line in fh:
    if(line.rstrip().split("\t")[0]=='<BioProject>'):
        BioProject_flag=1
    if(line.rstrip().split("\t")[0]=='</BioProject>'):
        BioProject_flag=0
    if(line.rstrip().split("\t")[0]=="<BioSample>"):
        BioSample_flag=1
    if(line.rstrip().split("\t")[0]=="</BioSample>"):
        BioSample_flag=0
    if(line.rstrip().split("\t")[0]=="<Reads>"):
        Reads_flag=1
    if(line.rstrip().split("\t")[0]=="</Reads>"):
        Reads_flag=0
    data_list=line.rstrip().split(",")
    if(BioProject_flag==1):
        bp.append(data_list)
    if(BioSample_flag):
        bs.append(data_list)
    if(Reads_flag):
        reads.append(data_list)
#print bp[1][0]
print "<Submission xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"../../xml/portal/submission.xsd\">"
print "\t<Description>"
print "\t\t<Comment> BP(1.0)+BS(1.0)+SRA</Comment>"
print "\t\t<Name>Icahn School of Medicine at Mount Sinai</Name>"
print "\t\t<Contact email=\"oliver.attie@mssm.edu\">"
print "\t\t<Name>"
print "\t\t\t<First>Oliver</First>"
print "\t\t\t<Last>Attie</Last>"
print "\t\t</Name>"
print "\t</Organization>"
print "\t<Hold release_date=\'"+bp[1][0]+"\'></Hold>"
print "\t</Description>"
print "<Action>"
print "\t<AddData target_db=\"BioProject\">"
print "\t\t<Data content_type=\"xml\">"
print "\t\t\t<XmlContent>"
print "\t\t\t<Project schema_version=\"2.0\">"
print "\t\t\t<ProjectID>"
print "\t\t\t<<SPUID spuid_namespace=\"ISMMS_PSP\">"+bp[1][1]+"</SPUID>"
print "\t\t\t</ProjectID>"
print "\t\t\t<Descriptor>"
print "\t\t\t<Title>"+bp[1][2]+"</Title>"
print "\t\t\t<Description><p>"+bp[1][2]+"</p></Description>"
print "\t\t\t<Relevance>"
print "\t\t\t<Medical>Yes</Medical>"
print "\t\t\t</Relevance>"
print "\t\t</Descriptor>"
print "\t\t<ProjectType>"
print "\t\t<ProjectTypeSubmission sample_scope=\"eMonoisolate\">"
print "\t\t\t<Organism>"
print "\t\t\t<OrganismName>"+bp[1][3]+"</OrganimName>"
print "\t\t</Organism>"
print "\t\t<IntendedDataTypeSet>"
print "\t\t\t<DataType>genome sequencing</DataType>"
print "\t\t</IntendedDataTypeSet>"
print "\t\t</ProjectTypeSubmission>"
print "\t</ProjectType>"
print "\t</Project>"
print "\t</XmlContent>"
print "</Data>"
print "<Identifier>"
print "\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bp[1][4]+"</SPUID>"
print "</Identifier>"
print "</AddData>"
print "</Action>"
print "<Action>"
print "\t<AddData target_db=\"BioSample\">"
print "\t\t<Data content_type=\"xml\">"
print "\t\t\t<XmlContent>"
print "\t\t\t<BioSample schema_version=\"2.0\">"
print "\t\t\t<SampleId>"
#print bs
print "\t\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[1][0]+"</SPUID>"
print "\t\t\t</SampleId>"
print "\t\t\t<Descriptor>"
print "\t\t\t<Title>"+bs[1][1]+"</Title>"
print "\t\t\t</Descriptor>"
print "\t\t\t<Organism>"
print "\t\t\t<OrganismName>"+bs[1][2]+"</OrganismName>"
print "\t\t\t</Organism>"
print "\t\t\t<Package>Pathogen.env.1.0</Package>"
print "\t\t\t<Attributes>"
print "\t\t\t<Attribute attribute_name=\"strain\">"+bs[1][0]+"</Attribute>"
print "\t\t\t<Attribute attribute_name=\"collected_by\">ISMMS</Attribute>"
print "\t\t\t<Attribute attribute_name=\"collected_date\">"+bs[1][3]+"</Attribute>"
print "\t\t\t<Attribute attribute_name=\"isolation_source\">"+bs[1][4]+"</Attribute>"
print "\t\t\t<Attribute attribute_name=\"geo_loc_name\">USA:New York City</Attribute>"
print "\t\t\t<Attribute attribute_name=\"lat_lon\">40N 73W</Attribute>"
print "\t\t</Attributes>"
print "\t\t</BioSample>"
print "\t</XmlContent>"
print "\t</Data>"
print "\t<Identifier>"
print "\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[1][0]+"</SPUID>"
print "\t</Identifier>"
print "</AddData>"
print "</Action>"
print "<Action>"
print "\t<AddFiles target_db=\"SRA\">"
#print reads
#print len(reads)
for i in range(1,len(reads)):
    if i%2==1:
        print "\t\t<File file_path=\""+reads[i][0]+"\">"
        print "\t\t<DataType>generic-data</DataType>"
        print "\t\t</File>"
        print "\t\t<Attribute name=\"instrument_model\">"+reads[i][1]+"</Attribute>"
        print "\t\t<Attribute name=\"library_name\">"+reads[i][2]+"</Attribute>"
        print "\t\t<Attribute name=\"library_strategy\">WGS</Attribute>"
        print "\t\t<Attribute name=\"library_source\">GENOMIC</Attribute>"
        print "\t\t<Attribute name=\"library_selection\">RANDOM</Attribute>"
        print "\t\t<Attribute name=\"library_layout\">FRAGMENTED</Attribute>"
        print "\t\t<AttributeRefId name=\"BioProject\">"
        print "\t\t<RefId>"
        print "\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[1][0]+"</SPUID>"
        print "\t\t</RefId>"
        print "\t</AttributeRefId>"
        print "\t<AttributeRefId name=\"BioSample\">"
        print "\t\t<RefId>"
        print "\t\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[1][0]+"</SPUID>"
        print "\t\t</RefId>"
        print "\t</AttributeRefId>"
        print "\t<Identifier>"
        print "\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+reads[i][3]+"</SPUID>"
        print "\t</Identifier>"
        print "</AddFiles>"
print "</Action>"
print "</Submission>"
