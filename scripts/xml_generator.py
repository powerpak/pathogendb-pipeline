import sys
bp=[]
bs=[]
reads=[]
reads_dict={}
bs_id=""
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
        if(data_list[0]!='<BioProject>'):
            bp.append(data_list)
    if(BioSample_flag):
        if(data_list[0]!='<BioSample>'):
            bs.append(data_list)
            bs_id=data_list[0]
            reads=[]
#            print bs_id
#        bs.append(data_list)
    if(Reads_flag):
        if(data_list[0]!='<Reads>'):
            reads.append(data_list)
    reads_dict[bs_id]=reads
#    print reads
#print bp[1][0]

print '''<Submission xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"../../xml/portal/submission.xsd\">
 <Description>
    <Comment> BP(1.0)+BS(1.0)+SRA</Comment>
     <Organization role=\"owner\" type=\"institute\">
      <Name>Icahn School of Medicine at Mount Sinai</Name>
       <Contact email=\"oliver.attie@mssm.edu\">
      <Name>
      <First>Oliver</First>
     <Last>Attie</Last>
     </Name>
     </Organization>
'''
print "\t<Hold release_date=\'"+bp[0][0]+"\'></Hold>"
print ''' </Description>
<Action>
 <AddData target_db=\"BioProject\">
   <Data content_type=\"xml\">
     <XmlContent>
      <Project schema_version=\"2.0\">
       <ProjectID>
'''
print "\t\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bp[0][1]+"</SPUID>"
print '''</ProjectID>
        <Descriptor>
'''
print "\t\t\t<Title>"+bp[0][2]+"</Title>"
print "\t\t\t<Description><p>"+bp[0][2]+"</p></Description>"
print '''    <Relevance>
             <Medical>Yes</Medical>
             </Relevance>
         </Descriptor>
         <ProjectType>
         <ProjectTypeSubmission sample_scope=\"eMonoisolate\">
         <Organism>
'''
print "\t\t\t<OrganismName>"+bp[0][3]+"</OrganimName>"
print '''  </Organism>
            <IntendedDataTypeSet>
             <DataType>genome sequencing</DataType>
            </IntendedDataTypeSet>
          </ProjectTypeSubmission>
      </ProjectType>
     </Project>
   </XmlContent>
</Data>
<Identifier>
'''
print "\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bp[0][4]+"</SPUID>"
print '''</Identifier>
            </AddData>
            </Action>
           <Action>
        <AddData target_db=\"BioSample\">
         <Data content_type=\"xml\">
          <XmlContent>
'''
#print len(bs)
for i in range(0,len(bs)):
    print '''     <BioSample schema_version=\"2.0\">
                      <SampleId>
              '''
#    print i
    print      "<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[i][0]+"</SPUID>"
    print '''         </SampleId>
                        <Descriptor>
        '''
    print "\t\t\t<Title>"+bs[i][1]+"</Title>"
    print '''       </Descriptor>
                     <Organism>'''
    print "\t\t\t<OrganismName>"+bs[i][2]+"</OrganismName>"
    print '''    </Organism>
                 <Package>Pathogen.env.1.0</Package>
                 <Attributes>'''
    print "\t\t\t<Attribute attribute_name=\"strain\">"+bs[i][0]+"</Attribute>"
    print '''<Attribute attribute_name=\"collected_by\">ISMMS</Attribute>
              '''
    print "\t\t\t<Attribute attribute_name=\"collected_date\">"+bs[i][3]+"</Attribute>"
    print "\t\t\t<Attribute attribute_name=\"isolation_source\">"+bs[i][4]+"</Attribute>"
    print '''    <Attribute attribute_name=\"geo_loc_name\">USA:New York City</Attribute>
                  <Attribute attribute_name=\"lat_lon\">40N 73W</Attribute>"
                  </Attributes>
                </BioSample>
              </XmlContent>
             </Data>
           <Identifier>
           '''
    print "\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[i][0]+"</SPUID>"
    print ''' </Identifier>
        </AddData>
        </Action>
        '''
    print '''<Action>
              <AddFiles target_db=\"SRA\">'''
#print BioSample_flag
#print reads
#print len(reads)
    for j in range(0,len(reads_dict[bs[i][0]])):
        print "\t\t<File file_path=\""+reads_dict[bs[i][0]][j][0]+"\">"
        print '''<DataType>generic-data</DataType>
        </File>'''

    print "\t\t<Attribute name=\"instrument_model\">"+reads_dict[bs[i][0]][0][1]+"</Attribute>"
    print "\t\t<Attribute name=\"library_name\">"+reads_dict[bs[i][0]][0][2]+"</Attribute>"
    print '''<Attribute name=\"library_strategy\">WGS</Attribute>
                 <Attribute name=\"library_source\">GENOMIC</Attribute>
                 <Attribute name=\"library_selection\">RANDOM</Attribute>
                 <Attribute name=\"library_layout\">FRAGMENTED</Attribute>
                 <AttributeRefId name=\"BioProject\">
                 <RefId>'''
    print "\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[i][0]+"</SPUID>"
    print '''</RefId>
                 </AttributeRefId>
                 <AttributeRefId name=\"BioSample\">
                 <RefId>'''
    print "\t\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+bs[i][0]+"</SPUID>"
    print '''</RefId>
                 </AttributeRefId>
                 <Identifier>'''
    print "\t\t<SPUID spuid_namespace=\"ISMMS_PSP\">"+reads_dict[bs[i][0]][0][3]+"</SPUID>"
    print '''</Identifier>
        </AddFiles>
     </Action>'''
print "</Submission>"
