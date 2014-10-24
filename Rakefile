require 'pp'
require_relative 'lib/colors'
include Colors

SOURCE_BASES = Rake::FileList.new("inputs/*_assembly.fa").map { |f| f.sub(/_assembly\.fa$/, '') }

def bwa_align_outputs
  Rake::FileList.new do |fl|
    SOURCE_BASES.pathmap("%{^inputs/,outputs/}p").each do |f|
      fl.add("#{f}_assembly_close.bam")
      fl.add("#{f}_assembly_close.bai")
      fl.add("#{f}_assembly_close.wig")
      fl.add("#{f}_assembly_close.bigwig")
    end
  end
end

def mpileup_outputs
  Rake::FileList.new do |fl|
    SOURCE_BASES.pathmap("%{^inputs/,outputs/}p").each do |f|
      fl.add("#{f}_assembly_close_ilm.vcf")
      fl.add("#{f}_assembly_close_ilm.fasta")
    end
  end
end

task :default => :mpileup

desc "Performs BWA Align on QCed Illumina reads + a closed PacBio assembly"
task :bwa_align => [:check] + bwa_align_outputs

desc "Performs mpileup on BWA Align output to produces VCF and FASTA"
task :mpileup => [:check] + mpileup_outputs

desc "Create some test files just to get started"
task :test do
  mkdir_p "inputs"
  touch "inputs/123_123_assembly.fa"
  touch "inputs/123_123_illumina.fastq"
end

directory "outputs"

desc "Checks that everything is good to go"
task :check do |t|
  unless Dir.exists? "inputs"
    fail red <<-FAIL.gsub(/^\s+/, '')
      You need to symlink or create inputs/ as a directory of Illumina FASTQ and PacBio FA inputs.
      As an example, run `rake test`, which will create this directory and some sample test files.
    FAIL
  end
end

# QC Filter
rule(/_QC\.fastq$/ => proc { |f| f.sub(/^outputs\//, 'inputs/').sub(/_QC\.fastq$/, '.fastq') }) do |t|
  mkdir_p t.name.pathmap("%d")
  cp t.source, t.name
end

# BWA Align
BWA_ALIGN_REQS = proc {|f| [f.ext("fa"), f.sub(/_assembly_close\.\w+$/, '_illumina_QC.fastq')] }
rule(/_assembly_close\.(bam)$/ => BWA_ALIGN_REQS) do |t|
  mkdir_p t.name.pathmap("%d")
  puts "Here, I'd make #{t.name} from: \n  - #{t.sources.join("\n  - ")}"
  touch t.name
end
rule(/_assembly_close\.(bai|wig|bigwig)/ => proc { |f| f.ext('.bam') }) do |t|
  mkdir_p t.name.pathmap("%d")
  puts "Here, I'd make #{t.name} from: \n  - #{t.sources.join("\n  - ")}"
  touch t.name
end

# mpileup recall consensus
MPILEUP_REQS = proc {|f| base = f.sub(/_ilm\.(vcf|fasta)$/, ''); ["#{base}.bam", "#{base}.bai"] }
rule(/_assembly_close_ilm\.(vcf|fasta)$/ => MPILEUP_REQS) do |t|
  mkdir_p t.name.pathmap("%d")
  puts "Here, I'd make #{t.name} from: \n  - #{t.sources.join("\n  - ")}"
  touch t.name
end

# Circularize + finish
CIRCULARIZE_REQS = proc {|f| f.sub(/^outputs\//, 'inputs/').sub(/_close\.fa$/, '.fa') }
rule(/_assembly_close\.fa$/ => CIRCULARIZE_REQS) do |t|
  mkdir_p t.name.pathmap("%d")
  cp t.source, t.name
end

def pull_down_raw_reads(jobid, dir1)
#     puts Dir.getwd
     Dir.chdir("#{dir1}")
#     puts Dir.getwd
     system("source /hpc/packages/minerva-mothra/smrtanalysis/2.2.0/ROOT/current/etc/setup.sh")
     system("python /hpc/users/attieo02/gitrepos/multiscale/bacterial_analysis/src/ccs_get.py --noprefix -e bax.h5 #{jobid} -i ")
     system("ls *bax.h5 > bash5.fofn")
end

def assemble_raw_reads(dir1)
    system("source /hpc/packages/minerva-mothra/smrtanalysis/2.2.0/ROOT/current/etc/setup.sh; fofnToSmrtpipeInput.py bash5.fofn > bash5.xml; cp /sc/orga/projects/InfectiousDisease/smrtpipe/example_params.xml \.; smrtpipe.py -D NPROC=16 -D CLUSTER=BASH -D MAX_THREADS=16 --params example_params.xml xml:bash5.xml -D TMP=/sc/orga/scratch/attieo02/; ")
end
 
def circularize_assembly
    system("gunzip data/polished_assembly.fasta.gz")
    system("export PERL5LIB=/usr/bin/perl5.10.1/;/sc/orga/work/attieo02/circularizeContig.pl data/polished_assembly.fasta")
end

def resequence_assembly(strain_name,dir1, species)
    system("mkdir circularized_sequence")
    system("source /hpc/packages/minerva-mothra/smrtanalysis/2.2.0/ROOT/current/etc/setup.sh; referenceUploader -c -p #{dir1}/circularized_sequence -n #{strain_name} -f data/polished_assembly_circularized.fasta")
system("cp /sc/orga/projects/InfectiousDisease/smrtpipe/resequence_example_params.xml \.")
system("export PERL5LIB=/usr/bin/perl5.10.1/; perl /sc/orga/work/attieo02/changeResequencingDirectory.pl resequence_example_params.xml #{dir1} circularized_sequence/#{strain_name} > resequence_params.xml")
    system("source /hpc/packages/minerva-mothra/smrtanalysis/2.2.0/ROOT/current/etc/setup.sh; samtools faidx circularized_sequence/#{strain_name}/sequence/#{strain_name}.fasta; smrtpipe.py -D NPROC=16 -D CLUSTER=BASH -D MAX_THREADS=16 --params resequence_params.xml xml:bash5.xml -D TMP=/sc/orga/scratch/attieo02/")
    system("gunzip data/consensus.fasta.gz")
    system("cp data/consensus.fasta data/#{strain_name}_consensus.fasta")
    rast_job=%x[export PERL5LIB=$PERL5LIB:/sc/orga/work/attieo02/sas/lib:/sc/orga/work/attieo02/sas/modules/lib;export PATH=$PATH:/sc/orga/work/attieo02/sas/bin;perl /sc/orga/work/attieo02/sas/plbin/svr_submit_status_retrieve.pl --user oattie --passwd sessiz_ev --fasta #{dir1}/data/#{strain_name}_consensus.fasta --domain Bacteria --bioname "#{species} #{strain_name}" --genetic_code 11 --gene_caller rast]
    system("export PERL5LIB=$PERL5LIB:/sc/orga/work/attieo02/sas/lib:/sc/orga/work/attieo02/sas/modules/lib;export PATH=$PATH:/sc/orga/work/attieo02/sas/bin; perl /sc/orga/work/attieo02/sas/test_server.pl oattie sessiz_ev genbank #{rast_job} ")
    sleep(120)
    system("export PERL5LIB=$PERL5LIB:/sc/orga/work/attieo02/sas/lib:/sc/orga/work/attieo02/sas/modules/lib;export PATH=$PATH:/sc/orga/work/attieo02/sas/bin; svr_retrieve_RAST_job oattie sessiz_ev #{rast_job} genbank > #{dir1}/data/#{strain_name}_consensus_rast.gbk")
    system("python /sc/orga/work/attieo02/genbank_to_fasta_v1.1/gb_to_fasta.py -i #{dir1}/data/#{strain_name}_consensus_rast.gbk -s \'aa\' -o #{dir1}/data/#{strain_name}_consensus_rast_aa.fa")
    system("python /sc/orga/work/attieo02/genbank_to_fasta_v1.1/gb_to_fasta.py -i #{dir1}/data/#{strain_name}_consensus_rast.gbk -s \'nt\' -o #{dir1}/data/#strain_name}_consensus_rast.fna")
end

task :resequence, :job_id, :dir1, :strain, :species do |t, args|
     job_id=args[:job_id]
     dir1=args[:dir1]
     strain_name=args[:strain]
     species=args[:species]
      pull_down_raw_reads(job_id,dir1)
      assemble_raw_reads(dir1)
      circularize_assembly
      resequence_assembly(strain_name, dir1, species)
#     resequence_assembly
end

task :circularize do
     circularize_assembly
end
 
#task :circularize=>:assemble do
#     circularize_assembly
#end

#task :assemble => :pull do
#     assemble_raw_reads
#end

#task :pull, :job_id, :dir1 do |t, args|
#     job_id=args[:job_id]
#     dir1=args[:dir1]
#     pull_down_raw_reads(job_id,dir1)
#end
