require 'pp'
require_relative 'lib/colors'
require_relative 'lib/lsf_client'
require 'shellwords'
include Colors

task :default => :check

LSF = LSFClient.new

REPO_DIR = File.dirname(__FILE__)
SAS_DIR = "#{REPO_DIR}/vendor/sas"
MUMMER_DIR = "#{REPO_DIR}/vendor/MUMmer3.23"

OUT = ENV['OUT'] || "#{REPO_DIR}/out"

#######
# Other environment variables that may be set by the user for specific tasks (see README.md)
#######
STRAIN_NAME = ENV['STRAIN_NAME']
SPECIES = ENV['SPECIES']
ILLUMINA_FASTQ = ENV['ILLUMINA_FASTQ']

#############################################################
#  IMPORTANT!
#  This Rakefile runs with the working directory set to OUT
#  All filenames from hereon are relative to that directory
#############################################################
Dir.chdir(OUT)

task :env do
  puts "Output directory: #{OUT}"
  mkdir_p OUT
  mkdir_p File.join(REPO_DIR, "vendor")
  
  sc_orga_scratch = "/sc/orga/scratch/#{ENV['USER']}"
  ENV['TMP'] ||= Dir.exists?(sc_orga_scratch) ? sc_orga_scratch : "/tmp"
  ENV['PERL5LIB'] ||= "/usr/bin/perl5.10.1"
end

file "#{REPO_DIR}/scripts/env.sh" => "#{REPO_DIR}/scripts/env.example.sh" do
  cp "#{REPO_DIR}/scripts/env.example.sh", "#{REPO_DIR}/scripts/env.sh"
end

ENV_ERROR = "Configure this in scripts/env.sh and run `source scripts/env.sh` before running rake."

desc "Checks environment variables and requirements before running tasks"
task :check => [:env, "#{REPO_DIR}/scripts/env.sh", :sas, :mummer] do
  unless `module avail 2>&1 | grep smrtpipe/2.2.0` != ''
    abort "FATAL: You must have the smrtpipe/2.2.0 module in your MODULEPATH."
  end
  unless ENV['SMRTPIPE'] && File.exists?("#{ENV['SMRTPIPE']}/example_params.xml")
    abort "FATAL: SMRTPIPE must be set to the directory containing example_params.xml for smrtpipe.py.\n#{ENV_ERROR}"
  end
  unless ENV['SMRTANALYSIS'] && File.exists?("#{ENV['SMRTANALYSIS']}/etc/setup.sh")
    abort <<-ERRMSG
      FATAL: SMRTANALYSIS must be set to the ROOT directory for the SMRT Analysis package, v2.3.0.
      This software can be downloaded from http://www.pacb.com/devnet/
      #{ENV_ERROR}"
    ERRMSG
  end
end

# pulls down http://blog.theseed.org/downloads/sas.tgz --> ./vendor/sas
#   then it adds SAS libs to PERL5LIB
#   then it adds SAS bins to PATH
task :sas => [:env, "#{SAS_DIR}/sas.tgz", "#{SAS_DIR}/modules/lib"] do
  ENV['PERL5LIB'] = "#{ENV['PERL5LIB']}:#{SAS_DIR}/lib:#{SAS_DIR}/modules/lib"
  ENV['PATH'] = "#{SAS_DIR}/bin:#{ENV['PATH']}"
  
  unless ENV['RAST_USER'] && ENV['RAST_USER'] != ''
    abort "FATAL: RAST_USER must be set to your username for http://rast.nmpdr.org/\n#{ENV_ERROR}"
  end
  unless ENV['RAST_PASSWORD'] && ENV['RAST_PASSWORD'] != ''
    abort "FATAL: RAST_PASSWORD must be set to your password for http://rast.nmpdr.org/\n#{ENV_ERROR}"
  end
end

directory SAS_DIR
file "#{SAS_DIR}/sas.tgz" => [SAS_DIR] do |t|
  Dir.chdir(File.dirname(t.name)) do
    system "curl -L -O 'http://blog.theseed.org/downloads/sas.tgz'" and
    system "tar xvzf sas.tgz"
  end
end

directory "#{SAS_DIR}/modules/lib"
file "#{SAS_DIR}/modules/lib" => ["#{SAS_DIR}/sas.tgz"] do |t|
  Dir.chdir("#{SAS_DIR}/modules") do
    system "./BUILD_MODULES"
  end
end

# pulls down and compiles MUMmer 3.23, which is used by scripts/circulizeContig.pl and others
# see http://mummer.sourceforge.net/
task :mummer => [:env, MUMMER_DIR, "#{MUMMER_DIR}/nucmer", "#{MUMMER_DIR}/show-coords"]
directory MUMMER_DIR
file "#{MUMMER_DIR}/show-coords" => "#{MUMMER_DIR}/nucmer"
file "#{MUMMER_DIR}/nucmer" do
  Dir.chdir(File.dirname(MUMMER_DIR)) do
    system <<-SH
      curl -L -o mummer.tar.gz 'http://sourceforge.net/projects/mummer/files/mummer/3.23/MUMmer3.23.tar.gz/download'
      tar xvzf mummer.tar.gz  # Creates MUMMER_DIR
    SH
  end
  Dir.chdir(MUMMER_DIR) { system "make install" }
end

file "pathogendb-pipeline.png" => [:graph]
desc "Generates a graph of tasks, intermediate files and their dependencies from this Rakefile"
task :graph do
  system <<-SH
    module load graphviz
    STRAIN_NAME=STRAIN_NAME rake -f #{Shellwords.escape(__FILE__)} -P \
        | #{REPO_DIR}/scripts/rake-prereqs-dot.rb --prune #{REPO_DIR} --replace-with REPO_DIR \
        | dot -Tpng -o pathogendb-pipeline.png
  SH
end


# =======================
# = pull_down_raw_reads =
# =======================

desc "Uses scripts/ccs_get.py to save raw reads from PacBio to OUT directory"
task :pull_down_raw_reads => [:check, "bash5.fofn"]  # <-- file(s) created by this task
file "bash5.fofn" do |t, args|                       # <-- implementation for generating each of these files
  job_id = ENV['SMRT_JOB_ID'] # an example that works is 019194
  abort "FATAL: Task pull_down_raw_reads requires specifying SMRT_JOB_ID" unless job_id
  
  system <<-SH
    module load python/2.7.6
    python #{REPO_DIR}/scripts/ccs_get.py --noprefix -e bax.h5 #{job_id} -i &&
    find #{OUT}/*bax.h5 > bash5.fofn
  SH
  # NOTE: we will change the above to not fetch the full sequence, but rather symlink to it on minerva, like so
  #       we could even skip straight to circularize_assembly if the polished_assembly.fasta is already there
  # cp /sc/orga/projects/pacbio/userdata_permanent/jobs/#{job_id[0..3]}/#{job_id}/input.fofn baxh5.fofn
  # mkdir_p "data"
  # ln -s /sc/orga/projects/pacbio/userdata_permanent/jobs/data/#{job_id[0..3]}/#{job_id}/polished_assembly.fasta \
  #       <input assembly file name>
end

# ======================
# = assemble_raw_reads =
# ======================

desc "Uses smrtpipe.py to assemble raw reads from PacBio within OUT directory"
task :assemble_raw_reads => [:check, "data/polished_assembly.fasta.gz"]
file "data/polished_assembly.fasta.gz" => "bash5.fofn" do |t|
  system <<-SH
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    fofnToSmrtpipeInput.py bash5.fofn > bash5.xml &&
    cp #{ENV['SMRTPIPE']}/example_params.xml \. &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=16 -D CLUSTER=LSF \
        -D MAX_THREADS=16 --distribute --params example_params.xml xml:bash5.xml 
  SH
end

# ========================
# = circularize_assembly =
# ========================

desc "Circularizes the PacBio assembly"
task :circularize_assembly => [:check, "data/polished_assembly_circularized.fasta"]
file "data/polished_assembly_circularized.fasta" => "data/polished_assembly.fasta.gz" do |t|
  system "gunzip -c data/polished_assembly.fasta.gz >data/polished_assembly.fasta" and
  system "#{REPO_DIR}/scripts/circularizeContig.pl data/polished_assembly.fasta"
end

# =======================
# = resequence_assembly =
# =======================

desc "Resequences the circularized assembly"
task :resequence_assembly => [:check, "data/#{STRAIN_NAME}_consensus.fasta"]
file "data/#{STRAIN_NAME}_consensus.fasta" => "data/polished_assembly_circularized.fasta" do |t|
  abort "FATAL: Task resequence_assembly requires specifying STRAIN_NAME" unless STRAIN_NAME 
  
  mkdir_p "circularized_sequence"
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    referenceUploader -c -p circularized_sequence -n #{STRAIN_NAME} -f data/polished_assembly_circularized.fasta
  SH
  cp "#{ENV['SMRTPIPE']}/resequence_example_params.xml", OUT
  system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl resequence_example_params.xml " +
      "#{OUT} circularized_sequence/#{STRAIN_NAME} > resequence_params.xml" and
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    samtools faidx circularized_sequence/#{STRAIN_NAME}/sequence/#{STRAIN_NAME}.fasta &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=16 -D CLUSTER=LSF \
        -D MAX_THREADS=16 --distribute --params resequence_params.xml xml:bash5.xml &&
    gunzip data/consensus.fasta.gz
  SH
  cp "data/consensus.fasta", "data/#{STRAIN_NAME}_consensus.fasta"
end


# =================
# = rast_annotate =
# =================

desc "Submits the circularized assembly to RAST for annotations"
task :rast_annotate => [:check, "data/#{STRAIN_NAME}_consensus_rast.fna", 
    "data/#{STRAIN_NAME}_consensus_rast.gbk", "data/#{STRAIN_NAME}_consensus_rast_aa.fa"]

file "data/#{STRAIN_NAME}_consensus_rast.gbk" => ["data/#{STRAIN_NAME}_consensus.fasta"] do |t|
  abort "FATAL: Task rast_annotate requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task rast_annotate requires specifying SPECIES" unless SPECIES 
  
  rast_job = %x[
    perl #{REPO_DIR}/scripts/svr_submit_status_retrieve.pl --user #{ENV['RAST_USER']} \
        --passwd #{ENV['RAST_PASSWORD']} --fasta data/#{STRAIN_NAME}_consensus.fasta --domain Bacteria \
        --bioname "#{SPECIES} #{STRAIN_NAME}" --genetic_code 11 --gene_caller rast
  ]
  system "perl #{REPO_DIR}/scripts/test_server.pl #{ENV['RAST_USER']} #{ENV['RAST_PASSWORD']} genbank #{rast_job}"
  sleep 120
  loop do
    success = system <<-SH
      svr_retrieve_RAST_job #{ENV['RAST_USER']} #{ENV['RAST_PASSWORD']} #{rast_job} genbank \
          > data/#{STRAIN_NAME}_consensus_rast.gbk
    SH
    break if success
    puts "RAST output not available yet, retrying..."
    sleep 60
  end
end

file "data/#{STRAIN_NAME}_consensus_rast_aa.fa" => "data/#{STRAIN_NAME}_consensus_rast.gbk" do |t|
  abort "FATAL: Task rast_annotate requires specifying STRAIN_NAME" unless STRAIN_NAME 
  system <<-SH
    module load python/2.7.6
    module load py_packages/2.7
    python #{REPO_DIR}/scripts/gb_to_fasta.py -i data/#{STRAIN_NAME}_consensus_rast.gbk -s aa \
        -o #{OUT}/data/#{STRAIN_NAME}_consensus_rast_aa.fa
  SH
end

file "data/#{STRAIN_NAME}_consensus_rast.fna" => "data/#{STRAIN_NAME}_consensus_rast.gbk" do |t|
  abort "FATAL: Task rast_annotate requires specifying STRAIN_NAME" unless STRAIN_NAME 
  system <<-SH
    module load python/2.7.6
    module load py_packages/2.7
    python #{REPO_DIR}/scripts/gb_to_fasta.py -i data/#{STRAIN_NAME}_consensus_rast.gbk -s nt \
        -o #{OUT}/data/#{STRAIN_NAME}_consensus_rast.fna
  SH
end


# ====================
# = recall_consensus =
# ====================

desc "Recalls a consensus for SNPs by piling Illumina reads onto a PacBio assembly"
task :recall_consensus => [:check, "data/#{STRAIN_NAME}_ref_flt.vcf"]

file "data/ref.aln.sam" => "data/#{STRAIN_NAME}_consensus.fasta" do |t|
  abort "FATAL: Task recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task recall_consensus requires specifying ILLUMINA_FASTQ" unless ILLUMINA_FASTQ
  
  LSF.set_out_err("log/recall_consensus.log", "log/recall_consensus.err.log")
  LSF.job_name "ref.aln.sam"
  LSF.bsub_interactive <<-SH
    module load bwa/0.7.8
    bwa index "data/#{STRAIN_NAME}_consensus.fasta"
    bwa mem "data/#{STRAIN_NAME}_consensus.fasta" #{Shellwords.escape(ILLUMINA_FASTQ)} > data/ref.aln.sam
  SH
end

file "data/ref.sort.bam" => "data/ref.aln.sam" do |t|
  # Use samtools to convert ref.aln.sam into a .bam file and sort the .bam file
  # This produces ref.sort.bam.
  LSF.set_out_err("log/recall_consensus.log", "log/recall_consensus.err.log")
  LSF.job_name "ref.sort.bam"
  LSF.bsub_interactive <<-SH
    module load samtools/1.1
    samtools view -bS data/ref.aln.sam > data/ref.aln.bam
    samtools sort data/ref.aln.bam data/ref.sort
  SH
end

file "data/#{STRAIN_NAME}_ref_raw.bcf" => "data/ref.sort.bam" do |t|
  abort "FATAL: Task recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
  # Use mpileup to do the consensus calling.
  LSF.set_out_err("log/recall_consensus.log", "log/recall_consensus.err.log")
  LSF.job_name "#{STRAIN_NAME}_ref_raw.bcf"
  LSF.bsub_interactive <<-SH
    module load samtools/1.1
    module load bcftools/1.1
    samtools mpileup -uf "data/#{STRAIN_NAME}_consensus.fasta" data/ref.sort.bam \
        | bcftools call -cv -Ob > "data/#{STRAIN_NAME}_ref_raw.bcf"
  SH
end

file "data/#{STRAIN_NAME}_ref_flt.vcf" => "data/#{STRAIN_NAME}_ref_raw.bcf" do |t|
  LSF.set_out_err("log/recall_consensus.log", "log/recall_consensus.err.log")
  LSF.job_name "#{STRAIN_NAME}_ref_flt.vcf"
  LSF.bsub_interactive <<-SH
    module load samtools/1.1
    module load bcftools/1.1
    bcftools view "data/#{STRAIN_NAME}_ref_raw.bcf" | vcfutils.pl varFilter > "data/#{STRAIN_NAME}_ref_flt.vcf"
  SH
end

desc "Fakes the prerequisites for the recall_consensus task"
task :recall_consensus_fake_prereqs do
  abort "FATAL: Task recall_consensus_fake_prereqs requires specifying STRAIN_NAME" unless STRAIN_NAME 
  mkdir_p "log"
  touch "bash5.fofn"                                  and sleep 1
  touch "data/polished_assembly.fasta.gz"             and sleep 1
  touch "data/polished_assembly_circularized.fasta"   and sleep 1
  touch "data/#{STRAIN_NAME}_consensus.fasta"
end