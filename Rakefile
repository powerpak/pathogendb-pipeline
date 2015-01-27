require 'pp'
require 'net/http'
require_relative 'lib/colors'
require_relative 'lib/lsf_client'
require_relative 'lib/subscreens'
require 'shellwords'
include Colors

task :default => :check

LSF = LSFClient.new

REPO_DIR = File.dirname(__FILE__)
SAS_DIR = "#{REPO_DIR}/vendor/sas"
MUMMER_DIR = "#{REPO_DIR}/vendor/MUMmer3.23"
BCFTOOLS_DIR = "#{REPO_DIR}/vendor/bcftools"
HTSLIB_DIR = "#{REPO_DIR}/vendor/htslib"

IGB_DIR = ENV['IGB_DIR'] || "#{ENV['HOME']}/www/igb"

OUT = File.expand_path(ENV['OUT'] || "#{REPO_DIR}/out")

#######
# Other environment variables that may be set by the user for specific tasks (see README.md)
#######
STRAIN_NAME = ENV['STRAIN_NAME']
SPECIES = ENV['SPECIES']
ILLUMINA_FASTQ = ENV['ILLUMINA_FASTQ']
TASK_FILE = ENV['TASK_FILE']

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
  # Always use our locally bundled (patched) perl modules over anything else
  ENV['PERL5LIB'] = "#{REPO_DIR}/lib/perl:#{ENV['PERL5LIB']}"
end

file "#{REPO_DIR}/scripts/env.sh" => "#{REPO_DIR}/scripts/example.env.sh" do
  cp "#{REPO_DIR}/scripts/example.env.sh", "#{REPO_DIR}/scripts/env.sh"
end

ENV_ERROR = "Configure this in scripts/env.sh and run `source scripts/env.sh` before running rake."

desc "Checks environment variables and requirements before running tasks"
task :check => [:env, "#{REPO_DIR}/scripts/env.sh", :sas, :mummer, :bcftools] do
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
  mkdir_p ENV['TMP'] or abort "FATAL: set TMP to a directory that can store scratch files"
  mkdir_p ENV['SHARED_DIR'] or abort "FATAL: set SHARED_DIR to a directory that can store scratch files"
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

# pulls down and compiles MUMmer 3.23, which is used by scripts/circularizeContig.pl and others
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

task :bcftools => [:env, "#{BCFTOOLS_DIR}/bcftools", HTSLIB_DIR, "#{HTSLIB_DIR}/bgzip", "#{HTSLIB_DIR}/tabix"]
directory BCFTOOLS_DIR
file "#{BCFTOOLS_DIR}/bcftools" do
  Dir.chdir(File.dirname(BCFTOOLS_DIR)) do
    system <<-SH
      git clone --branch=develop git://github.com/samtools/htslib.git
      git clone --branch=develop git://github.com/samtools/samtools.git
      git clone --branch=develop git://github.com/samtools/bcftools.git
    SH
  end
  Dir.chdir(BCFTOOLS_DIR) { system "make" }
  Dir.chdir("#{File.dirname(BCFTOOLS_DIR)}/samtools") { system "make" }
end
directory HTSLIB_DIR
file "#{HTSLIB_DIR}/bgzip" => "#{BCFTOOLS_DIR}/bcftools" do
  Dir.chdir(HTSLIB_DIR) { system "make" }
end
file "#{HTSLIB_DIR}/tabix" => "#{HTSLIB_DIR}/bgzip"

file "pathogendb-pipeline.png" => [:graph]
desc "Generates a graph of tasks, intermediate files and their dependencies from this Rakefile"
task :graph do
  system <<-SH
    module load graphviz
    STRAIN_NAME='${STRAIN_NAME}' SPECIES='${SPECIES}' rake -f #{Shellwords.escape(__FILE__)} -P \
        | #{REPO_DIR}/scripts/rake-prereqs-dot.rb --prune #{REPO_DIR} --replace-with '$REPO_DIR' \
        | unflatten -f -l5 -c 3 \
        | dot -Tpng -o pathogendb-pipeline.png
  SH
end


desc "Runs this pipeline in $n separate screens, doing all the tasks in $task_file"
task :multi, [:task_file, :n] => [:check] do |t, args|
  abort "FATAL: Task multi requires specifying task_file argument" unless args[:task_file]
  
  i = 0
  n = args[:n] && args[:n].to_i
  cmds = []
  File.open(args[:task_file], 'r') do |f|
    f.each_line do |line|
      next if line =~ /^\s*#/
      cmds[i] ||= "source scripts/env.sh"
      cmds[i] += "\nrake #{line.strip}"
      i = (n && n != 0) ? ((i + 1) % n) : (i + 1)
    end
  end
  puts "Creating screen session `nested` with #{cmds.size} windows."
  sleep 2
  Dir.chdir(REPO_DIR) do
    Subscreens.run(cmds)
  end
end


# =======================
# = pull_down_raw_reads =
# =======================

desc "Copies or downloads raw reads from a PacBio job to the OUT directory"
task :pull_down_raw_reads => [:check, "bash5.fofn"]  # <-- file(s) created by this task
file "bash5.fofn" do |t, args|                       # <-- implementation for generating each of these files
  job_id = ENV['SMRT_JOB_ID'] # Examples that work are: 019194, 020266
  abort "FATAL: Task pull_down_raw_reads requires specifying SMRT_JOB_ID" unless job_id
  job_id = job_id.rjust(6, '0')
  pacbio_job_dirs = ["/sc/orga/projects/pacbio/userdata_permanent/jobs/#{job_id[0..2]}/#{job_id}",
                     "/sc/orga/projects/InfectiousDisease/old_smrtportal_jobs/#{job_id}"]
  smrtpipe_log_url = "http://node1.1425mad.mssm.edu/pacbio/secondary/#{job_id[0..2]}/#{job_id}/log/smrtpipe.log"
  
  found_fofn_dir = pacbio_job_dirs.find {|dir| File.exist? "#{dir}/input.fofn" }
  if found_fofn_dir
    cp "#{found_fofn_dir}/input.fofn", "bash5.fofn"
    mkdir_p "data"
    if File.exist? "#{found_fofn_dir}/data/polished_assembly.fasta.gz"
      ln_s "#{found_fofn_dir}/data/polished_assembly.fasta.gz", "data/polished_assembly.fasta.gz"
    end
  else
    url = URI.parse(smrtpipe_log_url)
    req = Net::HTTP.new(url.host, url.port)
    res = req.request_head(url.path)
    unless res.code == "200"
      abort "FATAL: Task pull_down_raw_reads cannot find a way to fetch reads for that SMRT_JOB_ID"
    end
    
    puts "<< Fetching reads with ccs_get.py >>"
    system <<-SH
      module load python/2.7.6
      python #{REPO_DIR}/scripts/ccs_get.py --noprefix -e bax.h5 #{job_id} -i &&
      find #{OUT}/*bax.h5 > bash5.fofn
    SH
  end
end

# ======================
# = assemble_raw_reads =
# ======================

desc "Uses smrtpipe.py to assemble raw reads from PacBio within OUT directory"
task :assemble_raw_reads => [:check, "data/polished_assembly.fasta.gz"]
file "data/polished_assembly.fasta.gz" => "bash5.fofn" do |t|
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source "#{ENV['SMRTANALYSIS']}/etc/setup.sh" &&
    fofnToSmrtpipeInput.py bash5.fofn > bash5.xml
  SH
  
  lstat = File::lstat("data/polished_assembly.fasta.gz") rescue nil
  if lstat and lstat.symlink?
    puts "NOTICE: polished_assembly.fasta.gz is symlinked to an existing assembly, skipping assemble_raw_reads"
    next
  end
  
  cp "#{ENV['SMRTPIPE']}/example_params.xml", "."
  system <<-SH
    module load smrtpipe/2.2.0
    source "#{ENV['SMRTANALYSIS']}/etc/setup.sh" &&
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
  system "#{REPO_DIR}/scripts/circularizeContigs.pl -i data/polished_assembly.fasta"
end

# =======================
# = resequence_assembly =
# =======================

desc "Resequences the circularized assembly"
task :resequence_assembly => [:check, "data/#{STRAIN_NAME}_consensus.fasta"]
file "data/#{STRAIN_NAME}_consensus.fasta" => "data/polished_assembly_circularized.fasta" do |t|
  abort "FATAL: Task resequence_assembly requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: STRAIN_NAME can only contain letters, numbers, and underscores" unless STRAIN_NAME =~ /^[\w]+$/
  
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
    "data/#{STRAIN_NAME}_consensus_rast.gbk", "data/#{STRAIN_NAME}_consensus_rast_aa.fa",
    "data/rast_job_id"]

def submit_and_retrieve_rast(fasta, gbk_file, job_id_file="rast_job_id", task_name="rast_annotate") 
  abort "FATAL: Task #{task_name} requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task #{task_name} requires specifying SPECIES" unless SPECIES 
  
  if File.exist? "data/#{job_id_file}"
    rast_job = IO.read("data/#{job_id_file}").strip
  else
    rast_job = %x[
      perl #{REPO_DIR}/scripts/svr_submit_status_retrieve.pl --user #{Shellwords.escape ENV['RAST_USER']} \
          --passwd #{Shellwords.escape ENV['RAST_PASSWORD']} --fasta #{fasta} --domain Bacteria \
          --bioname "#{SPECIES} #{STRAIN_NAME}" --genetic_code 11 --gene_caller rast
    ]
    IO.write("data/#{job_id_file}", rast_job)
  end
  system <<-SH
    perl #{REPO_DIR}/scripts/test_server.pl #{Shellwords.escape ENV['RAST_USER']} \
        #{Shellwords.escape ENV['RAST_PASSWORD']} genbank #{rast_job}
  SH
  sleep 120
  loop do
    success = system <<-SH
      perl #{SAS_DIR}/plbin/svr_retrieve_RAST_job.pl #{Shellwords.escape ENV['RAST_USER']} \
          #{Shellwords.escape ENV['RAST_PASSWORD']} #{rast_job} genbank > #{gbk_file}
    SH
    break if success
    puts "RAST output not available yet, retrying..."
    sleep 60
  end
end

file "data/#{STRAIN_NAME}_consensus_rast.gbk" => ["data/#{STRAIN_NAME}_consensus.fasta"] do |t|
  submit_and_retrieve_rast("data/#{STRAIN_NAME}_consensus.fasta", t.name)
end
file "data/rast_job_id" => "data/#{STRAIN_NAME}_consensus_rast.gbk"

def gb_to_fasta(gb, fasta, seq_type=:nt, task_name="rast_annotate")
  abort "FATAL: Task #{task_name} requires specifying STRAIN_NAME" unless STRAIN_NAME
  abort "FATAL: gb_to_fasta called with invalid seq_type" unless [:nt, :aa].include? seq_type
  system <<-SH
    module load python/2.7.6
    module load py_packages/2.7
    python #{REPO_DIR}/scripts/gb_to_fasta.py -i #{gb} -s #{seq_type} -o #{fasta}
  SH
end

file "data/#{STRAIN_NAME}_consensus_rast_aa.fa" => "data/#{STRAIN_NAME}_consensus_rast.gbk" do |t|
  gb_to_fasta "data/#{STRAIN_NAME}_consensus_rast.gbk", "#{OUT}/#{t.name}", :aa
end

file "data/#{STRAIN_NAME}_consensus_rast.fna" => "data/#{STRAIN_NAME}_consensus_rast.gbk" do |t|
  gb_to_fasta "data/#{STRAIN_NAME}_consensus_rast.gbk", "#{OUT}/#{t.name}", :nt
end

# ===============
# = rast_to_igb =
# ===============

species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES
strain_igb_dir = "#{IGB_DIR}/#{species_clean}_#{STRAIN_NAME}"

desc "Creates an IGB Quickload-compatible directory for your genome in IGB_DIR"
task :rast_to_igb => [:check, strain_igb_dir]

directory IGB_DIR
file strain_igb_dir => [IGB_DIR, "data/rast_job_id"] do |t|
  abort "FATAL: Task rast_to_igb requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task rast_to_igb requires specifying SPECIES" unless SPECIES 
  rast_job = IO.read("data/rast_job_id").strip
  
  system <<-SH
    module load blat/3.0.5
    export SAS_DIR=#{SAS_DIR}
    perl #{REPO_DIR}/scripts/rast2igb.pl -u #{Shellwords.escape ENV['RAST_USER']} \
        -p #{Shellwords.escape ENV['RAST_PASSWORD']} -j #{rast_job} -g #{species_clean}_#{STRAIN_NAME} \
        -i #{IGB_DIR}
  SH
end


# ========================
# = recall_ilm_consensus =
# ========================

desc "Recalls a new consensus by piling Illumina reads onto a PacBio assembly"
task :recall_ilm_consensus => [:check, "data/#{STRAIN_NAME}_ref_flt.vcf", "data/#{STRAIN_NAME}_ilm_consensus.fasta"]

file "data/ref.sort.bam" => "data/#{STRAIN_NAME}_consensus.fasta" do |t|
  abort "FATAL: Task recall_ilm_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task recall_ilm_consensus requires specifying ILLUMINA_FASTQ" unless ILLUMINA_FASTQ
  
  LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
  LSF.job_name "ref.aln.sam"
  LSF.bsub_interactive <<-SH or abort
    module load bwa/0.7.8
    bwa index "data/#{STRAIN_NAME}_consensus.fasta"
    bwa mem "data/#{STRAIN_NAME}_consensus.fasta" #{Shellwords.escape(ILLUMINA_FASTQ)} > data/ref.aln.sam
    
    module load samtools/1.1
    samtools view -bS data/ref.aln.sam > data/ref.aln.bam
    samtools sort data/ref.aln.bam data/ref.sort
  SH
  
  # Can remove the SAM as it is huge and the .aln.bam should contain everything in it
  rm "data/ref.aln.sam"
end

file "data/#{STRAIN_NAME}_ref_raw.bcf" => "data/ref.sort.bam" do |t|
  abort "FATAL: Task recall_ilm_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
  # Use mpileup to do the consensus calling.
  LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
  LSF.job_name "#{STRAIN_NAME}_ref_raw.bcf"
  LSF.bsub_interactive <<-SH
    module load samtools/1.1
    module load bcftools/1.1
    samtools mpileup -uf "data/#{STRAIN_NAME}_consensus.fasta" data/ref.sort.bam \
        | bcftools call -cv -Ob > "data/#{STRAIN_NAME}_ref_raw.bcf"
  SH
end

file "data/#{STRAIN_NAME}_ref_flt.vcf" => "data/#{STRAIN_NAME}_ref_raw.bcf" do |t|
  LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
  LSF.job_name "#{STRAIN_NAME}_ref_flt.vcf"
  LSF.bsub_interactive <<-SH
    module load samtools/1.1
    module load bcftools/1.1
    bcftools view "data/#{STRAIN_NAME}_ref_raw.bcf" | vcfutils.pl varFilter > "data/#{STRAIN_NAME}_ref_flt.vcf"
  SH
end

file "data/#{STRAIN_NAME}_ilm_consensus.fasta" => 
    ["data/#{STRAIN_NAME}_ref_flt.vcf", "data/#{STRAIN_NAME}_consensus.fasta"] do |t|
  LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
  LSF.job_name "#{STRAIN_NAME}_ilm_consensus.fasta"
  system <<-SH
    module load vcftools/0.1.12b
    module load tabix/0.2.6 
    
    bgzip -c "data/#{STRAIN_NAME}_ref_flt.vcf" > "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
    tabix -p vcf "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
    cat "data/#{STRAIN_NAME}_consensus.fasta" | vcf-consensus "data/#{STRAIN_NAME}_ref_flt.vcf.gz" \
            > "data/#{STRAIN_NAME}_ilm_consensus.fasta"
    
    # New-style version of doing this with bcftools consensus, but it doesn't work (memory leak in bcftools)
    # #{HTSLIB_DIR}/bgzip -c "data/#{STRAIN_NAME}_ref_flt.vcf" > "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
    # #{HTSLIB_DIR}/tabix -p vcf "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
    # #{BCFTOOLS_DIR}/bcftools consensus -f "data/#{STRAIN_NAME}_consensus.fasta" "data/#{STRAIN_NAME}_ref_flt.vcf.gz" \
    #    > "data/#{STRAIN_NAME}_ilm_consensus.fasta"
  SH
end

desc "Fakes the prerequisites for the recall_ilm_consensus task"
task :recall_ilm_consensus_fake_prereqs do
  abort "FATAL: Task recall_ilm_consensus_fake_prereqs requires specifying STRAIN_NAME" unless STRAIN_NAME 
  mkdir_p "log"
  touch "bash5.fofn"                                  and sleep 1
  touch "data/polished_assembly.fasta.gz"             and sleep 1
  touch "data/polished_assembly_circularized.fasta"   and sleep 1
  touch "data/#{STRAIN_NAME}_consensus.fasta"
end


# =====================
# = rast_annotate_ilm =
# =====================

desc "Submits the Illumina-fixed consensus to RAST for annotations"
task :rast_annotate_ilm => [:check, "data/#{STRAIN_NAME}_ilm_consensus_rast.fna", 
    "data/#{STRAIN_NAME}_ilm_consensus_rast.gbk", "data/#{STRAIN_NAME}_ilm_consensus_rast_aa.fa",
    "data/ilm_rast_job_id"]

file "data/#{STRAIN_NAME}_ilm_consensus_rast.gbk" => ["data/#{STRAIN_NAME}_ilm_consensus.fasta"] do |t|
  fasta = "data/#{STRAIN_NAME}_ilm_consensus.fasta"
  submit_and_retrieve_rast(fasta, t.name, "data/ilm_rast_job_id", "rast_annotate_ilm")
end
file "data/ilm_rast_job_id" => "data/#{STRAIN_NAME}_consensus_rast.gbk"

file "data/#{STRAIN_NAME}_ilm_consensus_rast_aa.fa" => "data/#{STRAIN_NAME}_ilm_consensus_rast.gbk" do |t|
  gb_to_fasta "data/#{STRAIN_NAME}_ilm_consensus_rast.gbk", "#{OUT}/#{t.name}", :aa, "rast_annotate_ilm"
end

file "data/#{STRAIN_NAME}_ilm_consensus_rast.fna" => "data/#{STRAIN_NAME}_ilm_consensus_rast.gbk" do |t|
  gb_to_fasta "data/#{STRAIN_NAME}_ilm_consensus_rast.gbk", "#{OUT}/#{t.name}", :nt, "rast_annotate_ilm"
end
