require 'pp'
require 'net/http'
require_relative 'lib/colors'
require_relative 'lib/lsf_client'
require_relative 'lib/subscreens'
require 'shellwords'
require 'bundler/setup'
require 'rspec/core/rake_task'
include Colors

task :default => :check

LSF = LSFClient.new
LSF.disable! if ENV['LSF_DISABLED']  # Run everything locally if set (useful for debugging)

REPO_DIR = File.dirname(__FILE__)
SAS_DIR = "#{REPO_DIR}/vendor/sas"
MUMMER_DIR = "#{REPO_DIR}/vendor/MUMmer3.23"
ALIEN_DIR = "#{REPO_DIR}/vendor/alien_hunter-1.7"
BCFTOOLS_DIR = "#{REPO_DIR}/vendor/bcftools"
HTSLIB_DIR = "#{REPO_DIR}/vendor/htslib"

IGB_DIR = ENV['IGB_DIR'] || "#{ENV['HOME']}/www/igb"

OUT = File.expand_path(ENV['OUT'] || "#{REPO_DIR}/out")

#######
# Other environment variables that may be set by the user for specific tasks (see README.md)
#######
STRAIN_NAME = ENV['STRAIN_NAME']
SPECIES = ENV['SPECIES']
ILLUMINA_FASTQ = ENV['ILLUMINA_FASTQ'] && File.expand_path(ENV['ILLUMINA_FASTQ'])
ILLUMINA_FASTQ_2 = ENV['ILLUMINA_FASTQ_2'] && File.expand_path(ENV['ILLUMINA_FASTQ_2'])
ILLUMINA_REFERENCE = ENV['ILLUMINA_REFERENCE'] && File.expand_path(ENV['ILLUMINA_REFERENCE'])
TASK_FILE = ENV['TASK_FILE']
GENBANK_REFERENCES = ENV['GENBANK_REFERENCES'] && ENV['GENBANK_REFERENCES'].split(':')
CLUSTER = ENV['CLUSTER']
REPLACE_FASTA = ENV['REPLACE_FASTA'] && File.expand_path(ENV['REPLACE_FASTA'])
CURATED = ENV['CURATED']
PHAGE_DB = ENV['PHAGE_DB']

#############################################################
#  IMPORTANT!
#  This Rakefile runs with the working directory set to OUT
#  All filenames from hereon are relative to that directory
#############################################################
mkdir_p OUT
Dir.chdir(OUT)

task :env do
  puts "Output directory: #{OUT}"
  mkdir_p File.join(REPO_DIR, "vendor")
  
  sc_orga_scratch = "/sc/orga/scratch/#{ENV['USER']}"
  ENV['TMP'] ||= Dir.exists?(sc_orga_scratch) ? sc_orga_scratch : "/tmp"
  # Always use our locally bundled (patched) perl modules
  ENV['PERL5LIB'] = "#{REPO_DIR}/lib/perl"
end

file "#{REPO_DIR}/scripts/env.sh" => "#{REPO_DIR}/scripts/example.env.sh" do
  cp "#{REPO_DIR}/scripts/example.env.sh", "#{REPO_DIR}/scripts/env.sh"
end

ENV_ERROR = "Configure this in scripts/env.sh and run `source scripts/env.sh` before running rake."

desc "Checks environment variables and requirements before running tasks"
task :check => [:env, "#{REPO_DIR}/scripts/env.sh", :mummer, :bcftools, :alien_hunter] do
  unless `module avail 2>&1 | grep smrtpipe/2.2.0` != ''
    abort "FATAL: You must have the smrtpipe/2.2.0 module in your MODULEPATH."
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


# pulls down and compiles MUMmer 3.23, which is used by scripts/circularizeContigs.pl and others
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

task :alien_hunter => [:env, ALIEN_DIR, "#{ALIEN_DIR}/alien_hunter"]
directory ALIEN_DIR
file "#{ALIEN_DIR}/alien_hunter" do
  Dir.chdir("#{REPO_DIR}/vendor/") do
    system <<-SH
      curl -L -o alien_hunter.tar.gz 'ftp://ftp.sanger.ac.uk/pub/resources/software/alien_hunter/alien_hunter.tar.gz'
      tar xvzf alien_hunter.tar.gz  # Creates alien_hunter dir
    SH
  end
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

file "pathogendb-pipeline.png" => [:graph]
desc "Generates a graph of tasks, intermediate files and their dependencies from this Rakefile"
task :graph do
  system <<-SH
    module load graphviz
    STRAIN_NAME='${STRAIN_NAME}' SPECIES='${SPECIES}' SMRT_JOB_ID='${SMRT_JOB_ID}' rake -f \
        #{Shellwords.escape(__FILE__)} -P \
        | #{REPO_DIR}/scripts/rake-prereqs-dot.rb --prune #{REPO_DIR} --replace-with '$REPO_DIR' \
               --narrow-path old:check,check,default\
        | dot -Tpng -o pathogendb-pipeline.png
  SH
end

# Creates a special :spec task that runs all tests defined in spec/.
RSpec::Core::RakeTask.new(:spec, :speed) do |t, args|
  t.pattern = Dir.glob("#{REPO_DIR}/spec/**/*_spec.rb")
  t.rspec_opts = '--format documentation'
  t.rspec_opts << ' --fail-fast'
  t.rspec_opts << ' --color'
  if args[:speed]
    t.rspec_opts << " --tag speed:#{args[:speed]}" unless args[:speed] == 'all'
  else
    t.rspec_opts << ' --tag ~speed:slow'   # avoid slow tests by default
  end
end

desc "Clean all intermediate files from the OUT directory (and if $prereqs is set, all downloaded software in vendor/ too)"
task :clean, [:prereqs] do |t, args|
  rm_f "bash5.fofn"
  rm_f "pathogendb-pipeline.png"
  rm_rf "data"
  rm_rf Dir.glob("#{REPO_DIR}/vendor/*") if args[:prereqs]
end


# =======================
# = pull_down_raw_reads =
# =======================

desc "Copies or downloads raw reads from a PacBio job to the OUT directory"
task :pull_down_raw_reads => [:check, "bash5.fofn"]  # <-- file(s) created by this task
file "bash5.fofn" do |t, args|                       # <-- implementation for generating each of these files
  job_id = ENV['SMRT_JOB_ID']                        # Example SMRT_JOB_ID's that work are: 019194, 020266
  abort "FATAL: Task pull_down_raw_reads requires specifying SMRT_JOB_ID" unless job_id
  job_id = job_id.rjust(6, '0')
  pacbio_job_dirs = ["/sc/orga/projects/pacbio/userdata_permanent/jobs/#{job_id[0..2]}/#{job_id}",
                     "/sc/orga/projects/InfectiousDisease/old_smrtportal_jobs/#{job_id}","/sc/orga/scratch/attieo02/#{job_id}"]
  smrtpipe_log_url = "http://node1.1425mad.mssm.edu/pacbio/secondary/#{job_id[0..2]}/#{job_id}/log/smrtpipe.log"
  
  found_fofn_dir = pacbio_job_dirs.find {|dir| File.exist? "#{dir}/input.fofn" }
  if found_fofn_dir
    cp "#{found_fofn_dir}/input.fofn", "bash5.fofn"
    mkdir_p "data"
    if File.exist? "#{found_fofn_dir}/data/polished_assembly.fasta.gz"
      ln_s "#{found_fofn_dir}/data/polished_assembly.fasta.gz", "data/polished_assembly.fasta.gz"
      ln_s "#{found_fofn_dir}/data/corrected.fastq", "data/corrected.fastq"
      ln_s "#{found_fofn_dir}/data/celera-assembler.gkpStore", "data/celera-assembler.gkpStore"
      ln_s "#{found_fofn_dir}/data/celera-assembler.tigStore", "data/celera-assembler.tigStore"
      ln_s "#{found_fofn_dir}/data/4-unitigger/best.edges", "data/best.edges"
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
  
  # The optional REPLACE_FASTA parameter allows one to shunt a custom assembly into the rest of the pipeline 
  #     in place of the polished_assembly.fasta.gz built by SMRTPortal.
  # Note that the fastq, best edges data, etc., if they do not match the new fasta, may throw off QC analyses.
  if REPLACE_FASTA
    cp REPLACE_FASTA, "data/replaced_assembly.fasta"
    system "gzip data/replaced_assembly.fasta"  # creates data/replaced_assembly.fasta.gz
    if File.exist? "data/polished_assembly.fasta.gz"
        rm "data/polished_assembly.fasta.gz"
    end
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
  
  if REPLACE_FASTA
    puts "NOTICE: polished_assembly.fasta.gz has been replaced by user input, skipping assemble_raw_reads"
    cp "data/replaced_assembly.fasta.gz", "data/polished_assembly.fasta.gz"
    next
  end
  
  lstat = File::lstat("data/polished_assembly.fasta.gz") rescue nil
  if lstat and lstat.symlink?
    puts "NOTICE: polished_assembly.fasta.gz is symlinked to an existing assembly, skipping assemble_raw_reads"
    next
  end
  
  cp "#{REPO_DIR}/xml/example_params.xml", "."
  system <<-SH
    module load smrtpipe/2.2.0
    source "#{ENV['SMRTANALYSIS']}/etc/setup.sh" &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
        -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params example_params.xml xml:bash5.xml
  SH
end

# (Some steps were snipped out here because they have been replaced with circlator and prokka.)
load "#{REPO_DIR}/deprecated.rake"

# =================
# = run_circlator =
# =================

desc "Runs circlator on smrtpipe output."
task :run_circlator => [:check, "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta"]
file "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta" => "data/polished_assembly.fasta.gz" do |t|
  system <<-SH
    cp data/polished_assembly.fasta.gz data/circ_input.fasta.gz
    gunzip -f data/circ_input.fasta.gz
  SH
  puts "INFO: starting circlator."
  if CURATED
    system <<-SH or abort "FATAL: circlator failed to run to completion."
      module purge
      module load bwa/0.7.13
      module load prodigal/2.6.2
      module load samtools/1.1
      module load spades/3.6.0
      module load python/3.5.0  py_packages/3.5
      module load mummer/3.23
      mkdir -p data/#{STRAIN_NAME}_circlator
      circlator fixstart data/circ_input.fasta data/#{STRAIN_NAME}_circlator/06.fixstart
    SH
  else
    system <<-SH or abort "FATAL: circlator failed to run to completion."
      module purge
      module load bwa/0.7.13
      module load prodigal/2.6.2
      module load samtools/1.1
      module load spades/3.6.0
      module load python/3.5.0  py_packages/3.5
      module load mummer/3.23
      rm -rf data/#{STRAIN_NAME}_circlator
      circlator all data/circ_input.fasta data/corrected.fastq data/#{STRAIN_NAME}_circlator/
    SH
  end
end


# ==================
# = post_circlator =
# ==================

desc "Renames the reoriented assembly contigs using a shortened scheme. If reoriented too near the contig start, reorientate to middle of contig."
task :post_circlator => [:check, "data/#{STRAIN_NAME}_postcirc.fasta"]
file "data/#{STRAIN_NAME}_postcirc.fasta" => "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta" do |t|
  job_id = ENV['SMRT_JOB_ID']                        # Example SMRT_JOB_ID's that work are: 019194, 020266
  abort "FATAL: Task pull_down_raw_reads requires specifying SMRT_JOB_ID" unless job_id
  job_id = job_id.rjust(6, '0')
  system <<-SH
    # call script to rename contigs
    #{REPO_DIR}/scripts/post_circlator_contig_rename.py data/#{STRAIN_NAME}_circlator/ data/#{STRAIN_NAME}_postcirc.fasta data/#{STRAIN_NAME}_postcirc2.txt #{job_id}
  SH
end


# =======================
# = resequence_assembly =
# =======================

desc "Resequences the circularized assembly"
task :resequence_assembly => [:check, "data/#{STRAIN_NAME}_consensus_circ.fasta"]
file "data/#{STRAIN_NAME}_consensus_circ.fasta" => "data/#{STRAIN_NAME}_postcirc.fasta" do |t|
  abort "FATAL: Task resequence_assembly requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: STRAIN_NAME can only contain letters, numbers, and underscores" unless STRAIN_NAME =~ /^[\w]+$/
  
  rm_rf "circularized_sequence"
  mkdir_p "circularized_sequence"
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    referenceUploader -c -p circularized_sequence -n #{STRAIN_NAME} -f data/#{STRAIN_NAME}_postcirc.fasta
  SH
  # NOTE: sometimes referenceUploader auto-appends a timestamp to the reference name given by `-n` to avoid a conflict
  # Therefore, we must detect if it did this by checking how it named the directory
  reference_dir = Dir.glob("circularized_sequence/#{STRAIN_NAME}*").last
  cp "#{REPO_DIR}/xml/resequence_example_params.xml", OUT
  system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl resequence_example_params.xml " +
      "#{OUT} #{reference_dir} > resequence_params.xml" and
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    samtools faidx #{reference_dir}/sequence/#{STRAIN_NAME}.fasta &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
        -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params resequence_params.xml xml:bash5.xml &&
    gunzip -f data/consensus.fasta.gz
  SH
  cp "data/consensus.fasta", "data/#{STRAIN_NAME}_consensus_circ.fasta"
end


# ==============================
# = post_quiver_orient_correct =
# ==============================

desc "Corrects orientation post-resequencing and renames contigs using a shortened scheme"
task :post_quiver_orient_correct => [:check, "data/#{STRAIN_NAME}_prokka.fasta"]
file "data/#{STRAIN_NAME}_prokka.fasta" => "data/#{STRAIN_NAME}_consensus_circ.fasta" do |t|
  abort "FATAL: Task post_quiver_orient_correct requires specifying STRAIN_NAME" unless STRAIN_NAME 
  
  system <<-SH
    module load blast/2.2.26+
    #{REPO_DIR}/scripts/post_quiver_orient_correct.py data/#{STRAIN_NAME}_consensus_circ.fasta data/#{STRAIN_NAME}_postcirc2.txt data/#{STRAIN_NAME}_prokka.fasta data/pq_dir
  SH
end


# ===================
# = prokka_annotate =
# ===================

desc "Annotates the reoriented assembly with prokka"
task :prokka_annotate => [:check, "data/prokka/#{STRAIN_NAME}_prokka.gbk"]
file "data/prokka/#{STRAIN_NAME}_prokka.gbk" => "data/#{STRAIN_NAME}_prokka.fasta" do |t|
  abort "FATAL: Task prokka_annotate requires specifying STRAIN_NAME" unless STRAIN_NAME 
  
  system <<-SH
    module purge
    module load prokka/1.11  
    module load barrnap/0.6
    module unload rnammer/1.2
    module load minced/0.2.0
    module load signalp/4.1
        
    prokka --outdir data/prokka --force --prefix #{STRAIN_NAME}_prokka data/#{STRAIN_NAME}_prokka.fasta
  SH
end

# =====================
# = repeats_phage_pai =
# =====================

desc "Creates bedFile of repeats, phage and PAIs"
task :repeats_phage_pai => [:check, "data/www/wiggle/#{STRAIN_NAME}.rpi.phage.bed"]
file "data/www/wiggle/#{STRAIN_NAME}.rpi.phage.bed" => "data/#{STRAIN_NAME}_prokka.fasta" do |t|
  abort "FATAL: Task prokka_annotate requires specifying STRAIN_NAME" unless STRAIN_NAME

  system <<-SH
    module purge
    module load mummer
    module load blast
    module load python/2.7.6
    module load py_packages/2.7
    mkdir -p data/www/wiggle
    python #{REPO_DIR}/scripts/get_repeats_phage_pai.py -a #{ALIEN_DIR}/alien_hunter -d #{PHAGE_DB} -o data/www/wiggle/#{STRAIN_NAME}.rpi -f data/#{STRAIN_NAME}_prokka.fasta \
    --islands --repeats --phage
  SH
end


# =====================
# = create_QC_webpage =
# =====================

desc "Creates the QC webpage"
task :create_QC_webpage => [:check, "data/www/index.html"]
file "data/www/index.html" => "data/#{STRAIN_NAME}_prokka.fasta" do |t|
  job_id = ENV['SMRT_JOB_ID']
  abort "FATAL: Task create_QC_webpage requires specifying SMRT_JOB_ID" unless job_id
  abort "FATAL: Task create_QC_webpage requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task create_QC_webpage requires specifying SPECIES" unless SPECIES 
  species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES
  
  system <<-SH
    module purge
    module load blast/2.2.26+
    module load bwa/0.7.12
    module load celera/8.1
    module load python/2.7.6
    module load py_packages/2.7
    module load ucsc-utils/2015-04-07
    module load samtools/1.2
    #{REPO_DIR}/scripts/create_QC_webpage.py -o data/qc_wd -w data/www -f data/#{STRAIN_NAME}_prokka.fasta \
     -g data -r data/corrected.fastq -a #{species_clean}_#{STRAIN_NAME}_#{job_id}
  SH
end


# =================
# = prokka_QC_rpi =
# =================

desc "Run prokka and create the QC website"
task :prokka_QC_rpi => [:prokka_annotate, :create_QC_webpage, :repeats_phage_pai]




# =================
# = prokka_to_igb =
# =================

directory IGB_DIR

desc "Creates an IGB Quickload-compatible directory for your genome in IGB_DIR"
task :prokka_to_igb => [:check, :prokka_QC_rpi] do |t|
  job_id = ENV['SMRT_JOB_ID']
  abort "FATAL: Task prokka_to_igb requires specifying SMRT_JOB_ID" unless job_id
  abort "FATAL: Task prokka_to_igb requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task prokka_to_igb requires specifying SPECIES" unless SPECIES 
  abort "FATAL: Task prokka_to_igb requires specifying IGB_DIR" unless IGB_DIR 
  species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES
    
  system <<-SH
    module purge
    module load python/2.7.6
    module load py_packages/2.7
    module load blat
    module load bioperl
    export SAS_DIR=#{SAS_DIR}
    export REPO_DIR=#{REPO_DIR}
    perl #{REPO_DIR}/scripts/rast2igb.pl \
        -f data/prokka/#{STRAIN_NAME}_prokka.gbk \
        -g #{species_clean}_#{STRAIN_NAME}_#{job_id} \
        -q data/www/ \
        -w data/qc_wd/bigwig/ \
        -b data/qc_wd/alignment.sorted.bam \
        -i #{IGB_DIR} \
        -r #{REPO_DIR}
  SH
end


# =====================
# = igb_to_pathogendb =
# =====================

directory IGB_DIR

desc "Adds a new genome assembly to pathogendb from an IGB genome dir"
task :igb_to_pathogendb => [:check, :prokka_to_igb] do |t|
  job_id = ENV['SMRT_JOB_ID']
  abort "FATAL: Task igb_to_pathogendb requires specifying SMRT_JOB_ID" unless job_id
  abort "FATAL: Task igb_to_pathogendb requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task igb_to_pathogendb requires specifying SPECIES" unless SPECIES 
  abort "FATAL: Task igb_to_pathogendb requires specifying IGB_DIR" unless IGB_DIR 
  species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES
  
  system <<-SH
    export SAS_DIR=#{SAS_DIR}
    perl #{REPO_DIR}/scripts/igb2pathogendb.pl \
        -i #{IGB_DIR}/#{species_clean}_#{STRAIN_NAME}_#{job_id}
  SH
end


# ==================
# = motif_and_mods =
# ==================

desc "Reruns SMRTPipe for modification and motif analysis on the polished assembly"
task :motif_and_mods => [:check, "data/motif_summary.csv"]
file "data/motif_summary.csv" => "data/#{STRAIN_NAME}_prokka.fasta" do |t|
  abort "FATAL: Task motif_and_mods requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: STRAIN_NAME can only contain letters, numbers, and underscores" unless STRAIN_NAME =~ /^[\w]+$/
  
  rm_rf "polished_sequence"
  mkdir_p "polished_sequence"
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    referenceUploader -c -p polished_sequence -n #{STRAIN_NAME} -f data/#{STRAIN_NAME}_prokka.fasta
  SH
  # NOTE: sometimes referenceUploader auto-appends a timestamp to the reference name given by `-n` to avoid a conflict
  # Therefore, we must detect if it did this by checking how it named the directory
  reference_dir = Dir.glob("polished_sequence/#{STRAIN_NAME}*").last
  cp "#{REPO_DIR}/xml/motif_simple_example_params.xml", OUT
  system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl motif_simple_example_params.xml " +
      "#{OUT} #{reference_dir} > motif_simple_params.xml" and
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    samtools faidx #{reference_dir}/sequence/#{STRAIN_NAME}.fasta &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
        -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params motif_simple_params.xml xml:bash5.xml
  SH
end


# =======
# = all =
# =======

desc "Runs entire pipeline for PacBio-only data (up to prokka_to_IGB and motif_and_mods)"
task :all => [:prokka_to_igb, :motif_and_mods]



# ====================================================================================================
# = The following tasks are for assemblies where we want to incorporate Illumina reads to fix indels =
# ====================================================================================================

namespace :ilm do


  # ====================
  # = ilm:fake_prereqs =
  # ====================
  
  desc "Fakes the prerequisites for the Illumina tasks (for testing purposes)"
  task :fake_prereqs do |t|
    abort "FATAL: Task #{t.name} requires specifying STRAIN_NAME" unless STRAIN_NAME 
    mkdir_p "log"
    mkdir_p "data/#{STRAIN_NAME}_circlator"
    touch "bash5.fofn"                                        and sleep 1
    touch "data/polished_assembly.fasta.gz"                   and sleep 1
    if ILLUMINA_REFERENCE
      abort "FATAL: file '#{ILLUMINA_REFERENCE}' does not exist" unless File.exists?(ILLUMINA_REFERENCE)
      cp(ILLUMINA_REFERENCE, "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta") and sleep 1
      cp(ILLUMINA_REFERENCE, "data/#{STRAIN_NAME}_postcirc.fasta")              and sleep 1
      cp(ILLUMINA_REFERENCE, "data/#{STRAIN_NAME}_consensus_circ.fasta")        and sleep 1
      cp(ILLUMINA_REFERENCE, "data/#{STRAIN_NAME}_prokka.fasta")
    else
      touch "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta" and sleep 1
      touch "data/#{STRAIN_NAME}_postcirc.fasta"              and sleep 1
      touch "data/#{STRAIN_NAME}_consensus_circ.fasta"        and sleep 1
      touch "data/#{STRAIN_NAME}_prokka.fasta"
      puts "Replace #{OUT}/data/#{STRAIN_NAME}_prokka.fasta with the old reference sequence you are piling new reads onto."
    end
  end
  
  
  # ========================
  # = ilm:recall_consensus =
  # ========================

  desc "Recalls a new consensus by piling Illumina reads onto a circlator-finished PacBio assembly"
  task :recall_consensus => [:check, "data/#{STRAIN_NAME}_prokka_flt.vcf", "data/#{STRAIN_NAME}_ilm_prokka.fasta"]

  file "data/prokka.ref.sort.bam" => "data/#{STRAIN_NAME}_prokka.fasta" do |t|
    abort "FATAL: Task ilm:recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
    abort "FATAL: Task ilm:recall_consensus requires specifying ILLUMINA_FASTQ" unless ILLUMINA_FASTQ

    if ILLUMINA_FASTQ_2
      bwa_mem_args = <<-SH
        "data/#{STRAIN_NAME}_prokka.fasta" #{Shellwords.escape(ILLUMINA_FASTQ)} \
            #{Shellwords.escape(ILLUMINA_FASTQ_2)} > data/prokka.ref.aln.sam
      SH
    else
      bwa_mem_args = <<-SH
        "data/#{STRAIN_NAME}_prokka.fasta" #{Shellwords.escape(ILLUMINA_FASTQ)} > data/prokka.ref.aln.sam
      SH
    end
    
    LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
    LSF.job_name "prokka.ref.aln.sam"
    LSF.bsub_interactive <<-SH or abort
      module load bwa/0.7.12
      bwa index "data/#{STRAIN_NAME}_prokka.fasta"
      bwa mem #{bwa_mem_args}
      module load samtools/1.1
      samtools view -bS data/prokka.ref.aln.sam > data/prokka.ref.aln.bam
      samtools sort data/prokka.ref.aln.bam data/prokka.ref.sort
    SH

    # Can remove the SAM as it is huge and the .aln.bam should contain everything in it
    rm "data/prokka.ref.aln.sam"
  end

  file "data/#{STRAIN_NAME}_prokka_raw.bcf" => "data/prokka.ref.sort.bam" do |t|
    abort "FATAL: Task ilm:recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
    # Use mpileup to do the consensus calling.
    # Note: the -L and -d flags are important; they ensure samtools looks at up to 100k reads per base to call variants.
    # The default for -L is 250, which would turn off indel calling for deeply resequenced (depth >250) samples.
    LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
    LSF.job_name "#{STRAIN_NAME}_prokka_raw.bcf"
    LSF.bsub_interactive <<-SH
      module load samtools/1.1
      module load bcftools/1.1
      samtools mpileup -L100000 -d100000 -uf "data/#{STRAIN_NAME}_prokka.fasta" data/prokka.ref.sort.bam \
          | bcftools call -cv -Ob > "data/#{STRAIN_NAME}_prokka_raw.bcf"
    SH
  end

  file "data/#{STRAIN_NAME}_prokka_flt.vcf" => "data/#{STRAIN_NAME}_prokka_raw.bcf" do |t|
    LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
    LSF.job_name "#{STRAIN_NAME}_prokka_flt.vcf"
    LSF.bsub_interactive <<-SH
      module load samtools/1.1
      module load bcftools/1.1
      bcftools view "data/#{STRAIN_NAME}_prokka_raw.bcf" | vcfutils.pl varFilter > "data/#{STRAIN_NAME}_prokka_flt.vcf"
    SH
  end

  file "data/#{STRAIN_NAME}_ilm_fix.fasta" =>
      ["data/#{STRAIN_NAME}_prokka_flt.vcf", "data/#{STRAIN_NAME}_prokka.fasta"] do |t|
    system <<-SH
      module load vcftools/0.1.12b
      module load tabix/0.2.6 
  
      bgzip -c "data/#{STRAIN_NAME}_prokka_flt.vcf" > "data/#{STRAIN_NAME}_prokka_flt.vcf.gz"
      tabix -p vcf "data/#{STRAIN_NAME}_prokka_flt.vcf.gz"
      cat "data/#{STRAIN_NAME}_prokka.fasta" | vcf-consensus "data/#{STRAIN_NAME}_prokka_flt.vcf.gz" \
              > "data/#{STRAIN_NAME}_ilm_fix.fasta"
  
      # New-style version of doing this with bcftools consensus, but it doesn't work (memory leak in bcftools)
      # #{HTSLIB_DIR}/bgzip -c "data/#{STRAIN_NAME}_prokka_flt.vcf" > "data/#{STRAIN_NAME}_prokka_flt.vcf.gz"
      # #{HTSLIB_DIR}/tabix -p vcf "data/#{STRAIN_NAME}_prokka_flt.vcf.gz"
      # #{BCFTOOLS_DIR}/bcftools consensus -f "data/#{STRAIN_NAME}_prokka.fasta" "data/#{STRAIN_NAME}_prokka_flt.vcf.gz" \
      #    > "data/#{STRAIN_NAME}_ilm_fix.fasta"
    SH
  end

  file "data/#{STRAIN_NAME}_ilm_corrected.fasta" => "data/#{STRAIN_NAME}_ilm_fix.fasta" do |t|
    if ILLUMINA_FASTQ_2
      fix_repeats_ill = <<-SH
        #{REPO_DIR}/scripts/fix_repeats_ill.py -c data/ilm_coverage.cov -g data/#{STRAIN_NAME}_ilm_fix.fasta \
            -r #{Shellwords.escape(ILLUMINA_FASTQ)} -r2 #{Shellwords.escape(ILLUMINA_FASTQ_2)} -w data/ilm_fix \
            -o data/#{STRAIN_NAME}_ilm_corrected.fasta
      SH
    else
      fix_repeats_ill = <<-SH
        #{REPO_DIR}/scripts/fix_repeats_ill.py -c data/ilm_coverage.cov -g data/#{STRAIN_NAME}_ilm_fix.fasta \
            -r #{Shellwords.escape(ILLUMINA_FASTQ)} -w data/ilm_fix -o data/#{STRAIN_NAME}_ilm_corrected.fasta
      SH
    end
    system <<-SH or abort
      module purge
      module load bwa/0.7.12
      module load samtools/1.1
      module load bcftools/1.1
      module load tabix/0.2.6
      module load vcftools/0.1.12b
      module load bedtools/2.21.0
      
      cut -f1,2 "data/#{STRAIN_NAME}_prokka.fasta.fai" > "data/#{STRAIN_NAME}_prokka.chrom.sizes"
      genomeCoverageBed -d -ibam data/prokka.ref.sort.bam -g "data/#{STRAIN_NAME}_prokka.chrom.sizes" > data/ilm_coverage.cov
      
      module unload python
      module unload py_packages
      module load python/2.7.6
      module load py_packages/2.7
      #{fix_repeats_ill}
    SH
  end

  file "data/#{STRAIN_NAME}_ilm_prokka.fasta" => "data/#{STRAIN_NAME}_ilm_corrected.fasta" do |t|
    rm_rf "circularized_sequence"
    mkdir_p "circularized_sequence"
    system <<-SH or abort
      module load smrtpipe/2.2.0
      source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
      referenceUploader -c -p circularized_sequence -n #{STRAIN_NAME} -f data/#{STRAIN_NAME}_ilm_corrected.fasta
    SH
    # NOTE: sometimes referenceUploader auto-appends a timestamp to the reference name given by `-n` to avoid a conflict
    # Therefore, we must detect if it did this by checking how it named the directory
    reference_dir = Dir.glob("circularized_sequence/#{STRAIN_NAME}*").last
    cp "#{REPO_DIR}/xml/resequence_example_params.xml", OUT
    system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl resequence_example_params.xml " +
      "#{OUT} #{reference_dir} > resequence_params.xml" and
    system <<-SH or abort
      module load smrtpipe/2.2.0
      source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
      samtools faidx #{reference_dir}/sequence/#{STRAIN_NAME}.fasta &&
      smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
          -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params resequence_params.xml xml:bash5.xml &&
      gunzip -f data/consensus.fasta.gz
      python #{REPO_DIR}/scripts/smrt_vcf_to_consensus.py "data/#{STRAIN_NAME}_ilm_corrected.fasta" data/variants.vcf 40\
       > "data/#{STRAIN_NAME}_ilm_prokka.fasta"
    SH
  end
  

  # =======================
  # = ilm:prokka_annotate =
  # =======================

  desc "Annotates the Illumina-corrected assembly with prokka"
  task :prokka_annotate => [:check, "data/prokka/#{STRAIN_NAME}_ilm_prokka.gbk"]
  file "data/prokka/#{STRAIN_NAME}_ilm_prokka.gbk" => "data/#{STRAIN_NAME}_ilm_prokka.fasta" do |t|
    abort "FATAL: Task ilm:prokka_annotate requires specifying STRAIN_NAME" unless STRAIN_NAME 
  
    system <<-SH
      module load prokka  
      module load barrnap
      module unload rnammer
      module load minced
      module load signalp
        
      prokka --outdir data/prokka --force --prefix #{STRAIN_NAME}_ilm_prokka data/#{STRAIN_NAME}_ilm_prokka.fasta
    SH
  end




  # =========================
  # = ilm:create_QC_webpage =
  # =========================

  desc "Creates the QC webpage for the Illumina-corrected assembly"
  task :create_QC_webpage => [:check, "data/ilm_www/index.html"]
  file "data/ilm_www/index.html" => "data/#{STRAIN_NAME}_ilm_prokka.fasta" do |t|
    job_id = ENV['SMRT_JOB_ID']
    abort "FATAL: Task ilm:create_QC_webpage requires specifying SMRT_JOB_ID" unless job_id
    abort "FATAL: Task ilm:create_QC_webpage requires specifying STRAIN_NAME" unless STRAIN_NAME 
    abort "FATAL: Task ilm:create_QC_webpage requires specifying SPECIES" unless SPECIES 
    species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES
  
    system <<-SH
      module unload python
      module unload py_packages
      module load blast
      module load bwa/0.7.12
      module load celera
      module load python/2.7.6
      module load py_packages/2.7
      module load ucsc-utils
      module load samtools/1.2
      #{REPO_DIR}/scripts/create_QC_webpage.py -o data/ilm_qc_wd -w data/ilm_www -f data/#{STRAIN_NAME}_ilm_prokka.fasta \
       -g data -r data/corrected.fastq -a #{species_clean}_#{STRAIN_NAME}_#{job_id}
    SH
  end


  # =========================
  # = ilm:repeats_phage_pai =
  # =========================

  desc "Creates bedFile of repeats, phage and PAIs"
  task :repeats_phage_pai => [:check, "data/ilm_www/wiggle/#{STRAIN_NAME}.rpi.phage.bed"]
  file "data/ilm_www/wiggle/#{STRAIN_NAME}.rpi.phage.bed" => "data/#{STRAIN_NAME}_ilm_prokka.fasta" do |t|
    abort "FATAL: Task prokka_annotate requires specifying STRAIN_NAME" unless STRAIN_NAME

    system <<-SH
      module purge
      module load mummer
      module load blast
      module load python/2.7.6
      module load py_packages/2.7
      mkdir -p data/www/wiggle
      python #{REPO_DIR}/scripts/get_repeats_phage_pai.py -a #{ALIEN_DIR}/alien_hunter -d #{PHAGE_DB} -o data/ilm_www/wiggle/#{STRAIN_NAME}.rpi -f data/#{STRAIN_NAME}_ilm_prokka.fasta \
      --islands --repeats --phage
     SH
  end

  # =====================
  # = ilm:prokka_QC_rpi =
  # =====================

  desc "Run prokka and create the QC website for illumina correted data"
  task :prokka_QC_rpi => [ilm:prokka_annotate, ilm:create_QC_webpage, ilm:repeats_phage_pai]


  
  # =====================
  # = ilm:prokka_to_igb =
  # =====================

  desc "Creates an IGB Quickload-compatible directory for the Illumina-corrected assembly in IGB_DIR"
  task :prokka_to_igb => [:check, ilm:prokka_QC_rpi] do |t|
    job_id = ENV['SMRT_JOB_ID']
    abort "FATAL: Task ilm:prokka_to_igb requires specifying SMRT_JOB_ID" unless job_id
    abort "FATAL: Task ilm:prokka_to_igb requires specifying STRAIN_NAME" unless STRAIN_NAME 
    abort "FATAL: Task ilm:prokka_to_igb requires specifying SPECIES" unless SPECIES 
    abort "FATAL: Task ilm:prokka_to_igb requires specifying IGB_DIR" unless IGB_DIR 
    species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES
    
    system <<-SH
      module unload python
      module unload py_packages
      module load python/2.7.6
      module load py_packages/2.7
      module load blat
      module load bioperl
      export SAS_DIR=#{SAS_DIR}
      perl #{REPO_DIR}/scripts/rast2igb.pl \
          -f data/prokka/#{STRAIN_NAME}_ilm_prokka.gbk \
          -g #{species_clean}_#{STRAIN_NAME}_#{job_id} \
          -q data/ilm_www/ \
          -w data/ilm_qc_wd/bigwig/ \
          -b data/ilm_qc_wd/alignment.sorted.bam \
          -i #{IGB_DIR}
    SH
  end
  
end # namespace :ilm
