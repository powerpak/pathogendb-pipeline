require 'pp'
require 'net/http'
require_relative 'lib/colors'
require_relative 'lib/lsf_client'
require_relative 'lib/subscreens'
require 'shellwords'
include Colors

task :default => :check

LSF = LSFClient.new
LSF.disable! if ENV['LSF_DISABLED']  # Run everything locally if set (useful for debugging)

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
ILLUMINA_FASTQ = ENV['ILLUMINA_FASTQ'] && File.expand_path(ENV['ILLUMINA_FASTQ'])
ILLUMINA_REFERENCE = ENV['ILLUMINA_REFERENCE'] && File.expand_path(ENV['ILLUMINA_REFERENCE'])
TASK_FILE = ENV['TASK_FILE']
GENBANK_REFERENCES = ENV['GENBANK_REFERENCES'] && ENV['GENBANK_REFERENCES'].split(':')
CLUSTER = ENV['CLUSTER']

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
task :check => [:env, "#{REPO_DIR}/scripts/env.sh", :sas, :mummer, :bcftools] do
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
    STRAIN_NAME='${STRAIN_NAME}' SPECIES='${SPECIES}' SMRT_JOB_ID='${SMRT_JOB_ID}' rake -f \
        #{Shellwords.escape(__FILE__)} -P \
        | #{REPO_DIR}/scripts/rake-prereqs-dot.rb --prune #{REPO_DIR} --replace-with '$REPO_DIR' \
               --narrow-path check,default\
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
  job_id = ENV['SMRT_JOB_ID']                        # Example SMRT_JOB_ID's that work are: 019194, 020266
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
  
  cp "#{REPO_DIR}/xml/example_params.xml", "."
  system <<-SH
    module load smrtpipe/2.2.0
    source "#{ENV['SMRTANALYSIS']}/etc/setup.sh" &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
        -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params example_params.xml xml:bash5.xml
  SH
end

# ========================
# = circularize_assembly =
# ========================

desc "Circularizes the PacBio assembly"
task :circularize_assembly => [:check, "data/polished_assembly_circularized.fasta"]
file "data/polished_assembly_circularized.fasta" => "data/polished_assembly.fasta.gz" do |t|
  system "gunzip -c data/polished_assembly.fasta.gz >data/polished_assembly.fasta" and
  system "#{REPO_DIR}/scripts/circularizeContigs.pl -i data/polished_assembly.fasta -l 12000 2> polished_assembly_circularized.log"
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
  cp "#{REPO_DIR}/xml/resequence_example_params.xml", OUT
  system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl resequence_example_params.xml " +
      "#{OUT} circularized_sequence/#{STRAIN_NAME} > resequence_params.xml" and
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    samtools faidx circularized_sequence/#{STRAIN_NAME}/sequence/#{STRAIN_NAME}.fasta &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
        -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params resequence_params.xml xml:bash5.xml &&
    gunzip data/consensus.fasta.gz
  SH
  cp "data/consensus.fasta", "data/#{STRAIN_NAME}_consensus.fasta"
end


# =======================
# = reorient_assembly =
# =======================

desc "Reorients the circularized assembly to a given locus"
task :reorient_assembly => [:check, "data/#{STRAIN_NAME}_reorient.fasta"]
file "data/#{STRAIN_NAME}_reorient.fasta" => "data/#{STRAIN_NAME}_consensus.fasta" do |t|
  abort "FATAL: Task reorient_assembly requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: STRAIN_NAME can only contain letters, numbers, and underscores" unless STRAIN_NAME =~ /^[\w]+$/
  reorient_fasta = ENV['REORIENT_FASTA']
  reorient_fasta_doesnt_exist = reorient_fasta && !File.exist?(reorient_fasta)
  abort "FATAL: REORIENT_FASTA is nonempty but does not point to a file" if reorient_fasta_doesnt_exist
  reorient_flank = (ENV['REORIENT_FLANK'] || 25).to_i
  
  if reorient_fasta
    system <<-SH
      module load blat
      perl #{REPO_DIR}/scripts/fasta-orient-to-landmark.pl --key _circ --flank #{reorient_flank} \
          --landmark #{Shellwords.escape reorient_fasta} --type prot --matchlength 0.9 --orientsuffix reorient \
          --genome "data/#{STRAIN_NAME}_consensus.fasta" 2> "data/#{STRAIN_NAME}_reorient.log" >"data/#{STRAIN_NAME}_reorient.fasta"
    SH
    if File.size("data/#{STRAIN_NAME}_reorient.fasta") == 0
      rm "data/#{STRAIN_NAME}_reorient.fasta"
      abort "FATAL: Task reorient_assembly failed, exiting"
    end
  else
    # No reorient locus given, simply copy the assembly for the next step
    cp "data/#{STRAIN_NAME}_consensus.fasta", "data/#{STRAIN_NAME}_reorient.fasta"
  end
end

# =================
# = run_circlator =
# =================

desc "Runs circlator on smrtpipe output."
task :run_circlator => [:check, "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta"]
file "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta" => "data/polished_assembly.fasta.gz" do |t|
  system <<-SH
    module unload python
    module unload py_packages
    module load prodigal
    module load samtools
    module load spades
    module load python/3.5.0  py_packages/3.5
    module load gcc
    module load mummer
    module load bwa/0.7.12
    cp data/polished_assembly.fasta.gz data/circ_input.fasta.gz
    gunzip data/circ_input.fasta.gz
    circlator all data/circ_input.fasta data/corrected.fastq data/#{STRAIN_NAME}_circlator/
  SH
end


# ==================
# = post_circlator =
# ==================

desc "Renames the reoriented assembly contigs using a shortened scheme. If reorientated start too near contig start reorientate to middle of contig."
task :post_circlator => [:check, "data/#{STRAIN_NAME}_postcirc.fasta"]
file "data/#{STRAIN_NAME}_postcirc.fasta" => "data/#{STRAIN_NAME}_circlator/06.fixstart.fasta" do |t|
  job_id = ENV['SMRT_JOB_ID']                        # Example SMRT_JOB_ID's that work are: 019194, 020266
  abort "FATAL: Task pull_down_raw_reads requires specifying SMRT_JOB_ID" unless job_id
  job_id = job_id.rjust(6, '0')
  system <<-SH
    # call script to rename contigs
    #{REPO_DIR}/scripts/post_circlator_contig_rename.py data/#{STRAIN_NAME}_circlator/ data/#{STRAIN_NAME}_postcirc.fasta data/#{STRAIN_NAME}_postcirc2.fasta #{job_id}
  SH
end


# =========================
# = resequence_assembly_2 =
# =========================

desc "Resequences the circularized assembly"
task :resequence_assembly_2 => [:check, "data/#{STRAIN_NAME}_consensus_circ.fasta"]
file "data/#{STRAIN_NAME}_consensus_circ.fasta" => "data/#{STRAIN_NAME}_postcirc.fasta" do |t|
  abort "FATAL: Task resequence_assembly requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: STRAIN_NAME can only contain letters, numbers, and underscores" unless STRAIN_NAME =~ /^[\w]+$/
  
  mkdir_p "circularized_sequence"
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    referenceUploader -c -p circularized_sequence -n #{STRAIN_NAME} -f data/#{STRAIN_NAME}_postcirc.fasta
  SH
  cp "#{REPO_DIR}/xml/resequence_example_params.xml", OUT
  system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl resequence_example_params.xml " +
      "#{OUT} circularized_sequence/#{STRAIN_NAME} > resequence_params.xml" and
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    samtools faidx circularized_sequence/#{STRAIN_NAME}/sequence/#{STRAIN_NAME}.fasta &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
        -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params resequence_params.xml xml:bash5.xml &&
    gunzip data/consensus.fasta.gz
  SH
  cp "data/consensus.fasta", "data/#{STRAIN_NAME}_consensus_circ.fasta"
end

# ==============================
# = post_quiver_orient_correct =
# ==============================

desc "Renames the reoriented assembly contigs using a shortened scheme"
task :post_quiver_orient_correct => [:check, "data/#{STRAIN_NAME}_prokka.fasta"]
file "data/#{STRAIN_NAME}_prokka.fasta" => "data/#{STRAIN_NAME}_consensus_circ.fasta" do |t|
  system <<-SH
    module load blast
    #{REPO_DIR}/scripts/post_quiver_orient_correct.py data/#{STRAIN_NAME}_consensus_circ.fasta data/#{STRAIN_NAME}_postcirc2.fasta data/#{STRAIN_NAME}_prokka.fasta data/pq_dir
  SH
end



# ===================
# = prokka_annotate =
# ===================

desc "Annotates the reoriented assembly with prokka"
task :prokka_annotate => [:check, "data/prokka/#{STRAIN_NAME}_prokka.gbk"]
file "data/prokka/#{STRAIN_NAME}_prokka.gbk" => "data/#{STRAIN_NAME}_prokka.fasta" do |t|
  system <<-SH
    module load prokka  
    module load barrnap
    module unload rnammer
    module load minced
    module load signalp
        
    prokka --outdir data/prokka --force --prefix #{STRAIN_NAME}_prokka data/#{STRAIN_NAME}_prokka.fasta
  SH
end

# =====================
# = create_QC_webpage =
# =====================

desc "Creates the QC webpage"
task :create_QC_webpage => [:check, "data/www/index.html"]
file "data/www/index.html" => "data/#{STRAIN_NAME}_prokka.fasta" do |t|

  job_id = ENV['SMRT_JOB_ID']
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
    #{REPO_DIR}/scripts/create_QC_webpage.py -o data/qc_wd -w data/www -f data/#{STRAIN_NAME}_prokka.fasta \
     -g data -r data/corrected.fastq -a #{species_clean}_#{STRAIN_NAME}_#{job_id}
  SH
end

# =====================
# = Run QC and prokka =
# =====================

desc "Run prokka and create the QC website"
task :prokka_and_QC => [:prokka_annotate, :create_QC_webpage]
file "data/www/index.html" do |t|
end

# ==================
# = motif_and_mods =
# ==================

desc "Reruns SMRTPipe for modifcation and motif analysis on the reoriented assembly"
task :motif_and_mods => [:check, "data/motif_summary.csv"]
file "data/motif_summary.csv" => "data/#{STRAIN_NAME}_reorient.fasta" do |t|
  abort "FATAL: Task motif_and_mods requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: STRAIN_NAME can only contain letters, numbers, and underscores" unless STRAIN_NAME =~ /^[\w]+$/
  
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    referenceUploader -c -p reoriented_sequence -n #{STRAIN_NAME} -f data/#{STRAIN_NAME}_reorient.fasta
  SH
  cp "#{REPO_DIR}/xml/motif_simple_example_params.xml", OUT
  system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl motif_simple_example_params.xml " +
      "#{OUT} reoriented_sequence/#{STRAIN_NAME} > motif_simple_params.xml" and
  system <<-SH or abort
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh &&
    samtools faidx reoriented_sequence/#{STRAIN_NAME}/sequence/#{STRAIN_NAME}.fasta &&
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=12 -D CLUSTER=#{CLUSTER} \
        -D MAX_THREADS=16 #{CLUSTER != 'BASH' ? '--distribute' : ''} --params motif_simple_params.xml xml:bash5.xml
  SH
end


# =================
# = rast_annotate =
# =================

desc "Submits the circularized assembly to RAST for annotations"
task :rast_annotate => [:check, "data/#{STRAIN_NAME}_reorient_rast.fna", 
    "data/#{STRAIN_NAME}_reorient_rast.gbk", "data/#{STRAIN_NAME}_reorient_rast_aa.fa",
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

file "data/#{STRAIN_NAME}_reorient_rast.gbk" => ["data/#{STRAIN_NAME}_reorient.fasta"] do |t|
  submit_and_retrieve_rast("data/#{STRAIN_NAME}_reorient.fasta", t.name)
end
file "data/rast_job_id" => "data/#{STRAIN_NAME}_reorient_rast.gbk"

def gb_to_fasta(gb, fasta, seq_type=:nt, task_name="rast_annotate")
  abort "FATAL: Task #{task_name} requires specifying STRAIN_NAME" unless STRAIN_NAME
  abort "FATAL: gb_to_fasta called with invalid seq_type" unless [:nt, :aa].include? seq_type
  system <<-SH
    module load python/2.7.6
    module load py_packages/2.7
    python #{REPO_DIR}/scripts/gb_to_fasta.py -i #{gb} -s #{seq_type} -o #{fasta}
  SH
end

file "data/#{STRAIN_NAME}_reorient_rast_aa.fa" => "data/#{STRAIN_NAME}_reorient_rast.gbk" do |t|
  gb_to_fasta "data/#{STRAIN_NAME}_reorient_rast.gbk", "#{OUT}/#{t.name}", :aa
end

file "data/#{STRAIN_NAME}_reorient_rast.fna" => "data/#{STRAIN_NAME}_reorient_rast.gbk" do |t|
  gb_to_fasta "data/#{STRAIN_NAME}_reorient_rast.gbk", "#{OUT}/#{t.name}", :nt
end


# ================
# = improve_rast =
# ================

def improve_rast_genbank(input_gbk, output_gbk, task_name="improve_rast")
  abort "FATAL: Task #{task_name} requires specifying STRAIN_NAME" unless STRAIN_NAME 
  unless GENBANK_REFERENCES
    puts "WARN: No GenBank references for re-annotation provided, skipping this step"
    cp input_gbk, output_gbk
    return
  end
  
  # Once we want to integrate antibiotic resistance databases, we can start adding those in here
  # as improve_rask_gbk.rb supports them
  system <<-SH
    module load blast/2.2.26+
    ruby #{REPO_DIR}/scripts/improve_rast_gbk.rb \
        #{GENBANK_REFERENCES.map{|f| Shellwords.escape f }.join(' ')} #{Shellwords.escape input_gbk}\
        > #{Shellwords.escape output_gbk}
  SH
end

file "data/#{STRAIN_NAME}_reorient_rast_reannotate.gbk" => ["data/#{STRAIN_NAME}_reorient_rast.gbk",
    "data/#{STRAIN_NAME}_reorient_rast_aa.fa", "data/#{STRAIN_NAME}_reorient_rast.fna"] do |t|
  improve_rast_genbank("data/#{STRAIN_NAME}_reorient_rast.gbk", t.name)
end

desc "Improves GenBank output from RAST by re-annotating gene names from better references"
task :improve_rast => [:check, "data/#{STRAIN_NAME}_reorient_rast_reannotate.gbk"]


# ===============
# = rast_to_igb =
# ===============

job_id = ENV['SMRT_JOB_ID']
species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES

directory IGB_DIR

desc "Creates an IGB Quickload-compatible directory for your genome in IGB_DIR"
task :rast_to_igb => [:check, "data/#{STRAIN_NAME}_reorient_rast_reannotate.gbk"] do |t|
  abort "FATAL: Task rast_to_igb requires specifying SMRT_JOB_ID" unless job_id
  abort "FATAL: Task rast_to_igb requires specifying STRAIN_NAME" unless STRAIN_NAME 
  abort "FATAL: Task rast_to_igb requires specifying SPECIES" unless SPECIES 
  
  system <<-SH
    module load blat
    module load bioperl
    export SAS_DIR=#{SAS_DIR}
    perl #{REPO_DIR}/scripts/rast2igb.pl \
        -f data/#{STRAIN_NAME}_reorient_rast_reannotate.gbk \
        -g #{species_clean}_#{STRAIN_NAME}_#{job_id} \
        -i #{IGB_DIR}
  SH
end

# ===================
# = all =
# ===================

desc "Runs entire pipeline from top to bottom"
task :all => [:rast_to_igb, :motif_and_mods]
file "bash5.fofn" do |t|
end


# ====================================================================================================
# = The following tasks are for assemblies where we want to incorporate Illumina reads to fix indels =
# ====================================================================================================

namespace :ilm do

  # ====================
  # = ilm:fake_prereqs =
  # ====================
  desc "Fakes the prerequisites for the Illumina tasks"
  task :fake_prereqs do
    abort "FATAL: Task ilm:fake_prereqs requires specifying STRAIN_NAME" unless STRAIN_NAME 
    mkdir_p "log"
    mkdir_p "data"
    touch "bash5.fofn"                                  and sleep 1
    touch "data/polished_assembly.fasta.gz"             and sleep 1
    touch "data/polished_assembly_circularized.fasta"   and sleep 1
    if ILLUMINA_REFERENCE
      abort "FATAL: file '#{ILLUMINA_REFERENCE}' does not exist" unless File.exists?(ILLUMINA_REFERENCE)
      cp ILLUMINA_REFERENCE, "data/#{STRAIN_NAME}_consensus.fasta"
      cp ILLUMINA_REFERENCE, "data/#{STRAIN_NAME}_reorient.fasta"
    else
      touch "data/#{STRAIN_NAME}_consensus.fasta"         and sleep 1
      touch "data/#{STRAIN_NAME}_reorient.fasta"
      puts "Replace #{OUT}/data/#{STRAIN_NAME}_reorient.fasta with the old reference sequence you are piling new reads onto."
    end
  end
    

  # ========================
  # = ilm:recall_consensus =
  # ========================

  desc "Recalls a new consensus by piling Illumina reads onto a PacBio assembly"
  task :recall_consensus => [:check, "data/#{STRAIN_NAME}_ref_flt.vcf", "data/#{STRAIN_NAME}_ilm_reorient.fasta"]

  file "data/ref.sort.bam" => "data/#{STRAIN_NAME}_reorient.fasta" do |t|
    abort "FATAL: Task ilm:recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
    abort "FATAL: Task ilm:recall_consensus requires specifying ILLUMINA_FASTQ" unless ILLUMINA_FASTQ
  
    LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
    LSF.job_name "ref.aln.sam"
    LSF.bsub_interactive <<-SH or abort
      module load bwa/0.7.8
      bwa index "data/#{STRAIN_NAME}_reorient.fasta"
      bwa mem "data/#{STRAIN_NAME}_reorient.fasta" #{Shellwords.escape(ILLUMINA_FASTQ)} > data/ref.aln.sam
    
      module load samtools/1.1
      samtools view -bS data/ref.aln.sam > data/ref.aln.bam
      samtools sort data/ref.aln.bam data/ref.sort
    SH
  
    # Can remove the SAM as it is huge and the .aln.bam should contain everything in it
    rm "data/ref.aln.sam"
  end

  file "data/#{STRAIN_NAME}_ref_raw.bcf" => "data/ref.sort.bam" do |t|
    abort "FATAL: Task ilm:recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
    # Use mpileup to do the consensus calling.
    # Note: the -L and -d flags are important; they ensure samtools looks at up to 100k reads per base to call variants.
    # The default for -L is 250, which would turn off indel calling for deeply resequenced (depth >250) samples.
    LSF.set_out_err("log/recall_ilm_consensus.log", "log/recall_ilm_consensus.err.log")
    LSF.job_name "#{STRAIN_NAME}_ref_raw.bcf"
    LSF.bsub_interactive <<-SH
      module load samtools/1.1
      module load bcftools/1.1
      samtools mpileup -L100000 -d100000 -uf "data/#{STRAIN_NAME}_reorient.fasta" data/ref.sort.bam \
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

  file "data/#{STRAIN_NAME}_ilm_reorient.fasta" => 
      ["data/#{STRAIN_NAME}_ref_flt.vcf", "data/#{STRAIN_NAME}_reorient.fasta"] do |t|
    system <<-SH
      module load vcftools/0.1.12b
      module load tabix/0.2.6 
    
      bgzip -c "data/#{STRAIN_NAME}_ref_flt.vcf" > "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
      tabix -p vcf "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
      cat "data/#{STRAIN_NAME}_reorient.fasta" | vcf-consensus "data/#{STRAIN_NAME}_ref_flt.vcf.gz" \
              > "data/#{STRAIN_NAME}_ilm_reorient.fasta"
    
      # New-style version of doing this with bcftools consensus, but it doesn't work (memory leak in bcftools)
      # #{HTSLIB_DIR}/bgzip -c "data/#{STRAIN_NAME}_ref_flt.vcf" > "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
      # #{HTSLIB_DIR}/tabix -p vcf "data/#{STRAIN_NAME}_ref_flt.vcf.gz"
      # #{BCFTOOLS_DIR}/bcftools consensus -f "data/#{STRAIN_NAME}_reorient.fasta" "data/#{STRAIN_NAME}_ref_flt.vcf.gz" \
      #    > "data/#{STRAIN_NAME}_ilm_reorient.fasta"
    SH
  end

  # =====================
  # = ilm:rast_annotate =
  # =====================

  desc "Submits the Illumina-fixed consensus to RAST for annotations"
  task :rast_annotate => [:check, "data/#{STRAIN_NAME}_ilm_reorient_rast.fna", 
      "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", "data/#{STRAIN_NAME}_ilm_reorient_rast_aa.fa",
      "data/ilm_rast_job_id"]

  file "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk" => ["data/#{STRAIN_NAME}_ilm_reorient.fasta"] do |t|
    fasta = "data/#{STRAIN_NAME}_ilm_reorient.fasta"
    submit_and_retrieve_rast(fasta, t.name, "ilm_rast_job_id", "rast_annotate_ilm")
  end
  file "data/ilm_rast_job_id" => "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk"

  file "data/#{STRAIN_NAME}_ilm_reorient_rast_aa.fa" => "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk" do |t|
    gb_to_fasta "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", "#{OUT}/#{t.name}", :aa, "rast_annotate_ilm"
  end

  file "data/#{STRAIN_NAME}_ilm_reorient_rast.fna" => "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk" do |t|
    gb_to_fasta "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", "#{OUT}/#{t.name}", :nt, "rast_annotate_ilm"
  end


  # ====================
  # = ilm:improve_rast =
  # ====================

  desc "Improves GenBank output from RAST by re-annotating gene names from better references"
  task :improve_rast => [:check, "data/#{STRAIN_NAME}_ilm_reorient_rast_reannotate.gbk"]

  file "data/#{STRAIN_NAME}_ilm_reorient_rast_reannotate.gbk" => ["data/#{STRAIN_NAME}_ilm_reorient_rast.gbk"] do |t|
    improve_rast_genbank("data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", t.name)
  end


  # ===================
  # = ilm:rast_to_igb =
  # ===================

  desc "Creates an IGB Quickload-compatible directory for your Illumina-fixed genome in IGB_DIR"
  task :rast_to_igb => [:check, "data/#{STRAIN_NAME}_ilm_reorient_rast_reannotate.gbk"] do |t|
    abort "FATAL: Task ilm:rast_to_igb requires specifying SMRT_JOB_ID" unless job_id
    abort "FATAL: Task ilm:rast_to_igb requires specifying STRAIN_NAME" unless STRAIN_NAME 
    abort "FATAL: Task ilm:rast_to_igb requires specifying SPECIES" unless SPECIES 
  
    system <<-SH
      module load blat
      module load bioperl
      export SAS_DIR=#{SAS_DIR}
      perl #{REPO_DIR}/scripts/rast2igb.pl \
          -f data/#{STRAIN_NAME}_ilm_reorient_rast_reannotate.gbk \
          -g #{species_clean}_#{STRAIN_NAME}_#{job_id} \
          -i #{IGB_DIR}
    SH
  end

end # namespace :ilm
