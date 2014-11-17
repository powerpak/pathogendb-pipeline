require 'pp'
require_relative 'lib/colors'
include Colors

# DEPENDENCIES
# Centralized
#   /hpc/users/attieo02/gitrepos/multiscale/bacterial_analysis/src/ccs_get.py
#   /hpc/packages/minerva-mothra/smrtanalysis/2.2.0/ROOT/current/etc/setup.sh
#   /sc/orga/projects/InfectiousDisease/smrtpipe/example_params.xml
#   /sc/orga/projects/InfectiousDisease/smrtpipe/resequence_example_params.xml

# Oliver's
#   
#   All of the following come from http://blog.theseed.org/downloads/sas.tgz
#   /sc/orga/work/attieo02/sas/lib
#   /sc/orga/work/attieo02/sas/modules/lib
#   /sc/orga/work/attieo02/sas/plbin/svr_submit_status_retrieve.pl
#
#   Now in the repo
#   /sc/orga/work/attieo02/genbank_to_fasta_v1.1/gb_to_fasta.py

task :default => :check

REPO_DIR = File.dirname(__FILE__)
OUT = ENV['OUT'] || "#{REPO_DIR}/out"

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

desc "Checks environment variables and requirements before running tasks"
task :check => ["#{REPO_DIR}/scripts/env.sh", :env] do
  env_error = "Configure this in scripts/env.sh and run `source scripts/env.sh` before running rake."
  unless `module avail 2>&1 | grep smrtpipe/2.2.0` != ''
    abort "You must have the smrtpipe/2.2.0 module in your MODULEPATH."
  end
  unless ENV['SMRTPIPE'] && File.exists?("#{ENV['SMRTPIPE']}/example_params.xml")
    abort "SMRTPIPE must be set to the directory containing example_params.xml for smrtpipe.py.\n#{env_error}"
  end
  unless ENV['SMRTANALYSIS'] && File.exists?("#{ENV['SMRTANALYSIS']}/etc/setup.sh")
    abort <<-ERRMSG
      SMRTANALYSIS must be set to the ROOT directory for the SMRT Analysis package, v2.2.0.
      This software can be downloaded from http://www.pacb.com/devnet/
      #{env_error}"
    ERRMSG
  end
end

# pulls down http://blog.theseed.org/downloads/sas.tgz --> ./vendor/sas
#   then it adds SAS libs to PERL5LIB
#   then it adds SAS bins to PATH
task :sas => [:env, "#{REPO_DIR}/vendor/sas/sas.tgz"] do
  ENV['PERL5LIB'] = "#{ENV['PERL5LIB']}:#{REPO_DIR}/vendor/sas/lib:#{REPO_DIR}/vendor/sas/modules/lib"
  ENV['PATH'] = "#{REPO_DIR}/vendor/sas/bin:#{ENV['PATH']}"
end

directory "#{REPO_DIR}/vendor/sas"
file "#{REPO_DIR}/vendor/sas" => :env
file "#{REPO_DIR}/vendor/sas/sas.tgz" => [:env, "#{REPO_DIR}/vendor/sas"] do |t|
  Dir.chdir(File.dirname(t.name)) do
    system("curl -O 'http://blog.theseed.org/downloads/sas.tgz'")
    system("tar xvzf sas.tgz")
  end
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
    python #{REPO_DIR}/scripts/ccs_get.py --noprefix -e bax.h5 #{job_id} -i &&
    find #{OUT}/*bax.h5 > bash5.fofn
  SH
  # NOTE: we will change the above to not fetch the full sequence, but rather symlink to it on minerva, like so
  #       we could even skip straight to circularize_assembly if the polished_assembly.fasta is already there
  # cp /sc/orga/projects/pacbio/userdata_permanent/jobs/#{job_id[0..3]}/#{job_id}/input.fofn baxh5.fofn
  # mkdir_p "data"
  # ln -s /sc/orga/projects/pacbio/userdata_permanent/jobs/data/#{job_id[0..3]}/#{job_id}/polished_assembly.fasta <input assembly file name>

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
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=16 -D CLUSTER=LSF -D MAX_THREADS=16 --distribute --params example_params.xml xml:bash5.xml 
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
task :resequence_assembly => [:check, "data/consensus.fasta"]
file "data/consensus.fasta" => "data/polished_assembly_circularized.fasta" do |t|
  strain_name = ENV['STRAIN_NAME'] 
  abort "FATAL: Task resequence_assembly requires specifying STRAIN_NAME" unless strain_name 
  
  mkdir_p "circularized_sequence"
  system <<-SH
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh
    referenceUploader -c -p circularized_sequence -n #{strain_name} -f data/polished_assembly_circularized.fasta
  SH
  system "cp #{ENV['SMRTPIPE']}/resequence_example_params.xml \."
  system "perl #{REPO_DIR}/scripts/changeResequencingDirectory.pl resequence_example_params.xml " +
      "#{OUT} circularized_sequence/#{strain_name} > resequence_params.xml"
  system <<-SH
    module load smrtpipe/2.2.0
    source #{ENV['SMRTANALYSIS']}/etc/setup.sh
    samtools faidx circularized_sequence/#{strain_name}/sequence/#{strain_name}.fasta
    smrtpipe.py -D TMP=#{ENV['TMP']} -D SHARED_DIR=#{ENV['SHARED_DIR']} -D NPROC=16 -D CLUSTER=LSF -D MAX_THREADS=16 --distribute --params resequence_params.xml xml:bash5.xml
  SH
  system "gunzip data/consensus.fasta.gz"
  system "cp data/consensus.fasta data/#{strain_name}_consensus.fasta"
end


# =================
# = rast_annotate =
# =================


# Creates data/#{strain_name}_consensus.fasta
# Creates data/#{strain_name}_consensus_rast.gbk
# Creates data/#{strain_name}_consensus_rast_aa.fa
# Creates data/#{strain_name}_consensus_rast.fna
def rast_annotate(strain_name, dir1, species)
  rast_job=%x[export PERL5LIB=$PERL5LIB:/sc/orga/work/attieo02/sas/lib:/sc/orga/work/attieo02/sas/modules/lib;export PATH=$PATH:/sc/orga/work/attieo02/sas/bin;perl /sc/orga/work/attieo02/sas/plbin/svr_submit_status_retrieve.pl --user oattie --passwd sessiz_ev --fasta #{dir1}/data/#{strain_name}_consensus.fasta --domain Bacteria --bioname "#{species} #{strain_name}" --genetic_code 11 --gene_caller rast]
  system("export PERL5LIB=$PERL5LIB:/sc/orga/work/attieo02/sas/lib:/sc/orga/work/attieo02/sas/modules/lib;export PATH=$PATH:/sc/orga/work/attieo02/sas/bin; perl /sc/orga/work/attieo02/sas/test_server.pl oattie sessiz_ev genbank #{rast_job} ")
  sleep(120)
  system("export PERL5LIB=$PERL5LIB:/sc/orga/work/attieo02/sas/lib:/sc/orga/work/attieo02/sas/modules/lib;export PATH=$PATH:/sc/orga/work/attieo02/sas/bin; svr_retrieve_RAST_job oattie sessiz_ev #{rast_job} genbank > #{dir1}/data/#{strain_name}_consensus_rast.gbk")
  system("python /sc/orga/work/attieo02/genbank_to_fasta_v1.1/gb_to_fasta.py -i #{dir1}/data/#{strain_name}_consensus_rast.gbk -s \'aa\' -o #{dir1}/data/#{strain_name}_consensus_rast_aa.fa")
  system("python /sc/orga/work/attieo02/genbank_to_fasta_v1.1/gb_to_fasta.py -i #{dir1}/data/#{strain_name}_consensus_rast.gbk -s \'nt\' -o #{dir1}/data/#{strain_name}_consensus_rast.fna")
end

# task :resequence, :job_id, :dir1, :strain, :species do |t, args|
#   job_id=args[:job_id]
#   dir1=args[:dir1]
#   strain_name=args[:strain]
#   species=args[:species]
#   pull_down_raw_reads(job_id,dir1)
#   assemble_raw_reads(dir1)
#   circularize_assembly(dir1)
#   resequence_assembly(strain_name, dir1, species)
# end
#
# task :circularize do
#   circularize_assembly
# end
