require 'pp'
require 'net/http'
require_relative 'lib/colors'
require_relative 'lib/lsf_client'
require_relative 'lib/subscreens'
require 'shellwords'
include Colors

########
##
## The following are steps in the pipeline that we no longer need because they have been replaced with superior processes
## They have been moved into an "old" namespace
##
########

def gb_to_fasta(gb, fasta, seq_type=:nt, task_name="rast_annotate")
  abort "FATAL: Task #{task_name} requires specifying STRAIN_NAME" unless STRAIN_NAME
  abort "FATAL: gb_to_fasta called with invalid seq_type" unless [:nt, :aa].include? seq_type
  system <<-SH
    module load python/2.7.6
    module load py_packages/2.7
    python #{REPO_DIR}/scripts/gb_to_fasta.py -i #{gb} -s #{seq_type} -o #{fasta}
  SH
end

##
## Deprecated steps are now encapsulated in the "old" namespace
##

namespace :old do
  
  ##
  ## Circularization is now handled by Circlator (instead of custom scripts w/ nucmer)
  ##
  
  # ========================
  # = circularize_assembly =
  # ========================

  desc "deprecated - Circularizes the PacBio assembly"
  task :circularize_assembly => [:check, "data/polished_assembly_circularized.fasta"]
  file "data/polished_assembly_circularized.fasta" => "data/polished_assembly.fasta.gz" do |t|
    system "gunzip -c data/polished_assembly.fasta.gz >data/polished_assembly.fasta" and
    system "#{REPO_DIR}/scripts/circularizeContigs.pl -i data/polished_assembly.fasta -l 12000 2> polished_assembly_circularized.log"
  end

  # =======================
  # = resequence_assembly =
  # =======================

  desc "deprecated - Resequences the circularized assembly"
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

  desc "deprecated - Reorients the circularized assembly to a given locus"
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
  
  
  ##
  ## Annotation is now handled by prokka (instead of RAST)
  ##
  
  # =================
  # = rast_annotate =
  # =================

  desc "deprecated - Submits the circularized assembly to RAST for annotations"
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

  directory IGB_DIR

  desc "deprecated - Creates an IGB Quickload-compatible directory for your genome in IGB_DIR"
  task :rast_to_igb => [:check, "data/#{STRAIN_NAME}_reorient_rast_reannotate.gbk"] do |t|
    job_id = ENV['SMRT_JOB_ID']
    abort "FATAL: Task rast_to_igb requires specifying SMRT_JOB_ID" unless job_id
    abort "FATAL: Task rast_to_igb requires specifying STRAIN_NAME" unless STRAIN_NAME 
    abort "FATAL: Task rast_to_igb requires specifying SPECIES" unless SPECIES 
    species_clean = (SPECIES && SPECIES != '${SPECIES}') ? SPECIES.gsub(/[^a-z_]/i, "_") : SPECIES
  
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

end