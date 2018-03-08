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
  
  task :check => [:check, "old:sas"]
  
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
  
  ##
  ## Circularization is now handled by Circlator (instead of custom scripts w/ nucmer)
  ##
  
  # ============================
  # = old:circularize_assembly =
  # ============================

  desc "deprecated - Circularizes the PacBio assembly"
  task :circularize_assembly => ["old:check", "data/polished_assembly_circularized.fasta"]
  file "data/polished_assembly_circularized.fasta" => "data/polished_assembly.fasta.gz" do |t|
    system "gunzip -c data/polished_assembly.fasta.gz >data/polished_assembly.fasta" and
    system "#{REPO_DIR}/scripts/circularizeContigs.pl -i data/polished_assembly.fasta -l 12000 2> polished_assembly_circularized.log"
  end

  # ===========================
  # = old:resequence_assembly =
  # ===========================

  desc "deprecated - Resequences the circularized assembly"
  task :resequence_assembly => ["old:check", "data/#{STRAIN_NAME}_consensus.fasta"]
  file "data/#{STRAIN_NAME}_consensus.fasta" => "data/polished_assembly_circularized.fasta" do |t|
    abort "FATAL: Task old:resequence_assembly requires specifying STRAIN_NAME" unless STRAIN_NAME 
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


  # =========================
  # = old:reorient_assembly =
  # =========================

  desc "deprecated - Reorients the circularized assembly to a given locus"
  task :reorient_assembly => ["old:check", "data/#{STRAIN_NAME}_reorient.fasta"]
  file "data/#{STRAIN_NAME}_reorient.fasta" => "data/#{STRAIN_NAME}_consensus.fasta" do |t|
    abort "FATAL: Task old:reorient_assembly requires specifying STRAIN_NAME" unless STRAIN_NAME 
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
        abort "FATAL: Task old:reorient_assembly failed, exiting"
      end
    else
      # No reorient locus given, simply copy the assembly for the next step
      cp "data/#{STRAIN_NAME}_consensus.fasta", "data/#{STRAIN_NAME}_reorient.fasta"
    end
  end
  
  
  ##
  ## Annotation is now handled by prokka (instead of RAST)
  ##
  
  # =====================
  # = old:rast_annotate =
  # =====================

  desc "deprecated - Submits the circularized assembly to RAST for annotations"
  task :rast_annotate => ["old:check", "data/#{STRAIN_NAME}_reorient_rast.fna", 
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
    submit_and_retrieve_rast("data/#{STRAIN_NAME}_reorient.fasta", t.name, "rast_job_id", "old:rast_annotate")
  end
  file "data/rast_job_id" => "data/#{STRAIN_NAME}_reorient_rast.gbk"

  file "data/#{STRAIN_NAME}_reorient_rast_aa.fa" => "data/#{STRAIN_NAME}_reorient_rast.gbk" do |t|
    gb_to_fasta "data/#{STRAIN_NAME}_reorient_rast.gbk", "#{OUT}/#{t.name}", :aa, "old:rast_annotate"
  end

  file "data/#{STRAIN_NAME}_reorient_rast.fna" => "data/#{STRAIN_NAME}_reorient_rast.gbk" do |t|
    gb_to_fasta "data/#{STRAIN_NAME}_reorient_rast.gbk", "#{OUT}/#{t.name}", :nt, "old:rast_annotate"
  end


  # ====================
  # = old:improve_rast =
  # ====================

  def improve_rast_genbank(input_gbk, output_gbk, task_name="improve_rast")
    abort "FATAL: Task #{task_name} requires specifying STRAIN_NAME" unless STRAIN_NAME 
    unless GENBANK_REFERENCES
      puts "WARN: No GenBank references for re-annotation provided, skipping this step"
      cp input_gbk, output_gbk
      return
    end
  
    # Once we want to integrate antibiotic resistance databases, we can start adding those in here
    # as improve_rast_gbk.rb supports them
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

  desc "deprecated - Improves GenBank output from RAST by re-annotating gene names from better references"
  task :improve_rast => ["old:check", "data/#{STRAIN_NAME}_reorient_rast_reannotate.gbk"]


  # ===================
  # = old:rast_to_igb =
  # ===================

  directory IGB_DIR

  desc "deprecated - Creates an IGB Quickload-compatible directory for your genome in IGB_DIR"
  task :rast_to_igb => ["old:check", "data/#{STRAIN_NAME}_reorient_rast_reannotate.gbk"] do |t|
    job_id = ENV['SMRT_JOB_ID']
    abort "FATAL: Task #{t.name} requires specifying SMRT_JOB_ID" unless job_id
    abort "FATAL: Task #{t.name} requires specifying STRAIN_NAME" unless STRAIN_NAME 
    abort "FATAL: Task #{t.name} requires specifying SPECIES" unless SPECIES 
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
  
  
  namespace :ilm do
    
    ##
    ## Steps related to Illumina data that branch off of the old steps' outputs instead of the
    ## the circlator and prokka-annotated assembly
    ##
    
    # ========================
    # = old:ilm:fake_prereqs =
    # ========================
  
    desc "deprecated - Fakes the prerequisites for the old Illumina tasks"
    task :fake_prereqs do
      abort "FATAL: Task old:ilm:fake_prereqs requires specifying STRAIN_NAME" unless STRAIN_NAME 
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
  
  
    # ============================
    # = old:ilm:recall_consensus =
    # ============================

    desc "deprecated - Recalls a new consensus by piling Illumina reads onto a PacBio assembly"
    task :recall_consensus => ["old:check", "data/#{STRAIN_NAME}_ref_flt.vcf", "data/#{STRAIN_NAME}_ilm_reorient.fasta"]

    file "data/ref.sort.bam" => "data/#{STRAIN_NAME}_reorient.fasta" do |t|
      abort "FATAL: Task old:ilm:recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
      abort "FATAL: Task old:ilm:recall_consensus requires specifying ILLUMINA_FASTQ" unless ILLUMINA_FASTQ
  
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
      abort "FATAL: Task old:ilm:recall_consensus requires specifying STRAIN_NAME" unless STRAIN_NAME 
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

    # =========================
    # = old:ilm:rast_annotate =
    # =========================

    desc "deprecated - Submits the Illumina-fixed consensus to RAST for annotations"
    task :rast_annotate => ["old:check", "data/#{STRAIN_NAME}_ilm_reorient_rast.fna", 
        "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", "data/#{STRAIN_NAME}_ilm_reorient_rast_aa.fa",
        "data/ilm_rast_job_id"]

    file "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk" => ["data/#{STRAIN_NAME}_ilm_reorient.fasta"] do |t|
      fasta = "data/#{STRAIN_NAME}_ilm_reorient.fasta"
      submit_and_retrieve_rast(fasta, t.name, "ilm_rast_job_id", "old:ilm:rast_annotate")
    end
    file "data/ilm_rast_job_id" => "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk"

    file "data/#{STRAIN_NAME}_ilm_reorient_rast_aa.fa" => "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk" do |t|
      gb_to_fasta "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", "#{OUT}/#{t.name}", :aa, "old:ilm:rast_annotate"
    end

    file "data/#{STRAIN_NAME}_ilm_reorient_rast.fna" => "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk" do |t|
      gb_to_fasta "data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", "#{OUT}/#{t.name}", :nt, "old:ilm:rast_annotate"
    end


    # ========================
    # = old:ilm:improve_rast =
    # ========================

    desc "deprecated - Improves GenBank output from RAST by re-annotating gene names from better references"
    task :improve_rast => ["old:check", "data/#{STRAIN_NAME}_ilm_reorient_rast_reannotate.gbk"]

    file "data/#{STRAIN_NAME}_ilm_reorient_rast_reannotate.gbk" => ["data/#{STRAIN_NAME}_ilm_reorient_rast.gbk"] do |t|
      improve_rast_genbank("data/#{STRAIN_NAME}_ilm_reorient_rast.gbk", t.name)
    end


    # =======================
    # = old:ilm:rast_to_igb =
    # =======================

    desc "deprecated - Creates an IGB Quickload-compatible directory for your Illumina-fixed genome in IGB_DIR"
    task :rast_to_igb => ["old:check", "data/#{STRAIN_NAME}_ilm_reorient_rast_reannotate.gbk"] do |t|
      abort "FATAL: Task #{t.name} requires specifying SMRT_JOB_ID" unless job_id
      abort "FATAL: Task #{t.name} requires specifying STRAIN_NAME" unless STRAIN_NAME 
      abort "FATAL: Task #{t.name} requires specifying SPECIES" unless SPECIES 
  
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

end