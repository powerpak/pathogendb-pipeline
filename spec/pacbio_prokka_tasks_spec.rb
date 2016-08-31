require_relative './spec_helper'
require 'shellwords'
require 'fileutils'
require 'pp'

describe "pathogendb-pipeline" do
  
  around(:each) do |example|
    with_tempdir('OUT') do |dir|
      @out = dir
      example.run
      !(ENV['DEBUG'] && example.exception)
    end
  end
  
  context "when running tasks from scratch" do
  
    describe "task :pull_down_raw_reads" do
      it "pulls down a polished assembly for a job in SMRT Portal (023460)" do
        run "SMRT_JOB_ID=023460 rake --silent pull_down_raw_reads"
        expect(File).to exist("#{@out}/bash5.fofn")
        expect(File).to exist("#{@out}/data/polished_assembly.fasta.gz")
      end
      
      it "pulls down a polished assembly for a job in old_smrtportal_jobs (019203)" do
        run "SMRT_JOB_ID=019203 rake --silent pull_down_raw_reads"
        expect(File).to exist("#{@out}/bash5.fofn")
        expect(File).to exist("#{@out}/data/polished_assembly.fasta.gz")
      end
    end
    
    describe "task :assemble_raw_reads" do
      it "skips assembly for a job already assembled by SMRT Portal (023460)" do
        stdout = run "SMRT_JOB_ID=023460 rake --silent assemble_raw_reads"
        expect(File).to exist("#{@out}/data/polished_assembly.fasta.gz")
        expect(File).not_to exist("#{@out}/example_params.xml")
        expect(stdout).to match(/NOTICE: .* skipping assemble_raw_reads/)
      end
    end
    
    describe "task :run_circlator", :speed => 'slow' do
      it "produces a circularized FASTA file with md5 64ae02c8... as output (for 019203)" do
        strain = 'SM5478'
        run "SMRT_JOB_ID=019203 STRAIN_NAME=#{strain} rake --silent run_circlator"
        fasta_out = "#{@out}/data/#{strain}_circlator/06.fixstart.fasta"
        expect(File).to exist(fasta_out)
        expect(md5(fasta_out)).to eq('64ae02c8a8e68dc2134899b4281f93ca')
      end
    end
    
    describe "task :resequence_assembly", :speed => 'slow' do
      it "produces a consensus circularized FASTA file with md5 f6f1095c... as output (for 023460)" do
        strain = 'PS00098'
        run "SMRT_JOB_ID=023460 STRAIN_NAME=#{strain} CLUSTER=BASH rake --silent resequence_assembly"
        fasta_out = "#{@out}/data/#{strain}_consensus_circ.fasta"
        expect(File).to exist(fasta_out)
        expect(md5(fasta_out)).to eq('f6f1095c1848b7212450be0cc4937be9')
      end
    end
    
  end
  
  context "when running tasks on SM5478 (019203), after :run_circlator" do
    
    before(:each) do
      @strain = "SM5478"
      @env_vars = "SMRT_JOB_ID=019203 STRAIN_NAME=#{@strain}"
      copy_example "#{@strain}-run_circlator"
      touch_prereqs :run_circlator, "SMRT_JOB_ID=019203 STRAIN_NAME=#{@strain}"
    end
    
    describe "task :post_circlator", :speed => 'medium' do
      it "renames the first contig to u00000crxx_c_019203" do
        run "#{@env_vars} rake --silent post_circlator"
        fasta_out = "#{@out}/data/#{@strain}_postcirc.fasta"
        expect(File).to exist(fasta_out)
        first_line = File.open(fasta_out, &:readline).strip
        expect(first_line).to eq('>u00000crxx_c_019203')
      end
      
      context "when circlator input is absent" do
        before(:each) do
          FileUtils.rm "#{@out}/data/#{@strain}_circlator/00.input_assembly.fasta"
        end
      
        it "considers the assembly curated and renames the first contig to u00000crxm_c_019203" do
          run "#{@env_vars} rake --silent post_circlator"
          fasta_out = "#{@out}/data/#{@strain}_postcirc.fasta"
          expect(File).to exist(fasta_out)
          first_line = File.open(fasta_out, &:readline).strip
          expect(first_line).to eq('>u00000crxm_c_019203')
        end
      end
    end
    
  end
  
  context "when running tasks on PS00098 (023460), after :resequence_assembly" do
    
    before(:each) do
      @strain = "PS00098"
      copy_example "#{@strain}-resequence_assembly"
      touch_prereqs :resequence_assembly, "SMRT_JOB_ID=023460 STRAIN_NAME=#{@strain}"
    end
    
    describe "task :post_quiver_orient_correct", :speed => 'medium' do
      it "produces a reoriented FASTA file with md5 ee0ca0ee... as output" do
        run "SMRT_JOB_ID=019203 STRAIN_NAME=#{@strain} rake --silent post_quiver_orient_correct"
        fasta_out = "#{@out}/data/#{@strain}_prokka.fasta"
        expect(File).to exist(fasta_out)
        expect(md5(fasta_out)).to eq('ee0ca0eec43f0ca87b8a27c65af11640')
      end
    end
    
  end

  context "when running tasks on PS00098 (023460), after :post_quiver_orient_correct" do
    
    before(:each) do
      @strain = "PS00098"
      @env_vars = "SMRT_JOB_ID=019203 SPECIES='Serratia marcescens' STRAIN_NAME=#{@strain}"
      copy_example "#{@strain}-post_quiver_orient_correct"
      touch_prereqs :post_quiver_orient_correct, @env_vars
    end
    
    describe "task :prokka_annotate", :speed => 'medium' do
      it "produces an annotated GenBank file with md5 82490819... for everything but the first line" do
        run "#{@env_vars} rake --silent prokka_annotate"
        gbk_out = "#{@out}/data/prokka/#{@strain}_prokka.gbk"
        expect(File).to exist(gbk_out)
        # When calculating an MD5 for the GenBank output we have to exclude the 1st (header) line
        # because it contains the file creation date.
        gbk_tail = "#{@out}/data/prokka/#{@strain}_prokka.gbk.tail"
        `tail -n +2 #{Shellwords.escape gbk_out} >#{Shellwords.escape gbk_tail}`
        expect(md5(gbk_tail)).to eq('8249081908c7245f4a72fe9b06a8c9f5')
      end
    end
    
    # NOTE that for the purposes of this test we've truncated data/corrected.fastq
    # to a more reasonable size.
    describe "task :create_QC_webpage", :speed => 'medium' do
      it "produces a QC webpage with a dotplot, graphs, and BLAST results" do
        run "#{@env_vars} rake --silent create_QC_webpage"
        html_out = "#{@out}/data/www"
        expect(File).to exist("#{html_out}/index.html")
        expect(File).to exist("#{html_out}/qc_website/dot_plot.svg")
        expect(File).to exist("#{html_out}/qc_website/graph.svg")
        expect(File).to exist("#{html_out}/qc_website/blast/u00000_blast.html")
        expect(File).to exist("#{html_out}/qc_website/graphs/u00000_graphs.html")
      end
    end
    
  end
  
  context "when running tasks on PS00098 (023460), after :prokka_and_QC" do
    
    before(:each) do
      @strain = "PS00098"
      @env_vars = "SMRT_JOB_ID=023460 SPECIES='Serratia marcescens' STRAIN_NAME=#{@strain}"
      copy_example "#{@strain}-prokka_and_QC"
      touch_prereqs :prokka_and_QC, @env_vars
    end
    
    around(:each) do |example|
      with_tempdir('IGB_DIR') do |dir|
        @igb_dir = dir
        FileUtils.touch("#{@igb_dir}/contents.txt")
        example.run
        !(ENV['DEBUG'] && example.exception)
      end
    end
    
    describe "task :prokka_to_igb" do
      
      before(:each) do
        run "#{@env_vars} rake --silent prokka_to_igb"
        @genome_name = "Serratia_marcescens_PS00098_023460"
      end
      
      it "produces a correctly named directory in IGB_DIR for the genome" do
        expect(Dir).to exist(@igb_dir)
        expect(Dir).to exist("#{@igb_dir}/#{@genome_name}")
      end
      
      it "adds the genome to contents.txt in the IGB directory" do
        contents_txt = File.open("#{@igb_dir}/contents.txt", &:readline).strip
        expect(contents_txt).to eq("#{@genome_name}\t#{@genome_name}")
      end
      
      it "adds the contigs to genome.txt in the IGB genome's directory" do
        genome_txt = File.open("#{@igb_dir}/#{@genome_name}/genome.txt", &:readline).strip
        expect(genome_txt).to eq("u00000crpx_c_023460\t5400688")
      end
      
      it "creates a .2bit, .bed, and .fasta file in the IGB genome's directory" do
        expect(File).to exist("#{@igb_dir}/#{@genome_name}/#{@genome_name}.2bit")
        expect(File).to exist("#{@igb_dir}/#{@genome_name}/#{@genome_name}.bed")
        expect(File).to exist("#{@igb_dir}/#{@genome_name}/#{@genome_name}.fasta")
      end
      
      it "produces an annots.xml file for the genome that can be parsed as XML" do
        expect(File).to exist("#{@igb_dir}/#{@genome_name}/annots.xml")
        annots_txt = File.open("#{@igb_dir}/#{@genome_name}/annots.xml", &:read)
        require 'rexml/document'
        doc = REXML::Document.new(File.open("#{@igb_dir}/#{@genome_name}/annots.xml"))
        expect(doc.root.name).to eq("files")
      end
      
      it "includes the bed and the BAM files in the annots.xml file" do
        annots_txt = File.open("#{@igb_dir}/#{@genome_name}/annots.xml", &:read)
        expect(annots_txt).to match(/#{@genome_name}.bed/)
        expect(annots_txt).to match(/alignment.sorted.bam/)
      end
      
      it "copies the index.html from the QC website to the IGB genome's directory" do
        expect(md5("#{@out}/data/www/index.html")).to eq(md5("#{@igb_dir}/#{@genome_name}/index.html"))
      end
      
      it "copies the alignment.sorted.bam file into the IGB genome's directory" do
        expect(md5("#{@out}/data/qc_wd/alignment.sorted.bam")).to eq(
            md5("#{@igb_dir}/#{@genome_name}/alignment.sorted.bam"))
      end
      
    end
    
  end
  
end
