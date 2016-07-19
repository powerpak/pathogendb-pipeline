require_relative './spec_helper'
require 'pp'

describe "pathogendb-pipeline" do
  
  around(:each) do |example|
    setup_temporary_OUT
    example.run
    if ENV['DEBUG'] && example.exception
      $stderr.puts "DEBUG: Files left in #{$OUT}"
    else
      cleanup_temporary_OUT(example)
    end
  end
  
  context "when running tasks from scratch" do
  
    describe "task :pull_down_raw_reads" do
      it "pulls down a polished assembly for a job in SMRT Portal (023625)" do
        run "SMRT_JOB_ID=023625 rake --silent pull_down_raw_reads"
        expect(File).to exist("#{$OUT}/bash5.fofn")
        expect(File).to exist("#{$OUT}/data/polished_assembly.fasta.gz")
      end
      
      it "pulls down a polished assembly for a job in old_smrtportal_jobs (019203)" do
        run "SMRT_JOB_ID=019203 rake --silent pull_down_raw_reads"
        expect(File).to exist("#{$OUT}/bash5.fofn")
        expect(File).to exist("#{$OUT}/data/polished_assembly.fasta.gz")
      end
    end
    
    describe "task :assemble_raw_reads" do
      it "skips assembly for a job already assembled by SMRT Portal (023625)" do
        stdout = run "SMRT_JOB_ID=023625 rake --silent assemble_raw_reads"
        expect(File).to exist("#{$OUT}/data/polished_assembly.fasta.gz")
        expect(File).not_to exist("#{$OUT}/example_params.xml")
        expect(stdout).to match(/NOTICE: .* skipping assemble_raw_reads/)
      end
    end
    
    describe "task :run_circlator", :speed => 'slow' do
      it "produces a FASTA file with md5 64ae02c8... as output (for 019203)" do
        strain = 'SM5478'
        run "SMRT_JOB_ID=019203 STRAIN_NAME=#{strain} rake --silent run_circlator"
        fasta_out = "#{$OUT}/data/#{strain}_circlator/06.fixstart.fasta"
        expect(File).to exist(fasta_out)
        expect(md5(fasta_out)).to eq('64ae02c8a8e68dc2134899b4281f93ca')
      end
    end
    
  end
  
  context "when running tasks on SM5478 (019203), after run_circlator" do
    
    before(:each) do
      @strain = "SM5478"
      copy_example "#{@strain}-run_circlator"
      touch_prereqs :run_circlator, "SMRT_JOB_ID=019203 STRAIN_NAME=#{@strain}"
    end
    
    describe "task :post_circlator", :speed => 'medium' do
      it "renames the first contig to u00000crxx_c_019203" do
        run "SMRT_JOB_ID=019203 STRAIN_NAME=#{@strain} rake --silent post_circlator"
        fasta_out = "#{$OUT}/data/#{@strain}_postcirc.fasta"
        expect(File).to exist(fasta_out)
        first_line = File.open(fasta_out, &:readline).strip
        expect(first_line).to eq('>u00000crxx_c_019203')
      end
      
      it "renames the first contig to u00000crxm_c_019203 if circlator output is deleted" do
        FileUtils.rm "#{$OUT}/data/#{@strain}_circlator/00.input_assembly.fasta"
        run "SMRT_JOB_ID=019203 STRAIN_NAME=#{@strain} rake --silent post_circlator"
        fasta_out = "#{$OUT}/data/#{@strain}_postcirc.fasta"
        expect(File).to exist(fasta_out)
        first_line = File.open(fasta_out, &:readline).strip
        expect(first_line).to eq('>u00000crxm_c_019203')
      end
    end
    
    describe "task :resequence_assembly", :speed => 'slow' do
      it "produces a FASTA file with md5 ... as output" do
        run "SMRT_JOB_ID=019203 STRAIN_NAME=#{@strain} CLUSTER=BASH rake --silent resequence_assembly"
        fasta_out = "#{$OUT}/data/#{@strain}_consensus_circ.fasta"
        expect(File).to exist(fasta_out)
        puts md5(fasta_out)
        expect(md5(fasta_out)).to eq('64ae02c8a8e68dc2134899b4281f93ca')
      end
    end
    
  end
  
end
