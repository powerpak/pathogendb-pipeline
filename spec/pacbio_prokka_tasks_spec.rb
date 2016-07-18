require_relative './spec_helper'
require 'pp'

describe "pathogendb-pipeline" do
  
  context "when running prokka pacbio-only tasks" do
    
    before(:all) { setup_temporary_OUT }
    after(:all) { cleanup_temporary_OUT }
  
    describe "task :pull_down_raw_reads" do
      
      it "pulls down a polished assembly for a job in SMRT Portal" do
        %x{SMRT_JOB_ID=023625 rake --quiet pull_down_raw_reads}
        expect(File).to exist("#{$OUT}/bash5.fofn")
        expect(File).to exist("#{$OUT}/data/polished_assembly.fasta.gz")
      end
      
      it "pulls down a polished assembly for a job in old_smrtportal_jobs" do
        %x{SMRT_JOB_ID=019203 rake --quiet pull_down_raw_reads}
        expect(File).to exist("#{$OUT}/bash5.fofn")
        expect(File).to exist("#{$OUT}/data/polished_assembly.fasta.gz")
      end
      
    end
    
    describe "task :assemble_raw_reads" do
      
      it "skips assembly for a job already assembled by SMRT Portal" do
        stdout = %x{SMRT_JOB_ID=023625 rake --quiet assemble_raw_reads}
        expect(File).to exist("#{$OUT}/data/polished_assembly.fasta.gz")
        expect(File).to_not exist("#{$OUT}/example_params.xml")
        expect(stdout).to match(/NOTICE: .* skipping assemble_raw_reads/)
      end
      
    end
    
    describe "task :run_circlator" do
      
      # Warning, long-running test.
      it "produces a FASTA file as output" do
        %x{SMRT_JOB_ID=023625 STRAIN_NAME=CD01394 rake --quiet run_circlator}
        expect(File).to exist("#{OUT}/data/#{STRAIN_NAME}_circlator/06.fixstart.fasta")
      end
    
    end
    
  end
  
end