require_relative './spec_helper'

describe "pipeline" do
  
  context "running prokka pacbio-only tasks" do
    
    before(:all) { setup_temporary_OUT }
    after(:all) { cleanup_temporary_OUT }
  
    describe "pull_down_raw_reads" do
      
      it "pulls down reads for job 019194" do
        system "echo $OUT"
        expect("hi").to eq("nope")
        #system "SMRT_JOB_ID=019194 rake pull_down_raw_reads"
      end
      
    end
    
  end
  
end