require 'pp'
require_relative 'lib/colors'
include Colors

SOURCE_BASES = Rake::FileList.new("inputs/*_assembly.fa").map { |f| f.sub(/_assembly\.fa$/, '') }

def bwa_align_outputs
  Rake::FileList.new do |fl|
    SOURCE_BASES.pathmap("%{^inputs/,outputs/}p").each do |f|
      fl.add("#{f}_assembly_close.bam")
      fl.add("#{f}_assembly_close.bai")
      fl.add("#{f}_assembly_close.wig")
      fl.add("#{f}_assembly_close.bigwig")
    end
  end
end

def mpileup_outputs
  Rake::FileList.new do |fl|
    SOURCE_BASES.pathmap("%{^inputs/,outputs/}p").each do |f|
      fl.add("#{f}_assembly_close_ilm.vcf")
      fl.add("#{f}_assembly_close_ilm.fasta")
    end
  end
end

task :default => :mpileup

desc "Performs BWA Align on QCed Illumina reads + a closed PacBio assembly"
task :bwa_align => [:check] + bwa_align_outputs

desc "Performs mpileup on BWA Align output to produces VCF and FASTA"
task :mpileup => [:check] + mpileup_outputs

desc "Create some test files just to get started"
task :test do
  mkdir_p "inputs"
  touch "inputs/123_123_assembly.fa"
  touch "inputs/123_123_illumina.fastq"
end

directory "outputs"

desc "Checks that everything is good to go"
task :check do |t|
  unless Dir.exists? "inputs"
    fail red <<-FAIL.gsub(/^\s+/, '')
      You need to symlink or create inputs/ as a directory of Illumina FASTQ and PacBio FA inputs.
      As an example, run `rake test`, which will create this directory and some sample test files.
    FAIL
  end
end

# QC Filter
rule(/_QC\.fastq$/ => proc { |f| f.sub(/^outputs\//, 'inputs/').sub(/_QC\.fastq$/, '.fastq') }) do |t|
  mkdir_p t.name.pathmap("%d")
  cp t.source, t.name
end

# BWA Align
BWA_ALIGN_REQS = proc {|f| [f.ext("fa"), f.sub(/_assembly_close\.\w+$/, '_illumina_QC.fastq')] }
rule(/_assembly_close\.(bam)$/ => BWA_ALIGN_REQS) do |t|
  mkdir_p t.name.pathmap("%d")
  puts "Here, I'd make #{t.name} from: \n  - #{t.sources.join("\n  - ")}"
  touch t.name
end
rule(/_assembly_close\.(bai|wig|bigwig)/ => proc { |f| f.ext('.bam') }) do |t|
  mkdir_p t.name.pathmap("%d")
  puts "Here, I'd make #{t.name} from: \n  - #{t.sources.join("\n  - ")}"
  touch t.name
end

# mpileup recall consensus
MPILEUP_REQS = proc {|f| base = f.sub(/_ilm\.(vcf|fasta)$/, ''); ["#{base}.bam", "#{base}.bai"] }
rule(/_assembly_close_ilm\.(vcf|fasta)$/ => MPILEUP_REQS) do |t|
  mkdir_p t.name.pathmap("%d")
  puts "Here, I'd make #{t.name} from: \n  - #{t.sources.join("\n  - ")}"
  touch t.name
end

# Circularize + finish
CIRCULARIZE_REQS = proc {|f| f.sub(/^outputs\//, 'inputs/').sub(/_close\.fa$/, '.fa') }
rule(/_assembly_close\.fa$/ => CIRCULARIZE_REQS) do |t|
  mkdir_p t.name.pathmap("%d")
  cp t.source, t.name
end