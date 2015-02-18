#!/usr/bin/env ruby

require 'bio'
require 'optparse'
require 'fileutils'
require 'tempfile'
require 'pp'

class ImproveRastOutput

  def initialize(args)
    @options = {}

    opt_parser = OptionParser.new do |opts|
      opts.banner = <<-USAGE

Usage: improve_rast_gbk.rb ref1.gbk [ref2.gbk, ...] [options] rast_output.gbk

Improves default RAST genbank output by adding to annotated features as follows:
  1)  Other reference GenBank files are searched for homologous genes, and these
      names are added to RAST features as "/gene" qualifiers, with corresponding
      "/inference" qualifiers for the match evidence from BLAST.
  2)  Antibiotic resistance gene databases are searched for homologous genes,
      and phenotypic information is added as "/phenotype" qualifiers.  See the
      Specific Options below for how to include these databases for this purpose.
      USAGE

      opts.separator ""
      opts.separator "Specific options:"
      
      # See ftp://ftp.cbcb.umd.edu/pub/data/ARDB/
      # -> ftp://ftp.cbcb.umd.edu/pub/data/ARDB/ardbAnno1.0.tar.gz
      opts.on("-d", "--ardb [PATH]",
              "Path to the ardbAnno1.0 directory (from ardbAnno1.0.tar.gz)") do |path|
        @options[:ardb] = path
      end
      
      # See http://en.mediterranee-infection.com/article.php?laref=283&titre=arg-annot-
      # -> http://www.mediterranee-infection.com/arkotheque/client/ihumed/_depot_arko/articles/304/arg-annot-database_doc.zip
      opts.on("-a", "--arg-annot [PATH]",
              "Path to the arg-annot-database_doc dirctory (from arg-annot-database_doc.zip)") do |path|
        @options[:arg_annot] = path
      end
      
      # See http://arpcard.mcmaster.ca/?q=node/6626
      # -> http://arpcard.mcmaster.ca/blast/db/protein/AR-polypeptides.fa.gz
      # -> http://arpcard.mcmaster.ca/blast/db/nucleotide/AR-genes.fa.gz
      # -> http://arpcard.mcmaster.ca/obo-download/aro.obo.gz
      opts.on("-c", "--card [PATH]",
              "Path to the CARD directory containing AR-genes.fa and aro.obo") do |path|
        @options[:card] = path
      end
      
      # See https://cge.cbs.dtu.dk/services/data.php
      # -> POST to https://cge.cbs.dtu.dk/cge/download_data.php with ?folder=resfinder&filename=resfinder.zip
      opts.on("-r", "--res-finder [PATH]",
              "Path to the resfinder directory (from resfinder.zip)") do |path|
        @options[:res_finder] = path
      end

      opts.separator ""
      opts.separator "Common options:"

      # No argument, shows at tail.  This will print an options summary.
      # Try it and see!
      opts.on_tail("-h", "--help", "Show this message") do
        puts opts
        exit
      end
    end

    opt_parser.parse!(args)
    
    if args.size < 2
      puts opt_parser.help
      exit
    end
    
    open_files(args)
  end  # parse()
  
  def run!
    begin
      make_tmpdir
      make_blast_db
      print_features
    ensure
      remove_tmpdir
    end
  end

  protected
  
  def open_files(args)
    abort "FATAL: #{args.last} does not exist" unless File.exist? args.last
    @rast_output = Bio::FlatFile.auto(args.pop)
    abort "FATAL: #{@rast_output.path} is not a valid GenBank file" unless @rast_output.dbclass == Bio::GenBank
    
    @refs = args.map do |arg|
      abort "FATAL: #{arg} does not exist" unless File.exist? arg
      ref_file = Bio::FlatFile.auto(arg)
      abort "FATAL: #{arg} is not a valid GenBank file" unless ref_file.dbclass == Bio::GenBank
      ref_file
    end
  end
  
  def make_tmpdir
    @tmpdir = Dir.mktmpdir
  end
  
  def make_blast_db
    @refs.each do |ref|
      ref.each do |record|
        record.features.each do |feature|
          gene_qualifier = feature.qualifiers.find{|q| q.qualifier == 'gene' }
          translation_qualifier = feature.qualifiers.find{|q| q.qualifier == 'translation' }
          protein_id_qualifier = feature.qualifiers.find{|q| q.qualifier == 'protein_id' }
          if feature.feature == "CDS" && gene_qualifier
            protein_id = protein_id_qualifier && protein_id_qualifier.value
            if translation_qualifier
              puts ">#{gene_qualifier.value} #{File.basename ref.path} #{ protein_id }"
              puts translation_qualifier.value
            else
              STDERR.puts "WARN: No translation (AA sequence) given for #{gene_qualifier.value} in #{File.basename ref.path}"
            end
          end
        end
      end
    end
  end
  
  def print_features 
    @rast_output.each do |record|
      record.features.each do |feature|
        #pp feature
      end
    end
  end
  
  def remove_tmpdir
    FileUtils.remove_entry_secure(@tmpdir, true)
  end

end  # class ImproveRastOutput


# =====================================
# = If this script is called directly =
# =====================================

if __FILE__ == $0
  iro = ImproveRastOutput.new(ARGV)
  
  iro.run!
  
end