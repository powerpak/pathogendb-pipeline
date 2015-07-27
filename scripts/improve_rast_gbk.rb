#!/usr/bin/env ruby

require 'bio'
require 'optparse'
require 'fileutils'
require 'tempfile'
require 'pp'
require 'rexml/document'

class ImproveRastOutput

  REQUIRED_BINS = {"blastp"=>"BLAST+", "makeblastdb"=>"BLAST+"}

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
  end
  
  def run!
    begin
      make_tmpdir
      check_executables
      make_blast_db
      match_genes
      print_output
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
  
  def check_executables
    if missing = REQUIRED_BINS.keys.find{|b| `which #{b}`.strip.size == 0 }
      fail "FAIL: Could not find \`#{missing}\` in your $PATH; please ensure #{REQUIRED_BINS[missing]} is installed."
    end
  end
  
  def make_blast_db
    File.open("#{@tmpdir}/refgenes.fastp", "w") do |f|
      @refs.each do |ref|
        ref.each do |record|
          
          record.features.each do |feature|
            gene_qualifier = feature.qualifiers.find{|q| q.qualifier == 'gene' }
            
            next unless feature.feature == "CDS" && gene_qualifier
            
            translation_qualifier = feature.qualifiers.find{|q| q.qualifier == 'translation' }
            protein_id_qualifier = feature.qualifiers.find{|q| q.qualifier == 'protein_id' }
            
            unless translation_qualifier
              STDERR.puts "WARN: No translation (AA sequence) given for #{gene_qualifier.value} in #{File.basename ref.path}"
              # TODO: should be able to get this anyway by iterating on feature.locations, grabbing nt's, translating manually
              next
            end
            
            protein_id = protein_id_qualifier && protein_id_qualifier.value
            
            f.write ">#{gene_qualifier.value} #{ protein_id } #{File.basename ref.path}\n"
            f.write translation_qualifier.value
            f.write "\n"
          end
          
        end # ref.each
      end # @refs.each
    end
    Dir.chdir(@tmpdir) do
      `makeblastdb -in refgenes.fastp -dbtype prot`
    end
  end
  
  # The threshold for a blastp <Hit> becoming the annotated gene name for a CDS.
  # Currently we base this on E-value and <Hit_score> (normalized to product length) thresholds.
  def meets_threshold(hit, query)
    hit && hit[:evalue] < 1e-10 && hit[:score].to_f / query.size > 3.5
  end
  
  def match_genes
    counts = {:all => 0, :searched => 0, :matched => 0}
    @rast_output_records = []
    
    Dir.chdir(@tmpdir) do
      
      @rast_output.each do |record|
        record.features.each do |feature|
          counts[:all] += 1
          db_xref_qualifier = feature.qualifiers.find{|q| q.qualifier == 'db_xref' }
          translation_qualifier = feature.qualifiers.find{|q| q.qualifier == 'translation' }
          
          next unless feature.feature == "CDS" # TODO: should re-label tRNA and rRNA better too
          unless translation_qualifier
            STDERR.puts "WARN: No translation (AA sequence) given for #{db_xref_qualifier.value} in #{File.basename @rast_output.path}"
            # TODO: could use the nt's gathered from feature.locations, grabbing exons, and then blastx -query_gencode 11
            next
          end
          
          counts[:searched] += 1
          aa_query = translation_qualifier.value
          File.open("query.fastp", 'w') {|f| f.write(aa_query) }
          
          xml = `blastp -db refgenes.fastp -query query.fastp -outfmt 5`
          doc = REXML::Document.new(xml)
          
          hits = []
          doc.elements.each('BlastOutput/BlastOutput_iterations/Iteration/Iteration_hits/Hit') do |hit_el|
            hsps = []
            
            hit_el.elements.each('Hit_hsps/Hsp') do |hsp_el| 
              hsps << {
                :score => hsp_el.elements["Hsp_score"].text.to_i,
                :evalue=> hsp_el.elements["Hsp_evalue"].text.to_f,
                :midline => hsp_el.elements["Hsp_midline"].text
              }
            end
            top_hsp = hsps.max_by{|a| a[:score] }
            
            hit_def = hit_el.elements['Hit_def'].text.split(/ /, 3)
            hits << {:gene => hit_def[0], :protein_id => hit_def[1], :file => hit_def[2]}.merge(top_hsp)
          end
          top_hit = hits.min_by{|a| a[:evalue] }
          
          if top_hit && meets_threshold(top_hit, aa_query)
            feature.qualifiers << Bio::Feature::Qualifier.new("gene", top_hit[:gene])
            if top_hit[:protein_id]
              protein_id_type = top_hit[:protein_id][2] == '_' ? "RefSeq" : "GenBank"
              inference = "DESCRIPTION:similar to AA sequence (same species):blastp:#{protein_id_type}:#{top_hit[:protein_id]}"
            else
              inference = "DESCRIPTION:similar to AA sequence (same species):blastp:file:#{top_hit[:file]}"
            end
            feature.qualifiers << Bio::Feature::Qualifier.new("inference", inference)
            
            counts[:matched] += 1
          end
        end # record.features.each
        
        @rast_output_records << record
      end # @rast_output.each
      
    end # Dir.chdir
    
    STDERR.puts "INFO: #{counts[:all]} features processed, #{counts[:searched]} searched, #{counts[:matched]} matched"
  end
  
  def print_output
    @rast_output_records.each do |record|
      STDOUT.write record.to_biosequence.output(:genbank)
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