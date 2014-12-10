#!/usr/bin/env ruby
# Based on https://gist.github.com/mvidner/7911768

require 'optparse'

options = {}
OptionParser.new do |opts|
  opts.banner = <<-USAGE.gsub /^[ \t]+/, ""
    Usage: rake -P | rake-prereqs-dot.rb [options] | dotty -
           rake -P | rake-prereqs-dot.rb [options] | dot -Tpng -o graph.png
    
    This script reads the output of rake -P/--prereqs on STDIN and converts it
    to a graphviz graph.
    
    Options:
  USAGE

  opts.on("-p", "--prune PATH", 
      "Prune tasks corresponding to PATH and its ancestors (*not* its descendants)") do |p|
    options[:prune] = p
  end
  opts.on("-r", "--replace-with STRING", "Replaces the pruned filepath with STRING") do |r|
    options[:replace] = r
  end
end.parse!

def clean_task_name(task, options)
  if task and options[:prune]
    if task.size <= options[:prune].size and options[:prune].index(task) == 0
      task = nil # Prune (ignore) this task.
    elsif options[:replace] and task.index(options[:prune]) == 0
      task.sub!(options[:prune], options[:replace]) 
    end
  end
  task
end

puts "digraph g {"
puts 'size="20,20"'
puts 'ratio="compress"'

task = nil
STDIN.each_line do |line|
  line.chomp!
  if line =~ /^rake (.*)/
    task = $1
    file_task = $1.match(/\.|\//)
    task = clean_task_name(task, options)
    next unless task
    style = file_task ? "" : " [style=filled; fontname=\"Times bold\"]"
    puts "\"#{task}\"#{style};"
  else
    subtask = clean_task_name(line.strip, options)
    next unless subtask
    puts "\"#{task}\" -> \"#{subtask}\";"
  end
end

puts "}"