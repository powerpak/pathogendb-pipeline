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
  
  opts.on("-n", "--narrow-path TASK_NAME,PARENT_TASK", Array, 
      "Removes paths through TASK_NAME except the path coming from PARENT_TASK.") do |list|
    options[:narrow_path] = list
  end
  
  opts.on("-r", "--replace-with STRING", "Replaces the pruned filepath with STRING") do |r|
    options[:replace] = r
  end
  
  opts.on_tail("-h", "--help", "Show this message") do
    puts opts
    exit
  end
end.parse!

def clean_task_name(task, options, parent_task=nil)
  prune = options[:prune]
  narrow_path = options[:narrow_path]
  
  if task and prune
    if task.size <= prune.size and prune.index(task) == 0
      task = nil # Prune (ignore) this task, since it is an ancestor of the --prune argument
    elsif narrow_path and parent_task and task == narrow_path[0] and parent_task != narrow_path[1]
      task = nil # Prune (ignore) this task, since it is 
    elsif options[:replace] and task.index(prune) == 0
      task.sub!(options[:prune], options[:replace]) 
    end
  end
  task
end

puts "digraph g {"
puts 'size="30,12"'
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
    subtask = clean_task_name(line.strip, options, task)
    next unless subtask
    puts "\"#{task}\" -> \"#{subtask}\";"
  end
end

puts "}"