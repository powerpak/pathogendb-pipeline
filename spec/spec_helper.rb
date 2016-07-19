require 'tempfile'
require 'shellwords'
require 'fileutils'

require 'bundler/setup'

def run(cmd)
  stdout = `#{cmd}`
end

def setup_temporary_OUT
  ENV['OUT'] = $OUT = Dir.mktmpdir('rspec')
end

def cleanup_temporary_OUT(example)
  FileUtils.rm_r $OUT, :force => true, :secure => true
end

def copy_example(name)
  example_dir = File.expand_path("../data/#{name}", __FILE__)
  # trailing dot copies directory contents instead of directory itself
  FileUtils.cp_r "#{example_dir}/.", $OUT
end

def list_dependencies(task, env_string)
  stdout = `#{env_string} rake --silent -P`
  
  dependencies = Hash.new{|h, k| h[k] = [] }
  task = nil
  stdout.split("\n").each do |line|
    line.chomp!
    if line =~ /^rake (.*)/
      task = $1
    else
      subtask = line.strip
      next unless task
      dependencies[task] << subtask
    end
  end
  dependencies["check"] = []  # waste of time to touch these
  dependencies
end

def touch_prereqs(task, env_string, dependencies=nil)
  task = task.to_s
  dependencies = list_dependencies(task, env_string) unless dependencies
  
  dependencies[task].each do |subtask|
    if dependencies.has_key?(subtask)
      touch_prereqs(subtask, env_string, dependencies)
    end
    if subtask =~ /\.|\//  # All file tasks must contain one of these characters.
      $stderr.puts "DEBUG: touch #{subtask}" if ENV['DEBUG']
      FileUtils.touch(File.expand_path(subtask, $OUT))
      sleep 1
    end
  end
end

RSpec.configure do |config|
  config.expect_with :rspec do |c|
    c.syntax = :expect
  end
end