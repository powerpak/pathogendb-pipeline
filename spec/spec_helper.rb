require 'tempfile'
require 'shellwords'
require 'fileutils'

require 'bundler/setup'

# Runs a command, with logging and visible output if DEBUG=1.
def run(cmd)
  if ENV['DEBUG']
    log_to = File.expand_path("../../rake_spec.log", __FILE__)
    cmd += " 2>&1 | tee -a #{Shellwords.escape log_to}"
    system cmd
  else
    `#{cmd}`
  end
end

# Makes a temporary directory and adds it to the environment before yielding to a block.
# If the block returns false or nil, the directory isn't deleted afterward.
def with_tempdir(name)
  ENV[name] = tmpdir = Dir.mktmpdir("rspec-#{name}")
  if yield tmpdir
    FileUtils.rm_r tmpdir, :force => true, :secure => true
  else
    $stderr.puts "DEBUG: #{name} directory left at #{tmpdir}"
  end
end

# Calculates the MD5 of a file's contents
def md5(filename)
  `md5sum #{Shellwords.escape filename}`.split(/ /).first
end

# Copies example data to the OUT directory.
def copy_example(name)
  example_dir = File.expand_path("../data/#{name}", __FILE__)
  # trailing dot copies directory contents instead of directory itself
  FileUtils.cp_r "#{example_dir}/.", ENV['OUT']
end

# Figure out the dependencies for a particular task.
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

# Touch task dependencies that look like filenames in order, so that they are not regenerated.
def touch_prereqs(task, env_string, dependencies=nil, already_touched=nil)
  task = task.to_s
  dependencies = list_dependencies(task, env_string) unless dependencies
  already_touched = {} unless already_touched
  
  dependencies[task].each do |subtask|
    if dependencies.has_key?(subtask)
      touch_prereqs(subtask, env_string, dependencies, already_touched)
    end
    full_path = File.expand_path(subtask, ENV['OUT'])
    # All file tasks must contain one of these characters.
    if subtask =~ /\.|\// and !already_touched[full_path]
      $stderr.puts "DEBUG: touch #{subtask}" if ENV['DEBUG']
      FileUtils.touch(full_path)
      already_touched[full_path] = true
      sleep 1
    end
  end
end

RSpec.configure do |config|
  config.expect_with :rspec do |c|
    c.syntax = :expect
  end
end
