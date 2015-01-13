# Utility functions for running commands within GNU screen, particularly its split
# screen functionality
require 'shellwords'

module Subscreens
  extend self
  
  def create_screens(num)
    # If we aren't already in a screen session, re-run this command inside a screen session
    if ENV['STY'].nil?
      args = ARGV.map{|arg| Shellwords.escape(arg) }.join(' ')
      system "screen -c /dev/null -t run1 -S nested bash -c '#{Shellwords.escape($0)} #{args}'"
      return false
    end

    (num - 1).times do
      system "screen -X split"
    end
    true
  end
  
  # Run one cmd inside num_splits split screen sessions simultaneously
  def split(num_splits, cmd)
    return unless create_screens(num_splits)

    num_splits.times do |i|
      system "screen -X focus" if i != 0
      system "screen -X screen -t run#{i+1} #{i+1}"
      sleep 10
      system "screen -X -p #{i+1} stuff #{Shellwords.escape(cmd + "\n")}"
      sleep 10
      system "screen -X select #{i+1}" if i != 0
    end
  end
  
  # Run cmds inside their own screens simultaneously
  def run(cmds)    
    return unless create_screens(cmds.size)
    
    cmds.each_with_index do |cmd, i|
      system "screen -X focus" if i != 0
      system "screen -X screen -t run#{i+1} #{i+1}"
      sleep 10
      system "screen -X -p #{i+1} stuff #{Shellwords.escape(cmd + "\n")}"
      sleep 10
      system "screen -X select #{i+1}" if i != 0
    end
  end
end