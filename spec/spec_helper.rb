require 'tempfile'
require 'shellwords'
require 'fileutils'

require 'bundler/setup'

def setup_temporary_OUT
  ENV['OUT'] = $OUT = Dir.mktmpdir('rspec')
end

def cleanup_temporary_OUT
  FileUtils.rm_r $OUT, :force => true, :secure => true
end