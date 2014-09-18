# Prototype PathogenDB pipeline

## Requirements

Linux and Mac computers are supported.  You'll need Ruby ≥1.9.2, RubyGems, and some basic gems, which will be installed with Bundler.

On most Macs, you'll only need to run `sudo gem install bundler` because you already have Ruby 2.0.0 and RubyGems, as of Mavericks (10.9).

On Linux, I defer either to your package manager or your raging desire to build from source.

In both cases, [rbenv](https://github.com/sstephenson/rbenv) might help if your system Ruby is out of date.

## Usage

First, clone this repository to a directory and `cd` into it.  Then:

    $ bundle install

to install the required gems.  Bundler may prompt you for your password if it needs to write to a system directory.

Then, since we're just testing out bare functionality right now, use

    $ rake test

to create some test files in the `inputs` directory, and then

    $ rake

will run a fake pipeline on those inputs to produce a host of outputs.