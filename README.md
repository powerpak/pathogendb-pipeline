# Prototype PathogenDB pipeline

## Requirements

As of now, this only runs on [Minerva](http://hpc.mssm.edu) because it uses modules and software found on that cluster.  In due time, it might be made portable to other systems...

## Usage

First, clone this repository to a directory and `cd` into it.  You'll want to configure and setup your env first using an included script:

    $ cp scripts/env.example.sh scripts/env.sh
    $ $EDITOR scripts/env.sh
    $ source scripts/env.sh

At the moment, there's not much to edit there in the second step, as the defaults should work for any Minerva user.  Then:

    $ bundle install --deployment

to install the required gems locally into the `vendor/bundle` directory.  Then run `rake -T` to see the available tasks, and you can run them with `rake $TASK_NAME`.

## Environment variables

Certain tasks within the pipeline require you to specify some extra information as an environment variable.  You can do this by either editing them into `scripts/env.sh` and re-running `source scripts/env.sh`, or you can prepend them to the `rake` invocation, e.g.:

    $ SMRT_JOB_ID=019194 rake pull_down_raw_reads

If a required environment variable isn't present when a task is run, rake will abort with an error message.

Variable      | Required by           | Default | Purpose
--------------|-----------------------|---------|-----------------------------------
`OUT`         | all tasks             | ./out   | This is where your interim files are saved.
`SMRT_JOB_ID` | `pull_down_raw_reads` | (none)  | The ID of the job on the SMRT Portal with your reads.
`STRAIN_NAME` | `resequence_assembly` | (none)  | The strain name for your sample.
