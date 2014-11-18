# Prototype PathogenDB pipeline

## Requirements

As of now, this only runs on [Minerva](http://hpc.mssm.edu) because it uses modules and software found on that cluster.  In due time, it might be made portable to other systems.

## Usage

First, clone this repository to a directory and `cd` into it.  You'll want to configure your environment first using the included script:

    $ cp scripts/env.example.sh scripts/env.sh
    $ $EDITOR scripts/env.sh    

At a minimum, you will need to configure the following:

- `RAST_USER`: username for your [RAST][rast] account
- `RAST_PASSWORD`: password for your [RAST][rast] account

[rast]: http://rast.nmpdr.org/

For the rest of the variables, the defaults should work for any Minerva user.  Then, you can source the script into your shell and install required gems locally into the `vendor/bundle` directory as follows:

    $ source scripts/env.sh
    $ bundle install --deployment

When this is complete, you should be able to run rake to kick off the pipeline as follows. However, also read **[Environment variables](#environment-variables)** below, as certain tasks require more variables to be set before being invoked.

    $ rake -T                    # list the available tasks
    $ rake $TASK_NAME            # run the task named $TASK_NAME
    $ FOO="bar" rake $TASK_NAME  # run $TASK_NAME with FOO set to "bar"

When firing up the pipeline in a new shell, always remember to `source scripts/env.sh` before running `rake`.

### Environment variables

Certain tasks within the pipeline require you to specify some extra information as an environment variable.  You can do this by either editing them into `scripts/env.sh` and re-running `source scripts/env.sh`, or you can prepend them to the `rake` invocation, e.g.:

    $ SMRT_JOB_ID=019194 rake pull_down_raw_reads

If a required environment variable isn't present when a task is run and there is no default value, rake will abort with an error message.

Variable      | Required by                           | Default | Purpose
--------------|---------------------------------------|---------|-----------------------------------
`OUT`         | all tasks                             | ./out   | This is where your interim files are saved.
`SMRT_JOB_ID` | `pull_down_raw_reads`                 | (none)  | The ID of the job on the SMRT Portal with your reads.
`STRAIN_NAME` | `resequence_assembly` `rast_annotate` | (none)  | The strain name for your sample.
`SPECIES`     | `rast_annotate`                       | (none)  | The species for your sample.

### Dependency graph

This Rakefile is able to build a dependency graph of its intermediate files from itself.  Use the `rake graph` task for this; it will be generated at `$OUT/pathogendb-pipeline.png`.

![Dependency graph](https://pakt01.u.hpc.mssm.edu/pathogendb-pipeline.png)

## Other notes

This pipeline downloads and installs the [Network-based SEED API package](http://blog.theseed.org/servers/installation/distribution-of-the-seed-server-packages.html) into `vendor/sas`.  Documentation for some of the included executables and the Perl API are also on that page.
