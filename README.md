# PathogenDB assembly and annotation

## Requirements

As of now, this only runs on [Minerva](http://hpc.mssm.edu) because it uses modules and software found on that cluster.  In due time, it might be made portable to other systems.

Currently, you also need to be in the `pacbioUsers` group on Minerva and have access to the `premium` LSF queue and the `acc_PBG` LSF account.

To avoid hitting resource limits on the login nodes on Minerva, we recommend that you run the pipeline on the interactive1 or interactive2 nodes, which have no such limits.  To do so run `ssh interactive1` or `ssh interactive2` after logging into Minerva normally.

## Usage

First, clone this repository to a directory and `cd` into it.  You'll want to configure your environment first using the included script:

    $ cp scripts/example.env.sh scripts/env.sh
    $ $EDITOR scripts/env.sh    

At a minimum, you will need to configure the following:

- `RAST_USER`: username for your [RAST][rast] account
- `RAST_PASSWORD`: password for your [RAST][rast] account

[rast]: http://rast.nmpdr.org/

For the rest of the variables, the defaults should work for any Minerva user.  Then, you can source the script into your shell and install required gems locally into the `vendor/bundle` directory as follows:

    $ source scripts/env.sh
    $ bundle install --deployment

When this is complete, you should be able to run rake to kick off the pipeline as follows. However, first read **[Environment variables](#environment-variables)** below, as certain tasks require more variables to be set before being invoked.  A description of the typical sequence for assembling a genome is described below in **[Tasks](#tasks)**.

    $ rake -T                    # list the available tasks
    $ rake $TASK_NAME            # run the task named $TASK_NAME
    $ FOO="bar" rake $TASK_NAME  # run $TASK_NAME with FOO set to "bar"

When firing up the pipeline in a new shell, always remember to `source scripts/env.sh` before running `rake`.

### Environment variables

Certain tasks within the pipeline require you to specify some extra information as an environment variable.  You can do this by either editing them into `scripts/env.sh` and re-running `source scripts/env.sh`, or you can prepend them to the `rake` invocation, e.g.:

    $ SMRT_JOB_ID=019194 rake pull_down_raw_reads

If a required environment variable isn't present when a task is run and there is no default value, rake will abort with an error message.

Variable             | Required by                           | Default | Purpose
---------------------|---------------------------------------|---------|-----------------------------------
`OUT`                | all tasks                             | ./out   | This is where your interim and completed files are saved
`SMRT_JOB_ID`        | `pull_down_raw_reads`                 | (none)  | The ID of the job on the SMRT Portal with your reads.
`STRAIN_NAME`        | `resequence_assembly` `rast_annotate` `rast_annotate_ilm` `recall_ilm_consensus` `rast_to_igb` | (none)  | The strain name for your sample. **This cannot include anything but letters, numbers and underscores.**
`SPECIES`            | `rast_annotate` `rast_annotate_ilm` `rast_to_igb`         | (none)  | The species for your sample.
`ILLUMINA_FASTQ`     | `recall_ilm_consensus`                | (none)  | A path pointing to a FASTQ file containing the Illumina reads.

### Tasks

The typical series of tasks used to assemble a strain's genome from PacBio RS reads and then annotate with RAST are:

1. `pull_down_raw_reads`
2. `assemble_raw_reads`
3. `circularize_assembly`
4. `resequence_assembly`
5. `rast_annotate`
6. `rast_to_igb`

With some exceptions (for instance, if you need to manually edit interim files) you should be able to simply run `rake` with the last task you want to reach, and assuming you've specified all required [environment variables](#environment-variables), the pipeline will take care of running any necessary previous tasks, based on what's already present or missing from the `OUT` directory.

The final task, `rast_to_igb`, creates an [IGB](http://bioviz.org/igb/) Quickload-compatible directory so you can load the genome into IGB. By default, this occurs in `~/www/igb`, although you can override this by setting `IGB_DIR` in your `scripts/env.sh`. To view the genome in IGB, open IGB's preferences and add `https://YOUR_USERNAME.u.hpc.mssm.edu/igb/` as a Quickload data source (replacing `YOUR_USERNAME` with your Minerva username), and then you should be able to find your genome under the Species dropdown in the browser.

Optionally, if Illumina reads are also available the same isolate, they can be used to iron out small errors in the PacBio-produced assembly and then the new consensus can be re-annotated with these two extra steps:

6. `recall_ilm_consensus`
7. `rast_annotate_ilm`

### Multiple runs within `screen`

If you'd like to run the pipeline multiple times with different parameters placed into the environment variables, you might find it useful to try the special `rake multi[$TASK_FILE,$N]` task.

This takes one required parameter, `$TASK_FILE`, placed in the brackets. It should be a file that lists, one per line, the separate task names and environment variables you'd like to use.  Here's an example with two tasks:

    resequence_assembly OUT=$HOME/Steno/SM_278  SMRT_JOB_ID=017871 STRAIN_NAME=SM_278  SPECIES="Stenotrophomonas"
    rast_annotate       OUT=$HOME/Steno/SM_5478 SMRT_JOB_ID=019203 STRAIN_NAME=SM_5478 SPECIES="Stenotrophomonas"

When you run `rake multi[$TASK_FILE]`, with the filename of your task file in the brackets, a `screen` session will be created and split vertically into multiple windows, each of which will run `rake` with the various parameters you put on that line.

If you don't want as many splits, and would like some of the jobs to run in series, also specify a number `$N` that is smaller than the number of lines in your task file, like so:

    $ rake multi[/path/to/task/file,3]

You will almost certainly need to run `rake multi` on an interactive node or the process will hit resource limits.

### Dependency graph

This Rakefile is able to build a dependency graph of its intermediate files from itself.  Use the `rake graph` task for this; it will be generated at `$OUT/pathogendb-pipeline.png`.

![Dependency graph](https://pakt01.u.hpc.mssm.edu/pathogendb-pipeline.png)

## Other notes

This pipeline downloads and installs the [Network-based SEED API package](http://blog.theseed.org/servers/installation/distribution-of-the-seed-server-packages.html) into `vendor/sas`.  Documentation for some of the included executables and the Perl API are also on that page.  Within this package, we are overriding `RASTserver.pm` with our own version in `lib/perl`, because we need to support setting use of a proxy via the `HTTP_PROXY` environment variable.
