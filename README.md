# PathogenDB assembly and annotation

## Requirements

As of now, this only runs on [Minerva](http://hpc.mssm.edu) because it uses modules and software found on that cluster.  In due time, it might be made portable to other systems.

Currently, you also need to be in the `pacbioUsers` group on Minerva and have access to the `premium` LSF queue and the `acc_PBG` LSF account.

To avoid hitting resource limits on the login nodes on Minerva, we recommend that you run the pipeline on one of the interactive nodes, which have no such limits.  To do so run `ssh interactive1` or `ssh interactive2` ... up to `ssh interactive6` after logging into Minerva normally. You can also [submit pipeline runs to a bsub queue](#running-as-a-bsub-task), which will deduct resources against your LSF account.

## Usage

First, clone this repository to a directory and `cd` into it.  You'll want to configure your environment first using the included script:

    $ cp scripts/example.env.sh scripts/env.sh
    $ $EDITOR scripts/env.sh    

All the default variables should work for any Minerva user with appropriate permissions.  (You previously needed to configure your [RAST][rast] login details in this script, but this is no longer necessary now that we use [prokka][].)  Then, you can source the script into your shell and install required gems locally into the `vendor/bundle` directory as follows:

    $ source scripts/env.sh
    $ bundle install --deployment

When this is complete, you should be able to run `rake` to kick off the pipeline as follows. However, first read **[Environment variables](#required-environment-variables)** below, as certain tasks require more variables to be set before being invoked.  A description of the typical sequence for assembling a genome is described below in **[Tasks](#tasks)**.

    $ rake -T                    # list the available tasks
    $ rake $TASK_NAME            # run the task named $TASK_NAME
    $ FOO="bar" rake $TASK_NAME  # run $TASK_NAME with FOO set to "bar"

When firing up the pipeline in a new shell, **remember to always `source scripts/env.sh` before running `rake`.**

[rast]: http://rast.nmpdr.org/

### Required environment variables

Certain tasks within the pipeline require you to specify some extra information as an environment variable.  You can do this by either editing them into `scripts/env.sh` and re-running `source scripts/env.sh`, or you can prepend them to the `rake` invocation, e.g.:

    $ SMRT_JOB_ID=019194 rake pull_down_raw_reads

If a required environment variable isn't present when a task is run and there is no default value, rake will abort with an error message.

Variable             | Required by                                             | Default | Purpose
---------------------|---------------------------------------------------------|---------|-----------------------------------
`OUT`                | all tasks                                               | ./out   | This is where your interim and completed files are saved
`SMRT_JOB_ID`        | `pull_down_raw_reads` `prokka_to_igb` `post_circlator` `ilm:rast_to_igb`  | (none)  | The ID of the job on the SMRTPortal with your reads.
`STRAIN_NAME`        | `resequence_assembly` `prokka_annotate` `prokka_to_igb` `ilm:recall_consensus` `ilm:prokka_annotate` `ilm:prokka_to_igb`  | (none)  | The strain name for your sample. **This cannot include anything but letters, numbers and underscores.**
`SPECIES`            | `prokka_annotate` `prokka_to_igb` `ilm:prokka_annotate` `ilm:prokka_to_igb` | (none)  | The species for your sample.
`ILLUMINA_FASTQ`     | `ilm:recall_consensus`                                  | (none)  | A path pointing to a FASTQ file containing the Illumina unpaired reads.

### Optional environment variables

These variables may be provided to configure certain tasks within the pipeline, but are not required.

Variable             | Can be provided for                   | Default | Purpose
---------------------|---------------------------------------|---------|-----------------------------------
`LSF_DISABLED`       | all steps                             | (none)  | Set to a non-empty string to disable all submissions to LSF (everything only runs locally). Useful for debugging on an interactive node.
`REPLACE_FASTA`      | `pull_down_raw_reads`                 | (none)  | Shunts this FASTA file in place of the SMRTPortal-built `polished_assembly.fasta`. Use this to run the pipeline on a non-SMRTpipe (i.e., manually fixed) assembly.
`CURATED`            | `run_circlator` `post_circlator`      | (none)  | If set to a non-empty string, circlator will only reorient (not attempt to circularize) the assembly. Useful for manually fixed assemblies.
`CLUSTER`            | `assemble_raw_reads` `resequence_assembly` | LSF_PSP | Sets the `-D CLUSTER=` option for [`smrtpipe.py`][smrtpipe], which controls which job submission wrapper scripts are used. Set to `BASH` to disable job submissions by SMRT Pipe (all steps run local to the current node).
`ILLUMINA_REFERENCE` | `ilm:fake_prereqs`                    | (none)  | Path to the FASTA file containing the reference sequence that you want to shunt into the Illumina correction branch of the pipeline.
`REORIENT_FASTA`     | `old:reorient_assembly`               | (none)  | *This option is deprecated in favor of using [circlator][].* A path pointing to a FASTA file with a landmark that the assembly will be reoriented to.  If not given, reorientation will be skipped.
`REORIENT_FLANK`     | `old:reorient_assembly`               | 0       | *This option is deprecated in favor of using [circlator][].* This is the number of nt *before* the beginning of the landmark where the origin of the circular chromosome will be set.
`GENBANK_REFERENCES` | `old:improve_rast`                    | (none)  | *This option is deprecated in favor of using [prokka][].* Paths to to GenBank files that contain "good" gene names that will be lifted over to your RAST annotations.  Multiple paths should be separated with `:`, as with `PATH`.  If not given, `improve_rast` will be a no-op.

[smrtpipe]: http://www.pacb.com/wp-content/uploads/2015/09/SMRT-Pipe-Reference-Guide.pdf
[prokka]: http://www.vicbioinformatics.com/software.prokka.shtml

### Tasks

The typical series of tasks used to assemble a strain's genome from PacBio RS reads and then annotate with [prokka][] are:

1. `pull_down_raw_reads`
2. `assemble_raw_reads`
3. `run_circlator`
4. `post_circlator`
5. `resequence_assembly`
6. `post_quiver_orient_correct`
7. `prokka_annotate`
8. `create_QC_webpage`
9. `prokka_to_igb`

With some exceptions (for instance, if you need to manually edit interim files) you should be able to simply run `rake` with the last task you want to reach, and assuming you've specified all [required environment variables](#required-environment-variables), the pipeline will take care of running any necessary previous tasks, based on what's already present or missing from the `OUT` directory.

If you are trying to re-run the pipeline on a manually fixed and/or circularized assembly, you may find the `REPLACE_FASTA` and `SKIP_CIRCLATOR` [optional environment variables](#optional-environment-variables) useful.

The final task, `prokka_to_igb`, creates an [IGB](http://bioviz.org/igb/) Quickload-compatible directory so you can load the genome into IGB. By default, this occurs in `~/www/igb`, although you can override this by setting `IGB_DIR` in your `scripts/env.sh`. To view the genome in IGB, open IGB's preferences and add `https://YOUR_USERNAME.u.hpc.mssm.edu/igb/` as a Quickload data source (replacing `YOUR_USERNAME` with your Minerva username), and then you should be able to find your genome under the Species dropdown in the browser.

**TODO:** document `motif_and_mods` and `all`.

#### Illumina-specific tasks

Optionally, if unpaired Illumina reads are available for the same isolate (which you provide as the `ILLUMINA_FASTQ` parameter), they can be aligned to correct small (typically indel) errors in the PacBio assembly, and then this new consensus can be re-annotated and converted to a new IGB quickload directory with four extra steps.

1. `ilm:recall_consensus`
2. `ilm:prokka_annotate`
3. `ilm:create_QC_webpage`
4. `ilm:prokka_to_igb`

If you had discarded the intermediate files for the PacBio-only assembly, you can still shunt its final FASTA sequence into the Illumina branch of the pipeline by using the `ilm:fake_prereqs` task, which takes the `ILLUMINA_REFERENCE` parameter (the path to the FASTA file). This sets up a skeleton job directory so that you may run the four `ilm:` steps listed above on their own.

#### Deprecated tasks

Prior versions of the pipeline utilized custom scripts for circularization based on [MUMmer][] output, and [RAST][rast] for annotation of finished assemblies. These processes are now handled by the more mature tools [circlator][] and [prokka][], respectively, but we have preserved the old steps under the `old:` namespace within `deprecated.rake`.

[MUMmer]: (http://mummer.sourceforge.net/)
[circlator]: http://sanger-pathogens.github.io/circlator/

1.  `old:circularize_assembly`
2.  `old:resequence_assembly`
3.  `old:reorient_assembly`
4.  `old:rast_annotate`
5.  `old:improve_rast`
6.  `old:rast_to_igb`
7.  `old:ilm:recall_consensus`
8.  `old:ilm:rast_annotate`
9.  `old:ilm:improve_rast`
10. `old:ilm:rast_to_igb`

### Multiple runs within `screen`

If you'd like to run the pipeline multiple times with different parameters placed into the environment variables, you might find it useful to try the special `rake multi[$TASK_FILE,$N]` task.

This takes one required parameter, `$TASK_FILE`, placed in the brackets. It should be a file that lists, one per line, the separate task names and environment variables you'd like to use.  Here's an example with two tasks:

    resequence_assembly OUT=$HOME/SA_pt158_B  SMRT_JOB_ID=020044 STRAIN_NAME=SA_pt158_B SPECIES="Staphylococcus aureus"
    prokka_annotate     OUT=$HOME/SA_pt158_N  SMRT_JOB_ID=020095 STRAIN_NAME=SA_pt158_N SPECIES="Staphylococcus aureus"

When you run `rake multi[$TASK_FILE]`, with the filename of your task file in the brackets, a `screen` session will be created and split vertically into multiple windows, each of which will run `rake` with the various parameters you put on that line.

If you don't want as many splits, and would like some of the jobs to run in series, also specify a number `$N` that is smaller than the number of lines in your task file, like so:

    $ rake multi[/path/to/task/file,3]

You will almost certainly need to run `rake multi` on an interactive node or the process will exceed resource limits.

### Running as a `bsub` task

You may also want to run the pipeline as a non-interactive job on the cluster.  The benefit of this approach is that you can reserve specific resources in advance to decrease the likelihood of the job running out of memory or exceeding other system limits.  For this, the `scripts/example.post-assemble-pathogen` should be copied, modified as appropriate, and then can be submitted with `bsub` as in the following example:

    $ bsub -R 'rusage[mem=4000] span[hosts=1]' -m "bode mothra" -P acc_PBG -W "24:00" \
            -L /bin/bash -q premium -n 12 -J CD00246 \
            -o "%J.stdout" -eo "%J.stderr" \
        post-assemble-pathogen \
            SMRT_JOB_ID=020486 \
            STRAIN_NAME=CD00246 \
            SPECIES=Cdiff \
            OUT=scratch/out/CD00246_020486 \
            prokka_annotate

### Dependency graph

This Rakefile is able to build a dependency graph of its intermediate files from itself.  Use the `rake graph` task for this; it will be generated at `$OUT/pathogendb-pipeline.png`. Paths through the `check` task are filtered out of the diagram for clarity, since almost every task depends on it.

![Dependency graph](https://pakt01.u.hpc.mssm.edu/pathogendb-pipeline.png?)

## Regression tests

Can be run on any Minerva node:

    $ source scripts/env.sh
    $ bundle install --deployment
    $ rake spec

This runs all tests defined in `spec/*_spec.rb`, except by default, tests marked as `:speed => 'slow'` are skipped. `rake spec[all]` will run absolutely every test.

Prepend `DEBUG=1` to `rake spec` if you want the temporary directory created for a test saved if it fails. Its location will be printed to the terminal.

## Other notes

Although all [RAST][rast]-related tasks are now deprecated, this pipeline downloads and installs the [Network-based SEED API package](http://blog.theseed.org/servers/installation/distribution-of-the-seed-server-packages.html) into `vendor/sas`.  Documentation for some of the included executables and the Perl API are also on that page.  Within this package, we are overriding `RASTserver.pm` with our own version in `lib/perl`, because we need to support setting use of a proxy via the `HTTP_PROXY` environment variable.
