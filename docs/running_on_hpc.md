# Running BRIDGE on an HPC cluster

A step-by-step guide from a fresh copy of the code to a results table on an
HPC cluster with a PBS scheduler. Steps marked **CHECKPOINT** say what you
should see; if you do not see it, fix that before moving on.

Placeholders used below:

| Placeholder | Meaning |
|---|---|
| `<repo>` | where you put the BRIDGE code, e.g. `/home/<user>/bridge` |
| `<out>` | the results folder for this run, on a filesystem compute nodes can write |
| `<RUN>` | a sequencing run id, e.g. `RUN01`; it tags the baseline and every result |
| `<project>` | your PBS project / account code |

**How it fits together.** You run one command, `bridge run`, on a login node.
It submits one compute job per sample (the rMLST, 16S and AMR arms) and one
more job that waits for them all and writes the run tables (baseline,
reconciliation, outcomes). 

---

## 1. Build the Python environment (once)

BRIDGE needs Python 3.10 or newer. 

```bash
module load python3/3.12.1      # use whatever name `module avail python3` lists
python3 --version
```

Put that module name in `config.yaml` under `hpc: python_module:` so every job
loads the same Python.


```bash
cd <repo>
./setup_env.sh 
```

Alternatively, build a conda environment with every tool included
(`environment.yml`, see README) and set `tools: tools_env: null` in
`config.yaml`.

## 2. Configure

Edit `config.yaml`:

- `references:`: paths to the rMLST alleles, index and sequence definitions,
  the MIMt 16S FASTA and index, the rrnDB copy-number table and the ResFinder
  KMA index (see README, "Reference databases").
- `tools:`: binary names or absolute paths, and `tools_env` (the shell file
  each job sources).
- `hpc:`: the scheduler settings: `project`, `queue`, `storage` (if your site
  uses `-l storage`), `python_module`, the control-name prefixes and the
  resource tiers.

Then adapt `tools.env` to your cluster's module names (it is a template).

## 3. Activate and check (every session)

```bash
module purge
module load python3/3.12.1
source <repo>/.venv/bin/activate
source <repo>/tools.env
bridge check --out <out>
```

The `module purge` is worth the extra line: a shell can inherit the *list* of
loaded modules without the `PATH` entries that go with them, so `module load`
does nothing and the programs are missing. `tools.env` warns if that happens.

**CHECKPOINT.** The last line is `Everything needed is present.`

If a line says **MISSING**:

| Line | Fix |
|---|---|
| `cliquesnv` | it is a JAR, not a program on `PATH`: set `BRIDGE_CLIQUESNV` in `tools.env`, or put the JAR's full path under `tools: cliquesnv:` in `config.yaml` |
| another program | fix the `module load` line for it in `tools.env`, then repeat the lines above |
| a reference file | fix the path under `references:` in `config.yaml` |
| `output` | `<out>` is not writable |
| `qsub` shows as a note | you are not on a node that can submit jobs |

## 4. Prepare the sample list

One line per sample: the absolute path of its R1 file. The sample name is the
parent directory, and R2 is found by replacing `_R1_` with `_R2_` (or `_R1.`
with `_R2.`). Blank lines and lines starting with `#` are ignored.

```bash
ls /path/to/reads/*/*_R1*.fastq.gz > sample.list
```

Alternatively, `--reads-dir /path/to/reads` discovers
`<dir>/<SAMPLE>/*_R1*.fastq.gz` itself.

Include the run's negative control(s). Samples whose names start with one of
`hpc.controls_prefixes` (default `NEG`, `NTC`, `Blank`, case-insensitive) are
treated as controls; otherwise name them with `--control <name>`.

## 5. Dry run

```bash
bridge run --reads-list sample.list --out <out> --run-id <RUN> --dry-run
```

This writes every job script under `<out>/oe/scripts/` and submits nothing.

**CHECKPOINT.** Near the end, `dry_run : <number of samples>` and a line naming
the finalise script. `SKIP` lines name samples whose read files could not be
found.

## 6. Submit

```bash
bridge run --reads-list sample.list --out <out> --run-id <RUN>
```

It re-runs the check first and **refuses to submit anything** if a tool or
reference file is missing, so a broken environment costs one message rather
than one failed job per sample. It then builds the rMLST arm's shared BLAST
database once (`BLAST cache: ready: ...`) and prints one line per sample with
its job id.

**CHECKPOINT.** The tally says `submitted : <n>`, followed by
`Finalise job <id> waits for <n> sample job(s)`.

If some submissions fail (for example a per-user queue limit), run the same
command again later: it skips samples that are complete, queued or running, and
submits the rest. `--max-submit N` keeps each batch under a limit.

## 7. While it runs

```bash
qstat -u $USER
ls <out>/*/rmlst/species_calls.csv 2>/dev/null | wc -l
ls <out>/*/16s/marker16s_calls.csv 2>/dev/null | wc -l
ls <out>/*/amr/amr_hits.csv 2>/dev/null | wc -l
```

The counts climb towards the number of samples. Each job writes its log to
`<out>/oe/`. To find failures:

```bash
grep -l "Traceback\|ERROR" <out>/oe/*.o* 2>/dev/null
```

A failed sample is resubmitted by running the step 6 command again; the
finalise step then picks it up.

## 8. When it has finished

```bash
ls -l <out>/identification_summary_all.csv <out>/sample_outcomes_all.csv
```

**CHECKPOINT.** The run tables exist and are newer than the sample jobs.

If they are missing (the finalise job was rejected, or ran before a
resubmitted sample finished), run the step yourself on the login node:

```bash
bridge finalise --out <out> --run-id <RUN>
```

It prints the baseline summary, the number of samples reconciled and skipped
(controls), and the paths of the run tables. See README, "Outputs".

## 9. Changing a threshold later

Everything scientific is in `config.yaml`. After changing anything under
`rpm:`, `reconcile:`, `outcome:` or `taxonomy:`, re-run only
`bridge finalise`; after changing `amr:`, `bridge amr-report` is enough. No
sample job needs to run again unless you change something under `rmlst:` or
`marker16s:`, which the arms read; for those, re-run step 6 with `--force`
(add `--arms 16s` to redo only the 16S arm, which is quick).

## Troubleshooting

**`bridge: command not found`.** The venv is not active in this session
(step 3).

**`Unable to locate a modulefile` in a job log.** The module name under
`hpc: python_module:` or in `tools.env` is wrong for your site. Fix it and
resubmit with step 6.

**`Executable not found: <tool>` in a job log.** `tools.env` did not provide
the tool on the compute node. Fix it and resubmit.

**Every job fails within seconds with the same message.** Read one log in
`<out>/oe/`; the fault will be the same in all of them. Fix it once, try one
sample with `--max-submit 1`, then resubmit the rest.

**A job ends with a walltime error.** The sample was bigger than its resource
tier allowed. Raise `walltime` for that tier under `hpc: tiers:` and resubmit;
only unfinished samples go again.

**You want to start fresh.** Use a new `--out` folder; a new run never deletes
anything under an old one.

**Running without PBS.** `--scheduler local` runs every sample in the current
session (a workstation, or a cluster with a different scheduler inside one
interactive allocation), then finalises.
