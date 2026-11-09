# The rMLST arm

BRIDGE's rMLST arm calls one or more species from ribosomal MLST (rMLST) reads
using runtime allele clusters and cluster consensus sequences. It does not
require precomputed clusters.

It normally runs as part of `bridge run` / `bridge sample` (`--arms rmlst`),
which fills in the reference paths, tool paths, thresholds and the shared BLAST
cache from `config.yaml` and writes to `<out>/<SAMPLE>/rmlst/`. The arm can also
be run on its own, which is what this document describes:

```bash
python -m bridge.rmlst --help
```

## Requirements

The arm needs Python 3.10 or newer (installed with BRIDGE) and these programs on
`PATH`, or supplied with the corresponding command-line options:

- kallisto
- mash
- minimap2 (the default short-read aligner) or bwa-mem2
- samtools
- CliqueSNV (+ Java)
- ViralConsensus
- blastn and makeblastdb

Every normal run checks the Python and external dependencies required by its
selected route before processing data. The same check can be called directly:

```bash
python -m bridge.rmlst --check-deps
```

The check respects executable overrides, `--short-read-aligner`, and
`--viralconsensus-only`. It exits with status 1 and reports every missing or
outdated dependency when the environment is incomplete. `bridge check` runs the
equivalent check for every arm.

CliqueSNV may be supplied as either an executable or a JAR using `--cliquesnv`.
The `BRIDGE_CLIQUESNV` environment variable provides its default value when
set. A JAR is launched as `java` with the configured heap options followed by
`-jar`; it is never executed directly.

## Quick start

Run the arm on paired-end reads with an allele FASTA and an rMLST sequence
definition table:

```bash
python -m bridge.rmlst \
  --allele-target /path/to/rmlst_alleles.fasta \
  --rmlstseqdef /path/to/rmlstseqdef.tsv \
  --reads /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz \
  --sample-id sample \
  --output-folder results/sample \
  --threads 8 \
  --memory-limit-gb 32 \
  --cliquesnv-time-limit-seconds 600
```

The main species-level results are written to `results/sample/species_calls.csv`.
An editable PDF summary is written to `results/sample/called_species_pi.pdf`,
and `results/sample/rmlst_manifest.json` records the run configuration.

## Usage

The required arguments are:

- `--allele-target`: deduplicated rMLST allele sequences in FASTA format.
- `--rmlstseqdef`: allele-to-species and rST definition table in TSV format.
- `--reads READ1 READ2`: paired-end FASTQ files.
- `-o, --output-folder`: directory in which to write the run outputs.

Common options include:

- `--sample-id`: sample label used in output tables (default: `sample`).
- `-t, --threads`: total worker-thread budget (default: 4).
- `--memory-limit-gb`: overall memory budget used to size CliqueSNV jobs.
- `--kallisto-index`: reuse an existing full-reference kallisto index.
- `--short-read-aligner {bwa-mem2,minimap2}`: select the cluster-read aligner.
- `--cliquesnv-time-limit-seconds`: limit each CliqueSNV search while allowing
  CliqueSNV to return its best result when the limit is reached.
- `--viralconsensus-only`: disable subconsensus reconstruction and use one
  ViralConsensus sequence per runtime cluster.
- `--min-unambiguous-loci` (default 3) and `--require-accepted-consensus`: the
  reporting bar (see below; `rmlst:` in `config.yaml` under `bridge`).
- `--distinguish-subspecies`: use the `subspecies` column immediately after
  `species` in the sequence-definition table to call species/subspecies
  combinations separately. A blank subspecies is retained as a distinct
  `[no subspecies annotation]` category. Without this flag, all rows for a
  species remain pooled.
- `--blastdb-cache-dir`: directory holding the shared `.bridge_blastdb`
  cache. Defaults to the parent of `--output-folder`, which is right when that
  is `<run_root>/<SAMPLE>`. Pass the run root explicitly whenever the output
  folder sits deeper inside it (as `bridge` does for `<run_root>/<SAMPLE>/rmlst`),
  so concurrent samples keep sharing one database instead of each building
  their own.
- `--keeptmp`: retain `<output-folder>/work` and intermediate diagnostic
  outputs for troubleshooting (`--keep-temp` on `bridge run` / `bridge sample`).
- `--verbose`: show detailed pipeline messages, commands, and output from
  internally called tools in the terminal.

Executable paths can be overridden with `--kallisto`, `--mash`, `--bwa-mem2`,
`--minimap2`, `--samtools`, `--cliquesnv`, `--viral-consensus`, `--blastn`, and
`--makeblastdb`.

By default, terminal output is limited to main pipeline stages plus warnings and
errors. The complete detailed stream, including commands and output from
internally called tools, is always written to `<output-folder>/rmlst.log`. Use
`--verbose` to show that stream in the terminal as well.

## Workflow

1. kallisto quantifies the complete allele reference and emits 100 bootstrap
   abundance tables without constructing an enormous full-reference BAM header.
2. Candidate selection requests support in all 100 bootstraps independently
   for each locus. If that would leave a locus empty, only that locus relaxes
   to its highest observed bootstrap support; all alleles tied at the effective
   threshold are retained. Requested and effective thresholds are reported for
   every candidate.
3. A small candidate-only kallisto index is built and used to emit the
   pseudobAM required for read extraction.
4. Candidate alleles are clustered independently per sample and locus with
   Mash single linkage; the cluster medoid is the alignment representative.
5. Pseudoaligned fragments are assigned to every compatible runtime cluster
   (`--cluster-read-assignment compatible`) or only to an exclusive cluster.
6. Every runtime cluster with enough assigned read pairs is aligned with the
   selected short-read backend (minimap2 by default, or BWA-MEM2).
   Quality-filtered nucleotide counts are scanned for a possible minor allele.
   If the second-highest nucleotide count at every position is below 5,
   ViralConsensus produces a single consensus. Otherwise the cluster is passed
   to CliqueSNV, which may emit multiple linked haplotypes and estimated
   proportions. This routing gate uses neither strand balance nor allele
   frequency. Zero-read clusters are retained with an `insufficient_reads`
   status.
7. A candidate haplotype is retained only when its estimated proportion is
   strictly greater than 1%, its estimated support is at least 20 read pairs,
   its projected median depth is at least 10x, at least 80% of its locus reaches
   10x projected depth, and every callable distinguishing SNV is supported on
   both strands without being concentrated at read ends. Unsupported SNV calls
   are masked to `N` instead of causing the entire locus to be discarded. Short
   terminal truncations of at most 15 bp are padded with `N` only when the
   resulting partial sequence is compatible with a known runtime-cluster
   allele. An unambiguous exact member match can substitute for medoid-relative
   length and variant validation, but none of these corrections bypass
   abundance, support, depth, or breadth requirements. Projected depth is the
   observed cluster depth multiplied by the haplotype proportion, so the 1%
   threshold cannot override inadequate sequence evidence. Malformed CliqueSNV
   records are skipped individually rather than invalidating other haplotypes
   from the same cluster. A single 100% CliqueSNV haplotype is replaced by
   ViralConsensus because it provides no subconsensus information.
   ViralConsensus is also used as a fallback when CliqueSNV fails or none of
   its haplotypes pass these filters.
8. A cluster-member exact hit skips BLAST. Only unmatched complete consensuses
   and callable partial consensuses are combined into one multithreaded BLAST
   query against a cached reference database. A partial subconsensus is
   accepted when exactly one locus allele agrees at every callable base in a
   bounded, ungapped terminal alignment. Alleles may contain up to 15
   additional terminal bases. A shorter allele is allowed only when the
   discarded consensus overhang contains no callable bases. Multiple compatible
   alleles remain ambiguous.
9. By default, a fully callable consensus absent from its runtime cluster is
   checked position by position against same-length cluster members. If every
   consensus base is supplied by at least one member, it is marked
   `composite_from_members` and retained only as ambiguous evidence. Disable
   this filter with `--no-filter-composite-consensus`.
10. Runtime clusters producing the same consensus are merged. The cluster with
    the highest estimated count supplies the allele and species assignment; the
    duplicate clusters' estimated counts are reassigned to that winner.
    Distinct consensuses remain separate and can represent multiple species or
    strains.
11. Species-specific allele combinations are reported directly. Ambiguous and
    failed consensus evidence is passed to a mixture-aware EM fallback.

### The reporting bar

The caller does not collapse the result to a single organism: distinct
consensuses and multiple simultaneously supported species are retained. A
species is reported only when at least `--min-unambiguous-loci` loci (default 3)
support it **unambiguously**, that is, map to that species and no other.

Loci that are merely *compatible* with a species do not count toward that bar,
and the distinction matters. A profile-wide intersection collapses onto a
species as soon as one single locus is exclusive to it, because
`{X} ∩ {X, Y, Z} ∩ … = {X}`. With a bar of one unambiguous locus, a near
neighbour could be reported on one exclusive locus plus shared evidence
belonging to the organism actually present; on a synthetic mixture panel this
reported *Streptococcus pneumoniae* at around a thousandth of the sample, in
mixtures containing a mitis-group culture and no pneumococcus.

There is deliberately **no abundance floor** on a reported species. Abundance
does not separate a real minor organism from a near-neighbour artefact:
genuine organisms at 1% of a mixture were estimated below some of the spurious
calls, whereas evidence breadth separates them cleanly.

Two columns make the distinction visible:

- `exact_profile_unambiguous_loci`: loci mapping to this species alone.
- `exact_profile_loci`: loci merely compatible with it. Always the larger
  number, and for a near neighbour it is much larger.

### Demoted, not deleted

A species demoted by the bar is **kept**, with `reported=false`. Where it is a
second acceptable name for an organism that *was* reported (the usual case,
since mitis-group and other close relatives share allele evidence),
`alternative_to` names that organism. This separates "one organism with an
ambiguous name" from "two organisms", which two reported rows cannot express.

`alternative_to` is set only where **both** hold:

- the species held up as an allele profile in its own right and fell short
  *only* on exclusive loci. It is a near miss, not background; a species with
  no exclusive locus anywhere never named itself alone and is not a second name
  for anything.
- its genus has a reported species for it to be an alternative *to*. A demoted
  species standing alone in its genus is a weak call, and is left unlabelled.

On a synthetic mixture panel a shared genus alone matched thousands of
EM-fallback background rows for every real demotion; the evidence test
separates the two classes completely (every real demotion carried exactly one
exclusive locus, every EM-fallback row none).

NOTE!!! Two genuinely distinct members of the same genus in one
sample can still be labelled as one organism under two names, if the minor one
lands one or two exclusive loci and its profile intersection holds. Ruling that
out needs a test of whether the two species' evidence actually overlaps, which
this column does not attempt.

### Tools and resources

Supply an existing full-reference kallisto index with
`--kallisto-index /path/to/reference.kallisto.idx`. This takes precedence over
adjacent-index discovery and automatic index construction. The arm reads the
supplied index's k-mer size with `kallisto inspect` and uses the same value
when constructing the candidate-only index. The supplied index must have target
identifiers that correspond to `--allele-target`. Its path, source, and detected
k-mer size are recorded in `rmlst_manifest.json`.

By default, minimap2 performs runtime-cluster read alignment using the
parameters `-ax sr -k 17 -w 7 -m 20 -s 30` and a MAPQ threshold of 10 when
calculating variant evidence. Use `--short-read-aligner bwa-mem2` (or
`--aligner bwa-mem2`) to select BWA-MEM2 instead; it uses a MAPQ threshold of
20.

`--cliquesnv-time-limit-seconds` passes CliqueSNV's native `-tl` limit to
every CliqueSNV invocation. When the limit is reached, CliqueSNV stops its
haplotype search and returns a consensus haplotype rather than being killed.
If the option is omitted, CliqueSNV uses its own default of 10,800 seconds
(3 hours); `bridge` passes `hpc.cliquesnv_time_limit_seconds` from
`config.yaml` (600 s). `--cliquesnv-runtime-limit-seconds` is accepted as an
alias.

`--viralconsensus-only` disables CliqueSNV and always produces one
ViralConsensus sequence per runtime cluster, irrespective of detected minor
alleles. This bypasses subconsensus reconstruction but retains the normal
consensus depth, breadth, and allele-compatibility checks. CliqueSNV does not
need to be installed in this mode.

`--memory-limit-gb` sets the overall memory budget used for consensus work and
takes precedence over automatic physical-memory or cgroup-limit detection. The
arm reserves 20% for Python, the selected aligner, ViralConsensus, and native
Java overhead, then divides the remaining 80% evenly across the maximum number
of concurrent consensus tasks. Each CliqueSNV process receives that per-task
value as an explicit Java `-Xmx` limit. On a 16 GB, four-thread run, for
example, two concurrent tasks receive 6,553 MB of Java heap each.

During consensus generation, the single `--threads` budget is divided into two
threads per cluster, with up to `floor(threads / 2)` clusters processed
concurrently. The Java heap budget is also divided by that concurrent task
count. kallisto and batched BLAST can still use the full thread budget.

Important defaults:

```text
short-read aligner:          minimap2 (MAPQ >=10)
bootstrap support threshold: 100/100, with per-locus fallback
runtime Mash distance:       0.01
minimum assigned fragments:  1
subconsensus proportion:     >0.01
subconsensus support:        >=20 estimated read pairs
subconsensus depth:          >=10x projected median
subconsensus breadth:        >=0.80 at 10x projected depth
variant read-end distance:   >=5 bases
CliqueSNV routing depth:     second-highest base count >=5
CliqueSNV time limit:        native default 10800 seconds per invocation
minimum profile loci:        3
minimum unambiguous loci:    3
cluster read assignment:     compatible
composite consensus filter:  enabled
```

## Outputs

- `species_calls.csv`
- `consensus_groups.csv`
- `called_species_pi.pdf`
- `species_allele_fastas/<species>.fasta`
- `rmlst_manifest.json`
- `rmlst.log`

These are the only outputs retained after a normal successful run.

### Samples with nothing to find

A negative control, an extraction blank, or a sample whose pathogen load sits
below the assay floor will legitimately yield no candidate alleles and no
runtime clusters. That is a **result**, and the run ends cleanly: an empty
`species_calls.csv` and a manifest carrying

```json
"outcome": "insufficient_material",
"outcome_stage": "kallisto",
"outcome_reason": "No alleles passed the Kallisto candidate rule"
```

Because `species_calls.csv` exists, a batch driver that uses it as its
completion marker does not resubmit such a sample.

### Where a sample's time went

`rmlst_manifest.json` carries a `runtime_breakdown` with per-stage and per-tool
wall clock:

```json
"runtime_breakdown": {
  "stages": {"kallisto": 812.4, "cluster_consensus": 1503.9, ...},
  "tools":  {"cliquesnv": {"calls": 37, "seconds": 4211.6}, ...}
}
```

Stage totals alone cannot attribute cost: `cluster_consensus` wraps the
aligner, CliqueSNV, ViralConsensus, BLAST and a large amount of in-process
Python behind one number. Tool seconds are summed across calls, so with
concurrent clusters they can exceed the stage's own elapsed time; that is the
intended reading, being work done rather than time spent waiting.

### Diagnostic outputs

With `--keeptmp`, the arm also retains:

- `kallisto_candidates.csv`
- `runtime_clusters.csv`
- `runtime_cluster_abundance.csv`
- `runtime_cluster_read_counts.csv`
- `runtime_cluster_consensus.csv`
- `species_allele_fastas.csv`
- `em_fallback.csv`
- `cluster_consensus/<locus>/<cluster>/<cluster>_SC####.fasta`
- `cluster_consensus/<locus>/<cluster>/blast_<cluster>_SC####.tsv`
- `work/`

`kallisto_candidates.csv` records `requested_bootstrap_threshold`,
`effective_bootstrap_threshold`, and `bootstrap_threshold_relaxed` for every
retained allele, making each per-locus fallback auditable.

`runtime_cluster_consensus.csv` contains one row per retained subconsensus,
with `subconsensus_id`, `haplotype_proportion`, `supporting_fragments`,
`projected_median_depth`, `projected_breadth`, and `variant_evidence_passed`
fields. Its `consensus_method` field records `viralconsensus` or `cliquesnv`,
and `composite_from_members` records whether the default composite filter
demoted the sequence to ambiguous evidence. A cluster for which no candidate
passes has one `_SC0000` audit row with zero estimated abundance and a failure
status.

### `species_calls.csv`

Multiple rows may have `reported=true`. Direct calls use
`call_method=exact_allele_profile`; unresolved species use
`call_method=em_fallback`. Evidence locus counts are grouped by their source in
four columns:

- `allele_evidence_unambiguous_loci`: accepted consensus-allele evidence mapping
  to one species. This is the count compared against `--min-unambiguous-loci`.
- `allele_evidence_ambiguous_loci`: accepted consensus-allele evidence mapping
  to multiple species, including resolved genus placeholders.
- `em_evidence_unambiguous_loci`: all positive-abundance evidence supplied to
  EM that maps to one species, including accepted allele consensus evidence.
- `em_evidence_ambiguous_loci`: all positive-abundance evidence supplied to EM
  that maps to multiple species, including accepted allele consensus evidence.

Cluster fallback evidence can influence EM `pi` and the `em_evidence_*` fields,
but it does not become accepted consensus-allele evidence or qualify a direct
profile call.

`pi` is the EM estimate of each species' share of the sample. It ranks
correctly but compresses toward the middle (minor organisms are overestimated,
majority organisms underestimated), so use it as a rank, not a quantity.

### How well assembled the evidence was

Seven coverage columns summarise the assembly behind each call:

- `mean_callable_breadth`, `min_callable_breadth`: fraction of the allele length
  with callable bases. The mean says how complete the consensus assembly was;
  the minimum exposes the weakest locus supporting the call.
- `median_projected_depth`, `min_projected_depth`: per-locus projected read
  depth (a per-locus figure summarised across loci, not a total).
- `total_callable_bases`, `total_allele_length`: the same breadth information in
  absolute terms, so the ratio can be recomputed.
- `supporting_fragments`: read pairs supporting the winning consensus. This is
  a count observed during assembly, distinct from `unambiguous_est_counts`,
  which is kallisto's EM estimate.

All seven are aggregated over the loci that name the species **alone**, the
same loci `allele_evidence_unambiguous_loci` counts. Loci shared with a near
neighbour are excluded, or a well-covered relative would inflate the figure.
Where a locus resolved into several haplotypes, its best-assembled consensus
speaks for it, so a within-species mixture does not double-count a locus.

They are **blank, never zero**, for a species with no unambiguous evidence:
every `em_fallback` row, and every demoted near neighbour that fell short of the
bar. `consensus_groups.csv` holds the per-locus values these are computed from
and is retained on every run.

### A blank abundance is not a zero one

`pi`, `unambiguous_est_counts` and `ambiguous_est_counts` are written **empty**
for a species called from its allele profile whose consensus groups all carried
zero estimated abundance and so never reached EM. The species is still reported
(its locus evidence is real and is what the reporting bar reads), but no
abundance was estimated for it, and an empty cell says exactly that.

Writing `0.0` would assert "this organism was measured and there was none of
it", which is a different and much stronger claim, and a downstream read-count
gate would then silence a genuine call. Consumers should treat an empty cell as
"unknown", never as zero.

### rST matching and allele FASTAs

`matching_rSTs` lists up to 10 profiles compatible with all callable
species-supporting loci; `matching_rST_count` gives the full number before that
list is truncated. At every constrained locus, rST matching accepts any
compatible allele assigned to that species in the sequence-definition database
and searches all resulting multilocus combinations.

Each reported `exact_allele_profile` call also receives a multi-FASTA containing
the read-derived consensus sequence from every accepted consensus group
supporting that species. Partial consensuses retain their `N` bases; database
matches are described in each FASTA header. Unique exact and partial matches
name the matching allele; ambiguous partial matches report the number of
compatible alleles. In `species_calls.csv`, `allele_fasta` gives the supporting
multi-FASTA path and `allele_fasta_allele_count` gives its number of consensus
sequences. These fields are empty and zero, respectively, for EM-only calls.

## Interpretation and limitations

- A low-depth unique partial match is recorded separately from a complete exact
  match.
- Exactly 1% does not pass the relative threshold. A candidate over 1% can
  still fail the absolute fragment, projected-depth, projected-breadth, strand,
  or read-end filters. Rejection reasons are written to the run log.
- Haplotype estimated counts are the parent runtime-cluster count multiplied by
  the CliqueSNV proportion; evidence below the filters is not redistributed.
- Distinct consensuses at one locus may indicate different species, different
  strains of one species, or incorrect read partitioning.
- kallisto-compatible assignment can duplicate ambiguous fragments across
  clusters. Identical downstream consensuses are assigned to the strongest
  cluster and their counts are combined to prevent duplicate or conflicting EM
  evidence.
- A `Genus sp.` label is removed when definitive names occur in the same
  assignment. Generic-only evidence is shared across definitive species of
  that genus supported elsewhere in the sample; if none exist, the `sp.` name
  is retained.
- Partial matching requires the consensus and candidate allele to have the same
  length. Indel-containing partial consensuses remain unresolved rather than
  being forced to an allele.
- The default temporary directory is `<output-folder>/work` and is removed on
  success or failure. Use `--keeptmp` to retain it for debugging.
- Full-reference fallback queries are batched per sample. The reusable BLAST
  database is stored under `.bridge_blastdb/` and keyed by the reference
  SHA-256 digest. The database is built unlocked, so concurrent samples must
  share one directory and it should be built once up front; `bridge run` does
  this on the submitting node before any job starts.
- The allele reference defines what can be called: an organism outside the
  genera in the rMLST target set is invisible to this arm. BRIDGE's 16S arm is
  the fallback for that case.
