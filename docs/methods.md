# rmrodell/nanomodamp: Methods

> WP0 stub. WP6 completes this page (plan §7 WP6). Section headings list what it must cover.

## Preprocessing

### Orientation and adapter trimming (D4)
_To be written (WP6): both linked adapters required; legacy non-anchored behaviour (R-04)._

### Libraries without ONT adapters or UMIs (D29, D31)
The standard Nano-BID-Amp protocol adds ONT adapters on both ends of every amplicon and a
10-nt UMI at the 3′ end. The pipeline defaults (`orientation_adapters=ont`,
`ont_adapter_mode=linked`, `umi=true`) assume this and apply to **all** library types: endogenous,
MPRA in vitro and MPRA in cellulo. Reads are oriented by requiring both ONT adapters, the UMI is
extracted, and PCR duplicates are collapsed with UMICollapse.

Two opt-in exceptions reproduce how the paper processed **legacy libraries** that lacked some of
these elements. They are not chosen from the library type or from any sample label, and the
pipeline sets them only in `conf/test_golden.config` for the affected golden samples:

- **No ONT adapters and no UMIs** (`orientation_adapters=pool`, `umi=false`; R-26). The paper's
  in-cellulo MPRA samples from sequencing run 20241114 were sequenced from such a library. Reads
  are oriented with the MPRA pool adapters (sense `pool_adapter_5p…pool_adapter_3p`, antisense
  `pool_adapter_antisense_5p…pool_adapter_antisense_3p`, `trim1_*` length and overlap), the
  antisense reads are reverse-complemented and merged, and no UMI extraction or deduplication is
  done, so the counts include PCR duplicates.
- **No 5′ ONT adapter** (`ont_adapter_mode=three_prime_only`; R-29). The paper's endogenous
  samples from run 20250418 lack the 5′ adapter, and the paper trimmed only the 3′ adapter
  (`cutadapt -a`, untrimmed reads kept, no antisense pass). Under the default linked mode no read
  from that run would pass.

New data generated with the standard protocol, including new MPRA in-cellulo data, must use the
defaults.

### UMI extraction and deduplication (D7)
_To be written (WP6)._

### Mapping (D5)
_To be written (WP6): why `-ax sr` for MPRA and `-ax splice -uf` for endogenous; how presets change deletion representation; how and when to change them._

## Counting
_To be written (WP6): the denominator (D9), coverage filtering (D10, R-13 including the `>= 20` sweep wording), the BED convention (R-15)._

## Site calling
_To be written (WP6): treatment/TOST and factor models as equations; equivalence limitation at input = 0 (D25); replicate pairing (D17)._

## Golden test
_To be written (WP6): what `test_golden` checks and what is binding (`legacy_rerun`) vs reference (`published_rerun`, `paper_reference`; D30)._
