# Analyses config (plan §5.4)

YAML given with `--analyses`. Schema: `assets/schema_analyses.json`. At startup the pipeline
checks names (unique, `[A-Za-z0-9._-]+`), types, factor `levels`, and that every site set
(`of` and `quartiles_by.analysis`) refers to existing analyses; full validation happens in
site calling (WP4).

```yaml
analyses:
  - name: incellulo                 # treatment analysis
    type: treatment
    subset: {celltype: [HepG2, 293T]}      # optional; column -> allowed values
    random_effects: ""            # extra terms appended to "delrate ~ treat + (1|rep)"
    sesoi: 0.05
    fdr: 0.05
    plot_all_sites: true
    colors: {modified: "#c154c1", input: "#eee8aa"}
  - name: PUS7_dep_Both_WT_v_KD   # generic factor analysis
    type: factor
    subset: {celltype: [HepG2, 293T], vector: [P101, P102]}
    factor: vector
    levels: [P101, P102]          # [baseline, experimental]; dd = experimental - baseline
    random_effects: "(1|celltype)"
    direction: positive           # positive | both
    sesoi: 0.05
    fdr: 0.05
  - name: PUS7_dep_Both_WT_v_OE
    type: factor
    subset: {celltype: [HepG2, 293T], vector: [P4, P3]}
    factor: vector
    levels: [P4, P3]              # P4 = WT (baseline), P3 = OE (experimental)
    random_effects: "(1|celltype)"
site_sets:                        # optional
  - name: PUS7_dep_union
    op: union                     # union | intersection | difference
    of: [PUS7_dep_Both_WT_v_KD, PUS7_dep_Both_WT_v_OE]
    quartiles_by: {analysis: incellulo, column: delta_delrate}   # optional
```

Every name under `of` and `quartiles_by.analysis` must be an analysis defined above.

Defaults: `random_effects: ""`, `sesoi: 0.05`, `fdr: 0.05`, `direction: positive`,
`plot_all_sites: --plot_all_sites`, `p_adjust: BH` (BH across the tested sites of the analysis;
`legacy` reproduces the paper's per-site no-op, R-30). `assets/analyses_example.yaml` (WP4/WP5) reproduces the
paper's Figure 2 and Figure 3 analyses.

Example: [examples/analyses.yaml](examples/analyses.yaml).
