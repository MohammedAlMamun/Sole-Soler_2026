# Smc5 Project 2026

Analysis code for **“Smc5/6 couples ATPase-dependent translocation to sister chromatid disjunction during DNA replication.”**

## Abstract

The Smc5/6 complex is an ATP-dependent motor capable of topologically entrapping DNA segments and organizing transcription-induced supercoiled DNA. While Smc5/6 normally works during S phase to facilitate sister chromatid disjunction, its ATPase-dependent role on replicating DNA remains unexplored. Here, we show that inhibiting ATP hydrolysis traps Smc5/6 on DNA, restricting its movement away from its loading site. Using this defect as an in vivo reporter, we identify new Smc5/6 loading sites in yeast, including origins of replication and replication termination sites. Further, we show that following origin firing, Smc5/6 chases forks in an ATP hydrolysis-dependent manner. Moreover, we demonstrate that ATPase-dependent dynamics are essential for chromosome segregation, and loading of ATPase mutants severely interferes with sister chromatid disjunction. Our findings suggest that Smc5/6 couples ATP hydrolysis-dependent translocation on replicating DNA to the physical separation of sister chromatids behind forks.

<img width="468" height="210" alt="Smc5/6 project overview" src="https://github.com/user-attachments/assets/5312401b-02d2-46c7-8029-e666e09f5154" />

## Scripts

| Script | Main function(s) | Purpose |
|---|---|---|
| `Smc5_0826.R` | `ChIPseq_Primary_Analysis()`, `BrDUseq_Primary_Analysis()` and downstream analysis functions | Main analysis toolkit: aligns paired-end ChIP-seq and BrDU-seq reads, calculates strand-specific binned coverage and enrichment ratios, calls peaks, and produces origin-centred and pairwise comparison analyses. |
| `Smc5_0826_Run.R` | — | Example driver script showing how to run the primary ChIP/BrDU processing, early/late-origin analyses, combined figures, and origin-distance comparisons. |
| `Smc5_Plot_ChIP_Early_Late.R` | `Plot_ChIP_Early_Late()` | Generates early- and late-replication-origin profiles and Watson/Crick strand-bias summaries from previously processed ngsAnalyser results. |
| `Smc5_rDNA_analysis_script.R` | `SubRead_Multiple_Aligner_rDNA()`, `EqualDistRatio()`, `Calc_parallel()`, `PlotEnrichments()` | Performs the rDNA-specific multimapping workflow, calculates ChIP/Input enrichment, normalises rDNA profiles against genome-wide background, and plots the results. |
| `comparative_analysis_smc5.R` | `Calc()`, `PlotEnrichments()`, `ChIP_Plot_function()` | Compares four Smc5 conditions at consensus peaks, origins, tRNAs, Ty elements, and transcription-associated features, and generates average-enrichment and whole-genome plots. |
| `WholeChromosome_Plotter_0626_5c.R` | `WholeChromosome_ChIP_Plotter()`, `WholeGenome_Plotter()` | Prepares ChIP profiles, estimates non-peak background, integrates genomic annotations, and produces annotated chromosome- or genome-wide PDF plots. |
| `run_whole_chromosome_plotter_0626_5c.R` | — | Editable example configuration for running the whole-chromosome/whole-genome plotter on one sample. |
| `FourSample_Chromosome_Comparison_0626_5c.R` | `FourSample_Chromosome_Comparison_Plotter()` | Produces aligned four-sample ChIP profile comparisons for a chromosome or a selected genomic interval using shared annotations and plotting scales. |
| `run_four_sample_chromosome_comparison_0626_5c.R` | — | Editable example configuration for generating a four-sample chromosome comparison. |

The included BED, XLSX, image, and support files provide replication-origin coordinates, genomic annotations, and plotting resources used by the scripts. The example runner scripts contain project-specific absolute paths that must be adjusted for the local data and software installation.

## License

This repository is available under the [MIT License](LICENSE).
