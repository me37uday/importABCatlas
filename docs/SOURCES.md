# Sources and verification scope

Access/review date: 8 October 2026. Primary sources were consulted online. Their rendered example timings and file sizes were not reused as this package's measurements. Data schema documentation can change independently of cached file releases; record the actual selected manifest in every run.

S1. Allen Institute. ABC Atlas data access introduction, catalogue and release notes. Documents release 20260711 and SDK 1.4.0; includes earlier normalization and schema corrections.
https://alleninstitute.github.io/abc_atlas_access/intro.html

S2. Allen Institute. Official ABC Atlas access SDK source and installation.
https://github.com/AllenInstitute/abc_atlas_access
https://alleninstitute.github.io/abc_atlas_access/notebooks/getting_started.html

S3. Allen Institute. Measured Allen MERFISH resources, linked from the source catalogue; related primary atlas publication Yao et al. 2023.
https://alleninstitute.github.io/abc_atlas_access/intro.html
https://doi.org/10.1038/s41586-023-06812-z

S4. Allen Institute. Imputed MERFISH access notebook; Zhuang spatial atlas primary publication.
https://alleninstitute.github.io/abc_atlas_access/notebooks/merfish_imputed_genes_example.html
https://doi.org/10.1038/s41586-023-06808-9

S5. Allen Institute. PMDBS mapping to the WHB taxonomy; native expression overview linked in catalogue.
https://alleninstitute.github.io/abc_atlas_access/notebooks/asap_pmdbs_siletti_taxonomy.html

S6. Allen Institute. SEA-AD Multiregion clustering, identifiers, donor/library/disease metadata and taxonomy.
https://alleninstitute.github.io/abc_atlas_access/notebooks/sea-ad_multiregion_clustering_analysis_and_annotation.html

S7. Allen Institute. SEA-AD Multiregion expression access.
https://alleninstitute.github.io/abc_atlas_access/notebooks/sea-ad_multiregion_10X_RNASeq.html

S8. Allen Institute. SEA-AD CaH expression access.
https://alleninstitute.github.io/abc_atlas_access/notebooks/sea-ad_cah_10X_RNASeq.html

S9. Allen Institute. HMBA spatial expression and source panels.
https://alleninstitute.github.io/abc_atlas_access/notebooks/hmba_bg_spatial_genes.html

S10. Allen Institute. HMBA individual-species expression description. Native metadata discovery remains a live-validation gap.
https://alleninstitute.github.io/abc_atlas_access/descriptions/HMBA-10XMultiome-BG.html

S11. Allen Institute. HMBA Macaque Patch-seq notebook.
https://alleninstitute.github.io/abc_atlas_access/notebooks/hmba_bg_macaque_patchseq.html

S12. Allen Institute. Aging mouse expression notebook. Consensus and developing-mouse notebooks are linked in S1.
https://alleninstitute.github.io/abc_atlas_access/notebooks/Zeng_Aging_Mouse_10x_snRNASeq_tutorial.html

S13. SeuratObject. Object and assay constructors; transformed data must retain their transformation semantics.
https://satijalab.github.io/seurat-object/reference/CreateSeuratObject.html
https://satijalab.github.io/seurat-object/reference/CreateAssay5Object.html

S14. Official maintainer documentation for Git/GitHub, roxygen2, R and Bioconductor. These are reference destinations for the beginner guide, not evidence that the candidate has passed a check.
https://git-scm.com/docs/git-pull
https://git-scm.com/docs/git-push
https://docs.github.com/en/pull-requests
https://roxygen2.r-lib.org/articles/roxygen2.html
https://cran.r-project.org/
https://bioconductor.org/install/

S15. Bioinformatics author guidelines: Application Note length, availability, installation expectations, structured abstract and LLM policy. The lowercase and differently capitalised indexed versions can expose differing cached excerpts; confirm the live publisher policy immediately before submission.
https://academic.oup.com/bioinformatics/pages/author-guidelines
https://academic.oup.com/bioinformatics/pages/submission_online

S16. ISCB acceptable use of LLMs policy, updated 3 April 2025, applied by Bioinformatics. It disallows LLM-drafted paper sections/abstracts while allowing verified code/documentation assistance and background discovery. Merely disclosing or lightly editing a generated draft does not establish compliance.
https://www.iscb.org/about-iscb/policy-statements-bylaws-and-legal-documents/acceptable-use-of-large-language-models-policy

The package implementation and synthetic tests are new work in this delivery; their evidence is in validation/. Claims about runtime, operating-system support or all-dataset success require the supplied tests to run on actual target environments and source data. Data licences belong to their respective resources and are separate from the still-unresolved code licence.
