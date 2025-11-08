RNA-Seq Transcriptomics Analysis of Cutaneous Leishmaniasis 


This project presents analysis of Single-End Bulk RNA-seq data (https://www.ebi.ac.uk/ena/browser/view/PRJNA525604?show=reads) obtained from skin biopsies of Cutaneous Leishmaniasis (CL) patients and healthy controls. 

The goal was to identify Differentially Expressed Genes (DEGs), analyze Differential Transcript Usage (DTU)/Isoform Switching, and perform Pathway Enrichment Analysis to elucidate the host immune response and molecular mechanisms associated with chronic disease.


The analysis was executed using a structured pipeline consisting of Shell scripting for read processing and quantification and R scripts for DEGs, DTU, and GSEA.

The functional enrichment analysis revealed a distinct immune signature in the skin lesions of Cutaneous Leishmaniasis patients, suggesting a state of chronic, dysregulated inflammation:

Elevated Immune Activation: The disease cohort showed significant enrichment in pathways related to Immune Activation and B Cell-Mediated Immunity.

Molecular Example: A clear upregulation of Immunoglobulin genes (IGH, IGK, IGL) was observed, indicating a robust, and potentially hyperactive, B cell response at the site of infection.

Suppressed Regulatory Mechanisms: Conversely, pathways associated with Apoptosis (Programmed Cell Death) and Cell Adhesion were found to be significantly downregulated.



Biological Interpretation: This suppression suggests that inflamed or abnormal cells persist longer in the tissue (due to failed apoptosis) and that tissue regulation is impaired (reduced cell adhesion), contributing to the chronicity and persistence of the lesion.

