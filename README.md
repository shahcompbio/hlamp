# Dynamic evolution of oncogene amplification across 80,000 cancer cell genomes


This is the repository for the manuscript titled `Dynamic evolution of oncogene amplification across 80,000 cancer cell genomes`



**Jake June-Koo Lee<sup>\*1</sup>, Sohrab Salehi<sup>\*1,2,3</sup>, Matthew A. Myers<sup>1</sup>, Melissa Yao, Seongmin Choi<sup>1</sup>, Duaa H. Al-Rawi, Ignacio-Vazquez-Garcia, Eliyahu Havasov<sup>1</sup>, Michelle Wu<sup>1</sup>, Jin Lee, Fathema Uddin, Parvathy Manoj, Pedram Razavi, Samuel Aparicio, Natasha Rekhtman, Kenny K.H. Yu, Helena A. Yu, Charles M. Rudin, Andrea Ventura, Andrew McPherson<sup>1</sup>, Marc J. Williams<sup>1</sup>, Sohrab P. Shah<sup>1</sup>**


\* These authors contributed equally  


<sup>1</sup> The Halvorsen Center for Computational Oncology, Memorial Sloan Kettering Cancer Center, New York, NY, USA  
<sup>2</sup> Computational Oncology, Department of Epidemiology and Biostatistics, Memorial Sloan Kettering Cancer Center, New York, NY, USA  
<sup>3</sup> Irving Institute for Cancer Dynamics, Columbia University, New York, NY, USA  
<sup>5</sup> Human Oncology and Pathogenesis Program, Memorial Sloan Kettering Cancer Center, New York, NY, USA  

**Correspondence**: [leej39@mskcc.org](mailto:leej39@mskcc.org)


![plot](./figures/main.png)


## Abstract

Ddynamic evolution of oncogene amplification across 80,000in cancer cell genomes

Copy-number amplification is a major mechanism of oncogene activation and a therapeutic target, yet its evolution in human cancers is incompletely understood. We analyzed single-cell whole-genome sequencing (scWGS) data from >80,000 cancer cells in 100 individuals tumors from with major cancer types (ovary, breast, lung, and brain) and experimental models to resolve mechanisms and evolution of oncogene amplification. 
Copy-number distributions across cells revealed two patterns: narrow, uniform peaks from symmetric segregation, consistent with intrachromosomal amplifications (ICamps), and broad, heavy-tailed variationbility with extreme outliers (>100 copies/cell) indicative of extrachromosomal circular DNA (ecDNA). A probabilistic mixture model of these distributions classified 74 (15%) of 503 amplified regions as ecDNA, including all validated experimental models. These ecDNAs most frequently involved MYC, EGFR, and MDM2, and were enriched in glioblastoma and lung cancers, whereas high-grade serous ovarian and hormone receptor-negative breast cancers predominantly showed ICamps.
Notably, ICamp events showed significant subclonal specificitysubclonal variability, with 227 (53%) of 429 events displaying multiple modes in the copy-number peaksdistribution. These modes were congruent with other nuclear genomic features from phylogenetic analysis, suggesting lineage inheritance patterns of symmetric division and clonal expansion. Diversification arose via numeric mechanisms (aneuploidy or whole-genome doubling) and subclone-specific structural variants (e.g., ongoing breakage-fusion-bridge cycles or chromothripsis). Copy-number modulation was bidirectional, including loss of amplified derivative chromosome via subclone-specific aneuploidy.
Joint analysis of copy-number and structural variants at single-cell level revealedresolution uncovered several mechanisms of ecDNA-driven cancer evolution. First, ecDNA-positive tumors often harbored multiple species co-amplifying oncogenes. Remodeling of these ecDNAs through internal rearrangements and recombination between distinct species was often observed, leading to distinct patterns of oncogene co-evolution consistent with co-selected events. Second, a convergent evolution through recurrent acquisition of distinct EGFR-containing ecDNAs was observed in a glioblastoma, and suggesting the potential role of ecDNA loss in shaping subsequent evolution. Third, scWGS of isogenic cell lines clearly distinguished present ecDNAs from evolutionarily historical ecDNAs that underwent existent as homogeneously staining regions, revealing genomic rearrangements resulting in  events that led to chromosomal re-integration. Finally, we show that cComparison between circular genome graph-based prediction versus the copy-number distribution-based prediction of ecDNAs  we present hereof ecDNAs revealed substantial discrepancy with notable tissue-type specificity and potential over-calling of ecDNA by previous methods. 
In conclusion, the analysis of >80,000 single cell cancer genomes reported here reveals interpretable distributions of oncogene amplifications consistent with distinct ICamp and ecDNA generative processes. This approach, which refines ecDNA identification beyond standard genome graphs, further elucidates how distinct mechanisms of oncogene amplification diversify cancer cell populations in tumors. We posit new modes of amplification-driven tumor evolution, including distinguishing extant from vestigial ecDNAs in a tumor’s evolutionary history. scWGS enables accurate discrimination of ecDNA from ICamp, reveals their diversification and co-evolution, and refines ecDNA detection beyond bulk genome graphs, informingWe suggest these new insights will inform patient selection for emerging ecDNA-directed therapeutics.



## Organization

### Main Figures

For code to reproduce the main figures, see [Main Figures](src/main_figures).




