# neural_networks_in_biology

### *Neural Networks in Biology*

[ABodyBuilder2](https://github.com/oxpig/ImmuneBuilder/blob/0df4e2ad82a1aa60f37ea9dae335d1198159ef78/ImmuneBuilder/ABodyBuilder2.py#L86): predict the structure of antibodies

[AlphaFold3](https://github.com/google-deepmind/alphafold3): AF2 + ligands, DNA, RNA, and post-translational modifications 

[ATOMICA](https://github.com/mims-harvard/ATOMICA): universal representations of intermolecular interactions including protein-small molecule, protein-ion, small molecule-small molecule, protein-protein, protein-peptide, protein-RNA, protein-DNA, and nucleic acid-small molecule complexes

[BE-DICT](https://github.com/uzh-dqbm-cmi/crispr): predicts outcomes of CRISPR base editing 

[CellOT](https://github.com/bunnech/cellot?tab=readme-ov-file): neural optimal transport model predicting how cells transition across time, used in development and disease

[Chroma](https://github.com/generatebio/chroma): vector DB + retriever, use for integrating + querying protein design history, experimental metadata, model outputs 

[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/): tells you if a mutation causes a disease or not

[Cortex](https://github.com/prescient-design/cortex): modular architecture for deep learning systems, a general purpose scoring oracle for antibody optimization pipelines

[DBNN](https://www.researchgate.net/publication/344958092_A_Dynamical_Biomolecular_Neural_Network): biological molecular circuits in neural form instead of discrete, boolean gates; used for testing before wet lab verification  

[DeepAccNet](https://github.com/hiranumn/DeepAccNet): protein model accuracy evaluator 

[DeepAIR](https://github.com/TencentAILabHealthcare/DeepAIR): predicts TCR/BCR-antigen binding affinity and reactivity

[DeepRC](https://github.com/ml-jku/DeepRC): Hopfield network with attention that classifies disease 

[DeepVirFinder](https://github.com/jessieren/DeepVirFinder): identify viruses based on meta-genomic information 

[DeepZF](https://github.com/OrensteinLab/DeepZF): predicts ZF-DNA binding

[DiffDock](https://github.com/gcorso/DiffDock): diffusion model for ligand-protein docking, outperforms traditional scoring in benchmarks 

[dwJS](https://github.com/Genentech/walk-jump): walk-jump sampling, which combines score-based + energy-based models for antibody generation, used for exploring novel sequence variants far from known binders

[DyAb](https://www.biorxiv.org/content/10.1101/2025.01.28.635353v1): uses pair-wise representation to predict differences in protein properties, rather than absolute values; use when prioritizing/generating antibody mutations for property shifts 

[Enformer](https://www.notion.so/Neural-Networks-1db5b5e4de3180fe97ced9dec4c6012a?pvs=21): gene expression prediction over long range DNA sequences

[ESM models](https://github.com/facebookresearch/esm): protein folding models

[EvoDiff](https://github.com/microsoft/evodiff): generative protein design, $10^{30}$ drug-like molecules possible

[Evo2](https://github.com/ArcInstitute/evo2): predicts mutations, function, fitness, structure across the OpenGenome2 dataset

[FoldX](https://foldxsuite.crg.eu/documentation#manual): predicts how mutations affect a protein’s stability and interactions

[GEARS](https://github.com/snap-stanford/GEARS): predicts how gene expression will change when specific genes are edited (knocked out, over-expressed) (Shift Biosciences)

[GNINA](https://github.com/gnina/gnina): deep learning framework for molecular docking 

[GPN-MSA](https://www.biorxiv.org/content/10.1101/2023.10.10.561776v1.full): 3D models of RNA that predict mutations that affect RNA editing

[GraphBAN](https://www.nature.com/articles/s41467-025-57536-9): predict compound-protein interactions 

[ESM3](https://www.evolutionaryscale.ai/blog/esm3-release): protein foundation model, trained on 500M years of evolution

[mosGraphGPT](https://pmc.ncbi.nlm.nih.gov/articles/PMC11326168/): foundation model trained on multi-omic data to predict cell states and disease outcomes across species (transformer-based attention)

[LaMBO-2](https://github.com/ngruver/NOS?tab=readme-ov-file): multi-objective optimization in sequence space that is categorized, used when optimizing antibodies for analyzing multiple properties at once, like affinity and expression 

[LBSTER](https://github.com/prescient-design/lobster): base model you can use for fine-tuning or embedding antibodies, is learned representations from UniRef50 + antibody space 

[LigandMPNN](https://github.com/dauparas/LigandMPNN): deep learning model that allows explicit modeling of small molecule, nucleotide, metal, and other atomic contexts

[LungTCR](https://github.com/OpenGene/LungTCR): lung cancer prediction based on T-cell receptor data

[MMSeq2](https://github.com/soedinglab/MMseqs2): sequence clustering/deduplication, reducing redundancy in training sets; use when cleaning data before model training

[NEHVI](https://arxiv.org/abs/2105.08195): acquisition function that ranks candidates via expected Pareto improvement; use when selecting mathematically unique sequences under uncertainty 

[PanGenie](https://github.com/eblerjana/pangenie): Haplotype-aware genotyping using population graphs

[PepTune](https://huggingface.co/ChatterjeeLab/PepTune): de novo generation of new therapeutic peptides by gradually improving random sequences that have multiple design goals- like hitting a target, lasting in the body, and being easy to make

[ProtGPT2](https://huggingface.co/nferruz/ProtGPT2): de novo protein design, autoregressive transformer trained on UniRef50 for protein generation

[ProteinMPNN](https://github.com/dauparas/ProteinMPNN)**:** a sequence design neural network that takes a fixed 3D backbone structure ****and predicts amino acid sequences likely to fold into it, use it after generating a novel 3D structure, but before experimental expression

[PropEn](https://github.com/prescient-design/propen): transforms low-affinity to high-affinity sequences, use to find an average across strong and weak binders

[RFDiffusion2](https://www.nature.com/articles/s41586-023-06415-8): takes in atomic coordinates of scaffolding residues

[scGPT](https://github.com/bowang-lab/scGPT): identifies and classifies cells based on input RNA sequences (Shift Biosciences)

[scVI](https://github.com/scverse/scvi-tools): variational auto-encoder for analyzing single-cell RNA-seq, used in single-cell workflows 

[SeqVDM](https://arxiv.org/abs/2107.00630): variational diffusion model adapted for protein sequences, used for continuous control over diversity and stability of designs

[soNNia](https://github.com/statbiophys/soNNia): NN that predicts sequence-based TCR/BCR binding reactivity 

[StripedHyena](https://github.com/togethercomputer/stripedhyena): attention-based models for long sequences; use when modeling long protein sequences/genomes with limited compute 

[TCRAI](https://www.science.org/doi/10.1126/sciadv.abf5835): CNN + MIL model that classifies TCR-antigen binding and immune repertoire 

[Viral Mutation](https://github.com/brianhie/viral-mutation): language model for virus evolution 

### *Tools*

[AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina): ligand-protein docking and screening model (physics-based, not deep learning)

[DNAworks](https://github.com/davidhoover/DNAWorks): automatic oligonucleotide design for PCR-based gene synthesis 

[OpenGenome2](https://huggingface.co/datasets/arcinstitute/opengenome2): the 8.8T genome dataset Arc used to train Evo2

[geNomad](https://gnomad.broadinstitute.org/): database of human genetic variation from diverse populations- 807K genomes, 3+ petabytes raw, 35 TB variant summaries

[RifDock](https://github.com/rifdock/rifdock): rigid-body docking tool for initial binder-target orientations, assumes proteins don’t flex 

[Rosetta](https://rosettacommons.org/software/): protein structure energy scoring (uses Monte Carlo, not NN’s)

[Rosetta Fast Relax](https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/FastRelaxMover): structure refinement tool used to iteratively improve designs

[SELFIES](https://github.com/aspuru-guzik-group/selfies): generation of molecular graphs which are syntactically and semantically valid

[Savanna](https://github.com/Zymrael/savanna): pre-training infrastructure for multi-hybrid AI model architectures, like StripedHyena2
