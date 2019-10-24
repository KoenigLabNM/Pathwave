# Pathwave

## About PathWave

PathWave enables the identification of disease specific regulation patterns by combining gene expression data and network topology. PathWave is preferable to classical enrichment tests as it offers a much higher sensitivity for detecting such functionally related regulation patterns. Additionally, it is more precise to a permutation based enrichment method (Schramm et al 2010).

## Downloads

### PathWave version 2.1 (Piro et al. 2014):

The R package does already include all necessary data for BiGG/Human recon 1 (H. sapiens) and KEGG (H. sapiens, M. musculus, D. melanogaster, C. elegans, D. rerio, E. coli)

    R package PathWave v2.1.3 (Linux)
    R package PathWave v2.1.3 (Windows)
    PathWave v2.1 - User Manual
    PathWave v2.1 - Usage Example (Quick Guide)

### Additional software for PathWave version 2.1:

The following software is needed only to preprocess further pathway data sets from KEGG (KGML) or BiGG (SBML), not for applying PathWave to already available, preprocessed pathway data. Note: GridArranger is also included as external code in the PathWave 2.1 R package. We provide it here as an independent package because it may be useful for other purposes.

    GridArranger v1.0 - used to arrange pathways into compact 2D lattice grids (see Schramm et al 2010; Oswald et al., 2012; Piro et al 2014).
    GridArranger v1.0 README - basic information on installation and usage.
    ABACUS v2.4-alpha - required by GridArranger; this package is provided with kind permission of the original authors at the University of Cologne, Germany, because newer versions of ABACUS (available at the ABACUS webiste) are not compatible with GridArranger.

### PathWave version 1.0 (Schramm et al. 2010):

    Short manual
    Executable R file usePathWave.R
    R package PathWave v1.0 (Linux)
    R package PathWave v1.0 (Windows)
    PathWaveFilesFeb04_2009 (Linux)
    PathWaveFilesFeb04_2009 (Windows)
    Optimized grids for human (.RData file)

## Publications

### Network topology-based detection of differential gene regulation and regulatory switches in cell metabolism and signaling

BMC Systems Biology8:56, 2014.

Rosario M. Piro (1,2,#), Stefan Wiesberg (3), Gunnar Schramm (1,2), Nico Rebel (3), Marcus Oswald (1,4,5), Roland Eils (1,2), Gerhard Reinelt (3), and Rainer König (1,2,4,5*)

(1) Division of Theoretical Bioinformatics, German Cancer Research Center (Deutsches Krebsforschungszentrum, DKFZ), Heidelberg, Germany
(2) Department of Bioinformatics and Functional Genomics, Institute of Pharmacy and Molecular Biotechnology, BioQuant, University of Heidelberg, Germany
(3) Institute of Computer Science and Interdisciplinary Center for Scientific Computing, University of Heidelberg, Germany
(4) Center for Sepsis Control and Care, University Hospital Jena, Jena, Germany
(5) Hans-Knöll-Institute (HKI), Jena, Germany
(#) Current address: Dep. of Mathematics and Computer Science, Institute of Computer Science and Institute of Bioinformatics, Freie Universität Berlin, Germany
(*) corresponding author: rainer.koenig@uni-jena.de

-----

### PathWave: Discovering patterns of differentially regulated enzymes in metabolic pathways

Bioinformatics26(9):1225-1231, 2010.

Gunnar Schramm (1,2), Stefan Wiesberg (3,4), Nicolle Diessl (1), Anna-Lena Kranz (1,2), Vitalia Sagulenko (5), Marcus Oswald (3,4), Gerhard Reinelt (3), Frank Westermann (5), Roland Eils (1,2, *) and Rainer König (1,2, *)

(1) Department of Bioinformatics and Functional Genomics, Institute of Pharmacy and Molecular Biotechnology, and Bioquant, Im Neuenheimer Feld 267, 69120 Heidelberg, Germany
(2) Department of Theoretical Bioinformatics, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 280, 69120 Heidelberg, Germany
(3) Institute of Computer Science, University of Heidelberg, 69120 Heidelberg, Germany
(4) Interdisciplinary Center for Scientific Computing, University of Heidelberg, 69120 Heidelberg, Germany
(5) Department of Tumor Genetics, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 280, 69120 Heidelberg, Germany
(*) corresponding authors: r.eils@dkfz.de and rainer.koenig@uni-jena.de

-----

### Exact solution of the 2-dimensional grid arrangement problem

Discrete Optimization9:189-199, 2012.

Marcus Oswald (1,2), Gerhard Reinelt (1), Stefan Wiesberg (1,2)

(1) Institute of Computer Science, University of Heidelberg, 69120 Heidelberg, Germany
(2) Interdisciplinary Center for Scientific Computing, University of Heidelberg, 69120 Heidelberg, Germany

-----
Please cite: Piro et al (2014) for PathWave version 2.1, Schramm et al (2010) for PathWave version 1.0, and both papers (Schramm et al 2010; Piro et al 2014) for PathWave in general.
