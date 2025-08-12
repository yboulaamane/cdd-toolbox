# Awesome Drug Discovery [![Awesome](https://awesome.re/badge.svg)](https://awesome.re)

> Drug discovery is the process by which new candidate medications are identified, designed, and developed using experimental, computational, and informational techniques to address complex challenges in biology, chemistry, and medicine. — [Wikipedia](https://en.wikipedia.org/wiki/Drug_discovery)

---

## Contents
- [Databases and Chemical Libraries](#databases-and-chemical-libraries)  
  - [General Compound Libraries](#general-compound-libraries)  
  - [Natural Product Libraries](#natural-product-libraries)  
  - [Bioactivity Databases](#bioactivity-databases)  
- [Target and Protein Data](#target-and-protein-data)  
  - [Protein Structures](#protein-structures)  
  - [Binding Site and Pocket Detection](#binding-site-and-pocket-detection)  
  - [Protein Engineering and Modeling](#protein-engineering-and-modeling)  
- [Ligand Design and Optimization](#ligand-design-and-optimization)  
  - [Pharmacophore Modeling](#pharmacophore-modeling)  
  - [QSAR and Descriptor Tools](#qsar-and-descriptor-tools)  
  - [Descriptor and Featurization Tools](#descriptor-and-featurization-tools)
  - [Molecular Property Prediction](#molecular-property-prediction)  
  - [Fragment-Based Drug Design](#fragment-based-drug-design)  
- [Virtual Screening and Docking](#virtual-screening-and-docking)  
- [Interaction Analysis and Visualization](#interaction-analysis-and-visualization)  
- [Molecular Dynamics and Simulation](#molecular-dynamics-and-simulation)  
  - [Engines](#engines)  
  - [Topology and Force Field Tools](#topology-and-force-field-tools)  
  - [Analysis Tools](#analysis-tools)  
- [Synthesis and Retrosynthesis Planning](#synthesis-and-retrosynthesis-planning)  
- [Specialized Modalities](#specialized-modalities)  
  - [PROTACs and Ternary Complexes](#protacs-and-ternary-complexes)  
  - [Peptide Design](#peptide-design)  
- [Machine Learning and AI](#machine-learning-and-ai)  
  - [Core Libraries](#core-libraries)  
  - [Chemistry-focused ML Frameworks](#chemistry-focused-ml-frameworks)  
  - [Pretrained Models](#pretrained-models)  
  - [AutoML and Optimization](#automl-and-optimization)  
  - [Molecule Standardization](#molecule-standardization)  
- [Utility and Workflow Tools](#utility-and-workflow-tools)  
- [Learning Resources](#learning-resources)  
  - [Free Courses](#free-courses)  
  - [Blogs](#blogs)  
  - [Instructional Notebooks](#instructional-notebooks)  
- [License](#license)  

---

## Databases and Chemical Libraries

### General Compound Libraries
- [DrugBank](https://go.drugbank.com/) - Comprehensive data on approved and investigational drugs.
- [ZINC](https://zinc.docking.org/) - Free compounds for screening.  
- [ChemSpider](http://www.chemspider.com/) - Chemical structures and data.  
- [DrugSpaceX](https://drugspacex.simm.ac.cn/) - Chemical and biological spaces.  
- [Mcule](https://mcule.com/) - Virtual screening platform with purchasable compounds.  
- [Otava Chemicals](https://www.otavachemicals.com/) - Screening compounds and building blocks.  
- [Vitas-M Laboratory](https://vitasmlab.biz/) - Chemical libraries for HTS and lead discovery.  

### Natural Product Libraries
- [ZINC15 Natural Products](https://zinc15.docking.org/substances/subsets/natural-products/) - 200k+ natural compounds.  
- [COCONUT](https://coconut.naturalproducts.net/) - 400k+ natural products.  
- [LOTUS](https://lotus.naturalproducts.net/) - Annotated molecular data with sourcing organisms.  
- [NPASS](http://bidd.group/NPASS/index.php) - 94k activity-species links.  
- [AfroDB](http://african-compounds.org/about/afrodb/) - 4k+ African medicinal plant compounds.  
- [CMNPD](https://www.cmnpd.org/) - 31k+ marine natural products.  
- [SistematX](https://sistematx.ufpb.br/) - 8k+ secondary metabolites.  
- [Eximed](https://eximedlab.com/Screening-Compounds.html) - 5k+ natural product-like screening compounds.  
- [CoumarinDB](https://yboulaamane.github.io/CoumarinDB/) - 900 coumarins.  
- [ArtemisiaDB](https://yboulaamane.github.io/ArtemisiaDB/) - Artemisia genus compounds.  
- [OTAVA NP-like Library](https://otavachemicals.com/products/compound-libraries-for-hts/natural-product-like-library) - 1k+ NP-like compounds.  
- [BIAdb](https://webs.iiitd.edu.in/raghava/biadb/type.php?tp=natural) - Bioactive peptides and proteins.  
- [IMPPAT](https://cb.imsc.res.in/imppat/home) - Phytochemicals from Indian medicinal plants.  
- [NP-MRD](https://np-mrd.org/natural_products) - 280k+ NMR-based NP studies.  
- [IBS Natural Compounds](https://www.ibscreen.com/natural-compounds) - 60k+ compounds.  
- [Phytochemicals](https://www.phytochemicals.info/) - Comprehensive phytochemical info.  
- [NPACT](https://webs.iiitd.edu.in/raghava/npact/index.html) - Plant-based anticancer compounds.  
- [NaturAr](https://naturar.quimica.unlp.edu.ar/en/) - Argentine biodiversity compounds.  
- [DiaNat-DB](http://rdu.iquimica.unam.mx/handle/20.500.12214/1186) - Antidiabetic plant compounds.  
- [PhytoHub](https://phytohub.eu/) - Dietary phytochemicals and metabolites.  
- [Dr. Duke's Phytochemical DB](https://phytochem.nal.usda.gov/) - Plant compounds and uses.  
- [CyanoMetDB](https://zenodo.org/records/13854577) - Cyanobacterial metabolites.  
- [Seaweed Metabolite DB](https://www.swmd.co.in/) - Marine algae compounds.  
- [Arabidopsis.org](https://www.arabidopsis.org/) - Arabidopsis molecular biology.  

### Bioactivity Databases
- [ChEMBL](https://www.ebi.ac.uk/chembl/) - Bioactivity and ADMET data.  
- [SureChEMBL](https://www.surechembl.org/) - Patent chemistry search.  
- [BindingDB](https://www.bindingdb.org/) - Binding affinities for biomolecules.  
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) - Structures, properties, and bioassays.  
- [PDBbind](http://www.pdbbind.org.cn/index.php) - Protein-ligand affinity data.  
- [BRENDA](https://www.brenda-enzymes.org/) - Enzyme properties and functions.  
- [ExCAPE-DB](https://solr.ideaconsult.net/search/excape/) - Chemogenomics DB.  
- [Therapeutics Data Commons](https://tdcommons.ai/) - AI datasets for therapeutics.
- [Therapeutic Target Database (TTD)](https://idrblab.net/ttd/) - Drug targets with linked diseases and compounds.

---

## Target and Protein Data

### Protein Structures
- [RCSB PDB](https://www.rcsb.org/) - Repository for macromolecular structures.
- [PDBe](https://www.ebi.ac.uk/pdbe/) - European counterpart to RCSB PDB.
- [OPM](https://opm.phar.umich.edu/) - Orientation of proteins in membranes.
- [UniProt](https://www.uniprot.org/) - Protein sequences, structures, and functions.
- [InterPro](https://www.ebi.ac.uk/interpro/) - Protein classification and domain prediction.
- [AlphaFold DB](https://alphafold.ebi.ac.uk/) - Predicted structures from AlphaFold.
- [Proteopedia](https://proteopedia.org/wiki/index.php/Main_Page) - Interactive protein visualizations.

### Binding Site and Pocket Detection
- [ProteinsPlus](https://proteins.plus/) - Binding site analysis.
- [PrankWeb](https://prankweb.cz/) - Pocket prediction and analysis.
- [CASTp](http://sts.bioe.uic.edu/castp/index.html?2r7g) - Pocket geometry and volume analysis.
- [CavityPlus](http://www.pkumdl.cn:8000/cavityplus/index.php#/) - Pocket detection and druggability.
- [CaverWeb](https://loschmidt.chemi.muni.cz/caverweb/) - Tunnel and channel detection.
- [PASSer](https://passer.smu.edu/) - Allosteric site prediction.

### Protein Engineering and Modeling
- [DynaMut](https://biosig.lab.uq.edu.au/dynamut/) - Predicts mutation-induced stability changes.

---

## Ligand Design and Optimization

### Pharmacophore Modeling
- [ZINCPharmer](http://zincpharmer.csb.pitt.edu/) - Pharmacophore screening.
- [Pharmit](https://pharmit.csb.pitt.edu/) - Interactive pharmacophore modeling.
- [PharmMapper](https://www.lilab-ecust.cn/pharmmapper/) - Pharmacophore mapping.

### QSAR and Descriptor Tools
- [QSAR Toolbox](https://qsartoolbox.org/) - Hazard assessment and QSAR.
- [OCHEM](https://ochem.eu/home/show.do) - QSAR model building and prediction.
- [ChemMaster](https://crescent-silico.com/chemmaster/) - QSAR and cheminformatics suite.
- [3D-QSAR](https://www.3d-qsar.com/) - 3D-QSAR modeling resources.
- [QSAR-Co](https://sites.google.com/view/qsar-co/) - Robust multitarget QSAR modeling.
- [DataWarrior](https://openmolecules.org/datawarrior/) - Free software for chemical analysis, QSAR, and visualization.
- [KNIME](https://www.knime.com/) - Workflow platform for cheminformatics and ML integration.
  
### Descriptor and Featurization Tools
- [RDKit](https://www.rdkit.org/) - Open-source cheminformatics toolkit with descriptor, fingerprint, and molecular manipulation support.  
- [PaDEL-Descriptor](http://www.yapcwsoft.com/dd/padeldescriptor/) - Java tool for calculating molecular descriptors and fingerprints.  
- [Mordred](https://github.com/mordred-descriptor/mordred) - Python library with 1800+ molecular descriptors.  
- [CDK](https://cdk.github.io/) - Java cheminformatics library with descriptor calculators.  
- [alvaDesc](https://www.alvascience.com/alvadesc/) - Commercial software for molecular descriptors and fingerprints.  
- [MolFeat](https://molfeat.datamol.io/) - Python package for molecular featurization and embeddings.
- [DeepChem Featurizers](https://deepchem.io/) - Featurization for graph neural networks, CNNs, and sequence models.  
- [Dragon](https://www.talete.mi.it/products/dragon_description.htm) - Commercial molecular descriptor calculator (widely cited).  

### Molecular Property Prediction
- [SwissADME](http://www.swissadme.ch/) - Drug-likeness and PK.
- [pkCSM](https://biosig.lab.uq.edu.au/pkcsm/) - ADMET property prediction.
- [DeepPK](https://biosig.lab.uq.edu.au/deeppk/) - DL-based pharmacokinetics.
- [admetSAR 2.0](https://lmmd.ecust.edu.cn/admetsar2/) - Comprehensive ADMET.
- [ADMETlab 2.0](https://admetmesh.scbdd.com/) - PK, toxicity and drug-likeness.
- [ProTox-II](https://tox-new.charite.de/protox_II/) - Toxicity predictions.
- [PreADMET](https://preadmet.webservice.bmdrc.org/) - PK property predictions.
- [FAF-Drugs](https://bioserv.rpbs.univ-paris-diderot.fr/services.html) - ADMET filtering.
- [Admetboost](https://ai-druglab.smu.edu/admet) - ML-based ADMET prediction.

### Fragment-Based Drug Design

- [SwissSidechain](https://www.swisssidechain.ch/) - Fragment and linker library for small molecule design.  
- [FragBuilder](https://github.com/andersx/fragbuilder) - Python API for building peptide-like and small molecule fragments.  
- [SeeSAR](https://www.biosolveit.de/SeeSAR/) - Fragment growing and linking software (free academic version).
- [Enamine Fragment Libraries](https://enamine.net/compound-libraries/fragment-libraries) - Large curated collection of diverse fragments for FBDD.
  
---

## Virtual Screening and Docking
- [OpenBabel](https://openbabel.org/index.html) - Format conversion and ligand prep.
- [MGLTools](https://ccsb.scripps.edu/mgltools/) - Structure preparation.
- [AutoDockTools](https://autodocksuite.scripps.edu/adt/) - AutoDock GUI.
- [AutoDock Vina](https://vina.scripps.edu/) - Popular docking software.
- [EasyDockVina2](https://github.com/S3cr3t-SDN/EasyDockVina2) - Vina automation.
- [Webina](https://durrantlab.pitt.edu/webina/) - Web-based Vina.
- [Smina](https://github.com/mwojcikowski/smina) - Vina fork with extra features.
- [Gnina](https://github.com/gnina/gnina) - CNN-scoring docking.
- [EasyDock](https://github.com/ci-lab-cz/easydock) - Vina/Smina pipeline.
- [HADDOCK](https://wenmr.science.uu.nl/haddock2.4/) - Flexible docking suite.
- [PandaDock](https://github.com/pritampanda15/PandaDock) - Python docking tool.
- [ZDOCK](https://zdock.wenglab.org/) - Protein-protein docking.
- [ClusPro](https://cluspro.org/) - Protein-protein docking server.
- [pyDockWEB](https://life.bsc.es/pid/pydockweb/) - Electrostatics-based docking.
- [SwissDock](https://www.swissdock.ch/) - Web docking for beginners.
- [MzDOCK](https://github.com/Muzatheking12/MzDOCK) - GUI docking pipeline.
- [Uni-Mol Docking V2](https://www.bohrium.com/apps/unimoldockingv2/job?type=app) - AI-assisted docking.
- [Vina on Colab](https://autodock-vina.readthedocs.io/en/latest/colab_examples.html) - Run Vina in Google Colab.

---

## Interaction Analysis and Visualization
- [PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) - Protein-ligand interaction profiling.
- [LigPlot+](https://www.ebi.ac.uk/thornton-srv/software/LigPlus/) - 2D interaction diagrams.
- [Discovery Studio Visualizer](https://discover.3ds.com/discovery-studio-visualizer-download) - Advanced visualization.

---

## Molecular Dynamics and Simulation

### Engines
- [GROMACS](https://www.gromacs.org/) - Scalable MD engine.
- [LAMMPS](https://www.lammps.org/) - Parallel MD for materials and biomolecules.
- [NAMD](https://www.ks.uiuc.edu/Research/namd/) - High-performance biomolecular MD.
- [AMBER](https://ambermd.org/) - Suite for biomolecular simulations.
- [Desmond](https://www.deshawresearch.com/resources.html) - High-performance MD.

### Topology and Force Field Tools
- [CGenFF](https://cgenff.umaryland.edu/) - CHARMM force field parametrization.
- [SwissParam](https://www.swissparam.ch/) - Small molecule parameters.
- [ATB](https://atb.uq.edu.au/) - Automated topology builder.
- [CHARMM-GUI](https://www.charmm-gui.org/) - Input and topology preparation.
- [LigParGen](https://zarbi.chem.yale.edu/ligpargen/) - Ligand FF parameters.

### Analysis Tools
- [MD DaVis](https://md-davis.readthedocs.io/en/latest/index.html) - Interactive MD visualizations.
- [iMod](https://imods.iqfr.csic.es/) - Normal Mode Analysis.
- [MolAiCal](https://molaical.github.io/) - Binding free energy calculations.

---

## Synthesis and Retrosynthesis Planning
- [Spaya](https://spaya.ai/app/search) - AI retrosynthesis.
- [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder) - Monte Carlo retrosynthesis.
- [ASKCOS](https://askcos.mit.edu/) - Reaction planning with ML.
- [IBM RoboRXN](https://rxn.res.ibm.com/rxn/robo-rxn/welcome) - Automated reaction prediction.
- [MANIFOLD](https://app.postera.ai/manifold/) - Synthetic accessibility and search.

---

## Specialized Modalities

### PROTACs and Ternary Complexes
- [PROTAC-db](http://cadd.zju.edu.cn/protacdb/) - PROTAC data.
- [PROsettaC](https://prosettac.weizmann.ac.il/) - Ternary complex modeling.

### Peptide Design
- [PepDraw](https://pepdraw.com/) - Peptide visualization.
- [PepSite](http://pepsite2.russelllab.org/) - Predict peptide binding sites.
- [Peptimap](https://peptimap.bu.edu/) - Peptide mapping.
---

## Machine Learning and AI

### Core Libraries
- [scikit-learn](https://scikit-learn.org/) - General-purpose machine learning library for classification, regression, clustering, and model evaluation.
- [scikit-mol](https://github.com/datamol-io/scikit-mol) - Scikit-learn compatible cheminformatics extensions for molecular machine learning.
- [pytorch-geometric](https://pytorch-geometric.readthedocs.io/) - PyTorch-based library for graph neural networks (GNNs).
- [DGL-LifeSci](https://github.com/awslabs/dgl-lifesci) - Graph deep learning toolkit for life sciences using the Deep Graph Library.

### Chemistry-focused ML Frameworks
- [DeepChem](https://github.com/deepchem/deepchem) - Open-source deep learning framework for chemistry and biology.
- [Chemprop](https://github.com/chemprop/chemprop) - Directed message passing neural networks for molecular property prediction.
- [ChemML](https://github.com/hachmannlab/chemml) - Machine learning and informatics suite for analyzing, mining, and modeling chemical and materials data.
- [Oloren ChemEngine](https://github.com/Oloren-AI/olorenchemengine) - Unified API for molecular property prediction with uncertainty quantification, interpretability, and model tuning.
- [TorchDrug](https://torchdrug.ai/) - A machine learning library for drug discovery with support for GNNs and molecular datasets.

### Pretrained Models
- [MolBERT](https://github.com/BenevolentAI/MolBERT) - Transformer-based molecular representation learning.
- [ChemBERTa](https://huggingface.co/seyonec/ChemBERTa-zinc-base-v1) - Pretrained BERT-like models for molecules from SMILES.
- [Uni-Mol](https://github.com/dptech-corp/Uni-Mol) - 3D molecular representation learning framework.

### AutoML and Optimization
- [Auto-sklearn](https://automl.github.io/auto-sklearn/master/) - Automated machine learning for scikit-learn.
- [TPOT](https://epistasislab.github.io/tpot/) - Genetic programming-based AutoML for optimizing ML pipelines.
- [Optuna](https://optuna.org/) - Hyperparameter optimization framework for machine learning.

### Molecule Standardization
- [MolVS](https://github.com/mcs07/MolVS) - Molecule validation and standardization library based on RDKit.

---

## Utility and Workflow Tools
- [OPSIN](https://opsin.ch.cam.ac.uk) - IUPAC name to structure.
- [OSRA](https://cactus.nci.nih.gov/cgi-bin/osra/index.cgi) - Image to structure.
- [MetaPredict](http://metapredict.icoa.fr/) - Molecular property prediction.
- [ChemPlot](https://chemplot.streamlit.app/) - Chemical space visualization.
- [ChemDB](http://cdb.ics.uci.edu/) - Chemoinformatics portal.
- [BoBER](http://bober.insilab.org/) - Bioisosteric replacements.
- [Open Targets](https://platform.opentargets.org/) - Target identification.
- [Screening Explorer](http://stats.drugdesign.fr/) - Screening data analysis.
- [LigRMSD](https://ligrmsd.appsbio.utalca.cl/) - Ligand RMSD calculation.
- [NERDD](https://nerdd.univie.ac.at/) - Drug discovery resources.
- [MetaChemiBio](https://biochemia.uwm.edu.pl/metachemibio/) - Property prediction.
- [WenMR Portal](https://wenmr.science.uu.nl/) - Biomolecular interactions software.
- [LigBuilder3](http://www.pkumdl.cn:8080/ligbuilder3/) - Ligand design.
- [ChemMine Tools](https://chemminetools.ucr.edu/) - Cheminformatics tools.
- [MayaChemTools](http://www.mayachemtools.org/index.html) - Perl/Python scripts for cheminformatics.
- [SCBDD](http://www.scbdd.com/) - Cheminformatics and drug discovery software.
- [Click2Drug](https://www.click2drug.org/) - CADD software and DB directory.
- [Galaxy Europe](https://usegalaxy-eu.github.io/index-cheminformatics.html) - Galaxy instance for cheminformatics.
- [CADD Vault](https://drugbud-suite.github.io/CADD_Vault/) - CADD resources repository.
- [BioMoDes](https://abeebyekeen.com/biomodes-biomolecular-structure-prediction/) - Biomolecular modeling tools.
- [PlayMolecule](https://open.playmolecule.org/landing) - Molecular modeling simulations.
- [Venny 2.1](https://bioinfogp.cnb.csic.es/tools/venny/) - Venn diagram tool.

---

## Learning Resources

### Free Courses
- [TMP Chem Lectures](https://youtube.com/playlist?list=PLm8ZSArAXicIWTHEWgHG5mDr8YbrdcN1K) - Computational chemistry lectures.
- [Strasbourg Summer School in Chemoinformatics](https://youtube.com/playlist?list=PLhgURFExPmJsDuHevu5n8y0R41WsXfbnC) - Summer school lectures.
- [BIGCHEM](https://bigchem.eu/node/63) - Big data in chemistry course.
- [Geometric Deep Learning Course](https://geometricdeeplearning.com/lectures/) - GDL for ML in science.
- [Drug Discovery Course](https://www.stereoelectronics.org/webDD/DD_home.html) - Drug discovery fundamentals.
- [drugdesign.org](https://www.drugdesign.org/) - Drug design and cheminformatics courses.
- [Cheminformatics OLCC](https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics) - Cheminformatics theory and coding.
- [Python For Cheminformatics Docking](https://pdb101.rcsb.org/train/training-events/python4) - Python-based docking.
- [DDA CDD Workshop](https://wcair.dundee.ac.uk/training/training-resources/computational-drug-design/) - Generative drug design.

### Blogs
- [Practical Fragments](http://practicalfragments.blogspot.com/) - Fragment-based drug design.
- [avrilomics](https://avrilomics.blogspot.com/) - Genomics and bioinformatics.
- [Practical Cheminformatics](http://practicalcheminformatics.blogspot.com/) - Cheminformatics tools.
- [Cheminformania](https://www.cheminformania.com/) - Cheminformatics and deep learning.
- [Daily Dose of Data Science](https://www.blog.dailydoseofds.com/) - Data science insights.
- [Machine Learning Mastery](https://machinelearningmastery.com/) - ML tutorials.
- [Chem-Workflows](https://chem-workflows.com/index.html) - Jupyter chemistry tutorials.
- [Structural Bioinformatics](https://proteinstructures.com/) - Structure-based drug design guide.
- [Bioinformatics Answers](https://www.biostars.org/) - Bioinformatics Q&A.
- [McConnellsMedChem](https://mcconnellsmedchem.com/) - Medicinal chemistry blog.
- [DrugDiscovery.NET](http://www.drugdiscovery.net/) - AI in drug discovery.
- [MacinChem](https://macinchem.org/) - Comp chem on macOS.
- [Jeremy Monat](https://bertiewooster.github.io/) - Cheminformatics research.
- [Angelo Raymond Rossi](https://angeloraymondrossi.github.io/) - Computational chemistry research.

### Instructional Notebooks
- [TeachOpenCADD](https://projects.volkamerlab.org/teachopencadd/all_talktorials.html) - Jupyter tutorials for CADD.
- [intro_pharma_ai](https://github.com/kochgroup/intro_pharma_ai) - AI in pharma with notebooks.
- [AI/DL for Life Sciences](https://onlinelibrary.wiley.com/doi/10.1002/ardp.202200628) - Interactive AI/DL notebooks.

---

## Licensing
[![CC0](http://mirrors.creativecommons.org/presskit/buttons/88x31/svg/cc-zero.svg)](https://creativecommons.org/publicdomain/zero/1.0/)
