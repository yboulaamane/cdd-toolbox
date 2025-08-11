# Awesome Drug Discovery [![Awesome](https://awesome.re/badge.svg)](https://awesome.re) [![License: MIT](https://img.shields.io/badge/License-MIT-red.svg)](LICENSE)

> Drug discovery is the process by which new candidate medications are identified, designed, and developed, using experimental, computational, and informational techniques to address complex challenges in biology, chemistry, and medicine. — [Wikipedia](https://en.wikipedia.org/wiki/Drug_discovery)

---

## Table of Contents
1. [Databases & Chemical Libraries](#1-databases--chemical-libraries)  
   - [1.1 General Compound Libraries](#11-general-compound-libraries)  
   - [1.2 Natural Product Libraries](#12-natural-product-libraries)  
   - [1.3 Bioactivity Databases](#13-bioactivity-databases)  
2. [Target & Protein Data](#2-target--protein-data)  
   - [2.1 Protein Structures](#21-protein-structures)  
   - [2.2 Binding Site & Pocket Detection](#22-binding-site--pocket-detection)  
   - [2.3 Protein Engineering & Modeling](#23-protein-engineering--modeling)  
3. [Ligand Design & Optimization](#3-ligand-design--optimization)  
   - [3.1 Pharmacophore Modeling](#31-pharmacophore-modeling)  
   - [3.2 QSAR & Descriptor Tools](#32-qsar--descriptor-tools)  
   - [3.3 Molecular Property Prediction](#33-molecular-property-prediction)
   - [3.4 Fragment-Based Drug Design](#34-fragment-based-drug-design)
4. [Virtual Screening & Docking](#4-virtual-screening--docking)  
5. [Interaction Analysis & Visualization](#5-interaction-analysis--visualization)  
6. [Molecular Dynamics & Simulation](#6-molecular-dynamics--simulation)  
   - [6.1 Engines](#61-engines)  
   - [6.2 Topology & Force Field Tools](#62-topology--force-field-tools)  
   - [6.3 Analysis Tools](#63-analysis-tools)  
7. [Synthesis & Retrosynthesis Planning](#7-synthesis--retrosynthesis-planning)  
8. [Specialized Modalities](#8-specialized-modalities)  
   - [8.1 PROTACs & Ternary Complexes](#81-protacs--ternary-complexes)  
   - [8.2 Peptide Design](#82-peptide-design)  
9. [Machine Learning & AI for Drug Discovery](#9-machine-learning--ai-for-drug-discovery)  
10. [Utility & Workflow Tools](#10-utility--workflow-tools)  
11. [Learning Resources](#11-learning-resources)  
    - [11.1 Free Courses](#111-free-courses)  
    - [11.2 Blogs](#112-blogs)  
    - [11.3 Instructional Notebooks](#113-instructional-notebooks)  

---

## 1. Databases & Chemical Libraries

### 1.1 General Compound Libraries
- [DrugBank](https://go.drugbank.com/) – Comprehensive drug data (approved & investigational).  
- [ZINC](https://zinc.docking.org/) – Free compounds for screening.  
- [ZINC15 Natural Products](https://zinc15.docking.org/substances/subsets/natural-products/) – 200k+ natural compounds.  
- [ChemSpider](http://www.chemspider.com/) – Chemical structures and data.  
- [DrugSpaceX](https://drugspacex.simm.ac.cn/) – Chemical & biological spaces.  
- [Mcule](https://mcule.com/) – Virtual screening platform with purchasable compounds.  
- [Otava Chemicals](https://www.otavachemicals.com/) – Screening compounds & building blocks.  
- [Vitas-M Laboratory](https://vitasmlab.biz/) – Chemical libraries for HTS & lead discovery.  

### 1.2 Natural Product Libraries
- [COCONUT](https://coconut.naturalproducts.net/) – 400k+ natural products.  
- [LOTUS](https://lotus.naturalproducts.net/) – Annotated molecular data with sourcing organisms.  
- [NPASS](http://bidd.group/NPASS/index.php) – 94k activity–species links.  
- [AfroDB](http://african-compounds.org/about/afrodb/) – 4k+ African medicinal plant compounds.  
- [CMNPD](https://www.cmnpd.org/) – 31k+ marine natural products.  
- [SistematX](https://sistematx.ufpb.br/) – 8k+ secondary metabolites.  
- [Eximed](https://eximedlab.com/Screening-Compounds.html) – 5k+ natural product-like screening compounds.  
- [CoumarinDB](https://yboulaamane.github.io/CoumarinDB/) – 900 coumarins.  
- [ArtemisiaDB](https://yboulaamane.github.io/ArtemisiaDB/) – Artemisia genus compounds.  
- [OTAVA NP-like Library](https://otavachemicals.com/products/compound-libraries-for-hts/natural-product-like-library) – 1k+ NP-like compounds.  
- [BIAdb](https://webs.iiitd.edu.in/raghava/biadb/type.php?tp=natural) – Bioactive peptides & proteins.  
- [IMPPAT](https://cb.imsc.res.in/imppat/home) – Phytochemicals from Indian medicinal plants.  
- [NP-MRD](https://np-mrd.org/natural_products) – 280k+ NMR-based NP studies.  
- [IBS Natural Compounds](https://www.ibscreen.com/natural-compounds) – 60k+ compounds.  
- [Phytochemicals](https://www.phytochemicals.info/) – Comprehensive phytochemical info.  
- [NPACT](https://webs.iiitd.edu.in/raghava/npact/index.html) – Plant-based anticancer compounds.  
- [NaturAr](https://naturar.quimica.unlp.edu.ar/en/) – Argentine biodiversity compounds.  
- [DiaNat-DB](http://rdu.iquimica.unam.mx/handle/20.500.12214/1186) – Antidiabetic plant compounds.  
- [PhytoHub](https://phytohub.eu/) – Dietary phytochemicals & metabolites.  
- [Dr. Duke's Phytochemical DB](https://phytochem.nal.usda.gov/) – Plant compounds & uses.  
- [CyanoMetDB](https://zenodo.org/records/13854577) – Cyanobacterial metabolites.  
- [Seaweed Metabolite DB](https://www.swmd.co.in/) – Marine algae compounds.  
- [Arabidopsis.org](https://www.arabidopsis.org/) – Arabidopsis molecular biology.  

### 1.3 Bioactivity Databases
- [ChEMBL](https://www.ebi.ac.uk/chembl/) – Bioactivity & ADMET data.  
- [SureChEMBL](https://www.surechembl.org/) – Patent chemistry search.  
- [BindingDB](https://www.bindingdb.org/) – Binding affinities for biomolecules.  
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) – Structures, properties, and bioassays.  
- [PDBbind](http://www.pdbbind.org.cn/index.php) – Protein–ligand affinity data.  
- [BRENDA](https://www.brenda-enzymes.org/) – Enzyme properties & functions.  
- [ExCAPE-DB](https://solr.ideaconsult.net/search/excape/) – Chemogenomics DB.  
- [Therapeutics Data Commons](https://tdcommons.ai/) – AI datasets for therapeutics.
- [Therapeutic Target Database (TTD)](https://idrblab.net/ttd/) – Drug targets with linked diseases and compounds.

---

## 2. Target & Protein Data

### 2.1 Protein Structures
- [RCSB PDB](https://www.rcsb.org/) – Repository for macromolecular structures.
- [PDBe](https://www.ebi.ac.uk/pdbe/) – European counterpart to RCSB PDB.
- [OPM](https://opm.phar.umich.edu/) – Orientation of proteins in membranes.
- [UniProt](https://www.uniprot.org/) – Protein sequences, structures, and functions.
- [InterPro](https://www.ebi.ac.uk/interpro/) – Protein classification & domain prediction.
- [AlphaFold DB](https://alphafold.ebi.ac.uk/) – Predicted structures from AlphaFold.
- [Proteopedia](https://proteopedia.org/wiki/index.php/Main_Page) – Interactive protein visualizations.

### 2.2 Binding Site & Pocket Detection
- [ProteinsPlus](https://proteins.plus/) – Binding site analysis.
- [PrankWeb](https://prankweb.cz/) – Pocket prediction and analysis.
- [CASTp](http://sts.bioe.uic.edu/castp/index.html?2r7g) – Pocket geometry and volume analysis.
- [CavityPlus](http://www.pkumdl.cn:8000/cavityplus/index.php#/) – Pocket detection and druggability.
- [CaverWeb](https://loschmidt.chemi.muni.cz/caverweb/) – Tunnel and channel detection.
- [PASSer](https://passer.smu.edu/) – Allosteric site prediction.

### 2.3 Protein Engineering & Modeling
- [DynaMut](https://biosig.lab.uq.edu.au/dynamut/) – Predicts mutation-induced stability changes.

---

## 3. Ligand Design & Optimization

### 3.1 Pharmacophore Modeling
- [ZINCPharmer](http://zincpharmer.csb.pitt.edu/) – Pharmacophore screening.
- [Pharmit](https://pharmit.csb.pitt.edu/) – Interactive pharmacophore modeling.
- [PharmMapper](https://www.lilab-ecust.cn/pharmmapper/) – Pharmacophore mapping.

### 3.2 QSAR & Descriptor Tools
- [QSAR Toolbox](https://qsartoolbox.org/) – Hazard assessment & QSAR.
- [OCHEM](https://ochem.eu/home/show.do) – QSAR model building & prediction.
- [ChemMaster](https://crescent-silico.com/chemmaster/) – QSAR and cheminformatics suite.
- [3D-QSAR](https://www.3d-qsar.com/) – 3D-QSAR modeling resources.
- [QSAR-Co](https://sites.google.com/view/qsar-co/) – Robust multitarget QSAR modeling.
- [DataWarrior](https://openmolecules.org/datawarrior/) – Free software for chemical analysis, QSAR, and visualization.
- [KNIME](https://www.knime.com/) – Workflow platform for cheminformatics & ML integration.

### 3.3 Molecular Property Prediction
- [SwissADME](http://www.swissadme.ch/) – Drug-likeness and PK.
- [pkCSM](https://biosig.lab.uq.edu.au/pkcsm/) – ADMET property prediction.
- [DeepPK](https://biosig.lab.uq.edu.au/deeppk/) – DL-based pharmacokinetics.
- [admetSAR 2.0](https://lmmd.ecust.edu.cn/admetsar2/) – Comprehensive ADMET.
- [ADMETlab 2.0](https://admetmesh.scbdd.com/) – PK, toxicity & drug-likeness.
- [ProTox-II](https://tox-new.charite.de/protox_II/) – Toxicity predictions.
- [PreADMET](https://preadmet.webservice.bmdrc.org/) – PK property predictions.
- [FAF-Drugs](https://bioserv.rpbs.univ-paris-diderot.fr/services.html) – ADMET filtering.
- [Admetboost](https://ai-druglab.smu.edu/admet) – ML-based ADMET prediction.

### 3.4 Fragment-Based Drug Design

- [SwissSidechain](https://www.swisssidechain.ch/) – Fragment and linker library for small molecule design.  
- [FragBuilder](https://github.com/andersx/fragbuilder) – Python API for building peptide-like and small molecule fragments.  
- [SeeSAR](https://www.biosolveit.de/SeeSAR/) – Fragment growing and linking software (free academic version).
- [Enamine Fragment Libraries](https://enamine.net/compound-libraries/fragment-libraries) – Large curated collection of diverse fragments for FBDD.
  
---

## 4. Virtual Screening & Docking
- [OpenBabel](https://openbabel.org/index.html) – Format conversion & ligand prep.
- [MGLTools](https://ccsb.scripps.edu/mgltools/) – Structure preparation.
- [AutoDockTools](https://autodocksuite.scripps.edu/adt/) – AutoDock GUI.
- [AutoDock Vina](https://vina.scripps.edu/) – Popular docking software.
- [EasyDockVina2](https://github.com/S3cr3t-SDN/EasyDockVina2) – Vina automation.
- [Webina](https://durrantlab.pitt.edu/webina/) – Web-based Vina.
- [Smina](https://github.com/mwojcikowski/smina) – Vina fork with extra features.
- [Gnina](https://github.com/gnina/gnina) – CNN-scoring docking.
- [EasyDock](https://github.com/ci-lab-cz/easydock) – Vina/Smina pipeline.
- [HADDOCK](https://wenmr.science.uu.nl/haddock2.4/) – Flexible docking suite.
- [PandaDock](https://github.com/pritampanda15/PandaDock) – Python docking tool.
- [ZDOCK](https://zdock.wenglab.org/) – Protein–protein docking.
- [ClusPro](https://cluspro.org/) – Protein–protein docking server.
- [pyDockWEB](https://life.bsc.es/pid/pydockweb/) – Electrostatics-based docking.
- [SwissDock](https://www.swissdock.ch/) – Web docking for beginners.
- [MzDOCK](https://github.com/Muzatheking12/MzDOCK) – GUI docking pipeline.
- [Uni-Mol Docking V2](https://www.bohrium.com/apps/unimoldockingv2/job?type=app) – AI-assisted docking.
- [Vina on Colab](https://autodock-vina.readthedocs.io/en/latest/colab_examples.html) – Run Vina in Google Colab.

---

## 5. Interaction Analysis & Visualization
- [PLIP](https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index) – Protein–ligand interaction profiling.
- [LigPlot+](https://www.ebi.ac.uk/thornton-srv/software/LigPlus/) – 2D interaction diagrams.
- [Discovery Studio Visualizer](https://discover.3ds.com/discovery-studio-visualizer-download) – Advanced visualization.

---

## 6. Molecular Dynamics & Simulation

### 6.1 Engines
- [GROMACS](https://www.gromacs.org/) – Scalable MD engine.
- [LAMMPS](https://www.lammps.org/) – Parallel MD for materials and biomolecules.
- [NAMD](https://www.ks.uiuc.edu/Research/namd/) – High-performance biomolecular MD.
- [AMBER](https://ambermd.org/) – Suite for biomolecular simulations.
- [Desmond](https://www.deshawresearch.com/resources.html) – High-performance MD.

### 6.2 Topology & Force Field Tools
- [CGenFF](https://cgenff.umaryland.edu/) – CHARMM force field parametrization.
- [SwissParam](https://www.swissparam.ch/) – Small molecule parameters.
- [ATB](https://atb.uq.edu.au/) – Automated topology builder.
- [CHARMM-GUI](https://www.charmm-gui.org/) – Input & topology preparation.
- [LigParGen](https://zarbi.chem.yale.edu/ligpargen/) – Ligand FF parameters.

### 6.3 Analysis Tools
- [MD DaVis](https://md-davis.readthedocs.io/en/latest/index.html) – Interactive MD visualizations.
- [iMod](https://imods.iqfr.csic.es/) – Normal Mode Analysis.
- [MolAiCal](https://molaical.github.io/) – Binding free energy calculations.

---

## 7. Synthesis & Retrosynthesis Planning
- [Spaya](https://spaya.ai/app/search) – AI retrosynthesis.
- [AiZynthFinder](https://github.com/MolecularAI/aizynthfinder) – Monte Carlo retrosynthesis.
- [ASKCOS](https://askcos.mit.edu/) – Reaction planning with ML.
- [IBM RoboRXN](https://rxn.res.ibm.com/rxn/robo-rxn/welcome) – Automated reaction prediction.
- [MANIFOLD](https://app.postera.ai/manifold/) – Synthetic accessibility & search.

---

## 8. Specialized Modalities

### 8.1 PROTACs & Ternary Complexes
- [PROTAC-db](http://cadd.zju.edu.cn/protacdb/) – PROTAC data.
- [PROsettaC](https://prosettac.weizmann.ac.il/) – Ternary complex modeling.

### 8.2 Peptide Design
- [PepDraw](https://pepdraw.com/) – Peptide visualization.
- [PepSite](http://pepsite2.russelllab.org/) – Predict peptide binding sites.
- [Peptimap](https://peptimap.bu.edu/) – Peptide mapping.
---

## 9. Machine Learning & AI for Drug Discovery
- **Core Libraries**: [RDKit](https://www.rdkit.org/), [Pandas](https://pandas.pydata.org/), [NumPy](https://numpy.org/), [Scikit-learn](https://scikit-learn.org/stable/), [Matplotlib](https://matplotlib.org/), [Seaborn](https://seaborn.pydata.org/).
- **Deep Learning**: [Keras](https://keras.io/), [TensorFlow](https://www.tensorflow.org/), [PyTorch](https://pytorch.org/), [DeepChem](https://deepchem.io/), [TorchDrug](https://torchdrug.ai/), [DEEPScreen](https://github.com/cansyl/DEEPScreen), [GraphINVENT](https://github.com/MolecularAI/GraphINVENT).
- **Datasets & Platforms**: [MoleculeNet](https://moleculenet.org/), [Kaggle](https://www.kaggle.com/), [Hugging Face](https://huggingface.co/), [Code Ocean](https://codeocean.com/), [Zenodo](https://zenodo.org/), [ChemML](https://hachmannlab.github.io/chemml/index.html), [Datagrok](https://datagrok.ai/cheminformatics), [The Illustrated Machine Learning](https://illustrated-machine-learning.github.io/).

---

## 10. Utility & Workflow Tools
- [OPSIN](https://opsin.ch.cam.ac.uk) – IUPAC name to structure.
- [OSRA](https://cactus.nci.nih.gov/cgi-bin/osra/index.cgi) – Image to structure.
- [MetaPredict](http://metapredict.icoa.fr/) – Molecular property prediction.
- [ChemPlot](https://chemplot.streamlit.app/) – Chemical space visualization.
- [ChemDB](http://cdb.ics.uci.edu/) – Chemoinformatics portal.
- [BoBER](http://bober.insilab.org/) – Bioisosteric replacements.
- [Open Targets](https://platform.opentargets.org/) – Target identification.
- [Screening Explorer](http://stats.drugdesign.fr/) – Screening data analysis.
- [LigRMSD](https://ligrmsd.appsbio.utalca.cl/) – Ligand RMSD calculation.
- [NERDD](https://nerdd.univie.ac.at/) – Drug discovery resources.
- [MetaChemiBio](https://biochemia.uwm.edu.pl/metachemibio/) – Property prediction.
- [WenMR Portal](https://wenmr.science.uu.nl/) – Biomolecular interactions software.
- [LigBuilder3](http://www.pkumdl.cn:8080/ligbuilder3/) – Ligand design.
- [ChemMine Tools](https://chemminetools.ucr.edu/) – Cheminformatics tools.
- [MayaChemTools](http://www.mayachemtools.org/index.html) – Perl/Python scripts for cheminformatics.
- [SCBDD](http://www.scbdd.com/) – Cheminformatics & drug discovery software.
- [Click2Drug](https://www.click2drug.org/) – CADD software and DB directory.
- [Galaxy Europe](https://usegalaxy-eu.github.io/index-cheminformatics.html) – Galaxy instance for cheminformatics.
- [CADD Vault](https://drugbud-suite.github.io/CADD_Vault/) – CADD resources repository.
- [BioMoDes](https://abeebyekeen.com/biomodes-biomolecular-structure-prediction/) – Biomolecular modeling tools.
- [PlayMolecule](https://open.playmolecule.org/landing) – Molecular modeling simulations.
- [Venny 2.1](https://bioinfogp.cnb.csic.es/tools/venny/) – Venn diagram tool.

---

## 11. Learning Resources

### 11.1 Free Courses
- [TMP Chem Lectures](https://youtube.com/playlist?list=PLm8ZSArAXicIWTHEWgHG5mDr8YbrdcN1K) – Computational chemistry lectures.
- [Strasbourg Summer School in Chemoinformatics](https://youtube.com/playlist?list=PLhgURFExPmJsDuHevu5n8y0R41WsXfbnC) – Summer school lectures.
- [BIGCHEM](https://bigchem.eu/node/63) – Big data in chemistry course.
- [Geometric Deep Learning Course](https://geometricdeeplearning.com/lectures/) – GDL for ML in science.
- [Drug Discovery Course](https://www.stereoelectronics.org/webDD/DD_home.html) – Drug discovery fundamentals.
- [drugdesign.org](https://www.drugdesign.org/) – Drug design and cheminformatics courses.
- [Cheminformatics OLCC](https://chem.libretexts.org/Courses/Intercollegiate_Courses/Cheminformatics) – Cheminformatics theory & coding.
- [Python For Cheminformatics Docking](https://pdb101.rcsb.org/train/training-events/python4) – Python-based docking.
- [DDA CDD Workshop](https://wcair.dundee.ac.uk/training/training-resources/computational-drug-design/) – Generative drug design.

### 11.2 Blogs
- [Practical Fragments](http://practicalfragments.blogspot.com/) – Fragment-based drug design.
- [avrilomics](https://avrilomics.blogspot.com/) – Genomics & bioinformatics.
- [Practical Cheminformatics](http://practicalcheminformatics.blogspot.com/) – Cheminformatics tools.
- [Cheminformania](https://www.cheminformania.com/) – Cheminformatics & deep learning.
- [Daily Dose of Data Science](https://www.blog.dailydoseofds.com/) – Data science insights.
- [Machine Learning Mastery](https://machinelearningmastery.com/) – ML tutorials.
- [Chem-Workflows](https://chem-workflows.com/index.html) – Jupyter chemistry tutorials.
- [Structural Bioinformatics](https://proteinstructures.com/) – Structure-based drug design guide.
- [Bioinformatics Answers](https://www.biostars.org/) – Bioinformatics Q&A.
- [McConnellsMedChem](https://mcconnellsmedchem.com/) – Medicinal chemistry blog.
- [DrugDiscovery.NET](http://www.drugdiscovery.net/) – AI in drug discovery.
- [MacinChem](https://macinchem.org/) – Comp chem on macOS.
- [Jeremy Monat](https://bertiewooster.github.io/) – Cheminformatics research.
- [Angelo Raymond Rossi](https://angeloraymondrossi.github.io/) – Computational chemistry research.

### 11.3 Instructional Notebooks
- [TeachOpenCADD](https://projects.volkamerlab.org/teachopencadd/all_talktorials.html) – Jupyter tutorials for CADD.
- [intro_pharma_ai](https://github.com/kochgroup/intro_pharma_ai) – AI in pharma with notebooks.
- [AI/DL for Life Sciences](https://onlinelibrary.wiley.com/doi/10.1002/ardp.202200628) – Interactive AI/DL notebooks.

