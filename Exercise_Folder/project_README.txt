Single-Cell RNA-seq Workshop — Project Setup
============================================

This README explains how to organize the tutorial files on your own computer so
that everything runs cleanly and the intermediate results flow from one tutorial
to the next.


1. Create a project folder with a scripts/ subfolder
----------------------------------------------------

Make ONE project folder anywhere convenient (Desktop, Documents, …), and create
a "scripts" subfolder inside it. Save every tutorial you download into scripts/.

    my_scrnaseq_project/
    └── scripts/
        ├── Tutorial_00_Setup_RStudio_Packages.qmd
        ├── Tutorial_01_QC_Preprocessing.qmd
        ├── Tutorial_02_DimReduction_Clustering.qmd
        └── … (the rest as you reach them)

Tip: the website serves each tutorial as "Tutorial_XX_….qmd.txt". After
downloading, remove the trailing ".txt" so the name ends in ".qmd", then open it
in RStudio.


2. Run the tutorials FROM the scripts/ folder
---------------------------------------------

Open a tutorial's .qmd in RStudio from inside scripts/ and run the code chunks
(click the green ▶ on a chunk, or press Ctrl/Cmd+Enter). RStudio uses the folder
the .qmd lives in as the working directory, so the relative paths in the code
resolve correctly.

You can confirm where R is pointed with:

    getwd()        # should end in /my_scrnaseq_project/scripts


3. data/ and output/ are created for you
----------------------------------------

The first time you run a tutorial's setup chunk, it creates two folders as
SIBLINGS of scripts/ (i.e. at the project root):

    my_scrnaseq_project/
    ├── scripts/      ← your .qmd files (you created this)
    ├── data/         ← AUTO-CREATED: intermediate .rds objects (see below)
    └── output/       ← AUTO-CREATED: figures (.png) and tables (.csv), per module
        ├── Mod1/
        ├── Mod2/
        └── …

In the code these are written as "../data" and "../output/ModN" — the ".." means
"one level up from scripts/". Each setup chunk also prints the absolute path it
is writing to, e.g.:

    This module writes its figures & tables to: /…/my_scrnaseq_project/output/Mod1


4. data/ is the hand-off between tutorials
------------------------------------------

The core tutorials run as a chain. Each one reads the object the previous one
saved into data/, and writes its own:

    Tutorial 01  →  data/ifnb_preprocessed.rds
    Tutorial 02  →  data/ifnb_clustered.rds
    Tutorial 03  →  data/ifnb_annotated.rds
    Tutorial 04  →  data/ifnb_annotated_final.rds
    Tutorial 05  →  data/ifnb_integrated.rds
    Tutorial 06  →  data/ifnb_pseudobulk_de.csv
    Tutorial 08  →  data/ifnb_milo.rds

Because every tutorial shares the same project data/, KEEP ALL THE .qmd FILES
TOGETHER IN scripts/. If they are scattered across different folders, a later
tutorial won't find the file an earlier one saved and you'll see an error like
"cannot open the connection".

The bonus modules that use their own downloaded datasets (Module 13 WGCNA,
Module 15 scATAC) also read those files from data/ — the tutorial gives you the
download command; run it so the files land in your project's data/ folder.


5. Notes
--------

- The workshop dataset (ifnb) downloads automatically the first time Tutorial 01
  runs — nothing to place in data/ for the core tutorials.
- The published HTML pages show the code but do NOT run it (chunks are tagged
  eval: false). Running the chunks yourself (as above) is what creates the files.
- data/ and output/ can get large (~1.5 GB for Tutorials 01–08). They are fully
  reproducible from the tutorials, so it's safe to delete them and re-run.

Questions? See the Materials, Datasets, and FAQ pages on the workshop website.
