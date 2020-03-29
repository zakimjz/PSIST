# Protein Indexing (PSIST)

PSIST uses suffix trees to index protein 3D structure. It first converts the 3D structure into a structure-feature sequence over a new structural alphabet, which is then used to index protein structures. The PSIST index makes it very fast to query for a matching structural fragment.

Feng Gao and Mohammed J. Zaki. PSIST: indexing protein structures using suffix trees. In IEEE Computational Systems Bioinformatics Conference. August 2005.

Feng Gao and Mohammed J. Zaki. Psist: a scalable approach to indexing protein structures using suffix trees. Journal of Parallel and Distributed Computing, 68(1):55–63, January 2008. special issue on Parallel Techniques for Information Extraction. doi:10.1016/j.jpdc.2007.07.008.


Compiling the code:
-------------------

running the following command under psist directory:

        $ aclocal
        $ autoconf
        $ automake
        $ ./configure
        $ make


Running the code:
-----------------

1. PDB pre-processing:

   go to the psist/test/data subdirectory, extract the coordinates of the backbone atoms (Ca, C and N) for each pdb file in the pdblist:

       java Distance pdbchainlist pdbdir coorddir backbone

   you must download all the PDBs in the list file of the psist/test/scopdatabase.txt, which also includes all the query PDBs.

a test example:

       java Distance testpdb.list pdb/ coord/ backbone


2. Running:

   go to the psist/test direcory, you can get help by runing ./query

   test example:

       ./query -l testquery.txt -d data/coord/ -L testdatabase.txt -D data/coord/ -m 10 -e 0 -b 2 -w 3

   example in the paper:

       ./query -l scopquery.txt -d PDB_coord -L scopdatabase.txt -D PDB_coord -m 15 -e 0 -b 2 -w 3
