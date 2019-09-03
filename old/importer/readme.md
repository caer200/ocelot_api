#fomport
This is the main routine for generating data entries in FOMDB.
It consists of 3 parts:
1. Crawling CSD to get desired cif files of crystal structures,
2. Parsing the cif file into information stored in JSON,
3. Writing JSON into database. 

##crawl
- cifs20181102.zip, total counts 60209
    - at least 2 fused rings where for each ring, there're at least 2 unsaturated bonds
    - there's no metal/open-shell atom in the system
    - no polymer
    - one kind of molecule
##parser


##importer