# ncbi-acess

This script it's a toolbox to automatic recovery information of [NCBI]([https://www.ncbi.nlm.nih.gov/](https://www.ncbi.nlm.nih.gov/)). 

## Dependencies

This script was build on python 3.6.5+ and have only two dependencies:

- [argparse 1.1](https://docs.python.org/3/library/argparse.html)
- [Bio (biopython) 1.76](https://biopython.org/)

## Recomended lectures
- [biopython.Entrez](https://biopython.org/docs/1.74/api/Bio.Entrez.html): To understand the basic sintax.
- [ncbi Entrez ebook](https://www.ncbi.nlm.nih.gov/books/NBK25499/): To check the databases, output types and file formats that entrez can access.
## Usage
- To recovery genbank information from nucleotide sequences:
> python ncbi_seq_retrieve.py -in file_with_access_ids.txt -db nucleotide -ot gb

Or to recovery in xml format, just insert the parameter -tf xml.

- To recovery cds translated to aminoacids from nucleotide sequences:
> python ncbi_seq_retrieve.py -in file_with_acess_ids.txt -db nucleotide -ot fasta_cds_aa

Or to recovery cds not translated, just change fasta_cds_aa for fasta_cds_na

- To recovery nucleotide of aminoacid sequences
> python ncbi_seq_retrieve.py -in file_with_acess_ids.txt -db (nucleotide or protein) -ot fasta

Or to recovery in xml format, just insert the parameter -tf xml.

- To recovery taxonomy information of ncbi acess IDs
> python ncbi_seq_retrieve.py -in file_with_acess_ids.txt -db (nucleotide or protein) -ot gb -tx True
## Some considerations
If you have a file with IDs from nucleotide sequences, you can't use this file in a protein database, and vice-versa. If you call help function, a table with which text formats are allowed per output type, and which output types are allowed per database. 

## Disclaimer

- I'm not a computer engineer or some related profession, I'm just write this script to study python and to automatize some bioinformatics tasks. So fell free to commit changes that makes the code more efficient or more clean.
- This script will continue to be developed to englobe others functions, like recovery taxonomy information of sequences and features of sequences, for example.
