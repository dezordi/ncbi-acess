from Bio import Entrez
from Bio import SeqIO
import argparse, time, re
from pathlib import Path
from urllib.error import HTTPError
from socket import error as SocketError
parser = argparse.ArgumentParser()
parser.add_argument("-in", "--input", help="File with access IDs",  required=True)
args = parser.parse_args()
read_file = args.input
Entrez.email = str(input('Digite your e-mail to access NCBI:'))
sequences = []

with open(read_file,'r+') as lst_terms:
    print("Accesing list OK!")
    for i in lst_terms:
        out_name = i.rstrip('\n')
        with open(out_name+'.aa.cds.fasta','w+') as out_handle:
            try:
                handle = Entrez.efetch(db = 'sequences', id = i, rettype = 'fasta_cds_aa', retmode="text")
                out_handle.write(handle.read())
            #ver como reportar os IDs que não tem CDS...
            except RuntimeError:
                print("Runtime error in sequence acess: ",i)
                continue
            except HTTPError:
                print("HTTPError in sequence acess: ",i)
                time.sleep(0.35)
                continue
