from Bio import Entrez
from Bio import SeqIO
import argparse, time, re, sys
from pathlib import Path
from urllib.error import HTTPError
from socket import error as SocketError

parser = argparse.ArgumentParser(description = 'This toolkit automatize data recovery from NCBI',formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-in", "--input", help="File with access IDs",  required=True)
parser.add_argument("-db","--database",help="Select database that you want retrieve informations", required=True, choices=['nucleotide','protein'])
parser.add_argument("-ot","--outtype",help="Select the output format, folowing the rules:\n\n▛▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▜\n▎      Database       |       Options      |    format    ▐\n▎---------------------------------------------------------▐\n▎     Nucleotide      | fasta,gb,          |  text or xml ▐ \n▎                     | fasta_cds_(aa/na)  |     text     ▐\n▎---------------------------------------------------------▐\n▎      Protein        | fasta              |  text or xml ▐\n▙▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▟\n\n", default = 'fasta', choices=['fasta','fasta_cds_na','fasta_cds_aa','gb'])
parser.add_argument("-tf","--textformat",help="If you choose fasta or gb format output, you can choose the text format, between text our xml, default = text",default='text',choices=['text','xml'])

args = parser.parse_args()
read_file = args.input
entrez_database = args.database
out_type = args.outtype
text_format = args.textformat
protein_negative_list = ['gb','fasta_cds_na','fasta_cds_aa']

#any(x in a_string for x in matches)
if entrez_database == 'protein' and any(x in out_format for x in protein_negative_list):
    sys.exit("ERROR: Invalid options of output format for protein database")
if 'fasta_cds' in out_type and text_format == 'xml':
    sys.exit("ERROR: Invalid options of text format for output type")
    
Entrez.email = str(input('Digite your e-mail to access NCBI:'))

def efetch_function(var_data,var_id,var_rettype,var_retmode):
    """
    This function runs efetch functionalities
    
    Keyword arguments:
    var_data = entrez database
    var_id = access code
    var_rettype = type of data to retrieve
    var_retmode = text format to retrieve
    """
    
    handle = Entrez.efetch(db = var_data, id = var_id, rettype = var_rettype, retmode=var_retmode)
    out_handle.write(handle.read())
    return(handle)

with open(read_file,'r+') as lst_terms:
    for i in lst_terms:
        out_name = i.rstrip('\n')
        if text_format == 'xml':
            with open(out_name+'.'+out_type+'.xml','w+') as out_handle:
                try:
                    efetch_function(entrez_database,i,out_type,text_format)
                    #handle = Entrez.efetch(db = 'sequences', id = i, rettype = 'fasta_cds_aa', retmode="text")
                    #out_handle.write(handle.read())
                #ver como reportar os IDs que não tem CDS...
                except RuntimeError:
                    print("Runtime error in sequence acess: ",i)
                    continue
                except HTTPError:
                    print("HTTPError in sequence acess: ",i,"Possible reasons:\n1-No inthernet connection\n2-E-mail don't registred on NCBI\n3-The access IDs present on input file don't correspond to selected entrez database")
                    time.sleep(0.35)
                    continue
        else:
            with open(out_name+'.'+out_type,'w+') as out_handle:
                try:
                    efetch_function(entrez_database,i,out_type,text_format)
                    #handle = Entrez.efetch(db = 'sequences', id = i, rettype = 'fasta_cds_aa', retmode="text")
                    #out_handle.write(handle.read())
                #ver como reportar os IDs que não tem CDS...
                except RuntimeError:
                    print("Runtime error in sequence acess: ",i)
                    continue
                except HTTPError:
                    print("HTTPError in sequence acess: ",i,"Possible reasons:\n1-No inthernet connection\n2-E-mail don't registred on NCBI\n3-The access IDs present on input file don't correspond to selected entrez database")
                    time.sleep(0.35)
                    continue
