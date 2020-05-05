#!/usr/bin/python3
# -*- coding: utf-8 -*-

from Bio import Entrez
from Bio import SeqIO
import argparse, time, re, sys, csv
from pathlib import Path
from urllib.error import HTTPError
from socket import error as SocketError

parser = argparse.ArgumentParser(description = 'This toolkit automatize data recovery from NCBI',formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument("-in", "--input", help="File with access IDs",  required=True)
parser.add_argument("-db","--database",help="Select database that you want retrieve informations", required=True, choices=['nucleotide','protein'])
parser.add_argument("-ot","--outtype",help="Select the output format, folowing the rules:\n\n▛▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▜\n▎      Database       |       Options      |    format    ▐\n▎---------------------------------------------------------▐\n▎     Nucleotide      | fasta,gb,          |  text or xml ▐ \n▎                     | fasta_cds_(aa/na)  |     text     ▐\n▎---------------------------------------------------------▐\n▎      Protein        | fasta              |  text or xml ▐\n▙▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▟\n\n", default = 'fasta', choices=['fasta','fasta_cds_na','fasta_cds_aa','gb'])
parser.add_argument("-tf","--textformat",help="If you choose fasta or gb format output, you can choose the text format, between text our xml, default = text",default='text',choices=['text','xml'])
parser.add_argument("-tx","--taxonomy",help="This parameter provides an extra file with taxonomy of sequences", default='False', choices=['False','True'])
args = parser.parse_args()
read_file = args.input
entrez_database = args.database
out_type = args.outtype
text_format = args.textformat
out_taxonomy = args.taxonomy

protein_negative_list = ['fasta_cds_na','fasta_cds_aa']
out_tax_name = read_file.rstrip('\n')
if out_taxonomy == 'True':
    output_taxonomy = open(out_tax_name+'.tax.csv','w',newline='')
    writer_tax = csv.writer(output_taxonomy)
    writer_tax.writerow(["NCBI-Access", "Taxonomy-Access", "Superkingdom","Kingdom","Phylum","Class","Order","Family","Genus","Species"])

#any(x in a_string for x in matches)
if entrez_database == 'protein' and any(x in out_type for x in protein_negative_list):
    sys.exit("ERROR: Invalid options of output format for protein database")
if 'fasta_cds' in out_type and text_format == 'xml':
    sys.exit("ERROR: Invalid options of text format for output type")
#Entrez.email = str(input('Digite your e-mail to access NCBI:'))
Entrez.email = 'zimmer.filipe2@gmail.com'

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
    if out_taxonomy == 'True':
        list_of_lists = []
        list_of_tax = []
        list_of_tax.append(out_name)
        handle_tax = Entrez.efetch(db = var_data, id = var_id, rettype = var_rettype, retmode=var_retmode)
        handle_tax_var = str(re.findall(r'\/db_xref="taxon:.*',handle_tax.read()))
        handle_tax_var = re.sub(r"\[.*:", '', handle_tax_var)
        handle_tax_var = re.sub(r'".*\]', '', handle_tax_var)
        list_of_tax.append(handle_tax_var)
        handle_tax_data = Entrez.efetch(db="Taxonomy",id=handle_tax_var, retmode="xml")
        record_tax_data = Entrez.read(handle_tax_data)
        for x in record_tax_data[0]["LineageEx"]:
            if 'superkingdom' in x.values():
                superkingdom = x['ScientificName']
                list_of_tax.append(superkingdom)
            if 'kingdom' in x.values():
                kingdom = x['ScientificName']
            if 'phylum' in x.values():
                phylum = x['ScientificName']
                list_of_tax.append(phylum)
            if 'class' in x.values():
                classe = x['ScientificName']
                list_of_tax.append(classe)
            if 'order' in x.values():
                order = x['ScientificName']
                list_of_tax.append(order)
            if 'family' in x.values():
                family = x['ScientificName']
                list_of_tax.append(family)
            if 'genus' in x.values():
                genus = x['ScientificName']
                list_of_tax.append(genus)
            if 'species' in x.values():
                species = x['ScientificName']
                list_of_tax.append(species)
        writer_tax.writerow(list_of_tax)
    #ver uma forma de acelerar esse processo
    return()


#/db_xref="taxon:1608128"
#/db_xref="CDD:278997"


with open(read_file,'r+') as lst_terms:
    for i in lst_terms:
        out_name = i.rstrip('\n')
        if text_format == 'xml':
            with open(out_name+'.'+out_type+'.xml','w+') as out_handle:
                try:
                    efetch_function(entrez_database,i,out_type,text_format)
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
                except RuntimeError:
                    print("Runtime error in sequence acess: ",i)
                    continue
                except HTTPError:
                    print("HTTPError in sequence acess: ",i,"Possible reasons:\n1-No inthernet connection\n2-E-mail don't registred on NCBI\n3-The access IDs present on input file don't correspond to selected entrez database")
                    time.sleep(0.35)
                    continue
