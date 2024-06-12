import pandas as pd
import gzip
import time
import argparse
from pyfaidx import Fasta


parser = argparse.ArgumentParser(usage="GetLongestTransFromGTF_Mouse.py --gtffile Mus_musculus.GRCm39.110.chr.gtf.gz --genome Mus_musculus.GRCm39.dna.toplevel.fa --outfile longest_trans.fa",
                                description="Extract longest transcript based on gtf format annotation file in ensembl database. There will be three automatically generated output files: allCDS_info.csv, allTranscripts_info.csv, final_longest_transcripts_info.csv",
                                epilog="Thank your for your support.")

# read gtf file
parser.add_argument('-g','--gtffile', type=str,action="store",dest="gtffile",metavar="gtffile",
                    help='input your GTF file with ".gz" format.')
# genome fasta file
parser.add_argument('-fa','--genome',type=str,action="store",dest="genome",metavar="genome",
                    help='your genome fasta file matched with your GTF file with ".fa/.fasta" format.')
# output file
parser.add_argument('-o','--outfile', type=str,action="store",dest="longestfile",metavar="longestfile",
                    help='output your longest transcript file. (longest_trans.fa)')

args = parser.parse_args()


gtffile = args.gtffile
genomefile = args.genome
outfile = args.longestfile

# extract gene annotation info
def extract_gene_info(attributes):
    attribute_list = attributes.split("; ")
    attribute_dict = {}
    for attribute in attribute_list:
        attribute = attribute.split(" ")
        if attribute[0] == "transcript_support_level":
            pass
        else:
            attribute_dict[attribute[0]] = attribute[1].strip('"')
    return attribute_dict

def main():
    """
    Extract longest transcript from gtf format annotation file based on ensembl database.
    If the transcript have the longest CDS or the one with the longest CDS and longest transcript.
    
    """

    # main fuction
    print("Your job is running, please wait...")
    ######################################################################
    job_start = time.time()
    #######################
    info = {}
    cds_info = {}
    with gzip.open(gtffile,'rt') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            fields = line.split('\t')
            type = fields[2]
            if type == 'exon':
                gene_name = extract_gene_info(fields[8]).get("gene_name"," ")
                gene_id = extract_gene_info(fields[8])["gene_id"]
                trans_id = extract_gene_info(fields[8])["transcript_id"]
                biotype = extract_gene_info(fields[8])["gene_biotype"]
                key = ' '.join([trans_id,gene_id,biotype,gene_name])
                # length of exon
                start = int(fields[3])
                end = int(fields[4])
                length = end - start + 1
                info.setdefault(key,0)
                info[key] += length

            if type == "CDS":
                gene_name = extract_gene_info(fields[8]).get("gene_name"," ")
                gene_id = extract_gene_info(fields[8])["gene_id"]
                trans_id = extract_gene_info(fields[8])["transcript_id"]
                biotype = extract_gene_info(fields[8])["gene_biotype"]
                key = ' '.join([trans_id,gene_id,biotype,gene_name])
                start = int(fields[3])
                end = int(fields[4])
                length = end - start + 1
                cds_info.setdefault(key,0)
                cds_info[key] += length
    ######################################
    # convert to DataFrame
    res = pd.DataFrame(pd.Series(info), columns = ['transcript_length'])
    cds_res = pd.DataFrame(pd.Series(cds_info), columns = ['cds_length'])
    # save the DataFrames as CSV files
    res.to_csv(r'allTranscripts_info.csv', index=True)
    cds_res.to_csv(r'allCDS_info.csv', index=True)
    # add gene ID column
    res['gene_id'] = [line.split(sep=' ')[1] for line in list(res.index[:])]
    res['original_index'] = res.index
    new_order = ['original_index','gene_id','transcript_length']
    res_ordered = res[new_order]
    cds_res['gene_id'] = [line.split(sep=' ')[1] for line in list(cds_res.index[:])]
    cds_res['original_index'] = cds_res.index
    new_order_cds = ['original_index','gene_id','cds_length']
    cds_res_ordered = cds_res[new_order_cds]
    # sort by transcript length
    res_sorted = res_ordered.sort_values(by = ['gene_id','transcript_length'],ascending=False)

    longest_index = pd.DataFrame(columns=['original_index', 'gene_id', 'length'])

    for index, row in cds_res_ordered.iterrows():
        gene_id = row['gene_id']
        cds_len = int(row['cds_length'])
        gene_condition = longest_index['gene_id'] == gene_id
        current_row = longest_index.loc[gene_condition].iloc[0] if gene_condition.any() else None
        if current_row is not None:
            current_len = int(current_row['length'])
            current_index = current_row['original_index']

            if cds_len > current_len or (cds_len == current_len and res_ordered.loc[index, 'transcript_length'] > res_ordered.loc[current_index, 'transcript_length']):
                longest_index.loc[gene_condition, 'length'] = cds_len
                longest_index.loc[gene_condition, 'original_index'] = index

        else:
            longest_index = longest_index.append({'original_index': index, 'gene_id': gene_id, 'length': cds_len},ignore_index=True)

    # the raw longest transcript ID (without concidering CDS)
    longest_id = res_sorted.drop_duplicates(subset=['gene_id'],keep='first').index.values.tolist()
    # extract info
    longest_data = res_sorted.loc[res_sorted.index.isin(longest_id)]
    #index
    longest_data['original_index'] = longest_data.index
    # order columns
    longest_data = longest_data[['original_index','gene_id','transcript_length']]
    # create a new empty dataframe for final longst transcript info
    longest_t = pd.DataFrame(columns=res_sorted.columns)

# extracte transcript info corresponding to the Longest CDS
    longest_t_list = [longest_t]
    for value in longest_index['original_index'].values:
        longest_t_line = res_ordered.loc[res_ordered['original_index'] == value]
        longest_t_list.append(longest_t_line)
    longest_t = pd.concat(longest_t_list, ignore_index=True)

# extract the longest transcript information for all gene IDs that do not appear in the longest CDS.
    for index in longest_id:
        gene_id_value = res_sorted.loc[index, 'gene_id']
        if gene_id_value not in longest_t['gene_id'].values:
            selected_row = longest_data.loc[longest_data['gene_id'] == gene_id_value]
            longest_t = pd.concat([longest_t, selected_row], ignore_index=True)

#save the final longest transcript file
    longest_t.to_csv(r'final_longest_transcripts_info.csv', index=False)

    ###########################################################
    # ID of longest transcript 
    transid = {str(line).split(sep=' ')[0]:line for line in list(longest_t.original_index)}

    infolist = []
    with gzip.open(gtffile,'rt') as gtf:
        for line in gtf:
            # skip
            if line.startswith('#'):
                continue
            # split
            fields = line.split('\t')
            # feature type
            type = fields[2]
            if type == 'exon':
                # pos
                chr = fields[0]
                start = fields[3]
                end = fields[4]
                strand = fields[6]
                # name
                gene_name = extract_gene_info(fields[8]).get("gene_name"," ")
                gene_id = extract_gene_info(fields[8])["gene_id"]
                trans_id = extract_gene_info(fields[8])["transcript_id"]
                if trans_id in transid:
                    infolist.append([chr,start,end,strand,type,gene_name,gene_id,trans_id,transid[trans_id]])
                else:
                    pass
            else:
                pass

    # to dataframe
    dfinfo = pd.DataFrame(infolist,columns=['chr','start','end','strand','type','gene_name','gene_id','trans_id','id'])
    dfinfo_1_strand = dfinfo[dfinfo['strand'] == '+']

    # descrese coord by - strand gene
    dfinfo_2_strand = dfinfo[dfinfo['strand'] == '-']
    dfinfo_2_strand = dfinfo_2_strand.sort_values(by = ['trans_id','start','end'],ascending = False)

    # merge
    df_fianl = pd.concat([dfinfo_1_strand,dfinfo_2_strand],axis=0)

    ###########################################################
    # extact  sequnece from genome

    # load genome
    genome = Fasta(genomefile)

    # chrmosome info
    chrmosome_list = genome.keys()

    # save in dict
    res = {}
    for line in range(0,df_fianl.shape[0]):
        
        # chromosome strand
        fileds = df_fianl.iloc[line]
        chrom = fileds['chr']
        strand = fileds['strand']
        start = int(fileds['start'])
        end = int(fileds['end'])
        # key
        key = fileds['id']
        # filter chromoseome
        if chrom in chrmosome_list:
            # extarct sequence
            if strand == '+':
                seq = genome[chrom][(start-1):end].seq
            elif strand == '-':
                seq = genome[chrom][(start-1):end].complement.reverse.seq
            else:
                pass
            # save in dict
            res.setdefault(key,'')
            res[key] += seq
        else:
            pass
    
    ###########################################################
    # output file
    outputfile = open(outfile,'w')

    # set the length of sequence for each column
    my_length = 60
    ###########################################################

    for key,val in res.items():
        outputfile.write('>' + key + '\n')
        while len(val) > my_length:
            outputfile.write(val[0:my_length] + '\n')
            val = val[my_length:len(val)]
        outputfile.write(val + '\n')

    outputfile.close()

    ####################################################################################
    job_stop = time.time()
    print("Your job is done! ")
    print("Running with " + str(round(job_stop - job_start,2)) + " seconds!")

if __name__=="__main__":
    main()
