import sys

def merge_files(info_link_file, mm_tpm_file, nmr_tpm_file, output_file):
    try:
        # open input and output files
        with open(info_link_file, 'r') as info_link, \
             open(mm_tpm_file, 'r') as MM_TPM, \
             open(nmr_tpm_file, 'r') as NMR_TPM, \
             open(output_file, 'w') as final:
            
            
            mm = {}
            nmr = {}
            
            for i, line in enumerate(MM_TPM):
                fields = line.strip().split('\t')
                if i == 0:  #write header
                    final.write('\t'.join(fields) + '\t')
                else:
                    mm[fields[0]] = line.strip()  
            
            
            for i, line in enumerate(NMR_TPM):
                fields = line.strip().split('\t')
                if i == 0:  
                    final.write('\t'.join(fields) + '\n')
                else:
                    nmr[fields[0]] = line.strip()  # 
            
            # process info_link file
            for line in info_link:
                fields = line.strip().split()
                if len(fields) >= 6:
                    mmid = fields[2].split('.')[0]
                    nmrid = fields[5].split('.')[0]
                    if mmid in mm and nmrid in nmr:
                        final.write(mm[mmid] + '\t' + nmr[nmrid] + '\n')
    
    except FileNotFoundError as e:
        print(f"Error: {e}. Please check if all input files exist.")
        sys.exit(1)
    except Exception as e:
        print(f"Error occurred: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py info_link_file mm_tpm_file nmr_tpm_file output_file")
        sys.exit(1)
    
    info_link_file = sys.argv[1]
    mm_tpm_file = sys.argv[2]
    nmr_tpm_file = sys.argv[3]
    output_file = sys.argv[4]
    
    merge_files(info_link_file, mm_tpm_file, nmr_tpm_file, output_file)

