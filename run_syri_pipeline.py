import subprocess
import pandas as pd
from pathlib import Path

# minimap2 -ax asm5 -t 4 --eqx CGMH058_GCA_013378355.1_genomic.fna SK36_GCA_000014205.1_genomic.fna | samtools sort -O BAM - > CGMH058_SK36.bam
# syri -c CGMH058_SK36.bam -r CGMH058_GCA_013378355.1_genomic.fna -q SK36_GCA_000014205.1_genomic.fna -F B --prefix CGMH058_SK36
# plotsr --sr CGMH058_SK36syri.out --genomes genomes.txt -o output.png -H 1 -W 5

def check_files(manifest_file, file_folder=None):
    df = pd.read_csv(manifest_file, sep='\t')
    if file_folder:
        df['ref_fasta'] = df['ref_fasta'].apply(lambda x: str(Path(file_folder).joinpath(x)))
        df['query_fasta'] = df['query_fasta'].apply(lambda x: str(Path(file_folder).joinpath(x)))

    for index, row in df.iterrows():
        ref_fasta_path = Path(row['ref_fasta'])
        query_fasta_path = Path(row['query_fasta'])
        
        if not ref_fasta_path.exists():
            raise FileNotFoundError(f"Reference fasta file not found: {ref_fasta_path}")
            
        if not query_fasta_path.exists():
            raise FileNotFoundError(f"Query fasta file not found: {query_fasta_path}")
    return df


def cp_genomes(ref_fasta, query_fasta, output):
    # cp ref_fasta query_fastas output_folder
    subprocess.run(['cp', ref_fasta, query_fasta, output],check=True)


def make_genome_file(bam_output, ref_fasta, ref_name, query_fasta, query_name):
    genome_file =  str(bam_output).strip('.bam') + '_genomes.txt'
    with open(genome_file, 'w') as f:
        f.write('#file\tname\n')
        f.write(ref_fasta + '\t' + ref_name +'\n')
        f.write(query_fasta + '\t' + query_name +'\n')
    
    print('\n')
    print(f'Genome file created: {genome_file}')


def plotsr(syri_out_folder):
    # plotsr --sr CGMH058_SK36syri.out --genomes genomes.txt -o output.png -H 1 -W 5
    #subprocess.run(['plotsr', '--sr', syri_out, '--genomes', genomes_file, '-o', output.png, '-H', '1', '-W', '5'],check=True)
    pass
def main(df):
    for i, row in df.iterrows():

        align_name = row['bam_output'].strip('.bam')
        Path(align_name).mkdir(parents=True, exist_ok=True)
        bam_output = Path(align_name).joinpath(row['bam_output'])

        cp_genomes(row['ref_fasta'],row['query_fasta'], Path(align_name))
        
        ######## minimap2 and samtools command ########
        minimap2_cmd = ["minimap2", "-ax", "asm5", "-t", "4", "--eqx", row['ref_fasta'], row['query_fasta']]
        samtools_cmd = ["samtools", "sort", "-O", "BAM", "-"]
        minimap2_process = subprocess.Popen(minimap2_cmd,stdout=subprocess.PIPE,text=True)

        with open(bam_output, "wb") as bam_file:
            subprocess.run( samtools_cmd, stdin=minimap2_process.stdout, stdout=bam_file, check=True )
            minimap2_process.wait()
        print('\n')
        print(f"Minimap2 and samtools command executed for: {row['bam_output']}")
        ################################################


        ################ syri command ##################
        syri_cmd = ["syri", "-c", bam_output, "-r", row['ref_fasta'], "-q", row['query_fasta'], "-F", "B", "--dir", Path(align_name)]
        subprocess.run(syri_cmd, check=True)
        print(f"SyRI command executed for: {row['bam_output']}")
        ################################################

        make_genome_file(bam_output, row['ref_fasta'], row['ref_name'], row['query_fasta'], row['query_name'])

        # break
        # '''



if __name__ == '__main__':
    manifest = 'mainfest_test.txt'
    a = check_files(manifest,r'20250211_S.sanguinis_comaprison_genome')
    # make_genome_file(a)
    main(a)
    # subprocess.run(["cat"],stdin=None)
