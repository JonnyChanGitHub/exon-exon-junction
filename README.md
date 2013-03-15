## Exon-Exon Junction Database Creation

# Summary 

Updates to the Putative Exon-Exon junction database creation code found here:
http://www.zcni.zju.edu.cn/en/download.htm


# Prerequistes

1) Install bioperl. 

2) Install ensembl-api.

Download ensembl-api from here http://www.ensembl.org/info/docs/api/index.html.

Extract the contents and add the resulting ensemble/modules directory
to your PERL5LIB path (adjust based on where you extract ensemble-api.tar.gz).

tar xzvf ensembl-api.tar.gz
export PERL5LIB=$PERL5LIB:$HOME/ensembl/modules 

## Running

    > git clone https://github.com/jmchilton/exon-exon-junction.git
    > cd exon-exon-junction
    # Build human junction-junction database
    > perl Build_Novel_Compatible_PEEJ_Peptide_Database_V1.pl 
    # Build pig junction-junction database
    > perl Build_Novel_Compatible_PEEJ_Peptide_Database_V1.pl 1 18 sus_scrofa_core_70_102

This will create the database Novel_Compatible_EEJ_peptide.fa in your
current directory.
