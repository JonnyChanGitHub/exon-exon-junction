##########################################################################################################################################
# Script for constructing the putative compatible exon-exon junction database.                                                           #
# Copyright (c) 2008, Fan Mo, Xu Hong, System Biology Division, Zhejiang-California Nanosystems Institute(ZCNI) of Zhejiang University.  #
# Released under the terms of the GNU General Public License(version 2.1 or later).                                                      #
# Contact Email: mofan.hz@gmail.com, hongxu@zju.edu.cn                                                                                   #
##########################################################################################################################################


#!/usr/bin/perl -w
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my %nucleohash=("GCA","A","GCC","A","GCG","A","GCU","A","AGA","R","AGG","R","CGA","R","CGC","R","CGG","R","CGU","R","AAC","N","AAU","N","GAC","D","GAU","D","UGC","C","UGU","C","GAA","E","GAG","E","CAA","Q","CAG","Q","GGA","G","GGC","G","GGG","G","GGU","G","CAC","H","CAU","H","AUA","I","AUC","I","AUU","I","UUA","L","UUG","L","CUA","L","CUC","L","CUG","L","CUU","L","AAA","K","AAG","K","AUG","M","UUC","F","UUU","F","CCA","P","CCC","P","CCG","P","CCU","P","AGC","S","AGU","S","UCA","S","UCC","S","UCG","S","UCU","S","ACA","T","ACC","T","ACG","T","ACU","T","UGG","W","UAC","Y","UAU","Y","GUA","V","GUC","V","GUG","V","GUU","V","UAA","*","UAG","*","UGA","*",);
my $count=0;

# function to get the reverse complement sequence
sub revComp 
{
  my $r=reverse($_[0]);
  $r=~tr/ACGTacgt/TGCAtgca/;
  return $r;
}

# function to translate the nuclear acid sequence into amono acid sequence
sub trans
{
      my($str,$conj_pos)=@_;
      $str=~tr/T/U/;
      my $peptide="";
      my $i=0;
      my($nuc,$aa,$index);
      for($i=0;$i<int(length($str)/3);$i++)
      {
        $nuc=substr($str,$i*3,3);

        if(defined($nucleohash{$nuc}))
        {
			$aa=$nucleohash{$nuc};
            $peptide=$peptide.$aa;
        }
      }
      $index=index($peptide,"*");
      if(($index>=0) && ($index<=$conj_pos-1))
      {
        $peptide="";
      } 
      elsif(($index>$conj_pos-1) && ($index < length($peptide)))
      {
        $peptide=substr($peptide,0,$index);
      } 
      return $peptide;
}

my $start=1;
my $end=100000000;
my $output="Novel_Compatible_EEJ_peptide.fa";
open (OUT_P,">$output") || die "can not open $output for writing\n";
$loop_start = $ARGV[0] || 1;
print $loop_start, "\n";
$num_chroms = $ARGV[1] || 25;
$db_name = $ARGV[2] || 'homo_sapiens_core_70_37';
for ($loop=$loop_start;$loop<$num_chroms;$loop++)
{
	if ($loop == ($num_chroms-2))
	{
		$chrom="X";
	}
	elsif ($loop == ($num_chroms-1))
	{
		$chrom="Y";
	}
    else
    {
        $chrom=$loop;
    }

	# Bioperl mysql database interface
	my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => 'useastdb.ensembl.org',
									            -user => 'anonymous',
									            -dbname => $db_name,
                              -port => '5306',
									            -pass => '');
	$slice_adaptor = $db->get_SliceAdaptor();
	$slice = $slice_adaptor->fetch_by_region('chromosome',$chrom, $start, $end);

	# Main function to construct the putative compatible exon-exon junction database.
	foreach my $gene (@{$slice->get_all_Genes()}) 
	{
		#trans_exon hash here is used to exclude those previously described exon-exon junction sequences
		 my %trans_exon;    
		 my %exons;
		 foreach my $trans (@{$gene->get_all_Transcripts()}) 
		 {
			my $pre_exon="";
			foreach my $exon (@{$trans->get_all_Exons()}) 
			{
				if(exists($trans_exon{$exon->stable_id}))
				{
					if ($trans_exon{$exon->stable_id}=~/$pre_exon/)
					{
						$pre_exon=$exon->stable_id;
						next;
					}
					else
					{
						$trans_exon{$exon->stable_id}=$trans_exon{$exon->stable_id}."|$pre_exon";

					}
				}
				else
				{
					$trans_exon{$exon->stable_id}=$pre_exon;
				}
				$pre_exon=$exon->stable_id;
			}
		}

		my ($i,$j,$pos);
		my $gene_id=$gene->stable_id();
		my $gene_strand=$gene->strand();
		my $seq_latter;
		my $seq_former;
		my $seq_conj;
		my $seq_re_com;
		my %exon_phase;
		my %exon_end_phase;
		my %exon_start;
		my %exon_end;
		my %exon_strand;
		my %exon_length;

		foreach my $exon (@{$gene->get_all_Exons()}) 
		{
			$exon_id=$exon->stable_id;
			$exon_strand=$exon->strand;
			$exon_start=$exon->start;
			$exon_end=$exon->end;
			$exon_phase=$exon->phase;
			$exon_end_phase=$exon->end_phase;
			$exon_length=$exon_end-$exon_start+1;
			if (exists($exon_phase{$exon_id}))
			{
				next;
			}
			else
			{
				$exon_phase{$exon_id}=$exon_phase;
				$exon_end_phase{$exon_id}=$exon_end_phase;
				$exon_strand{$exon_id}=$exon_strand;
				$exon_start{$exon_id}=$exon_start;
				$exon_end{$exon_id}=$exon_end;
				$exon_length{$exon_id}=$exon_length;
			}
		}

		# If the gene is located in forward strand
		if($gene_strand==1)          
		{
			foreach $id (keys(%exon_start))
			{
				# In the case the exon could not produce 25 amino acids, we took as many amino acids as possible from the exon for the junction sequence
				if($exon_length{$id}<76)
				{   
					$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$id}, $exon_end{$id});
				}
				else
				{
					$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$id}, $exon_start{$id}+75);
				}

				#get sequence of acceptor exon
				$seq_acceptor=$seq_conj->seq();          
				
				foreach $sub_id (keys(%exon_start)) 
				{
					if ($exon_start{$id} > $exon_end{$sub_id} && $trans_exon{$id}!~$sub_id && $trans_exon{$sub_id}!~$id && $exon_phase{$id} eq $exon_end_phase{$sub_id} && $exon_phase{$id} ne "-1")
					{
						if ($exon_end_phase{$sub_id} eq "0")
						{
							if ($exon_length{$sub_id}<75)
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_end{$sub_id});
								$shift=$exon_length{$sub_id}%3;
								$pos=int($exon_length{$sub_id}/3)+1;
							}
							else
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_end{$sub_id}-74, $exon_end{$sub_id});
								$pos=26;
								$shift=0;
							}
						}
						elsif ($exon_end_phase{$sub_id} eq "1")
						{
							if ($exon_length{$sub_id}<76)
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_end{$sub_id});
								$shift=($exon_length{$sub_id}-1)%3;
								$pos=int(($exon_length{$sub_id}-1)/3)+1;
							}
							else
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_end{$sub_id}-75, $exon_end{$sub_id});
								$pos=26;
								$shift=0;
							}
						}
						else
						{
							if ($exon_length{$sub_id}<74)
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_end{$sub_id});
								$shift=($exon_length{$sub_id}-2)%3;
								$pos=int(($exon_length{$sub_id}-2)/3)+1;
							}
							else
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_end{$sub_id}-73, $exon_end{$sub_id});
								$pos=26;
								$shift=0;
							}
						}

						#get sequence of donor exons
						my $seq=$seq_conj->seq();   
						
						#get sequence of exon-exon junctions
						$seq.=$seq_acceptor;
						$jun_seq=substr($seq,$shift);

						# get peptide sequence according to phase of exons, write into fasta file
						$seq_peptide=trans($jun_seq,$pos);
						if ($seq_peptide ne "")
						{
							 print OUT_P ">$gene_id|$chrom|$gene_strand|".$exon_end_phase{$sub_id}."|$sub_id:".$exon_length{$sub_id}."-$id:".$exon_length{$id}."|$pos\n";
							 print OUT_P "$seq_peptide\n";
						}
					}
					else
					{
						next;
					}
				}
			}
		}

		#for gene strand eq -1
		else             
		{
			foreach $id (keys(%exon_start))
			{
				# In case the exon could not produce 25 amino acids, we took as many amino acids as possible from the exon for the junction sequence
				if($exon_length{$id}<76)
				{   
					$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$id}, $exon_end{$id});
				}
				else
				{
					$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_end{$id}-75, $exon_end{$id});
				}

				#get sequence of acceptor exon
				$seq_acceptor=$seq_conj->seq(); 
				foreach $sub_id (keys(%exon_start))
				{
					if ($exon_start{$sub_id} > $exon_end{$id} && $trans_exon{$id}!~$sub_id && $trans_exon{$sub_id}!~$id && $exon_phase{$id} eq $exon_end_phase{$sub_id} && $exon_phase{$id} ne "-1")
					{
						if ($exon_phase{$sub_id} eq "0")
						{
							if ($exon_length{$sub_id}<75)
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_end{$sub_id});
								$shift=$exon_length{$sub_id}%3;
								$pos=int($exon_length{$sub_id}/3)+1;
							}
							else
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_start{$sub_id}+74);
								$pos=26;
								$shift=0;
							}
						}
						elsif ($exon_end_phase{$sub_id} eq "1")
						{
							if ($exon_length{$sub_id}<76)
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_end{$sub_id});
								$shift=($exon_length{$sub_id}-1)%3;
								$pos=int(($exon_length{$sub_id}-1)/3)+1;
							}
							else
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_start{$sub_id}+75);
								$pos=26;
								$shift=0;
							}
						}
						else
						{
							if ($exon_length{$sub_id}<74)
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_end{$sub_id});
								$shift=($exon_length{$sub_id}-2)%3;
								$pos=int(($exon_length{$sub_id}-2)/3)+1;
							}
							else
							{
								$seq_conj= $slice_adaptor->fetch_by_region('chromosome',$chrom, $exon_start{$sub_id}, $exon_start{$sub_id}+73);
								$pos=26;
								$shift=0;
							}
						}

						#get sequence of donor exon
						my $seq=$seq_conj->seq();    
						
						#get sequence of exon-exon junctions						
						$seq.=$seq_acceptor.$seq;
						$seq_re_com=revComp($seq);
						$jun_seq=substr($seq,$shift);

						# get peptide sequence according to phase of exons, write into fasta file
						$seq_peptide=trans($jun_seq,$pos);
						if ($seq_peptide ne "")
						{
							 print OUT_P ">$gene_id|$chrom|$gene_strand|".$exon_end_phase{$sub_id}."|$sub_id:".$exon_length{$sub_id}."-$id:".$exon_length{$id}."|$pos\n";
							 print OUT_P "$seq_peptide\n";
						}
					}
				}
			}
		}			
	}
}
close OUT_P;

