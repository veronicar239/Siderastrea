# Siderastrea holobiont de novo transcriptome assembly
## Workflow - summary of files

Original directory:
> /cm/shared/courses/dbarshis/15AdvBioinf/classdata/Siderastrea_radians/QCFastqs/nofilter/newpaired/trinity_out_PE-150_SE-150_control-pop/

Moved to my directory (current location):
> /cm/shared/courses/dbarshis/barshislab/VRad/taxons/Siderastrea_siderea/assembly_denovo/


- Original Trinity output:

Trinity.fasta

- Renamed:

834637_Trinity_Sid_PE150-SE150_Control.fasta

- Longest isoform:

690120_Trinity_Sid_Control_LongestIsoform.fasta

- Edited to remove >> that was introduced:

690120_Trinity_Sid_Control_LongestIsoform_edited.fasta

- 500 length threshold:

146542_Trinity_Control_LongestIsoform500lnThresh.fasta

- Minus SILVA rRNA:

146151_Trinity_Control_LongIso-500_MinusSilva.fasta

- Edited to remove >> that was introduced:

146151_Trinity_Control_LongIso-500_MinusSilva_edited.fasta

- Renamed:

146151_Trinity_Control_LongIso-500_Minus_Silva.fasta


## Host

- Blasted against Dirty Coral database, minus Clean Symbiont database:

26352_Trinity_Control_LongIso500_MinusSilva_DCminusCS.fasta

- Dirty coral database checked:

23398_Sid_DCcleaned.fasta 

- Cross-checked against Clean Coral database:

23842_Sid_DCcleaned_CCchecked.fasta

- Still in Blast format:

19272_Sid_GoodCoral.fasta

- Edited to contig name and sequence (removed extra Blast info):

19272_Sid_GoodCoral_Final.fasta

- Renamed to contig names that are more informative:

19272_Sid_GoodCoral_Final_renamed.fasta

- Rechecked 500 length threshold because sequences from DC/CC databases may have reintroduced smaller contigs (removed small contigs):

19222_Sid_GoodCoral_500lnThresh_Final.fasta



## Symbiont

- Blast against Dirty Symbiont database:

7301_SidSymb_DScleaned.txt

- Clean Symbiont database checked:

7651_SidSymb_DScleaned_CSchecked.fasta

- 3081 contigs:

SidSymb_GoodSymb.fasta

- Removed extra blast stuff:

3081_Sid_GoodSymb_Final.fasta

- Renamed to useful contig names:

3081_Sid_GoodSymb_Final_renamed.fasta
