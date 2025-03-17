#! /bin/bash
cd /public/agis/changyuxiao_group/chenhaotian/mywork/2024.11.14-rice_RIL/data_sequencing_test/RIL-ELF6/lib_data_split
for m in $(cat Lib_name.list)
do 
	cd ./${m}
	for n in $(cat ../enzyme_name.list)
	do bsub -J ${n}_grep8bp -o %J_${n}.out -e %J_${n}.err -n 12 -q arm  "bash /public/agis/changyuxiao_group/chenhaotian/script/seqkit_8bp_split_without_index.sh /public/agis/changyuxiao_group/chenhaotian/mywork/2024.11.14-rice_RIL/rice_RIL_12enzyme_digest/250-500bp/8bp_bed/8bp_fasta/without_pos/unique_barcode/${n}.8bp.fa ${m}_R1.fq.gz ${m}_R2.fq.gz ${n}"
	done
	cd /public/agis/changyuxiao_group/chenhaotian/mywork/2024.11.14-rice_RIL/data_sequencing_test/RIL-ELF6/lib_data_split
done
