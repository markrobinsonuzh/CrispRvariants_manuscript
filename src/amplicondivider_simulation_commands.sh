cd ../ampliconDIVider-master
source ampliconDIV_minimal.sh
samtools view -hb ../simulation/nlgn4a_0mut_300wt_100offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_0mut_300wt_100offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn4a_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_0mut_300wt_100offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_0mut_300wt_100offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap4_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_0mut_300wt_100offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_0mut_300wt_100offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn1_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_0mut_300wt_100offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_0mut_300wt_100offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_0mut_300wt_100offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_0mut_300wt_100offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_0mut_300wt_100offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_0mut_300wt_100offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_0mut_300wt_100offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_0mut_300wt_100offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_0mut_300wt_100offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_0mut_300wt_100offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_0mut_300wt_100offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_0mut_300wt_100offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn3a_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_0mut_300wt_100offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_0mut_300wt_100offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/cnsta_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_0mut_300wt_100offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_0mut_300wt_100offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_0mut_300wt_100offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_0mut_300wt_100offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap2b_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_0mut_300wt_100offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_0mut_300wt_100offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_0mut_300wt_100offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_0mut_300wt_100offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_0mut_300wt_100offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_0mut_300wt_100offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1a_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_0mut_300wt_100offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_0mut_300wt_100offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_0mut_300wt_100offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_0mut_300wt_100offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_0mut_300wt_100offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_0mut_300wt_100offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1b_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_0mut_300wt_100offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_0mut_300wt_100offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/gjd1a_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_0mut_300wt_100offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_0mut_300wt_100offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_0mut_300wt_100offtarget_200readlen ../simulation/amplicondivider/cntn3b_0mut_300wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_100mut_200wt_100offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_100mut_200wt_100offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn4a_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_100mut_200wt_100offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_100mut_200wt_100offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap4_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_100mut_200wt_100offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_100mut_200wt_100offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn1_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_100mut_200wt_100offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_100mut_200wt_100offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_100mut_200wt_100offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_100mut_200wt_100offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_100mut_200wt_100offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_100mut_200wt_100offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_100mut_200wt_100offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_100mut_200wt_100offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_100mut_200wt_100offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_100mut_200wt_100offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_100mut_200wt_100offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_100mut_200wt_100offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn3a_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_100mut_200wt_100offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_100mut_200wt_100offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/cnsta_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_100mut_200wt_100offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_100mut_200wt_100offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_100mut_200wt_100offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_100mut_200wt_100offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap2b_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_100mut_200wt_100offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_100mut_200wt_100offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_100mut_200wt_100offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_100mut_200wt_100offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_100mut_200wt_100offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_100mut_200wt_100offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1a_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_100mut_200wt_100offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_100mut_200wt_100offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_100mut_200wt_100offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_100mut_200wt_100offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_100mut_200wt_100offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_100mut_200wt_100offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1b_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_100mut_200wt_100offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_100mut_200wt_100offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/gjd1a_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_100mut_200wt_100offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_100mut_200wt_100offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_100mut_200wt_100offtarget_200readlen ../simulation/amplicondivider/cntn3b_100mut_200wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_200mut_100wt_100offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_200mut_100wt_100offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn4a_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_200mut_100wt_100offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_200mut_100wt_100offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap4_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_200mut_100wt_100offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_200mut_100wt_100offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn1_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_200mut_100wt_100offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_200mut_100wt_100offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_200mut_100wt_100offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_200mut_100wt_100offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_200mut_100wt_100offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_200mut_100wt_100offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_200mut_100wt_100offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_200mut_100wt_100offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_200mut_100wt_100offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_200mut_100wt_100offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_200mut_100wt_100offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_200mut_100wt_100offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn3a_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_200mut_100wt_100offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_200mut_100wt_100offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/cnsta_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_200mut_100wt_100offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_200mut_100wt_100offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_200mut_100wt_100offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_200mut_100wt_100offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap2b_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_200mut_100wt_100offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_200mut_100wt_100offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_200mut_100wt_100offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_200mut_100wt_100offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_200mut_100wt_100offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_200mut_100wt_100offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1a_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_200mut_100wt_100offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_200mut_100wt_100offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_200mut_100wt_100offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_200mut_100wt_100offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_200mut_100wt_100offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_200mut_100wt_100offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1b_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_200mut_100wt_100offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_200mut_100wt_100offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/gjd1a_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_200mut_100wt_100offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_200mut_100wt_100offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_200mut_100wt_100offtarget_200readlen ../simulation/amplicondivider/cntn3b_200mut_100wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_270mut_30wt_100offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_270mut_30wt_100offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn4a_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_270mut_30wt_100offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_270mut_30wt_100offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap4_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_270mut_30wt_100offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_270mut_30wt_100offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn1_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_270mut_30wt_100offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_270mut_30wt_100offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_270mut_30wt_100offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_270mut_30wt_100offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_270mut_30wt_100offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_270mut_30wt_100offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_270mut_30wt_100offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_270mut_30wt_100offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_270mut_30wt_100offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_270mut_30wt_100offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_270mut_30wt_100offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_270mut_30wt_100offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nlgn3a_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_270mut_30wt_100offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_270mut_30wt_100offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/cnsta_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_270mut_30wt_100offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_270mut_30wt_100offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_270mut_30wt_100offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_270mut_30wt_100offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/cntnap2b_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_270mut_30wt_100offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_270mut_30wt_100offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_270mut_30wt_100offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_270mut_30wt_100offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_270mut_30wt_100offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_270mut_30wt_100offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1a_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_270mut_30wt_100offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_270mut_30wt_100offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_270mut_30wt_100offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_270mut_30wt_100offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_270mut_30wt_100offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_270mut_30wt_100offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/cntn1b_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_270mut_30wt_100offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_270mut_30wt_100offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/gjd1a_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_270mut_30wt_100offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_270mut_30wt_100offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_270mut_30wt_100offtarget_200readlen ../simulation/amplicondivider/cntn3b_270mut_30wt_100offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_0mut_300wt_33offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_0mut_300wt_33offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn4a_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_0mut_300wt_33offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_0mut_300wt_33offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap4_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_0mut_300wt_33offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_0mut_300wt_33offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn1_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_0mut_300wt_33offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_0mut_300wt_33offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_0mut_300wt_33offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_0mut_300wt_33offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_0mut_300wt_33offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_0mut_300wt_33offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_0mut_300wt_33offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_0mut_300wt_33offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_0mut_300wt_33offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_0mut_300wt_33offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_0mut_300wt_33offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_0mut_300wt_33offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn3a_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_0mut_300wt_33offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_0mut_300wt_33offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/cnsta_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_0mut_300wt_33offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_0mut_300wt_33offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_0mut_300wt_33offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_0mut_300wt_33offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap2b_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_0mut_300wt_33offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_0mut_300wt_33offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_0mut_300wt_33offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_0mut_300wt_33offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_0mut_300wt_33offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_0mut_300wt_33offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1a_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_0mut_300wt_33offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_0mut_300wt_33offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_0mut_300wt_33offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_0mut_300wt_33offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_0mut_300wt_33offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_0mut_300wt_33offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1b_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_0mut_300wt_33offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_0mut_300wt_33offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/gjd1a_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_0mut_300wt_33offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_0mut_300wt_33offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_0mut_300wt_33offtarget_200readlen ../simulation/amplicondivider/cntn3b_0mut_300wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_100mut_200wt_33offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_100mut_200wt_33offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn4a_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_100mut_200wt_33offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_100mut_200wt_33offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap4_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_100mut_200wt_33offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_100mut_200wt_33offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn1_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_100mut_200wt_33offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_100mut_200wt_33offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_100mut_200wt_33offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_100mut_200wt_33offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_100mut_200wt_33offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_100mut_200wt_33offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_100mut_200wt_33offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_100mut_200wt_33offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_100mut_200wt_33offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_100mut_200wt_33offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_100mut_200wt_33offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_100mut_200wt_33offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn3a_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_100mut_200wt_33offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_100mut_200wt_33offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/cnsta_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_100mut_200wt_33offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_100mut_200wt_33offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_100mut_200wt_33offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_100mut_200wt_33offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap2b_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_100mut_200wt_33offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_100mut_200wt_33offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_100mut_200wt_33offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_100mut_200wt_33offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_100mut_200wt_33offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_100mut_200wt_33offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1a_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_100mut_200wt_33offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_100mut_200wt_33offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_100mut_200wt_33offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_100mut_200wt_33offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_100mut_200wt_33offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_100mut_200wt_33offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1b_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_100mut_200wt_33offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_100mut_200wt_33offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/gjd1a_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_100mut_200wt_33offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_100mut_200wt_33offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_100mut_200wt_33offtarget_200readlen ../simulation/amplicondivider/cntn3b_100mut_200wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_200mut_100wt_33offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_200mut_100wt_33offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn4a_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_200mut_100wt_33offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_200mut_100wt_33offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap4_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_200mut_100wt_33offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_200mut_100wt_33offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn1_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_200mut_100wt_33offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_200mut_100wt_33offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_200mut_100wt_33offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_200mut_100wt_33offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_200mut_100wt_33offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_200mut_100wt_33offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_200mut_100wt_33offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_200mut_100wt_33offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_200mut_100wt_33offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_200mut_100wt_33offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_200mut_100wt_33offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_200mut_100wt_33offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn3a_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_200mut_100wt_33offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_200mut_100wt_33offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/cnsta_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_200mut_100wt_33offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_200mut_100wt_33offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_200mut_100wt_33offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_200mut_100wt_33offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap2b_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_200mut_100wt_33offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_200mut_100wt_33offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_200mut_100wt_33offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_200mut_100wt_33offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_200mut_100wt_33offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_200mut_100wt_33offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1a_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_200mut_100wt_33offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_200mut_100wt_33offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_200mut_100wt_33offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_200mut_100wt_33offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_200mut_100wt_33offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_200mut_100wt_33offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1b_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_200mut_100wt_33offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_200mut_100wt_33offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/gjd1a_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_200mut_100wt_33offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_200mut_100wt_33offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_200mut_100wt_33offtarget_200readlen ../simulation/amplicondivider/cntn3b_200mut_100wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_270mut_30wt_33offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_270mut_30wt_33offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn4a_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_270mut_30wt_33offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_270mut_30wt_33offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap4_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_270mut_30wt_33offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_270mut_30wt_33offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn1_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_270mut_30wt_33offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_270mut_30wt_33offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_270mut_30wt_33offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_270mut_30wt_33offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_270mut_30wt_33offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_270mut_30wt_33offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_270mut_30wt_33offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_270mut_30wt_33offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_270mut_30wt_33offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_270mut_30wt_33offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_270mut_30wt_33offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_270mut_30wt_33offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nlgn3a_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_270mut_30wt_33offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_270mut_30wt_33offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/cnsta_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_270mut_30wt_33offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_270mut_30wt_33offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_270mut_30wt_33offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_270mut_30wt_33offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/cntnap2b_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_270mut_30wt_33offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_270mut_30wt_33offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_270mut_30wt_33offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_270mut_30wt_33offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_270mut_30wt_33offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_270mut_30wt_33offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1a_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_270mut_30wt_33offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_270mut_30wt_33offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_270mut_30wt_33offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_270mut_30wt_33offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_270mut_30wt_33offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_270mut_30wt_33offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/cntn1b_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_270mut_30wt_33offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_270mut_30wt_33offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/gjd1a_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_270mut_30wt_33offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_270mut_30wt_33offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_270mut_30wt_33offtarget_200readlen ../simulation/amplicondivider/cntn3b_270mut_30wt_33offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_0mut_300wt_0offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_0mut_300wt_0offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn4a_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_0mut_300wt_0offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_0mut_300wt_0offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap4_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_0mut_300wt_0offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_0mut_300wt_0offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn1_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_0mut_300wt_0offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_0mut_300wt_0offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_0mut_300wt_0offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_0mut_300wt_0offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_0mut_300wt_0offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_0mut_300wt_0offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_0mut_300wt_0offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_0mut_300wt_0offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_0mut_300wt_0offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_0mut_300wt_0offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_0mut_300wt_0offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_0mut_300wt_0offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn3a_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_0mut_300wt_0offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_0mut_300wt_0offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/cnsta_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_0mut_300wt_0offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_0mut_300wt_0offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_0mut_300wt_0offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_0mut_300wt_0offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap2b_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_0mut_300wt_0offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_0mut_300wt_0offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_0mut_300wt_0offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_0mut_300wt_0offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_0mut_300wt_0offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_0mut_300wt_0offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1a_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_0mut_300wt_0offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_0mut_300wt_0offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_0mut_300wt_0offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_0mut_300wt_0offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_0mut_300wt_0offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_0mut_300wt_0offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1b_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_0mut_300wt_0offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_0mut_300wt_0offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/gjd1a_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_0mut_300wt_0offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_0mut_300wt_0offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_0mut_300wt_0offtarget_200readlen ../simulation/amplicondivider/cntn3b_0mut_300wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_100mut_200wt_0offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_100mut_200wt_0offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn4a_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_100mut_200wt_0offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_100mut_200wt_0offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap4_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_100mut_200wt_0offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_100mut_200wt_0offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn1_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_100mut_200wt_0offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_100mut_200wt_0offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_100mut_200wt_0offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_100mut_200wt_0offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_100mut_200wt_0offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_100mut_200wt_0offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_100mut_200wt_0offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_100mut_200wt_0offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_100mut_200wt_0offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_100mut_200wt_0offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_100mut_200wt_0offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_100mut_200wt_0offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn3a_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_100mut_200wt_0offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_100mut_200wt_0offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/cnsta_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_100mut_200wt_0offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_100mut_200wt_0offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_100mut_200wt_0offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_100mut_200wt_0offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap2b_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_100mut_200wt_0offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_100mut_200wt_0offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_100mut_200wt_0offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_100mut_200wt_0offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_100mut_200wt_0offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_100mut_200wt_0offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1a_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_100mut_200wt_0offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_100mut_200wt_0offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_100mut_200wt_0offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_100mut_200wt_0offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_100mut_200wt_0offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_100mut_200wt_0offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1b_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_100mut_200wt_0offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_100mut_200wt_0offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/gjd1a_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_100mut_200wt_0offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_100mut_200wt_0offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_100mut_200wt_0offtarget_200readlen ../simulation/amplicondivider/cntn3b_100mut_200wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_200mut_100wt_0offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_200mut_100wt_0offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn4a_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_200mut_100wt_0offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_200mut_100wt_0offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap4_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_200mut_100wt_0offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_200mut_100wt_0offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn1_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_200mut_100wt_0offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_200mut_100wt_0offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_200mut_100wt_0offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_200mut_100wt_0offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_200mut_100wt_0offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_200mut_100wt_0offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_200mut_100wt_0offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_200mut_100wt_0offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_200mut_100wt_0offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_200mut_100wt_0offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_200mut_100wt_0offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_200mut_100wt_0offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn3a_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_200mut_100wt_0offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_200mut_100wt_0offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/cnsta_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_200mut_100wt_0offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_200mut_100wt_0offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_200mut_100wt_0offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_200mut_100wt_0offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap2b_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_200mut_100wt_0offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_200mut_100wt_0offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_200mut_100wt_0offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_200mut_100wt_0offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_200mut_100wt_0offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_200mut_100wt_0offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1a_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_200mut_100wt_0offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_200mut_100wt_0offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_200mut_100wt_0offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_200mut_100wt_0offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_200mut_100wt_0offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_200mut_100wt_0offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1b_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_200mut_100wt_0offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_200mut_100wt_0offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/gjd1a_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_200mut_100wt_0offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_200mut_100wt_0offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_200mut_100wt_0offtarget_200readlen ../simulation/amplicondivider/cntn3b_200mut_100wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn4a_270mut_30wt_0offtarget_200readlen_merged.bam chr1:32080517-32080949 > temp.bam
parseBam temp.bam nlgn4a_270mut_30wt_0offtarget_200readlen 32080712 32080754; rm temp.bam
mv frameshift_summary_nlgn4a_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn4a_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap4_270mut_30wt_0offtarget_200readlen_merged.bam chr10:2281271-2281703 > temp.bam
parseBam temp.bam cntnap4_270mut_30wt_0offtarget_200readlen 2281466 2281508; rm temp.bam
mv frameshift_summary_cntnap4_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap4_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn1_270mut_30wt_0offtarget_200readlen_merged.bam chr11:10037445-10037877 > temp.bam
parseBam temp.bam nlgn1_270mut_30wt_0offtarget_200readlen 10037640 10037682; rm temp.bam
mv frameshift_summary_nlgn1_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn1_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/CNTNAP5_270mut_30wt_0offtarget_200readlen_merged.bam chr11:34313638-34314070 > temp.bam
parseBam temp.bam CNTNAP5_270mut_30wt_0offtarget_200readlen 34313833 34313875; rm temp.bam
mv frameshift_summary_CNTNAP5_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/CNTNAP5_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_L_270mut_30wt_0offtarget_200readlen_merged.bam chr12:25855514-25855946 > temp.bam
parseBam temp.bam nrxn1a_L_270mut_30wt_0offtarget_200readlen 25855709 25855751; rm temp.bam
mv frameshift_summary_nrxn1a_L_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_L_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1a_S_270mut_30wt_0offtarget_200readlen_merged.bam chr12:26069046-26069478 > temp.bam
parseBam temp.bam nrxn1a_S_270mut_30wt_0offtarget_200readlen 26069241 26069283; rm temp.bam
mv frameshift_summary_nrxn1a_S_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1a_S_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_L_270mut_30wt_0offtarget_200readlen_merged.bam chr13:605215-605647 > temp.bam
parseBam temp.bam nrxn1b_L_270mut_30wt_0offtarget_200readlen 605410 605452; rm temp.bam
mv frameshift_summary_nrxn1b_L_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_L_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn1b_S_270mut_30wt_0offtarget_200readlen_merged.bam chr13:645809-646241 > temp.bam
parseBam temp.bam nrxn1b_S_270mut_30wt_0offtarget_200readlen 646004 646046; rm temp.bam
mv frameshift_summary_nrxn1b_S_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn1b_S_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nlgn3a_270mut_30wt_0offtarget_200readlen_merged.bam chr14:18215057-18215489 > temp.bam
parseBam temp.bam nlgn3a_270mut_30wt_0offtarget_200readlen 18215252 18215294; rm temp.bam
mv frameshift_summary_nlgn3a_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nlgn3a_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cnsta_270mut_30wt_0offtarget_200readlen_merged.bam chr17:11951995-11952427 > temp.bam
parseBam temp.bam cnsta_270mut_30wt_0offtarget_200readlen 11952190 11952232; rm temp.bam
mv frameshift_summary_cnsta_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/cnsta_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3a_L_270mut_30wt_0offtarget_200readlen_merged.bam chr17:17271310-17271742 > temp.bam
parseBam temp.bam nrxn3a_L_270mut_30wt_0offtarget_200readlen 17271505 17271547; rm temp.bam
mv frameshift_summary_nrxn3a_L_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3a_L_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntnap2b_270mut_30wt_0offtarget_200readlen_merged.bam chr2:50866239-50866671 > temp.bam
parseBam temp.bam cntnap2b_270mut_30wt_0offtarget_200readlen 50866434 50866476; rm temp.bam
mv frameshift_summary_cntnap2b_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/cntnap2b_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_L_270mut_30wt_0offtarget_200readlen_merged.bam chr20:5595801-5596233 > temp.bam
parseBam temp.bam nrxn3b_L_270mut_30wt_0offtarget_200readlen 5595996 5596038; rm temp.bam
mv frameshift_summary_nrxn3b_L_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_L_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/nrxn3b_S_270mut_30wt_0offtarget_200readlen_merged.bam chr20:5854358-5854790 > temp.bam
parseBam temp.bam nrxn3b_S_270mut_30wt_0offtarget_200readlen 5854553 5854595; rm temp.bam
mv frameshift_summary_nrxn3b_S_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/nrxn3b_S_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1a_270mut_30wt_0offtarget_200readlen_merged.bam chr25:13719-14151 > temp.bam
parseBam temp.bam cntn1a_270mut_30wt_0offtarget_200readlen 13914 13956; rm temp.bam
mv frameshift_summary_cntn1a_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1a_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_L_270mut_30wt_0offtarget_200readlen_merged.bam chr25:33186887-33187319 > temp.bam
parseBam temp.bam tjp1b_L_270mut_30wt_0offtarget_200readlen 33187082 33187124; rm temp.bam
mv frameshift_summary_tjp1b_L_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_L_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/tjp1b_B_270mut_30wt_0offtarget_200readlen_merged.bam chr25:33212989-33213421 > temp.bam
parseBam temp.bam tjp1b_B_270mut_30wt_0offtarget_200readlen 33213184 33213226; rm temp.bam
mv frameshift_summary_tjp1b_B_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/tjp1b_B_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn1b_270mut_30wt_0offtarget_200readlen_merged.bam chr4:12806870-12807302 > temp.bam
parseBam temp.bam cntn1b_270mut_30wt_0offtarget_200readlen 12807065 12807107; rm temp.bam
mv frameshift_summary_cntn1b_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/cntn1b_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/gjd1a_270mut_30wt_0offtarget_200readlen_merged.bam chr5:38574931-38575363 > temp.bam
parseBam temp.bam gjd1a_270mut_30wt_0offtarget_200readlen 38575126 38575168; rm temp.bam
mv frameshift_summary_gjd1a_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/gjd1a_270mut_30wt_0offtarget_200readlen.txt

samtools view -hb ../simulation/cntn3b_270mut_30wt_0offtarget_200readlen_merged.bam chr6:44589143-44589575 > temp.bam
parseBam temp.bam cntn3b_270mut_30wt_0offtarget_200readlen 44589338 44589380; rm temp.bam
mv frameshift_summary_cntn3b_270mut_30wt_0offtarget_200readlen ../simulation/amplicondivider/cntn3b_270mut_30wt_0offtarget_200readlen.txt

