---
title: |
  ![](../aux_files/forKnitting/logo_blk.png){width=5in}  
  QC report for sequencing run 2022-09-08_COVIDSeq
author: "Philadelphia Public Health Labs"
date: \today
geometry: margin=3cm
output: 
    pdf_document:
        template: ../aux_files/forKnitting/toc_after.tex
        keep_tex: false
        toc: true
        toc_depth: 3
        keep_md: true
        includes:
            in_header: ../aux_files/forKnitting/TeX_packages_commands.sty

---

\tableofcontents



<!-- ========================================================== -->
<!--   Beginning of Preamble : Preamble seldom requires change  -->
<!-- ========================================================== -->

































































\newpage

# Introduction

This report summarizes the results of 51 COVID samples sequenced on 2022-09-08 and analyzed with BaseSpace DRAGEN Lineage App version 3.5.11 on 2022-09-09. The results section contains lineage information using the pangolin software version 4.1.2 (https://cov-lineages.org/resources/pangolin.html) and Nextclade version 2.5.0 (https://clades.nextstrain.org/)

# Number of samples sequenced

\begin{table}[H]
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{lr}
\toprule
\textcolor[HTML]{7C0A02}{\textbf{Sample Type}} & \textcolor[HTML]{7C0A02}{\textbf{n}}\\
\midrule
Empty well control & 1\\
\cmidrule{1-2}
Mock DNA positive control & 1\\
\cmidrule{1-2}
Nasal swab & 27\\
\cmidrule{1-2}
Nasal swab - No Variant Calls & 22\\
\bottomrule
\end{tabular}}
\end{table}

\newpage

# FASTQ generation, demultiplexing, and quality control

## GenerateFASTQ

The MiSeq machine uses the GenerateFASTQ workflow to demultiplex the samples and generate FASTQ files from BCL files (raw data files that contain sequence information and quality scores). GenerateFASTQ version 2.5.56.9 was used to generate the FASTQ files

## Number of read pairs per sample after demultiplexing

Samples were sequenced on the MiSeq as a 76:76 length run and demultiplexed on BaseSpace by the GenerateFASTQ module. The demultiplexing step involves matching the IDT barcode sequence associated with each sample to the index sequence of each read

### Number of samples and their pass filter percentages

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/percent of reads that passed filter-1.pdf)<!-- --> 

### Distributions of reads across all samples

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/reads_histogram-1.pdf)<!-- --> 

### Number of reads per sample type

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/reads per sample type-1.pdf)<!-- --> 

### Percentage of reads that make up the lane

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/percent reads of lane-1.pdf)<!-- --> 

## Top 10 unassigned indices

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/unassigned indices-1.pdf)<!-- --> 

These are the adapters: \newline
adapter: CTGTCTCTTATACACATCT

\blandscape

## Read counts and the DNA concentration/viral values of each sample



![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/reads and RLU-1.pdf)<!-- --> 

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/reads and DNA-1.pdf)<!-- --> 



![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/RLU and DNA-1.pdf)<!-- --> 

\newpage

## DNA concentration and read counts values per plate

The samples are processed on a 96-well plate. These figures illustrate the viral values, DNA concentration, and read counts of each sample on a well. Wells without any samples are white while samples with low values begin with dark blue. 



![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/RLU on plate-1.pdf)<!-- --> 

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/DNA on plate-1.pdf)<!-- --> 

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/reads on plate-1.pdf)<!-- --> 

\elandscape

# Results section

DRAGEN COVID Lineage App (on BaseSpace cloud) version 3.5.11 \newline
Lineage assignment software and version PUSHER-v1.14 \newline
Scorpio software version 0.3.17 \newline
Constellation database version 0.1.10

## Coverage

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/coverage results-1.pdf)<!-- --> 

\newpage

## Kmer results (number of unique SARS-CoV-2 fragments detected)

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/kmer numbers-1.pdf)<!-- --> 

## Kmer results (percentage of reference SARS-CoV-2 fragments detected)

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/COVID detection-1.pdf)<!-- --> 

## Pangolin Lineage Results

### Variant assignments from pangolin

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/variant assignments-1.pdf)<!-- --> 

\newpage

### Figure showing different pUSHER placements

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/multiple variant assignments-1.pdf)<!-- --> 

\newpage

### Figure showing conflicting variant assignments from pangolin, scorpio, and nextclade
Stars in boxes mean that the variants/sub-variants from different software agree

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/conflicting assignments-1.pdf)<!-- --> 







## Nextclade

### Coverage and QC scores

![](C:\Users\vincent.tu\OneDrive - City of Philadelphia\Sequencing_results\COVIDSeq_Runs_VT1\2022-09-08_COVIDSeq\output\2022-09-08_COVIDSeq.QC.report_gen.on.2022-09-09_files/figure-latex/nextclade coverage-1.pdf)<!-- --> 

## Summary of variants
Percentages of variants only in samples that meet criteria (Pangolin results)

\begin{table}[H]
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{llll}
\toprule
\textcolor[HTML]{7C0A02}{\textbf{Sample Type}} & \textcolor[HTML]{7C0A02}{\textbf{WHO}} & \textcolor[HTML]{7C0A02}{\textbf{Variant}} & \textcolor[HTML]{7C0A02}{\textbf{\makecell[l]{Variant\\percentage}}}\\
\midrule
 &  & BA.5.2.1 & 31.25\%\\
\cmidrule{3-4}
 &  & BA.5.2 & 18.75\%\\
\cmidrule{3-4}
 &  & BA.5.1 & 12.5\%\\
\cmidrule{3-4}
 &  & BE.3 & 12.5\%\\
\cmidrule{3-4}
 &  & BA.4.6 & 6.25\%\\
\cmidrule{3-4}
 &  & BA.5.1.2 & 6.25\%\\
\cmidrule{3-4}
 &  & BA.5.1.5 & 6.25\%\\
\cmidrule{3-4}
\multirow[t]{-8}{*}{\raggedright\arraybackslash Nasal swab} & \multirow[t]{-8}{*}{\raggedright\arraybackslash Omicron} & BF.5 & 6.25\%\\
\bottomrule
\end{tabular}}
\end{table}




\newpage

# Appendix

## Number of reads before and after filtering
Red samples refer to results that will not be reported (controls or no lineage results) \newline
Orange samples also refer to results that will not be reported (bad or mediocre NextClade QC with <80% 30X coverage) \newline
The total number of samples with acceptable results for reporting is 16 out of 49 nasal samples sequenced. \newline
Number of samples to send to epidemiologists: 16 \newline
Number of samples to upload to GISAID database: 16 \newline
Number of samples to upload to SRA and number of BioSamples to create: 16 and 16\newline

\begingroup\fontsize{10}{12}\selectfont

\begin{longtable}{llrr}
\toprule
\textcolor[HTML]{7C0A02}{\textbf{sample\_id}} & \textcolor[HTML]{7C0A02}{\textbf{sample\_type}} & \textcolor[HTML]{7C0A02}{\textbf{total\_raw\_reads}} & \textcolor[HTML]{7C0A02}{\textbf{read\_counts}}\\
\midrule
\endfirsthead
\multicolumn{4}{@{}l}{\textit{(continued)}}\\
\toprule
\textcolor[HTML]{7C0A02}{\textbf{sample\_id}} & \textcolor[HTML]{7C0A02}{\textbf{sample\_type}} & \textcolor[HTML]{7C0A02}{\textbf{total\_raw\_reads}} & \textcolor[HTML]{7C0A02}{\textbf{read\_counts}}\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\textcolor{black}{PHL2-A-A05-20220908} & Nasal swab & 1164874 & 1085095\\
\textcolor{black}{PHL2-A-B03-20220908} & Nasal swab & 978255 & 925715\\
\textcolor{red}{Undetermined} & Unassigned reads & 1454095 & 891127\\
\textcolor{black}{PHL2-A-H03-20220908} & Nasal swab & 918796 & 862937\\
\textcolor{black}{PHL2-A-H05-20220908} & Nasal swab & 865022 & 813917\\
\addlinespace
\textcolor{black}{PHL2-A-E04-20220908} & Nasal swab & 851382 & 782027\\
\textcolor{red}{PHL2-A-A01-20220908} & Nasal swab - No Variant Calls & 831913 & 780357\\
\textcolor{black}{PHL2-A-G04-20220908} & Nasal swab & 834772 & 756713\\
\textcolor{black}{PHL2-A-A06-20220908} & Nasal swab & 789827 & 740897\\
\textcolor{black}{PHL2-A-G03-20220908} & Nasal swab & 770550 & 721545\\
\addlinespace
\textcolor{black}{PHL2-A-C05-20220908} & Nasal swab & 730723 & 689622\\
\textcolor{red}{PHL2-A-C07-20220908} & Mock DNA positive control & 728519 & 682722\\
\textcolor{black}{PHL2-A-C06-20220908} & Nasal swab & 640871 & 582877\\
\textcolor{orange}{PHL2-A-C04-20220908} & Nasal swab & 621970 & 581416\\
\textcolor{black}{PHL2-A-D06-20220908} & Nasal swab & 640866 & 573002\\
\addlinespace
\textcolor{black}{PHL2-A-E05-20220908} & Nasal swab & 605526 & 557386\\
\textcolor{black}{PHL2-A-H06-20220908} & Nasal swab & 584061 & 549577\\
\textcolor{black}{PHL2-A-G05-20220908} & Nasal swab & 471246 & 443454\\
\textcolor{black}{PHL2-A-D02-20220908} & Nasal swab & 471037 & 432709\\
\textcolor{black}{PHL2-A-F02-20220908} & Nasal swab & 452315 & 421010\\
\addlinespace
\textcolor{red}{PHL2-A-E06-20220908} & Nasal swab - No Variant Calls & 450458 & 414410\\
\textcolor{orange}{PHL2-A-F06-20220908} & Nasal swab & 415403 & 387193\\
\textcolor{red}{PHL2-A-G06-20220908} & Nasal swab - No Variant Calls & 405810 & 374827\\
\textcolor{red}{PHL2-A-D05-20220908} & Nasal swab - No Variant Calls & 400610 & 372292\\
\textcolor{orange}{PHL2-A-F04-20220908} & Nasal swab & 398561 & 370833\\
\addlinespace
\textcolor{orange}{PHL2-A-A02-20220908} & Nasal swab & 395389 & 369785\\
\textcolor{orange}{PHL2-A-B04-20220908} & Nasal swab & 387143 & 359290\\
\textcolor{red}{PHL2-A-A07-20220908} & Nasal swab - No Variant Calls & 364437 & 341935\\
\textcolor{red}{PHL2-A-E02-20220908} & Nasal swab & 345556 & 313788\\
\textcolor{orange}{PHL2-A-F01-20220908} & Nasal swab & 297604 & 272000\\
\addlinespace
\textcolor{red}{PHL2-A-A03-20220908} & Nasal swab - No Variant Calls & 246787 & 227474\\
\textcolor{red}{PHL2-A-B01-20220908} & Nasal swab - No Variant Calls & 237760 & 222593\\
\textcolor{orange}{PHL2-A-B02-20220908} & Nasal swab & 234490 & 220107\\
\textcolor{red}{PHL2-A-H01-20220908} & Nasal swab - No Variant Calls & 220755 & 205585\\
\textcolor{red}{PHL2-A-E03-20220908} & Nasal swab - No Variant Calls & 228069 & 205547\\
\addlinespace
\textcolor{red}{PHL2-A-E01-20220908} & Nasal swab - No Variant Calls & 189590 & 176009\\
\textcolor{red}{PHL2-A-B06-20220908} & Nasal swab - No Variant Calls & 188554 & 174513\\
\textcolor{red}{PHL2-A-F05-20220908} & Nasal swab - No Variant Calls & 183459 & 165719\\
\textcolor{red}{PHL2-A-H04-20220908} & Nasal swab - No Variant Calls & 176234 & 164739\\
\textcolor{red}{PHL2-A-C01-20220908} & Nasal swab - No Variant Calls & 166918 & 153996\\
\addlinespace
\textcolor{orange}{PHL2-A-C03-20220908} & Nasal swab & 155375 & 142075\\
\textcolor{red}{PHL2-A-C02-20220908} & Nasal swab - No Variant Calls & 141681 & 129440\\
\textcolor{red}{PHL2-A-F03-20220908} & Nasal swab - No Variant Calls & 119205 & 111603\\
\textcolor{red}{PHL2-A-D03-20220908} & Nasal swab & 119910 & 106242\\
\textcolor{red}{PHL2-A-D04-20220908} & Nasal swab - No Variant Calls & 98436 & 90467\\
\addlinespace
\textcolor{red}{PHL2-A-A04-20220908} & Nasal swab & 92425 & 84866\\
\textcolor{red}{PHL2-A-G02-20220908} & Nasal swab - No Variant Calls & 88807 & 81505\\
\textcolor{red}{PHL2-A-G01-20220908} & Nasal swab - No Variant Calls & 69942 & 64474\\
\textcolor{red}{PHL2-A-D01-20220908} & Nasal swab - No Variant Calls & 70349 & 63880\\
\textcolor{red}{PHL2-A-B05-20220908} & Nasal swab - No Variant Calls & 48851 & 44845\\
\addlinespace
\textcolor{red}{PHL2-A-B07-20220908} & Empty well control & 41223 & 37426\\
\textcolor{red}{PHL2-A-H02-20220908} & Nasal swab - No Variant Calls & 36429 & 32976\\*
\end{longtable}
\endgroup{}

\blandscape

# Review and Approval

\begin{table}[H]
\centering
\resizebox{\linewidth}{!}{
\begin{tabular}{l>{}ll}
\toprule
\textcolor[HTML]{7C0A02}{\textbf{Role and Name}} & \textcolor[HTML]{7C0A02}{\textbf{Signature}} & \textcolor[HTML]{7C0A02}{\textbf{Date}}\\
\midrule
\makecell[l]{Bioinformatician:\\Vincent Tu, PhD} & \includegraphics[width=1.67in, height=0.6in]{C:/Users/vincent.tu/OneDrive - City of Philadelphia/Sequencing_results/COVIDSeq_Runs_VT1/aux_files/signature/VT.png} & 2022-09-09\\
\cmidrule{1-3}
\makecell[l]{Sequencing Manager:\\Mazen Sid Ahmed, PhD} & \includegraphics[width=1.67in, height=0.6in]{C:/Users/vincent.tu/OneDrive - City of Philadelphia/Sequencing_results/COVIDSeq_Runs_VT1/aux_files/signature/blank.png} & \\
\cmidrule{1-3}
\makecell[l]{Laboratory Directory:\\Bernadette Matthis, MSBA} & \includegraphics[width=1.67in, height=0.6in]{C:/Users/vincent.tu/OneDrive - City of Philadelphia/Sequencing_results/COVIDSeq_Runs_VT1/aux_files/signature/blank.png} & \\
\bottomrule
\end{tabular}}
\end{table}

\elandscape




