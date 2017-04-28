# all values are not filtered if at that value

INFO:
  exon: true
  exon_flag:
  - 'splicing'
  - 'exonic'
  - 'exonic,splicing'
  max_AC: 30
  max_1KG: 0.01
  max_gnomAD: 0.01
  max_ExAC: 0.01
  excluded_gene:
  - MUC
  - HLA
  excluded_chrom:
  - GL
  max_seqdup: 0.95
  #min_map: 1
  #rmsk: true

FILTER:
  VQSRSNP: 99.8
  VQSRINDEL: 99.7

READS:
  # proband filter
  min_proband_AD: 5
  min_proband_PL: 60
  min_proband_alt_freq: 0.1 #0.2 # AD < 10
  # parents filter
  #min_parents_DP: 10
  #min_parents_GQ: 30
  #min_parents_ref_freq: 0.95

SNP:
  max_FS: 25
  min_QD: 2

INDEL:
  max_FS: 25
  min_QD: 1
  min_ReadPosRankSum: -3

