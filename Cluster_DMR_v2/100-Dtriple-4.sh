mkdir  -p   100-Dtriple-4
danpos.py     dtriple       \
2-BED/4_P1M-E113-M-Mother-H3K4me3_Rep1.bed,2-BED/4_P2M-E113-M-Mother-H3K27me3_Rep1.bed,\
2-BED/4_P3M-E113-M-Mother-H3K27ac_Rep1.bed,2-BED/4_P4M-E113-M-Mother-H3K4me1_Rep1.bed     \
--bg   2-BED/4_Input-P5M-E113-M-Mother_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-4   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-4/100-Dtriple-4.runLog  2>&1



