mkdir  -p   100-Dtriple-1
danpos.py     dtriple       \
2-BED/1_P1C-E113-C-Boy-H3K4me3_Rep1.bed,2-BED/1_P2C-E113-C-Boy-H3K27me3_Rep1.bed,\
2-BED/1_P3C-E113-C-Boy-H3K27ac_Rep1.bed,2-BED/1_P4C-E113-C-Boy-H3K4me1_Rep1.bed     \
--bg   2-BED/1_Input-P5C-E113-C-Boy_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-1   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-1/100-Dtriple-1.runLog  2>&1



