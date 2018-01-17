mkdir  -p   100-Dtriple-3
danpos.py     dtriple       \
2-BED/3_P1F-E113-F-Father-H3K4me3_Rep1.bed,2-BED/3_P2F-E113-F-Father-H3K27me3_Rep1.bed,\
2-BED/3_P3F-E113-F-Father-H3K27ac_Rep1.bed,2-BED/3_P4F-E113-F-Father-H3K4me1_Rep1.bed     \
--bg   2-BED/3_Input-P5F-E113-F-Father_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-3   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-3/100-Dtriple-3.runLog  2>&1



