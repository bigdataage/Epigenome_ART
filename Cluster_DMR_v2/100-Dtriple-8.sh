mkdir  -p   100-Dtriple-8
danpos.py     dtriple       \
2-BED/8_Y1M-E24-M-Mother-H3K4me1_Rep1.bed,2-BED/8_Y2M-E24-M-Mother-H3K4me3_Rep1.bed,\
2-BED/8_Y3M-E24-M-Mother-H3K27me3_Rep1.bed,2-BED/8_Y4M-E24-M-Mother-H3K27ac_Rep1.bed     \
--bg   2-BED/8_Input-Y5M-E24-M-Mother_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-8   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-8/100-Dtriple-8.runLog  2>&1



