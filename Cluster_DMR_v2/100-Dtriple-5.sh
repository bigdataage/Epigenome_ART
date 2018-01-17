mkdir  -p   100-Dtriple-5
danpos.py     dtriple       \
2-BED/5_Y1C-E24-C-Boy-H3K4me1_Rep1.bed,2-BED/5_Y2C-E24-C-Boy-H3K4me3_Rep1.bed,\
2-BED/5_Y3C-E24-C-Boy-H3K27me3_Rep1.bed,2-BED/5_Y4C-E24-C-Boy-H3K27ac_Rep1.bed     \
--bg   2-BED/5_Input-Y5C-E24-C-Boy_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-5   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-5/100-Dtriple-5.runLog  2>&1



