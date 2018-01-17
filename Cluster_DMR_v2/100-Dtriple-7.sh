mkdir  -p   100-Dtriple-7
danpos.py     dtriple       \
2-BED/7_Y1F-E24-F-Father-H3K4me1_Rep1.bed,2-BED/7_Y3F-E24-F-Father-H3K27me3_Rep1.bed,\
2-BED/7_Y2F-E24-F-Father-H3K4me3_Rep1.bed,2-BED/7_Y4F-E24-F-Father-H3K27ac_Rep1.bed     \
--bg   2-BED/7_Input-Y5F-E24-F-Father_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-7   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-7/100-Dtriple-7.runLog  2>&1



