mkdir  -p   100-Dtriple-6
danpos.py     dtriple       \
2-BED/6_Y1D-E24-D-Boy-H3K4me1_Rep1.bed,2-BED/6_Y2D-E24-D-Boy-H3K4me3_Rep1.bed,\
2-BED/6_Y3D-E24-D-Boy-H3K27me3_Rep1.bed,2-BED/6_Y4D-E24-D-Boy-H3K27ac_Rep1.bed     \
--bg   2-BED/6_Input-Y5D-E24-D-Boy_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-6   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-6/100-Dtriple-6.runLog  2>&1



