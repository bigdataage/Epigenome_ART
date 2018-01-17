mkdir  -p   100-Dtriple-2
danpos.py     dtriple       \
2-BED/2_P1D-E113-D-Boy-H3K4me3_Rep1.bed,2-BED/2_P2D-E113-D-Boy-H3K27me3_Rep1.bed,\
2-BED/2_P3D-E113-D-Boy-H3K27ac_Rep1.bed,2-BED/2_P4D-E113-D-Boy-H3K4me1_Rep1.bed     \
--bg   2-BED/2_Input-P5D-E113-D-Boy_Rep1.bed      \
--paired 0   --pheight 0.01  --height 0   --testcut 0    --out 100-Dtriple-2   --save 1     --edge 1     \
--count 10000000    --span 20    --clonalcut 1e-10       >>  100-Dtriple-2/100-Dtriple-2.runLog  2>&1



