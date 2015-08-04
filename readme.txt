==========================================================================================
This software can be used freely for research purposes.
Published reports of research  using this code (or a modified version) 
should cite the article that describes the algorithm: 

X. Dong, J. Shen, L. Shao, and M. H. Yang, Interactive co-segmentation 
using global and local energy optimization, IEEE Trans. on Image Processing, 
24(11):3966-3977, 2015 

This code has been compiled and tested using matlab R2011b 

version 1.0, 
rewritten and updated by Xingping Dong,  Aug. 2, 2015  

Email:  shenjianbing@bit.edu.cn / dongxingping@bit.edu.cn
============================================================================================
How to use
============================================================================================

1. There is one demo function. You can directly run it. 
 
  demo_lsr_iCoseg_hist.m
       
2. The results are shown in '.\results'.

3. You can add some scribbles into the images in '.\Datasets\scribbles\cow'to 		
   improve the performance. 

4. You can also add new images group into '.\Datasets\images' and corresponding 		
   scribbled images into '.\Datasets\scribbles'. And then you should modify the code: 
   dataset = 'cow\'; in demo_lsr_iCoseg_hist.m


============================================================================================
X. Dong, J. Shen, L. Shao, and M. H. Yang, Interactive co-segmentation 
using global and local energy optimization, IEEE Trans. on Image Processing, 
24(11):3966-3977, 2015 