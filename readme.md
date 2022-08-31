# The introduction of CTV-RPCA

## The simulation experiments are placed in folder "Simulation"
    
    ### Run "test_demo.m" to get metrics
      The mean value of the 50 results
          rank_ratio  sparse_ratio  error of CTV-RPCA     error of RPCA
            0.05          0.05          0.001142            0.002047
            0.10          0.10          0.000593            0.000893
            0.30          0.10          0.011400            0.018900
            0.10          0.30          0.008380            0.021556
    
    ### run "lambda_selection.m" to get visualizations
        we can verify that lambda= 3/sqrt(M) is reasonale, 
         where M,N are size of data, and M > N. 

## Real experiments

   ### Run "Demo_of_HSI_denoising.m" to get metrics and visualizations
        one result of small Pure DC mall dataset is:

         gaussian_level sparse_level  psnr/ssim of CTV   psnr/ssim of RPCA
              0.10             0.10       32.18/0.9697     30.44/0.9691 
              0.00             0.20       48.46/0.9923     43.54/0.9983
              0.20             0.00       27.82/0.9317     26.72/0.9335

   
   ### Run "Demo_of_Video_Extraction.m" to get metrics and visualizations
         result of airport dataset is:

                    auc of CTV   auc of RPCA
                      0.9156         0.8721

## If you have any questions, please contact the author, 
      the author's email address is：andrew.pengjj@gmail.com.


        
