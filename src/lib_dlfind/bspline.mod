  �  A   k820309    &          9.1         &<V                                                                                                           
       dlf_util_d.F90 BSPLINE              gen@SPLINE_INIT gen@SPLINE_CREATE gen@SPLINE_GET gen@SPLINE_DESTROY                                             
       RK                      @                      
       ALLOCATE DEALLOCATE          @  "   X                            u #SPLINE_INIT             @  "   X                            u #SPLINE_CREATE             @  "   X                            u #SPLINE_GET             @  "   X                           u #SPLINE_DESTROY                   �                              u #ALLOCATE_R1    #ALLOCATE_R2 
   #ALLOCATE_R3    #ALLOCATE_I1    #ALLOCATE_I2    #ALLOCATE_L1    #         @     @                                      #ARRAY    #UBOUND 	             @                                        
               &                                                     
   @                      	           #         @     @                  
                   #ARRAY    #UBOUND1    #UBOUND2              @                                        
               &                   &                                                     
   @                                           
   @                                 #         @     @                                      #ARRAY    #UBOUND1    #UBOUND2    #UBOUND3              @                                        
               &                   &                   &                                                     
   @                                           
   @                                           
   @                                 #         @     @                                      #ARRAY    #UBOUND              @                                          	              &                                                     
   @                                 #         @     @                                      #ARRAY    #UBOUND1    #UBOUND2              @                                          
              &                   &                                                     
   @                                           
   @                                 #         @     @                                     #ARRAY    #UBOUND              @                                                        &                                                     
   @                                                �                              u #DEALLOCATE_R1    #DEALLOCATE_R2     #DEALLOCATE_R3 #   #DEALLOCATE_I1 &   #DEALLOCATE_I2 )   #DEALLOCATE_L1 ,   #         @     @                                     #DEALLOCATE_R1%SIZE    #ARRAY                  @                         SIZE          @                                        
               &                                           #         @     @                                     #DEALLOCATE_R2%SIZE !   #ARRAY "                 @                    !     SIZE          @                     "                   
               &                   &                                           #         @     @                   #                  #DEALLOCATE_R3%SIZE $   #ARRAY %                 @                    $     SIZE          @                     %                   
               &                   &                   &                                           #         @     @                   &                  #DEALLOCATE_I1%SIZE '   #ARRAY (                 @                    '     SIZE          @                      (                                  &                                           #         @     @                   )                  #DEALLOCATE_I2%SIZE *   #ARRAY +                 @                    *     SIZE          @                      +                                  &                   &                                           #         @     @                  ,                  #DEALLOCATE_L1%SIZE -   #ARRAY .                 @                    -     SIZE          @                      .                                  &                                           #         @  "   X                                      #SPLINE_INIT%ALLOCATED /   #LENGTH_IN 0   #NFUNC_IN 1                 @                    /     ALLOCATED           
   @                      0                     
   @                      1           #         @  "   X                                       #IFUNC 2   #X_IN 3   #Y_IN 5             
   @                      2                    
   @                     3                    
    p          5 r 4       5 r 4                              
   @                     5                    
    p          5 r 4       5 r 4                     #         @  "   X                                       #IFUNC 6   #XVAL 7   #YVAL 8   #DYVAL 9             
   @                      6                     
   @                     7     
                D  @                     8     
                 D  @                     9     
       #         @  "   X                                                  @ @@                      4               �         fn#fn    �   T   b   uapp(BSPLINE %     ;   J  DLF_PARAMETER_MODULE    N  L   J  DLF_ALLOCATE     �  I       gen@SPLINE_INIT "   �  K       gen@SPLINE_CREATE    .  H       gen@SPLINE_GET #   v  L       gen@SPLINE_DESTROY *   �  �       gen@ALLOCATE+DLF_ALLOCATE )   `  W      ALLOCATE_R1+DLF_ALLOCATE /   �  �   e   ALLOCATE_R1%ARRAY+DLF_ALLOCATE 0   ;  8   e   ALLOCATE_R1%UBOUND+DLF_ALLOCATE )   s  e      ALLOCATE_R2+DLF_ALLOCATE /   �  �   e   ALLOCATE_R2%ARRAY+DLF_ALLOCATE 1   t  8   e   ALLOCATE_R2%UBOUND1+DLF_ALLOCATE 1   �  8   e   ALLOCATE_R2%UBOUND2+DLF_ALLOCATE )   �  r      ALLOCATE_R3+DLF_ALLOCATE /   V  �   e   ALLOCATE_R3%ARRAY+DLF_ALLOCATE 1   
  8   e   ALLOCATE_R3%UBOUND1+DLF_ALLOCATE 1   B  8   e   ALLOCATE_R3%UBOUND2+DLF_ALLOCATE 1   z  8   e   ALLOCATE_R3%UBOUND3+DLF_ALLOCATE )   �  W      ALLOCATE_I1+DLF_ALLOCATE /   	  �   e   ALLOCATE_I1%ARRAY+DLF_ALLOCATE 0   �  8   e   ALLOCATE_I1%UBOUND+DLF_ALLOCATE )   �  e      ALLOCATE_I2+DLF_ALLOCATE /   *	  �   e   ALLOCATE_I2%ARRAY+DLF_ALLOCATE 1   �	  8   e   ALLOCATE_I2%UBOUND1+DLF_ALLOCATE 1   �	  8   e   ALLOCATE_I2%UBOUND2+DLF_ALLOCATE )   6
  W      ALLOCATE_L1+DLF_ALLOCATE /   �
  �   e   ALLOCATE_L1%ARRAY+DLF_ALLOCATE 0     8   e   ALLOCATE_L1%UBOUND+DLF_ALLOCATE ,   I  �       gen@DEALLOCATE+DLF_ALLOCATE +   �  c      DEALLOCATE_R1+DLF_ALLOCATE 5   V  5      DEALLOCATE_R1%SIZE+DLF_ALLOCATE=SIZE 1   �  �   e   DEALLOCATE_R1%ARRAY+DLF_ALLOCATE +     c      DEALLOCATE_R2+DLF_ALLOCATE 5   r  5      DEALLOCATE_R2%SIZE+DLF_ALLOCATE=SIZE 1   �  �   e   DEALLOCATE_R2%ARRAY+DLF_ALLOCATE +   C  c      DEALLOCATE_R3+DLF_ALLOCATE 5   �  5      DEALLOCATE_R3%SIZE+DLF_ALLOCATE=SIZE 1   �  �   e   DEALLOCATE_R3%ARRAY+DLF_ALLOCATE +   �  c      DEALLOCATE_I1+DLF_ALLOCATE 5   �  5      DEALLOCATE_I1%SIZE+DLF_ALLOCATE=SIZE 1   '  �   e   DEALLOCATE_I1%ARRAY+DLF_ALLOCATE +   �  c      DEALLOCATE_I2+DLF_ALLOCATE 5     5      DEALLOCATE_I2%SIZE+DLF_ALLOCATE=SIZE 1   C  �   e   DEALLOCATE_I2%ARRAY+DLF_ALLOCATE +   �  c      DEALLOCATE_L1+DLF_ALLOCATE 5   B  5      DEALLOCATE_L1%SIZE+DLF_ALLOCATE=SIZE 1   w  �   e   DEALLOCATE_L1%ARRAY+DLF_ALLOCATE    �  x       SPLINE_INIT &   s  :      SPLINE_INIT%ALLOCATED &   �  8   a   SPLINE_INIT%LENGTH_IN %   �  8   a   SPLINE_INIT%NFUNC_IN      _       SPLINE_CREATE $   |  8   a   SPLINE_CREATE%IFUNC #   �  �   a   SPLINE_CREATE%X_IN #   @  �   a   SPLINE_CREATE%Y_IN    �  j       SPLINE_GET !   6  8   a   SPLINE_GET%IFUNC     n  8   a   SPLINE_GET%XVAL     �  8   a   SPLINE_GET%YVAL !   �  8   a   SPLINE_GET%DYVAL      @       SPLINE_DESTROY    V  8      LENGTH 