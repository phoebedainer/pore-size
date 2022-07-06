clc
clear all
close all

%THIS PROGRAM COMPARES 2 MATRICES. Returns 1 for each value of a third
%matrix if each number is within a certain number from the other.

%PRECONDITIONS: Vertical, 1 column any number rows. Each matrix is same
%size. 

%VARIABLES TO CHANGE
m1dilation = 1;
m2dilation = (1.0e+03); 

within = 500; %# u want both to be within of the other

%insert matrix1 here
matrix1 = [           8
         162
         394
         518
         406
         698
         676
         703
         743
         116
         651
         659
         669
         938
         468
         137
         258
         277
         584
         574
         358
         371
         632
         253
         554
         550
         345
         500
         215
          14
         182
         392
         519
         690
         542
         178
         570
        1205
         193
         413
         306
         561
         875
         684
          79
          63
         224
          22
         312
         170
         221
         137
         376
          74
         560
         783
         362
         134
         413
         180
         553
        1442
         162
         414
         251
          83
         438
         709
         441
         202
         931
         727
         567
         383
         378
         346
         264
         578
         542
         489
         398
         553
         271
         460
         637
          30
         847
         873
         156
         506
          84
         855
         253
         183
          59
          43
         300
         448
         549
         261
          35
         767
          51
         493
         305
         490
         961
         574
         145
         177
         292
        1046
          41
         429
         478
         238
         397
         391
          84
          25
         325];

%insert matrix2 here
matrix2 = [         0
    0.2600
    0.4483
    0.7465
    0.4325
    0.5193
    0.7386
    0.3590
    0.7895
    0.0031
    0.6079
    0.7180
    0.9343
    1.0467
    0.4172
         0
    0.2833
    0.1874
    0.3538
    0.5942
    0.3707
    0.6749
    0.4365
    0.2969
    0.6621
    0.6126
    0.0222
    0.5024
    0.0798
         0
    0.1613
    0.4524
    0.4640
    0.7358
    0.4669
    0.2406
    0.6067
    1.0642
    0.1570
    0.3821
    0.2833
    0.4255
    0.7832
    0.6907
         0
         0
    0.2003
         0
    0.2026
    0.1342
    0.0645
         0
    0.4411
         0
    0.3988
    0.5614
    0.3811
         0
    0.2714
         0
    0.6348
    1.6664
    0.1871
    0.5402
    0.2415
         0
    0.4956
    0.7701
    0.4222
    0.2318
    0.8238
    1.0284
    0.6091
    0.6061
    0.3971
    0.3180
    0.2801
    0.5283
    0.2482
    0.6151
    0.4531
    0.5757
    0.2938
    0.3398
    0.6468
         0
    0.8537
    0.7162
    0.2691
    0.8622
         0
    0.8404
    0.0967
         0
         0
         0
    0.4552
    0.5229
    0.5941
    0.2097
         0
    1.0332
         0
    0.5489
    0.3746
    0.4122
    0.9472
    0.6163
         0
    0.0387
    0.2299
    1.1862
         0
    0.3432
    0.4613
    0.1103
    0.4796
    0.3051
         0
         0
    0.2950];

matrix1 = matrix1 .* m1dilation;
matrix2 = matrix2 .* m2dilation;

sizem1 = size(matrix1,1);
sizem2 = size(matrix2,1);

matrix3 = zeros(sizem1, 1);

if(sizem1 == sizem2)
    
    for i=1:sizem1
        
        minus = matrix1(i,1) - matrix2(i,1); 
        
        if(minus <= within && minus >= -within)
            matrix3(i,1) = 1;
        else
            matrix3(i,1) = 0;
        end
    
    end

end

disp(matrix3);