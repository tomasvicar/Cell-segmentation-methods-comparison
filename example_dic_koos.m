clc;clear all;close all;
addpath('detection','evaluation','foreground_segmentation','final_segmentation','example_data')
addpath(genpath('util'))

%Composed by Tomas Vicar 24/05/2018, 
% Department of Biomedical Engineering, Brno University of Technology  
% vicar@feec.vutbr.cz



I=double(imread('example_data/DIC_02_rDIC-Koos.tif'));%load example data
gt_mask=imread('example_data/DIC_02_groundtruth_mask.png')==1;

rmin=20;%estimated celll radius range and area
rmax=56;
area_min=230;
min_and_max_of_img1=[0 2.7198]; % minimal and maximal intensity of first image (for normalization)








%% simple treshold
t=0.03;
I_norm=mat2gray(I,min_and_max_of_img1);
res_simple_treshold=I_norm>t;
figure;imshow(I,[]);hold on;
visboundaries(res_simple_treshold)
title('simple treshold')



%% otsu treshold 
I_norm=mat2gray(I);
t = graythresh(I_norm);
res_otsu_treshold=I_norm>t;
figure;imshow(I,[]);hold on;
visboundaries(res_otsu_treshold)
title('otsu treshold')



%% poisson treshold
I_norm=imadjust(I);%works better after elimination of extreme values
t=poisson_tresh(I_norm);
res_poisson_treshold=I_norm>t;
figure;imshow(I,[]);hold on;
visboundaries(res_poisson_treshold)
title('poisson treshold')



%% EGT
min_object_size=area_min;
min_hole_size=322;
treshold_finetune=11.1;
res_EGT = EGT_Segmentation(I, min_object_size,min_hole_size,treshold_finetune);
figure;imshow(I,[]);hold on;
visboundaries(res_EGT)
title('EGT')


%% sJuneau
window_size=3;
t=0.01;
min_object_size=area_min;
I_norm=mat2gray(I,min_and_max_of_img1);
[res_juneau] = juneau(I_norm,window_size,t, min_object_size);
figure;imshow(I,[]);hold on;
visboundaries(res_juneau)
title('sJuneau')




%% level-set Chaselles
smoothness=0.0011111;
aditional_force=0;
I_norm=mat2gray(I,min_and_max_of_img1);
initialization=imdilate(res_EGT,strel('disk',5));
res_chaselles=activecontour(I_norm,initialization,500,'edge','ContractionBias',-aditional_force,'SmoothFactor',smoothness);
figure;imshow(I,[]);hold on;
visboundaries(res_chaselles)
title('sLS-Chasellses')



%% level-set Chan-Vese
smoothness=0.2;
aditional_force=0.17;
I_norm=mat2gray(I,min_and_max_of_img1);
initialization=imdilate(res_EGT,strel('disk',5));
res_chanvese=activecontour(I_norm,initialization,500,'Chan-Vese','ContractionBias',-aditional_force,'SmoothFactor',smoothness);
figure;imshow(I,[]);hold on;
visboundaries(res_chanvese)
title('sLS-Chan-Vese')






fg_mask=res_EGT;%choose foreground segmentation




%% dLoGm-Peng
sigma=rmin/sqrt(2):rmax/sqrt(2);
gamma=0.37;
t=-0.0022;
I_norm=mat2gray(I,min_and_max_of_img1);
res_log_peng=mlog_peng(I_norm,fg_mask,sigma,gamma,t);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_log_peng>0,ones(5)))
title('dLoGm-Peng')




%% dLoGm-Kong
sigma=rmin/sqrt(2):rmax/sqrt(2);
gamma=0.1;
max_overlap=0.05;
t=0;
I_norm=mat2gray(I,min_and_max_of_img1);
res_mlog_kong=mlog_kong(I_norm,fg_mask,sigma,gamma,t,max_overlap);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_mlog_kong>0,ones(5)))
title('dLoGm-Kong')




%% dLoGh-Zhang
gamma=3;
I_norm=mat2gray(I,min_and_max_of_img1);
res_hlog_zhang=hlog_zhang(I_norm,gamma,fg_mask);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_hlog_zhang>0,ones(5)))
title('dLoGh-Zhang')








%% dLoGg-Kong
sigma=rmin;
alpha=2;
I_norm=mat2gray(I,min_and_max_of_img1);
res_glog_kong=glog_kong(I_norm,fg_mask,sigma,alpha);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_glog_kong>0,ones(5)))
title('dLoGg-Kong')




%% dLoGg-Xu
sigma=rmin/sqrt(2):rmax/sqrt(2);
mean_shift_window=20;
I_norm=mat2gray(I,min_and_max_of_img1);
res_glog_xu=glog_xu(I_norm,fg_mask,sigma,mean_shift_window);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_glog_xu>0,ones(5)))
title('dLoGg-Xu')



%% dFRST
r=rmin:rmax;
alpha=2.6778;
t=4.10e-03;
k=6.25;
I_norm=mat2gray(I,min_and_max_of_img1);
res_frst=frst(I_norm,fg_mask,r,t,k,alpha);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_frst>0,ones(5)))
title('dFRST')


%% dGRST
r=rmin:rmax;
alpha=2.6778;
t=0.003;
k=5;
I_norm=mat2gray(I,min_and_max_of_img1);
res_gfrst=gfrst(I_norm,fg_mask,r,t,k,alpha);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_gfrst>0,ones(5)))
title('dGRST')


%% dRV-Qi
r=rmin:rmax;
sigma1=2.275;
w=19;
sigma2=5.333;
I_norm=mat2gray(I,min_and_max_of_img1);
res_rv_qi=rv_qi(I_norm,fg_mask,r,sigma1,w,sigma2);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_rv_qi>0,ones(5)))
title('dRV-Qi')



%% DT
t=0.069;
min_object_size=180;
min_hole_size=227.5;
h=4.66666;
I_norm=mat2gray(I,min_and_max_of_img1);
res_dt=dt(I_norm,fg_mask,min_object_size,min_hole_size,h,t);
figure;imshow(I,[]);hold on;
visboundaries(imdilate(res_dt>0,ones(5)))
title('DT')



%% MSER
% VLFeat package needed
% http://www.vlfeat.org/download.html

% delta=2.166666667;
% allowed_area_change=0.75;
% min_object_size=360;
% I_norm=mat2gray(I,min_and_max_of_img1);
% res_mser=mser_det(I_norm,fg_mask,delta,allowed_area_change,min_object_size);
% figure;imshow(I,[]);hold on;
% visboundaries(imdilate(res_mser>0,ones(5)))
% title('MSER')
% 



detection=res_dt;%choose detection



%% seeded_watershed_dt
res_seeded_watershed_dt=seeded_watershed_dt(I,detection,fg_mask);
figure;imshow(I,[]);hold on;
visboundaries(res_seeded_watershed_dt>0)
visboundaries(imdilate(detection>0,ones(5)),'Color','b')




%% seeded_watershed
res_seeded_watershed=seeded_watershed(I,detection,fg_mask);
figure;imshow(I,[]);hold on;
visboundaries(res_seeded_watershed>0)
visboundaries(imdilate(detection>0,ones(5)),'Color','b')




dice=dice_foreground(gt_mask,res_EGT);
dice
dice=dice_points(gt_mask,res_dt);
dice
dice=dice_final_segmentation(gt_mask,res_seeded_watershed);
dice




