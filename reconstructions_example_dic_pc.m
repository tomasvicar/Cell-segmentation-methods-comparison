clc;clear all;close all;
addpath('detection','evaluation','foreground_segmentation','final_segmentation','example_data')
addpath(genpath('util'))

% For functionality the third party libraries have to be installed:
% http://www.highcontentanalysis.org/
% https://www.albany.edu/celltracking/downloads.html


I_dic=imread('DIC_02_raw.tif');

I_pc=imread('PC_02_raw.tif');
I_pc=mat2gray(I_pc,[1484 16383]);%normalization based on the minimal and maximal intensity of first image



%% reconstruction DIC Koos
%[1] Koos, K., Moln´ar, J., Kelemen, L., Tam´as, G., Horvath, P.: DIC image reconstruction using an energy
% minimization framework to visualize optical path length distribution. Scientific reports 6, 30420 (2016).doi:10.1038/srep30420

ws=0.015;
wa=0.02;
I_dic_rek_Koos=dicEMM(I_dic,'direction',-135,'numiter',20000,'wSmooth',ws,'wAccept',wa);%takes ~30min
figure;
imshow(I_dic_rek_Koos,[])


%% reconstruction DIC Yin
%[1] Yin, Z., Ker, D.F.E., Kanade, T.: Restoring DIC microscopy images from multiple shear directions. In: Lecture
% Notes in Computer Science (including Subseries Lecture Notes in Artificial Intelligence and Lecture Notes in
% Bioinformatics), vol. 6801 LNCS, pp. 384–397 (2011). doi:10.1007/978-3-642-22092-0

ws=50;
wr=0.005;
I_dic_rek_Yin=dicZhaozhengYin(I_dic,'direction',-135, 's', ws, 'r', 0.005);
figure;
imshow(I_dic_rek_Yin,[])


%% reconstruction PC Yin
% [1]Yin, Z., Kanade, T., Chen, M.: Understanding the phase contrast optics to restore artifact-free microscopy
% images for segmentation. Medical Image Analysis 16(5), 1047–1062 (2012). doi:10.1016/j.media.2011.12.006.NIHMS150003

kernparas=struct('R',4000,'W',800,'radius',2);
optparas= struct('w_smooth_spatio',1,'w_sparsity',0.5,'epsilon',100,'gamma',3,'m_scale',1,'maxiter',100,'tol',eps);
I_pc_rek_Yin=precondition_linear_model(I_pc,optparas,kernparas,0);%takes ~30min
figure;
imshow(I_pc_rek_Yin,[])


%% reconstruction TopHat
r=17;
I_pc_rek_tophat=imtophat(I_pc,strel('disk',r));
figure;
imshow(I_pc_rek_tophat,[])

