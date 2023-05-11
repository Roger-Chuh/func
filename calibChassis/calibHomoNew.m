function [HH, camPoseVec,metricList1, ptUndist, intrMatHomo,pixList] = calibHomoNew(inputDir, paraDir)
%	The expected input properties are:
%
%        'cluster_min_size' : Minimum number of pixels in the clusters of approximately
%                             collinear feature pixels. The default value is 10.
%
%   'cluster_min_deviation' : Minimum accepted distance between a feature pixel and the
%                             line segment defined by the end points of its cluster.
%                             The default value is 2.
%
%                   'delta' : Discretization step for the parameter space. The default
%                             value is 0.5.
%
%       'kernel_min_height' : Minimum height for a kernel pass the culling operation.
%                             This property is restricted to the [0,1] range. The
%                             default value is 0.002.
%
%                'n_sigmas' : Number of standard deviations used by the Gaussian kernel
%                             The default value is 2.

% % close all
use_fisheye = 1;
global cfg

cbSize = 70;

cfg.cb_size = cbSize;
cfg.cb_row = 9;
cfg.cb_col = 14;
cfg.check_depth = 1;

calibFuncDir = 'E:\bk_20180627\SLAM\slam_algo\nanhuCalib_bak_20200921\func';

try
    if ~exist(fullfile(inputDir,'homo.mat'))
        if 0
            load(fullfile(inputDir,'CC.mat'),'CC');
            try
                load(fullfile(inputDir,'DD.mat'),'DD');
                CC = [CC;DD];
            catch
                wdvnlkj = 1;
            end
            len = 300;600;
            
            len1 = 600; len2 = 600/2;
        end
        
        
        load(fullfile(paraDir, 'calib.mat'));
        %     load('D:\Temp\20180625\calibCalTag1\calib_use\Omni_Calib_Results.mat');
        %     U_same = ocam_undistort_map(calib_data.ocam_model,'OutputView','same');
        %     U_full = ocam_undistort_map(calib_data.ocam_model,'OutputView','full');
        %     intrMat_same = U_same.K';
        %     intrMat_full = U_full.K';
        
        
        
        
        
        k1 = camParam.kc(1);
        k2 = camParam.kc(2);
        k3 = camParam.kc(5);
        p1 = camParam.kc(3);
        p2 = camParam.kc(4);
        if use_fisheye == 0
            intrMatRight = [camParam.foc(1),0,camParam.cen(1); 0,camParam.foc(2),camParam.cen(2); 0,0,1];
        else
            % %         load('D:\Temp\20180625\calibCalTag1\calib_use\oCamModel.mat');
            % %         U_same = ocam_undistort_map(oCamModel,'OutputView','same');
            % %         U_full = ocam_undistort_map(oCamModel,'OutputView','full');
            
            % %             load('E:\bk_20180627\pc\Temp\20180724\calibFishEye1\oCamModel.mat');
            
            load(fullfile(paraDir,'oCamModel.mat'));
            
            U_same = ocam_undistort_map(oCamModel,'OutputView','same');
            U_full = ocam_undistort_map(oCamModel,'OutputView','full');
            intrMat_same = U_same.K';
            intrMat_full = U_full.K';
            intrMatRight = intrMat_same;
            
            % %         intrMatRight = intrMat_full;
            % %         U_same = U_full;
        end
        
        dirInfo = dir(fullfile(inputDir,'*.png'));
        if length(dirInfo) == 0
            dirInfo = dir(fullfile(inputDir,'*.jpg'));
        end
        I = imread(fullfile(inputDir,dirInfo(1).name));
        %     I = imread('D:\Temp\20180702\homo1\calib\calib_frame_000001.png');
        %     I = imread('D:\Temp\20180706\calibHomo1\calib\calib_frame_000046.png');
        
        if 0
            if 0
                pixLine{1,1} = CC(1:7,:);
                metricLine{1,1} = [0 0 0; 1 0 0; 2 0 0; 6 0 0; 7 0 0; 8 0 0; 9 0 0];
                pixLine{2,1} = CC(8:17,:);
                metricLine{2,1} = [0 1 0; 1 1 0; 2 1 0; 3 1 0; 4 1 0; 5 1 0;6 1 0; 7 1 0; 8 1 0; 9 1 0];
                pixLine{3,1} = CC(18:28,:);
                metricLine{3,1} = [-1 2 0; 0 2 0; 1 2 0; 2 2 0; 3 2 0; 4 2 0; 5 2 0; 6 2 0; 7 2 0; 8 2 0; 9 2 0];
                pixLine{4,1} = CC(29:39,:);
                metricLine{4,1} = [-1 3 0; 0 3 0; 1 3 0; 2 3 0; 3 3 0; 4 3 0; 5 3 0; 6 3 0; 7 3 0; 8 3 0; 9 3 0];
                pixLine{5,1} = CC(40:51,:);
                metricLine{5,1} = [-2 4 0; -1 4 0; 0 4 0; 1 4 0; 2 4 0; 3 4 0; 4 4 0; 5 4 0; 6 4 0; 7 4 0; 8 4 0; 9 4 0];
                pixLine{6,1} = CC(52:58,:);
                metricLine{6,1} = [-2 5 0; -1 5 0; 0 5 0; 1 5 0; 7 5 0; 8 5 0; 9 5 0];
            elseif 0
                pixLine{1,1} = CC(1:6,:);
                metricLine{1,1} = [0 0 0; 1 0 0; 6 0 0; 7 0 0; 8 0 0; 9 0 0];
                pixLine{2,1} = CC(7:16,:);
                metricLine{2,1} = [0 1 0; 1 1 0; 2 1 0; 3 1 0; 4 1 0; 5 1 0;6 1 0; 7 1 0; 8 1 0; 9 1 0];
                pixLine{3,1} = CC(17:27,:);
                metricLine{3,1} = [-1 2 0; 0 2 0; 1 2 0; 2 2 0; 3 2 0; 4 2 0; 5 2 0; 6 2 0; 7 2 0; 8 2 0; 9 2 0];
                pixLine{4,1} = CC(28:39,:);
                metricLine{4,1} = [-2 3 0; -1 3 0; 0 3 0; 1 3 0; 2 3 0; 3 3 0; 4 3 0; 5 3 0; 6 3 0; 7 3 0; 8 3 0; 9 3 0];
                pixLine{5,1} = CC(40:51,:);
                metricLine{5,1} = [-2 4 0; -1 4 0; 0 4 0; 1 4 0; 2 4 0; 3 4 0; 4 4 0; 5 4 0; 6 4 0; 7 4 0; 8 4 0; 9 4 0];
                pixLine{6,1} = CC(52:57,:);
                metricLine{6,1} = [-2 5 0; -1 5 0; 0 5 0; 1 5 0; 7 5 0; 8 5 0];
            elseif 0
                pixLine{1,1} = CC(1:10,:);
                metricLine{1,1} = [0 1 0; 1 1 0; 2 1 0; 3 1 0; 4 1 0; 5 1 0;6 1 0; 7 1 0; 8 1 0; 9 1 0];
                pixLine{2,1} = CC(11:20,:);
                metricLine{2,1} = [0 2 0; 1 2 0; 2 2 0; 3 2 0; 4 2 0; 5 2 0;6 2 0; 7 2 0; 8 2 0; 9 2 0];
                pixLine{3,1} = CC(21:29,:);
                metricLine{3,1} = [0 3 0; 1 3 0; 2 3 0; 3 3 0; 4 3 0; 5 3 0; 6 3 0; 7 3 0; 8 3 0];
                pixLine{4,1} = CC(30:38,:);
                metricLine{4,1} = [0 4 0; 1 4 0; 2 4 0; 3 4 0; 4 4 0; 5 4 0; 6 4 0; 7 4 0; 8 4 0];
                % %     pixLine{5,1} = CC(30:38,:);
                % %     metricLine{5,1} = [-2 4 0; -1 4 0; 0 4 0; 1 4 0; 2 4 0; 3 4 0; 4 4 0; 5 4 0; 6 4 0; 7 4 0; 8 4 0; 9 4 0];
                % %     pixLine{6,1} = CC(52:57,:);
                % %     metricLine{6,1} = [-2 5 0; -1 5 0; 0 5 0; 1 5 0; 7 5 0; 8 5 0];
            elseif 0
                pixLine{1,1} = CC(1:10,:);
                metricLine{1,1} = [0 1 0; 1 1 0; 2 1 0; 3 1 0; 4 1 0; 5 1 0; 6 1 0; 7 1 0; 8 1 0; 9 1 0];
                pixLine{2,1} = CC(11:20,:);
                metricLine{2,1} = [0 2 0; 1 2 0; 2 2 0; 3 2 0; 4 2 0; 5 2 0; 6 2 0; 7 2 0; 8 2 0; 9 2 0];
                pixLine{3,1} = CC(21:30,:);
                metricLine{3,1} = [0 3 0; 1 3 0; 2 3 0; 3 3 0; 4 3 0; 5 3 0; 6 3 0; 7 3 0; 8 3 0; 9 3 0];
                pixLine{4,1} = CC(31:40,:);
                metricLine{4,1} = [0 4 0; 1 4 0; 2 4 0; 3 4 0; 4 4 0; 5 4 0; 6 4 0; 7 4 0; 8 4 0; 9 4 0];
                try
                    pixLine{5,1} = DD(1:9,:);
                    metricLine{5,1} = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 5 0 0; 6 0 0; 7 0 0; 8 0 0; 9 0 0];
                    pixLine{6,1} = DD(10:12,:);
                    metricLine{6,1} = [0 5 0; 7 5 0; 8 5 0];
                catch
                    ake = 1;
                end
                
            elseif 0
                
                pixLine{1,1} = CC(1:10,:);
                metricLine{1,1} = [0 0 0; 1 0 0; 2 0 0; 3 0 0; 4 0 0; 5 0 0; 6 0 0; 7 0 0; 8 0 0; 9 0 0];
                pixLine{2,1} = CC(11:20,:);
                metricLine{2,1} = [-1 1 0; 0 1 0; 1 1 0; 2 1 0; 3 1 0; 4 1 0; 5 1 0; 6 1 0; 7 1 0; 8 1 0];
                pixLine{3,1} = CC(21:30,:);
                metricLine{3,1} = [-1 2 0; 0 2 0; 1 2 0; 2 2 0; 3 2 0; 4 2 0; 5 2 0; 6 2 0; 7 2 0; 8 2 0];
                pixLine{4,1} = CC(31:40,:);
                metricLine{4,1} = [-1 3 0; 0 3 0; 1 3 0; 2 3 0; 3 3 0; 4 3 0; 5 3 0; 6 3 0; 7 3 0; 8 3 0];
                pixLine{5,1} = CC(41:45,:);
                metricLine{5,1} = [-1 4 0; 0 4 0; 1 4 0; 7 4 0; 8 4 0];
                pixLine{6,1} = CC(46:47,:);
                metricLine{6,1} = [7 -1 0; 8 -1 0];
                
            elseif 0
                
                pixLine{1,1} = CC(1:18,:);
                pixLine{2,1} = CC(19:27,:);
                pixLine{3,1} = CC(28:45,:);
                pixLine{4,1} = CC(46:54,:);
                pixLine{5,1} = CC(55:71,:);
                pixLine{6,1} = CC(72:80,:);
                pixLine{7,1} = CC(81:97,:);
                pixLine{8,1} = CC(98:106,:);
                pixLine{9,1} = CC(107:117,:);
                pixLine{10,1} = CC(118:119,:);
                pixLine{11,1} = CC(120,:);
                
                % %         pixLine{2,1} = CC(18:26,:);
                % %         pixLine{4,1} = CC(44:52,:);
                % %         pixLine{6,1} = CC(69:77,:);
                % %         pixLine{8,1} = CC(95:103,:);
                
                metricLine{1,1} = len2.*[[1 13:17]' repmat(-2,6,1) repmat(0,6,1)];
                metricLine{3,1} = len2.*[[0:17]' repmat(0,18,1) repmat(0,18,1)];
                metricLine{5,1} = len2.*[[0:17]' repmat(2,18,1) repmat(0,18,1)];
                metricLine{7,1} = len2.*[[0:17]' repmat(4,18,1) repmat(0,18,1)];
                metricLine{9,1} = len2.*[[0:17]' repmat(6,18,1) repmat(0,18,1)];
                metricLine{11,1} = len2.*[[0:3 13:17]' repmat(8,9,1) repmat(0,9,1)];
                
                
                metricLine{2,1} = len2.*[[0:2:16]' repmat(-1,9,1) repmat(0,9,1)];
                metricLine{4,1} = len2.*[[0:2:16]' repmat(1,9,1) repmat(0,9,1)];
                metricLine{6,1} = len2.*[[0:2:16]' repmat(3,9,1) repmat(0,9,1)];
                metricLine{8,1} = len2.*[[0:2:16]' repmat(5,9,1) repmat(0,9,1)];
                metricLine{10,1} = len2.*[[0:2:16]' repmat(7,9,1) repmat(0,9,1)];
            elseif 0
                
                pixLine{1,1} = CC(1:10,:);
                pixLine{2,1} = CC(11:28,:);
                pixLine{3,1} = CC(29:38,:);
                pixLine{4,1} = CC(39:56,:);
                pixLine{5,1} = CC(57:66,:);
                pixLine{6,1} = CC(67:85,:);
                pixLine{7,1} = CC(86:95,:);
                pixLine{8,1} = CC(96:114,:);
                pixLine{9,1} = CC(115:124,:);
                pixLine{10,1} = CC(125:139,:);
                pixLine{11,1} = CC(140:142,:);
                pixLine{12,1} = CC(143:147,:);
                pixLine{13,1} = CC(148,:);
                
                % %         pixLine{2,1} = CC(18:26,:);
                % %         pixLine{4,1} = CC(44:52,:);
                % %         pixLine{6,1} = CC(69:77,:);
                % %         pixLine{8,1} = CC(95:103,:);
                
                
                
                metricLine{1,1} = len2.*[[0:2:18]' repmat(-1,10,1) repmat(0,10,1)];
                metricLine{3,1} = len2.*[[0:2:18]' repmat(1,10,1) repmat(0,10,1)];
                metricLine{5,1} = len2.*[[0:2:18]' repmat(3,10,1) repmat(0,10,1)];
                metricLine{7,1} = len2.*[[0:2:18]' repmat(5,10,1) repmat(0,10,1)];
                metricLine{9,1} = len2.*[[0:2:18]' repmat(7,10,1) repmat(0,10,1)];
                metricLine{11,1} = len2.*[[0 16 18]' repmat(9,3,1) repmat(0,3,1)];
                metricLine{13,1} = len2.*[[18]' repmat(-3,1,1) repmat(0,1,1)];
                
                
                metricLine{2,1} = len2.*[[1:18]' repmat(0,18,1) repmat(0,18,1)];
                metricLine{4,1} = len2.*[[0 [2:18]]' repmat(2,18,1) repmat(0,18,1)];
                metricLine{6,1} = len2.*[[0:18]' repmat(4,19,1) repmat(0,19,1)];
                metricLine{8,1} = len2.*[[0:18]' repmat(6,19,1) repmat(0,19,1)];
                metricLine{10,1} = len2.*[[[0:5] [10:18]]' repmat(8,15,1) repmat(0,15,1)];
                metricLine{12,1} = len2.*[[14:18]' repmat(-2,5,1) repmat(0,5,1)];
                
                
            elseif 0
                
                pixLine{1,1} = CC(1:18,:);
                pixLine{2,1} = CC(19:27,:);
                pixLine{3,1} = CC(28:45,:);
                pixLine{4,1} = CC(46:54,:);
                pixLine{5,1} = CC(55:71,:);
                pixLine{6,1} = CC(72:80,:);
                pixLine{7,1} = CC(81:97,:);
                pixLine{8,1} = CC(98:106,:);
                pixLine{9,1} = CC(107:117,:);
                pixLine{10,1} = CC(118:119,:);
                pixLine{11,1} = CC(120,:);
                
                % %         pixLine{2,1} = CC(18:26,:);
                % %         pixLine{4,1} = CC(44:52,:);
                % %         pixLine{6,1} = CC(69:77,:);
                % %         pixLine{8,1} = CC(95:103,:);
                
                
                
                metricLine{1,1} = len2.*[[1:18]' repmat(-6,18,1) repmat(0,18,1)];
                metricLine{3,1} = len2.*[[1:18]' repmat(-4,18,1) repmat(0,18,1)];
                metricLine{5,1} = len2.*[[2:18]' repmat(-2,17,1) repmat(0,17,1)];
                metricLine{7,1} = len2.*[[2:18]' repmat(0,17,1) repmat(0,17,1)];
                metricLine{9,1} = len2.*[[[2 3] [10:18]]' repmat(2,11,1) repmat(0,11,1)];
                metricLine{11,1} = len2.*[[18]' repmat(4,1,1) repmat(0,1,1)];
                
                
                metricLine{2,1} = len2.*[[2:2:18]' repmat(-5,9,1) repmat(0,9,1)];
                metricLine{4,1} = len2.*[[2:2:18]' repmat(-3,9,1) repmat(0,9,1)];
                metricLine{6,1} = len2.*[[2:2:18]' repmat(-1,9,1) repmat(0,9,1)];
                metricLine{8,1} = len2.*[[2:2:18]' repmat(1,9,1) repmat(0,9,1)];
                metricLine{10,1} = len2.*[[16 18]' repmat(3,2,1) repmat(0,2,1)];
                
            elseif 1
                pixLine{1,1} = CC(1:5,:);
                pixLine{2,1} = CC(6:14,:);
                pixLine{3,1} = CC(15:33,:);
                pixLine{4,1} = CC(34:43,:);
                pixLine{5,1} = CC(44:62,:);
                pixLine{6,1} = CC(63:72,:);
                pixLine{7,1} = CC(73:91,:);
                pixLine{8,1} = CC(92:101,:);
                pixLine{9,1} = CC(102:120,:);
                pixLine{10,1} = CC(121:130,:);
                pixLine{11,1} = CC(131:149,:);
                pixLine{12,1} = CC(150:153,:);
                pixLine{13,1} = CC(154:end,:);
                
                % %         pixLine{2,1} = CC(18:26,:);
                % %         pixLine{4,1} = CC(44:52,:);
                % %         pixLine{6,1} = CC(69:77,:);
                % %         pixLine{8,1} = CC(95:103,:);
                
                
                
                metricLine{1,1} = len2.*[[0 1 16:18]' repmat(-2,5,1) repmat(0,5,1)];
                metricLine{3,1} = len2.*[[0:18]' repmat(0,19,1) repmat(0,19,1)];
                metricLine{5,1} = len2.*[[0:18]' repmat(2,19,1) repmat(0,19,1)];
                metricLine{7,1} = len2.*[[0:18]' repmat(4,19,1) repmat(0,19,1)];
                metricLine{9,1} = len2.*[[0:18]' repmat(6,19,1) repmat(0,19,1)];
                metricLine{11,1} = len2.*[[[-1:5] [7:18]]' repmat(8,19,1) repmat(0,19,1)];
                metricLine{13,1} = len2.*[[[17:18]]' repmat(10,2,1) repmat(0,2,1)];
                
                metricLine{2,1} = len2.*[[[0:2:6] [10:2:18]]' repmat(-1,9,1) repmat(0,9,1)];
                metricLine{4,1} = len2.*[[0:2:18]' repmat(1,10,1) repmat(0,10,1)];
                metricLine{6,1} = len2.*[[0:2:18]' repmat(3,10,1) repmat(0,10,1)];
                metricLine{8,1} = len2.*[[0:2:18]' repmat(5,10,1) repmat(0,10,1)];
                metricLine{10,1} = len2.*[[0:2:18]' repmat(7,10,1) repmat(0,10,1)];
                metricLine{12,1} = len2.*[[-2 14 16 18]' repmat(9,4,1) repmat(0,4,1)];
            elseif 0
                pixLine{1,1} = CC(1:2,:);
                pixLine{2,1} = CC(3:16,:);
                pixLine{3,1} = CC(17:27,:);
                pixLine{4,1} = CC(28:47,:);
                pixLine{5,1} = CC(48:57,:);
                pixLine{6,1} = CC(58:76,:);
                pixLine{7,1} = CC(77:86,:);
                pixLine{8,1} = CC(87:105,:);
                pixLine{9,1} = CC(106:115,:);
                pixLine{10,1} = CC(116:134,:);
                pixLine{11,1} = CC(135:144,:);
                pixLine{12,1} = CC(145:154,:);
                pixLine{13,1} = CC(155,:);
                
                % %         pixLine{2,1} = CC(18:26,:);
                % %         pixLine{4,1} = CC(44:52,:);
                % %         pixLine{6,1} = CC(69:77,:);
                % %         pixLine{8,1} = CC(95:103,:);
                
                
                
                metricLine{2,1} = len2.*[[[-1:4] [11:18]]' repmat(4,14,1) repmat(0,14,1)];
                metricLine{4,1} = len2.*[[-1:18]' repmat(6,20,1) repmat(0,20,1)];
                metricLine{6,1} = len2.*[[0:18]' repmat(8,19,1) repmat(0,19,1)];
                metricLine{8,1} = len2.*[[0:18]' repmat(10,19,1) repmat(0,19,1)];
                metricLine{10,1} = len2.*[[0:18]' repmat(12,19,1) repmat(0,19,1)];
                metricLine{12,1} = len2.*[[[1 2] [11:18]]' repmat(14,10,1) repmat(0,10,1)];
                
                
                metricLine{1,1} = len2.*[[16 18]' repmat(3,2,1) repmat(0,2,1)];
                metricLine{3,1} = len2.*[[-2:2:18]' repmat(5,11,1) repmat(0,11,1)];
                metricLine{5,1} = len2.*[[0:2:18]' repmat(7,10,1) repmat(0,10,1)];
                metricLine{7,1} = len2.*[[0:2:18]' repmat(9,10,1) repmat(0,10,1)];
                metricLine{9,1} = len2.*[[0:2:18]' repmat(11,10,1) repmat(0,10,1)];
                metricLine{11,1} = len2.*[[0:2:18]' repmat(13,10,1) repmat(0,10,1)];
                metricLine{13,1} = len2.*[[16]' repmat(15,1,1) repmat(0,1,1)];
            elseif 0
                
                pixLine{1,1} = CC(1,:);
                pixLine{2,1} = CC(2:4,:);
                pixLine{3,1} = CC(5:18,:);
                pixLine{4,1} = CC(19:27,:);
                pixLine{5,1} = CC(28:46,:);
                pixLine{6,1} = CC(47:55,:);
                pixLine{7,1} = CC(56:73,:);
                pixLine{8,1} = CC(74:82,:);
                pixLine{9,1} = CC(83:99,:);
                pixLine{10,1} = CC(100:108,:);
                pixLine{11,1} = CC(109:126,:);
                pixLine{12,1} = CC(127:132,:);
                pixLine{13,1} = CC(133:136,:);
                
                % %         pixLine{2,1} = CC(18:26,:);
                % %         pixLine{4,1} = CC(44:52,:);
                % %         pixLine{6,1} = CC(69:77,:);
                % %         pixLine{8,1} = CC(95:103,:);
                
                
                
                metricLine{1,1} = len2.*[[-9]' repmat(4,1,1) repmat(0,1,1)];
                metricLine{3,1} = len2.*[[[-9:-2] [5:10]]' repmat(6,14,1) repmat(0,14,1)];
                metricLine{5,1} = len2.*[[-9:9]' repmat(8,19,1) repmat(0,19,1)];
                metricLine{7,1} = len2.*[[-9:8]' repmat(10,18,1) repmat(0,18,1)];
                metricLine{9,1} = len2.*[[-9:7]' repmat(12,17,1) repmat(0,17,1)];
                metricLine{11,1} = len2.*[[-9:8]' repmat(14,18,1) repmat(0,18,1)];
                metricLine{13,1} = len2.*[[[6:9]]' repmat(16,4,1) repmat(0,4,1)];
                
                metricLine{2,1} = len2.*[[[-8 -6 10] ]' repmat(5,3,1) repmat(0,3,1)];
                metricLine{4,1} = len2.*[[-8:2:8]' repmat(7,9,1) repmat(0,9,1)];
                metricLine{6,1} = len2.*[[-8:2:8]' repmat(9,9,1) repmat(0,9,1)];
                metricLine{8,1} = len2.*[[-8:2:8]' repmat(11,9,1) repmat(0,9,1)];
                metricLine{10,1} = len2.*[[-8:2:8]' repmat(13,9,1) repmat(0,9,1)];
                metricLine{12,1} = len2.*[[-8 -6 2 4 6 8]' repmat(15,6,1) repmat(0,6,1)];
                
            else
                
                pixLine{1,1} = CC(1:2,:);
                pixLine{2,1} = CC(3:5,:);
                pixLine{3,1} = CC(6:24,:);
                pixLine{4,1} = CC(25:33,:);
                pixLine{5,1} = CC(34:51,:);
                pixLine{6,1} = CC(52:60,:);
                pixLine{7,1} = CC(61:79,:);
                pixLine{8,1} = CC(80:84,:);
                pixLine{9,1} = CC(85:94,:);
                pixLine{10,1} = CC(95:99,:);
                % %         pixLine{2,1} = CC(18:26,:);
                % %         pixLine{4,1} = CC(44:52,:);
                % %         pixLine{6,1} = CC(69:77,:);
                % %         pixLine{8,1} = CC(95:103,:);
                
                
                
                metricLine{1,1} = len2.*[[-9 -8]' repmat(10,2,1) repmat(0,2,1)];
                metricLine{3,1} = len2.*[[-9:9]' repmat(12,19,1) repmat(0,19,1)];
                metricLine{5,1} = len2.*[[[-9:7] [9]]' repmat(14,18,1) repmat(0,18,1)];
                metricLine{7,1} = len2.*[[-9:9]' repmat(16,19,1) repmat(0,19,1)];
                metricLine{9,1} = len2.*[[-9:0]' repmat(18,10,1) repmat(0,10,1)];
                
                metricLine{2,1} = len2.*[[[-8 -6 10] ]' repmat(11,3,1) repmat(0,3,1)];
                metricLine{4,1} = len2.*[[-8:2:8]' repmat(13,9,1) repmat(0,9,1)];
                metricLine{6,1} = len2.*[[-8:2:8]' repmat(15,9,1) repmat(0,9,1)];
                metricLine{8,1} = len2.*[[-8:2:0]' repmat(17,5,1) repmat(0,5,1)];
                metricLine{10,1} = len2.*[[-8:2:0]' repmat(19,5,1) repmat(0,5,1)];
                
            end
            
            
            
            
            pixLine = pixLine([1:end]);
            metricLine = metricLine([1:end]);
            
            
            
            
            
            % pixLine{7,1} = [50 981; 138 1019; 253 1062; 1744 1042];
            % metricLine{7,1} = [-2 5 0; -1 5 0; 0 5 0; 8 5 0];
            pixList = cell2mat(pixLine);
            metricList = cell2mat(metricLine);  %.*len;
            % %
            % %     pixList = pixList([[[2;3;4;8;9;10;11;12;13;14;15;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;116;117;118;119;120;121;122;123;126;127;128;130;131]]],:);
            % %     metricList = metricList([[2;3;4;8;9;10;11;12;13;14;15;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;116;117;118;119;120;121;122;123;126;127;128;130;131]],:);
            % %
            % %     pixList = pixList([[[1;2;3;4;5;6;7;8;9;10;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118]]],:);
            % %     metricList = metricList([[1;2;3;4;5;6;7;8;9;10;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118]],:);
            
            
            % % %             pixList = pixList([[[2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;24;25;26;27;28;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;116;117;118;119;120;121;122;123;124;127;128;130;131;132;133;134;135;136;137;138;141;143;144;145]]],:);
            % % %             metricList = metricList([[[2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;24;25;26;27;28;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;116;117;118;119;120;121;122;123;124;127;128;130;131;132;133;134;135;136;137;138;141;143;144;145]]],:);
            % % %         % % % %
            % % %             pixList = pixList([[[1;2;3;4;5;6;7;8;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;55;56;57;58;59;60;61;62;63;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;132;133;134]]],:);
            % % %             metricList = metricList([[[1;2;3;4;5;6;7;8;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;55;56;57;58;59;60;61;62;63;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;132;133;134]]],:);
            
            
            
            if 1
                pixList = pixList([[[7;8;9;10;11;12;13;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;34;35;36;37;38;39;40;41;42;44;45;46;47;48;49;50;51;52;53;54;55;56;57;59;60;61;62;63;64;65;66;67;68;69;71;72;73;74;76;77;78;79;80;81;82;83;85;89;90;91;92;95;96;97;98;100;101;102;107;108;109;110;111;112;113;114;115;118;119;120;121;124;125;126;127;128;129;130;133;138;139;140;141;142;143;144;145;146;147;148;149;151;152]]],:);
                metricList = metricList([[[7;8;9;10;11;12;13;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;34;35;36;37;38;39;40;41;42;44;45;46;47;48;49;50;51;52;53;54;55;56;57;59;60;61;62;63;64;65;66;67;68;69;71;72;73;74;76;77;78;79;80;81;82;83;85;89;90;91;92;95;96;97;98;100;101;102;107;108;109;110;111;112;113;114;115;118;119;120;121;124;125;126;127;128;129;130;133;138;139;140;141;142;143;144;145;146;147;148;149;151;152]]],:);
                % %
                % %             pixList = pixList([[[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;114]]],:);
                % %             metricList = metricList([[[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;114]]],:);
                % %
                % %             pixList = pixList([[[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;104;105;106;107;108;109;110;111;112]]],:);
                % %             metricList = metricList([[[1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;104;105;106;107;108;109;110;111;112]]],:);
            end
            
            kjvn = 1;
            
            metricList1 = metricList;
            %     camPoseVec = posest(ptUndist', metricList, 0.95, intrMatRight, 'repr_err');
            
            metricList(:,end) = 1;
            
        else
            
            [pt3dX, pt3dY] = meshgrid([0:cfg.cb_size:cfg.cb_size*(cfg.cb_col-1)], [0:cfg.cb_size:cfg.cb_size*(cfg.cb_row-1)]);
            metricList = [pt3dX(:) pt3dY(:)];
            metricList(:,3) = 1;
            
            metricList1 = metricList;
            metricList1(:,end) = 0;
            
            modeE = 2;
            try
                imgLst1 = getImgLst(inputDir,'png');
            catch
                try
                    imgLst1 = getImgLst(inputDir,'bmp');
                catch
                    imgLst1 = getImgLst(inputDir,'jpg');
                end
            end
            config.dX = cbSize; % 23.7; 24.1; 70; 18.4;
            config.dY = cbSize; % 23.7; 24.1; 70; 18.4;
            config.estDistortion = [1;1;1;1;0];
            [camParamL, cbcXYL, cbGridL, config, ~, cbcXYTmpL, cbGridTmpL, goodIdL] = CbCalibSingleUse4(imgLst1, [], calibFuncDir, modeE, config);
            pixList = cbcXYL{1}';
        end
        
        
        figure,subplot(1,2,1);imshow(I);hold on;plot(pixList(:,1),pixList(:,2),'-r');hold on;plot(pixList(1,1),pixList(1,2),'*b');axis equal;subplot(1,2,2);plot3(metricList(:,1),metricList(:,2),metricList(:,3),'-r');axis equal;hold on;plot3(metricList(1,1),metricList(1,2),metricList(1,3),'*b');axis equal
        %     ptUndistRight = normalize_pixel(pixList',camParam.foc,camParam.cen,camParam.kc,0);
        if use_fisheye == 0
            ptUndistRight = normalize_pixel(pixList',camParam.foc,camParam.cen,camParam.kc,0);
            ptUndist3D = [ptUndistRight; ones(1,size(ptUndistRight,2))]';
            ptUndist = intrMatRight*[ptUndistRight;ones(1,size(ptUndistRight,2))];
            ptUndist = [ptUndist(1,:)./ptUndist(3,:);ptUndist(2,:)./ptUndist(3,:)];
            camPoseVec = posest(ptUndist', metricList1, 0.95, intrMatRight, 'repr_err');
            % %         ptUndistRight1 = normalize_pixel(pixList',camParam.foc,camParam.cen,camParam.kc,0);
            % %         ptUndist3D1 = [ptUndistRight1; ones(1,size(ptUndistRight1,2))]';
            % %         ptUndist1 = intrMatRight*[ptUndistRight1;ones(1,size(ptUndistRight1,2))];
            % %         ptUndist1 = [ptUndist1(1,:)./ptUndist1(3,:);ptUndist1(2,:)./ptUndist1(3,:)];
            
            
            
        else
            [ptUndist, ptUndist3D] = remapFishEyePix(pixList', intrMatRight, oCamModel);
            ptUndist3D = ptUndist3D';
            camPoseVec = posest(ptUndist', metricList1, 0.95, intrMatRight, 'repr_err');
            
            
        end
        %     ptUndist = [ptUndist(1,:)./ptUndist(3,:);ptUndist(2,:)./ptUndist(3,:)];
        %     ptUndist3D = [ptUndistRight; ones(1,size(ptUndistRight,2))]';
        ptUndistIn = ptUndist;
        ptUndistIn = ptUndistIn(:,ptUndistIn(1,:) > 0 & ptUndistIn(1,:) < size(I,2) & ptUndistIn(2,:) > 0 & ptUndistIn(2,:) < size(I,1));
        if use_fisheye == 0
            nim = undistortimage(I, [intrMatRight(1,1) intrMatRight(2,2)]', intrMatRight(1,3), intrMatRight(2,3), k1, k2, k3, p1, p2);
        else
            nim = ocam_undistort(I,U_same);
        end
        % % lines = kht(nim, 'cluster_min_size', 10, 'cluster_min_deviation', 2, 'delta', 0.5, 'kernel_min_height', 0.002, 'n_sigmas', 2);
        % lines2 = kht(nim, 'cluster_min_size', 10, 'cluster_min_deviation',2, 'delta', 0.5, 'kernel_min_height', 0.05, 'n_sigmas', 2);
        if 0
            lines2 = kht(nim, 'cluster_min_size', 10, 'cluster_min_deviation',2, 'delta', 0.5, 'kernel_min_height', 0.04, 'n_sigmas', 2);
            lines2 = kht(nim, 'cluster_min_size', 10, 'cluster_min_deviation',2, 'delta', 0.5, 'kernel_min_height', 0.004, 'n_sigmas', 2);
            % % figure, linePara2 = plotLinePara(lines2, nim);
            lines = lines2([4 5 7 8 9 10 13 19 20 21 23],:);
            lines = lines([2 4 3 1 8 10 9 6 5 7 11],:);
            % % figure, linePara = plotLinePara(lines,nim);
            x = repmat([0 : 7],4,1);
            y = repmat([0 : 3]',1,8);
            xyz = [x(:) y(:) ones(size(x,1)*size(x,2),1)];
            rowList = [1 2 3 4]; colList = [1 2 3 5 6 7 8];
            cnt = 1;
            for k = 1 : length(rowList)
                xyzTmp = [x(rowList(k),:);y(rowList(k),:);ones(1,length(x(rowList(k),:)))];
            end
            
            for k = 1 : length(colList)
                xyzTmp = [x(:,colList(k)) y(:,colList(k)) ones(length(x(:,colList(k))),1)]';
            end
        end
        
        [tform,inlierPtsDistorted,inlierPtsOriginal] = estimateGeometricTransform(metricList(:,1:2), ptUndist3D(:,1:2), 'projective');
        [H,Hnorm,inv_Hnorm] = compute_homography(metricList', ptUndist3D');
        if 1 % all pixels are considered as inliers
            % [H2, inliers2] = ransacfithomography(metricList', ptUndist3D', 0.01);
            [H2, inliers2] = ransacfithomography(metricList', ptUndist3D', 0.0002);
        else % only those who are within the image are considered as inliers
            [H2, inliers2] = ransacfithomography(metricList(ptUndist(1,:) > 0 & ptUndist(1,:) < size(I,2) & ptUndist(2,:) > 0 & ptUndist(2,:) < size(I,1),:)', ptUndist3D(ptUndist(1,:) > 0 & ptUndist(1,:) < size(I,2) & ptUndist(2,:) > 0 & ptUndist(2,:) < size(I,1),:)', 0.01);
        end
        flag = ptUndist(1,:) > 0 & ptUndist(1,:) < size(I,2) & ptUndist(2,:) > 0 & ptUndist(2,:) < size(I,1);
        % % [H2, inliers2] = ransacfithomography(metricList(ptUndist(1,:) > 495 & ptUndist(1,:) < 1610 & ptUndist(2,:) > 328 & ptUndist(2,:) < 970,:)', ptUndist3D(ptUndist(1,:) > 495 & ptUndist(1,:) < 1610 & ptUndist(2,:) > 328 & ptUndist(2,:) < 970,:)', 0.01);
        % % flag = ptUndist(1,:) > 495 & ptUndist(1,:) < 1610 & ptUndist(2,:) > 328 & ptUndist(2,:) < 970;
        
        
        % % Hline = estHomoWithLines(img1,img2,lmatch);
        % % [tform,inlierPtsDistorted,inlierPtsOriginal] = estHomo(img2,img1,lineCoeff1(:,1:2),lineCoeff2(:,1:2));
        % % Hline = tform.T;
        
        
        
        
        
        
        
        if 1
            mappedPt = inv(H)*metricList';
            mappedPtMetric = H*ptUndist3D';
            HH = H;
        else
            mappedPt = H2*metricList';
            mappedPtMetric = inv(H2)*ptUndist3D';
            HH = inv(H2);
        end
        
        mappedPtMetric = [mappedPtMetric(1,:)./mappedPtMetric(3,:);mappedPtMetric(2,:)./mappedPtMetric(3,:);ones(1,size(mappedPtMetric,2))];
        [~,lenMetric] = NormalizeVector(mappedPtMetric' - metricList);
        mappedPt = intrMatRight*mappedPt;
        mappedPt = [mappedPt(1,:)./mappedPt(3,:); mappedPt(2,:)./mappedPt(3,:)]';
        figure,imshow(I);hold on;plot(pixList(:,1),pixList(:,2),'or');
        figure,hold on;axis equal;plot(ptUndist(1,:),ptUndist(2,:),'-r');plot(ptUndist(1,:),ptUndist(2,:),'or'); plot(mappedPt(:,1),mappedPt(:,2),'.g');plot(mappedPt(inliers2,1),mappedPt(inliers2,2),'*b');
        
        figure,imshow(nim);hold on;plot(ptUndistIn(1,:),ptUndistIn(2,:),'-r');plot(ptUndistIn(1,:),ptUndistIn(2,:),'or'); plot(mappedPt(:,1),mappedPt(:,2),'.g');
        
        figure;plot(ptUndist(1,:),ptUndist(2,:),'or');axis equal
        
        inliers22 = 1:size(mappedPt,1);
        figure,subplot(1,2,1);plot(mappedPt(inliers22,1)-ptUndist(1,inliers22)',mappedPt(inliers22,2)-ptUndist(2,inliers22)','+r');axis equal;title('pixel');
        subplot(1,2,2);plot(metricList(inliers22,1)-mappedPtMetric(1,inliers22)',metricList(inliers22,2)-mappedPtMetric(2,inliers22)','+r');axis equal;title('metric');
        
        
        
        figure,subplot(1,2,1);plot(mappedPt(flag,1)-ptUndist(1,flag)',mappedPt(flag,2)-ptUndist(2,flag)','+r');axis equal; title('pixel');
        subplot(1,2,2);plot(metricList(flag,1)-mappedPtMetric(1,flag)',metricList(flag,2)-mappedPtMetric(2,flag)','+r');axis equal; title('metric');
        intrMatHomo = intrMatRight;
        
        
        save(fullfile(inputDir,'homo.mat'),'HH','camPoseVec','metricList1', 'ptUndist', 'intrMatHomo','pixList');
        
        inliers222 = 1:size(mappedPt,1);
        errr = [mappedPt(inliers222,1)-ptUndist(1,inliers222)',mappedPt(inliers222,2)-ptUndist(2,inliers222)'];
        errr_ = [metricList(inliers222,1)-mappedPtMetric(1,inliers222)',metricList(inliers222,2)-mappedPtMetric(2,inliers222)'];
        [~,err2] = NormalizeVector(errr);
        [~,err2_] = NormalizeVector(errr_);
        inliers3 = find(err2<3);
        inliers3_ = find(err2_<7);
        
        
        figure,imshow(I);hold on;
        
        for i = 1 : size(mappedPt,1)
            
            plot(pixList(i,1),pixList(i,2),'*r');
            text(pixList(i,1),pixList(i,2),num2str(round(errr_(i))), 'Color',[1 0 0],'FontSize',15,'FontWeight','bold');
        end
        
        
        
        
        
        sagb = 1;
        
        
    else
        load(fullfile(inputDir,'homo.mat'));
    end
catch
    load(fullfile(inputDir,'homo.mat'));
end

if ~exist('pixList','var')
    pixList = [];
end


end
