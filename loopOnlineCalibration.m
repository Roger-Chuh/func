function loopOnlineCalibration()

% close all



inputDir = 'G:\matlab\data\vio\with_epiplane';
% inputDir = 'G:\matlab\data\vio\with_epiplane_good';
% inputDir = 'G:\matlab\data\vio\with_epiplane_back1';
% inputDir = 'G:\matlab\data\vio\with_epiplane_back2';
inputDir = 'G:\matlab\data\vio\with_epiplane_1';
inputDir = 'G:\matlab\data\vio\with_epiplane_back0';
inputDir = 'G:\matlab\data\vio\with_epiplane_back1';
inputDir = 'G:\matlab\data\vio\with_epiplane_calib';
inputDir = 'G:\matlab\data\vio\with_epiplane_calib\reproj_epiplane';
inputDir = 'G:\matlab\data\vio\with_epiplane_calib\reproj_epiplane_non_fix_tbc0';
inputDir = 'G:\matlab\data\vio\with_epiplane_1\reproj_epiplane_non_fix_tbc0';
inputDir = 'G:\matlab\data\vio\with_epiplane_back3\reproj_epiplane_non_fix_tbc0';
inputDir = 'G:\matlab\data\vio\with_epiplane_back3\epiplane_fix_tbc0';
inputDir = 'G:\matlab\data\vio\with_epiplane_back0\reproj_epiplane_non_fix_tbc0'; 
inputDir = 'G:\matlab\data\vio\with_epiplane_back0\epiplane_fix_tbc0';
inputDir = 'G:\matlab\data\vio\gt_2\reproj_epiplane_non_fix_tbc0';
% inputDir = 'G:\matlab\data\vio\with_epiplane_1\reproj_epiplane_non_fix_tbc0';
inputDir = 'G:\matlab\data\vio\gt_2\epiplane_fix_tbc0';
inputDir = 'G:\matlab\data\vio\gt_6\epiplane_fix_tbc0';
inputDir = 'G:\matlab\data\vio\gt_6\with_prior\0.1';
% inputDir = 'G:\matlab\data\vio\with_epiplane_back0\with_prior\0.1';
inputDir = 'G:\matlab\data\vio\with_epiplane_back2\with_prior\0.1';
% inputDir = 'G:\matlab\data\vio\with_epiplane_back1\with_prior\0.1';
inputDir = 'G:\matlab\data\vio\gt_2\0.1\with_prior';
inputDir = 'G:\matlab\data\vio\gt_6\with_prior\fast';
inputDir = 'G:\matlab\data\vio\ci\008\on1';
inputDir = 'G:\matlab\data\vio\ci\001\on1';
inputDir = 'G:\matlab\data\vio\ci\009\on1';
inputDir = 'G:\matlab\data\vio\ci\009\on2';
inputDir = 'G:\matlab\data\vio\ci\004\on2';
inputDir = 'G:\matlab\data\vio\ci\009\on3';
inputDir = 'G:\matlab\data\vio\ci\009\on4';
inputDir = 'G:\matlab\data\vio\gt_1\with_prior\on1';

inputDir = 'G:\matlab\data\vio\ci\009\on5';
inputDir = 'G:\matlab\data\vio\ci\002\on1';

% inputDir = 'G:\matlab\data\vio\gt_6\reproj_epiplane_non_fix_tbc0';
% inputDir = 'G:\matlab\data\vio\gt_6\non_overlapping';


% inputDir = 'G:\matlab\data\vio\gt_3\epiplane_fix_tbc0';
% inputDir = 'G:\matlab\data\vio\gt_5\epiplane_fix_tbc0';



% inputDir = 'G:\matlab\data\vio\with_epiplane_5';
% inputDir = 'G:\matlab\data\vio\with_epiplane_good2';
dirInfo = dir(fullfile(inputDir,'oc_time_*'));

Trans = [];
for i = 1 : size(dirInfo)
    folder = fullfile(inputDir,dirInfo(i).name) ;
    try
        [trans, trans_cat, rot, rot_cat] = testOnlineCalibration(folder, 0);
        
        if 1
            figure(1),subplot(2,4,1);hold on;plot(trans(:,1,1),'-r');plot(trans(:,2,1),'-g');plot(trans(:,3,1),'-b');title('cam0 trans');grid on;
            subplot(2,4,2);hold on;plot(trans(:,1,2),'-r');plot(trans(:,2,2),'-g');plot(trans(:,3,2),'-b');title('cam1 trans');grid on;
            subplot(2,4,3);hold on;plot(trans(:,1,3),'-r');plot(trans(:,2,3),'-g');plot(trans(:,3,3),'-b');title('cam2 trans');grid on;
            subplot(2,4,4);hold on;plot(trans(:,1,4),'-r');plot(trans(:,2,4),'-g');plot(trans(:,3,4),'-b');title('cam3 trans');grid on;
            subplot(2,4,5);hold on;plot(rot(:,1,1),'-r');plot(rot(:,2,1),'-g');plot(rot(:,3,1),'-b');title('cam0 rot');grid on;
            subplot(2,4,6);hold on;plot(rot(:,1,2),'-r');plot(rot(:,2,2),'-g');plot(rot(:,3,2),'-b');title('cam1 rot');grid on;
            subplot(2,4,7);hold on;plot(rot(:,1,3),'-r');plot(rot(:,2,3),'-g');plot(rot(:,3,3),'-b');title('cam2 rot');grid on;
            subplot(2,4,8);hold on;plot(rot(:,1,4),'-r');plot(rot(:,2,4),'-g');plot(rot(:,3,4),'-b');title('cam3 rot');grid on;
        else
            figure(1),subplot(1,2,1),hold on; grid on;plot(trans_cat(:,1),'-r');plot(trans_cat(:,2),'-g');plot(trans_cat(:,3),'-b'); title('trans');
                      subplot(1,2,2),hold on; grid on;plot(rot_cat(:,1),'-r');plot(rot_cat(:,2),'-g');plot(rot_cat(:,3),'-b'); title('rot');
        end
    catch
        sdfjgh = 1;
    end
end

if 0
    
    inputDir = 'G:\matlab\data\vio\without_epiplane';
    % inputDir = 'G:\matlab\data\vio\without_epiplane_good';
    % inputDir = 'G:\matlab\data\vio\without_epiplane_back1';
    % inputDir = 'G:\matlab\data\vio\without_epiplane_back2';
    % inputDir = 'G:\matlab\data\vio\without_epiplane_1';
    % inputDir = 'G:\matlab\data\vio\without_epiplane_5';
    % inputDir = 'G:\matlab\data\vio\without_epiplane_good2';
    inputDir = 'G:\matlab\data\vio\with_epiplane_calib/non_overlapping';
    inputDir = 'G:\matlab\data\vio\with_epiplane_calib\reproj_epiplane';
    inputDir = 'G:\matlab\data\vio\gt_2\epiplane_fix_tbc0';
    inputDir = 'G:\matlab\data\vio\gt_6\epiplane_fix_tbc0';
    inputDir = 'G:\matlab\data\vio\gt_2\epiplane_fix_tbc0';
    
    
    inputDir = 'G:\matlab\data\vio\with_epiplane_back2\with_prior\0.1';
    inputDir = 'G:\matlab\data\vio\with_epiplane_back0\with_prior\0.1';
    inputDir = 'G:\matlab\data\vio\gt_5\with_prior\0.1';
%     inputDir = 'G:\matlab\data\vio\gt_6\with_prior\0.1';
    
    inputDir = 'G:\matlab\data\vio\ci\004\on1';
    inputDir = 'G:\matlab\data\vio\gt_6\with_prior\on1';
    
    
    dirInfo = dir(fullfile(inputDir,'oc_time_*'));
    
    Trans = [];
    for i = 1 : size(dirInfo)
        folder = fullfile(inputDir,dirInfo(i).name) ;
        try
            [trans, trans_cat, rot, rot_cat] = testOnlineCalibration(folder, 0);
            
            if 0
                figure(2),subplot(2,4,1);hold on;plot(trans(:,1,1),'-r');plot(trans(:,2,1),'-g');plot(trans(:,3,1),'-b');title('cam0 trans');grid on;
                subplot(2,4,2);hold on;plot(trans(:,1,2),'-r');plot(trans(:,2,2),'-g');plot(trans(:,3,2),'-b');title('cam1 trans');grid on;
                subplot(2,4,3);hold on;plot(trans(:,1,3),'-r');plot(trans(:,2,3),'-g');plot(trans(:,3,3),'-b');title('cam2 trans');grid on;
                subplot(2,4,4);hold on;plot(trans(:,1,4),'-r');plot(trans(:,2,4),'-g');plot(trans(:,3,4),'-b');title('cam3 trans');grid on;
                subplot(2,4,5);hold on;plot(rot(:,1,1),'-r');plot(rot(:,2,1),'-g');plot(rot(:,3,1),'-b');title('cam0 rot');grid on;
                subplot(2,4,6);hold on;plot(rot(:,1,2),'-r');plot(rot(:,2,2),'-g');plot(rot(:,3,2),'-b');title('cam1 rot');grid on;
                subplot(2,4,7);hold on;plot(rot(:,1,3),'-r');plot(rot(:,2,3),'-g');plot(rot(:,3,3),'-b');title('cam2 rot');grid on;
                subplot(2,4,8);hold on;plot(rot(:,1,4),'-r');plot(rot(:,2,4),'-g');plot(rot(:,3,4),'-b');title('cam3 rot');grid on;
            else
                figure(2),subplot(1,2,1);hold on; grid on;plot(trans_cat(:,1),'-r');plot(trans_cat(:,2),'-g');plot(trans_cat(:,3),'-b');title('trans');
                          subplot(1,2,2);hold on; grid on;plot(rot_cat(:,1),'-r');plot(rot_cat(:,2),'-g');plot(rot_cat(:,3),'-b');title('rot');
            end
        catch
            safjg  =2;
        end
    end
    
end

end