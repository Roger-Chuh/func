function checkIntr()
% close all;

aprilSize = 6;


use_van_plane = 0;1; 0;1;0; 1;0;


inputDir = 'D:\yvr\data\intr\outlier_img';
% inputDir = 'D:\yvr\data\intr\outlier_img1';
% inputDir = 'D:\yvr\data\intr\outlier_img2';
% inputDir = 'D:\yvr\data\intr\outlier_img3';
inputDir = 'D:\yvr\data\intr\outlier_img3_1';
inputDir = 'D:\yvr\data\intr\outlier_img4';
inputDir = 'D:\yvr\data\intr\outlier_img5';
inputDir = 'D:\yvr\data\intr\outlier_img5_1';
inputDir = 'D:\yvr\data\intr\outlier_img5_2';
% inputDir = 'D:\yvr\data\intr\outlier_img6';
% inputDir = 'D:\yvr\data\intr\outlier_img7';
% inputDir = 'D:\yvr\data\intr\outlier_img7_1';
% inputDir = 'D:\yvr\data\intr\outlier_img8';
% inputDir = 'D:\yvr\data\intr\outlier_img8_1';
% inputDir = 'D:\yvr\data\intr\outlier_img8_1_opt';
% inputDir = 'D:\yvr\data\intr\outlier_img9_1_opt';
% inputDir = 'D:\yvr\data\intr\outlier_img_9';
% % inputDir = 'D:\yvr\data\intr\outlier_img10';
% inputDir = 'D:\yvr\data\intr\outlier_img11';
inputDir = 'D:\yvr\data\intr\outlier_img12';
inputDir = 'D:\yvr\data\intr\outlier_img13';
inputDir = 'D:\yvr\data\intr\outlier_img14';
inputDir = 'D:\yvr\data\intr\outlier_img15';
inputDir = 'D:\yvr\data\intr\outlier_img16';
inputDir = 'D:\yvr\data\intr\outlier_img17';
inputDir = 'D:\yvr\data\intr\outlier_img18';
inputDir = 'D:\yvr\data\intr\outlier_img19';
inputDir = 'D:\yvr\data\intr\outlier_img22';
inputDir = 'D:\yvr\data\intr\outlier_img23';




inputDir = 'G:\matlab\data\outlier_img1';
inputDir = 'G:\matlab\data\outlier_img1_2';
% inputDir = 'G:\matlab\data\outlier_img2';
inputDir = 'G:\matlab\data\outlier_img2_2';
% inputDir = 'G:\matlab\data\outlier_img3';
inputDir = 'G:\matlab\data\outlier_img3_2';
inputDir = 'G:\matlab\data\outlier_img3_3';


inputDir = 'G:\matlab\data\outlier_img_rep';
inputDir = 'G:\matlab\data\outlier_img_dir';

inputDir = 'G:\matlab\data\outlier_img1_3';


inputDir = 'G:\matlab\data\outlier_img1_4_oneself';
inputDir = 'G:\matlab\data\outlier_img1_5_module';
inputDir = 'G:\matlab\data\outlier_img1_6_module';
inputDir = 'G:\matlab\data\outlier_img1_7_module';
inputDir = 'G:\matlab\data\outlier_img1_8_mode3';
inputDir = 'G:\matlab\data\outlier_img1_9_mode2';
inputDir = 'G:\matlab\data\outlier_img1_10_mode1';

inputDir = 'G:\matlab\data\outlier_img1_11_module_mode0';
inputDir = 'G:\matlab\data\outlier_img1_12_module_mode1';
inputDir = 'G:\matlab\data\outlier_img1_13_module_mode2';
inputDir = 'G:\matlab\data\outlier_img1_14_module_mode3';

inputDir = 'G:\matlab\data\outlier_img1_15_module_mode3_weight';
% inputDir = 'G:\matlab\data\outlier_img1_16_module_mode2_weight';
% inputDir = 'G:\matlab\data\outlier_img1_17_module_mode1_weight';
% inputDir = 'G:\matlab\data\outlier_img1_18_module_mode3_weight_rep';
% inputDir = 'G:\matlab\data\outlier_img1_19_module_mode3_weight2_rep';
inputDir = 'G:\matlab\data\outlier_img1_20_module_mode3_weight2_dir';
inputDir = 'G:\matlab\data\outlier_img1_21_module_mode3_weightsq_dir';
inputDir = 'G:\matlab\data\outlier_img1_22_oneself_mode3_weightsq_dir';
inputDir = 'G:\matlab\data\outlier_img1_23_module_mode2_weightsq_dir';
inputDir = 'G:\matlab\data\outlier_img1_24_MulOne';
inputDir = 'G:\matlab\data\outlier_img1_25';
inputDir = 'G:\matlab\data\outlier_img1_26';
% inputDir = 'G:\matlab\data\outlier_img27';
inputDir = 'G:\matlab\data\outlier_img28';
inputDir = 'G:\matlab\data\outlier_img29';
inputDir = 'G:\matlab\data\outlier_img30';
inputDir = 'G:\matlab\data\outlier_img31';
inputDir = 'G:\matlab\data\outlier_img32';
inputDir = 'G:\matlab\data\outlier_img33';
inputDir = 'G:\matlab\data\outlier_img34';
inputDir = 'G:\matlab\data\outlier_img35';
inputDir = 'G:\matlab\data\outlier_img36';
inputDir = 'G:\matlab\data\outlier_img37';
inputDir = 'G:\matlab\data\outlier_img38';
inputDir = 'G:\matlab\data\outlier_img39';
inputDir = 'G:\matlab\data\outlier_img40';
inputDir = 'G:\matlab\data\outlier_img41';
inputDir = 'G:\matlab\data\outlier_img42';
inputDir = 'G:\matlab\data\outlier_img43';
% inputDir = 'G:\matlab\data\outlier_img44';
% inputDir = 'G:\matlab\data\outlier_img45';
% inputDir = 'G:\matlab\data\outlier_img46';
% inputDir = 'G:\matlab\data\outlier_img47';
% inputDir = 'G:\matlab\data\outlier_img48';
% inputDir = 'G:\matlab\data\outlier_img49';
inputDir = 'G:\matlab\data\outlier_img50';
inputDir = 'G:\matlab\data\outlier_img51';
inputDir = 'G:\matlab\data\outlier_img52';
inputDir = 'G:\matlab\data\outlier_img53';
inputDir = 'G:\matlab\data\outlier_img54';

inputDir = 'G:\matlab\data\outlier_img55';
inputDir = 'G:\matlab\data\outlier_img56';
inputDir = 'G:\matlab\data\outlier_img57';
inputDir = 'G:\matlab\data\outlier_img58';
inputDir = 'G:\matlab\data\outlier_img59';


inputDir = 'G:\matlab\data\outlier_img60';
% inputDir = 'G:\matlab\data\outlier_img61';
% inputDir = 'G:\matlab\data\outlier_img62';
% inputDir = 'G:\matlab\data\outlier_img63';
% inputDir = 'G:\matlab\data\outlier_img64';

inputDir = 'G:\matlab\data\outlier_img65';
inputDir = 'G:\matlab\data\outlier_img66';
inputDir = 'G:\matlab\data\outlier_img67';
inputDir = 'G:\matlab\data\outlier_img68';
inputDir = 'G:\matlab\data\outlier_img69';

inputDir = 'G:\matlab\data\dump\008\dump1\outlier_img_005';





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputDir = 'G:\matlab\data\fac\outlier_img';







%%%%%%%%%%%%%%%%%%%%% NEW TRACE %%%%%%%%%%%%%%%%
inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_1'; %default
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_2'; %use weight
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_3'; %use weight opt intr
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_4'; %opt intr
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_5'; %use weight new opt intr
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_6'; %use weight new
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_7'; %use weight new2

% inputDir = 'G:\matlab\data\new_trace\outlier_img_ce_1'; %use weight new2
% % inputDir = 'G:\matlab\data\new_trace\outlier_img_ce_2'; %use weight new2 opt intr
% % inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_8'; %use weight new2 opt intr
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ce_3'; %use weight new
% % inputDir = 'G:\matlab\data\new_trace\outlier_img_ce_4'; %default
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ce_5'; %use weight new0.5
% % inputDir = 'G:\matlab\data\new_trace\outlier_img_ce_6'; %use weight new0.25
% 
% inputDir = 'G:\matlab\data\new_trace\outlier_img_ju_9'; %use weight new0.5




inputDir = 'G:\matlab\data\outlier_img_d1';  % six edge intr opt no weight
inputDir = 'G:\matlab\data\new_trace\outlier_img_hao';  % six edge intr opt no weight

inputDir = 'G:\matlab\data\new_trace\outlier_img_hao_2';  % six edge intr opt no weight new MulPlateOneCam
inputDir = 'G:\matlab\data\new_trace\outlier_img_hao_3';  % six edge intr opt no weight new MulPlateOneCam close scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_hao_4';  % four edge intr opt no weight new MulPlateOneCam close scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_hao_5';  % six edge intr opt no weight new MulPlateOneCam far scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_hao_6';  % four edge intr opt no weight new MulPlateOneCam far scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_hao_7';  % four edge intr opt with weight new MulPlateOneCam far scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_fail1_1'; % four edge intr opt no weight new MulPlateOneCam old scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_fail1_2'; % six edge intr opt no weight new MulPlateOneCam old scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_fail2_1'; % four edge intr opt no weight new MulPlateOneCam old scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_fail2_2'; % six edge intr opt no weight new MulPlateOneCam old scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_fail2_3'; % four edge no intr opt no weight new MulPlateOneCam old scene
inputDir = 'G:\matlab\data\new_trace\outlier_img_fail2_4'; % six edge no intr opt no weight new MulPlateOneCam old scene
inputDir = 'G:\matlab\data\new_trace\hao\outlier_img_hao_3'; % six edge with intr opt no weight new MulPlateOneCam handheld scene
% inputDir = 'G:\matlab\data\new_trace\hao\outlier_img_hao_4_rep'; % six edge with intr opt no weight new MulPlateOneCam handheld scene

inputDir = 'G:\matlab\data\new_trace\outlier_img_close_1';
inputDir = 'G:\matlab\data\close\ning\outlier_img_ning_1';
inputDir = 'G:\matlab\data\close\hao\Device_D2HD232328D9000859_2022_08_01_19_10_57\controller_2022-08-01-19-08-13\results\outlier_img';
inputDir = 'G:\matlab\data\close\hao\Device_D2HD232328D9000859_2022_08_01_19_10_57\controller_2022-08-01-19-08-13\results\outlier_img_4';
inputDir = 'G:\matlab\data\close\hao\outlier_img_four_edge_no_eq';
inputDir = 'G:\matlab\data\close\hao\outlier_img_six_edge_no_eq';

inputDir = 'G:\matlab\data\close\new_tag\hui\Device_D2HD432713D9000609_2022_08_16_17_07_46\controller_2022-08-15-20-09-14\results\outlier_img';

inputDir = 'G:\matlab\data\outlier_img_dvt5';
inputDir = 'G:\matlab\data\outlier_img_dvt5_2';
inputDir = 'G:\matlab\data\outlier_img_dvt5_3';




inputDir = 'G:\matlab\data\outlier_img_bug_camera2';




%%
inputDir = 'G:\matlab\data\intr\outlier_img_check_extend_1';
% inputDir = 'G:\matlab\data\intr\outlier_img_check_1';

inputDir = 'G:\matlab\data\intr\outlier_img_no_epi_check_extend_1';
% inputDir = 'G:\matlab\data\intr\outlier_img_no_epi_check_1';


aprilSize = 10;




dirInfo = dir(fullfile(inputDir, 'matlab_corners_cam_*.txt'));

dataNum = length(dirInfo);
ids = 0 : aprilSize*aprilSize*4-1;
idsMat = reshape(ids, aprilSize*2,[]);

fitErr=[];
errMat = zeros(2*480, 2*640);
cntMat = zeros(2*480, 2*640);
eachCam = {[0] [];[0] [];[0] [];[0] []};
eachK = {};
kk = 1;
for k = 1 : dataNum
    
    data = load(fullfile(inputDir, dirInfo(k).name));
    camId = str2num(dirInfo(k).name(21))+1;
    if 0
        if camId == 1
            eachCam{1,1} = eachCam{1,1} + 1;
            eachCam{1,2} = [eachCam{1,2}; k];
        end
        if camId == 2
            eachCam{2,1} = eachCam{2,1} + 1;
            eachCam{2,2} = [eachCam{2,2}; k];
        end
        if camId == 3
            eachCam{3,1} = eachCam{3,1} + 1;
            eachCam{3,2} = [eachCam{3,2}; k];
        end
        if camId == 4
            eachCam{4,1} = eachCam{4,1} + 1;
            eachCam{4,2} = [eachCam{4,2}; k];
        end
    end
%     data = [data data(:,1:2) data(:,9:12)];
    K = eye(3);
    K(1,1) = data(1,9);
    K(2,2) = data(1,10);
    K(1,3) = data(1,11);
    K(2,3) = data(1,12);
    
    K_before = eye(3);
    K_before(1,1) = data(1,15);
    K_before(2,2) = data(1,16);
    K_before(1,3) = data(1,17);
    K_before(2,3) = data(1,18);
    eachK{camId,1} = K;
    eachK{camId,2} = K_before;
    horiLine = {};
    horiCnt = 1;
    vertLine = {};
    vertCnt = 1;
    hori_line_fit_err = [];
    vert_line_fit_err = [];
    hori_paras = [];
    vert_paras = [];
    for i = 1:size(idsMat,1)
        hori_all = idsMat(i,:);
        hori_candi = [];
        
        for j = 1 : size(data,1)
            if ismember(data(j,6),hori_all)
                hori_candi = [hori_candi; data(j,[1:2 7:8 13:14])];
            end
        end
        if size(hori_candi,1) > 8
            if use_van_plane
                [hori_err, hori_para] = fitline_plane(hori_candi(:,1:2), K);
                [hori_err_before, hori_para_before] = fitline_plane(hori_candi(:,5:6), K_before);
            else
                [hori_err, hori_para] = fitline(hori_candi(:,1:2));
                [hori_err_before, hori_para_before] = fitline(hori_candi(:,5:6));
            end
            hori_paras(horiCnt,:) = [hori_para hori_para_before];
            hori_line_fit_err(horiCnt,:) = [mean(abs(hori_err)) mean(abs(hori_err_before))];
            hori_ind = sub2ind(size(errMat), round(hori_candi(:,4)), round(hori_candi(:,3)));
            errMat(hori_ind) = errMat(hori_ind)+abs(hori_err)';
            cntMat(hori_ind) = cntMat(hori_ind) + 1;
            horiLine{horiCnt,1} = hori_candi;
            horiCnt = horiCnt + 1;
        end
    end
    fitErr = [fitErr; hori_line_fit_err];
    for i = 1:size(idsMat,2)
        vert_all = idsMat(:,i);
        vert_candi = [];
        for j = 1 : size(data,1)
            if ismember(data(j,6),vert_all)
                vert_candi = [vert_candi; data(j,[1:2 7:8 13:14])];
            end
        end
        if size(vert_candi,1) > 8
            if use_van_plane
                [vert_err, vert_para] = fitline_plane(vert_candi(:,1:2), K);
                [vert_err_before, vert_para_before] = fitline_plane(vert_candi(:,5:6), K_before);
            else
                [vert_err, vert_para] = fitline(vert_candi(:,1:2));
                [vert_err_before, vert_para_before] = fitline(vert_candi(:,5:6));
            end
            vert_paras(vertCnt,:) = [vert_para vert_para_before];
            vert_line_fit_err(vertCnt,:) = [mean(abs(vert_err)) mean(abs(vert_err_before))];
            vert_ind = sub2ind(size(errMat), round(vert_candi(:,4)), round(vert_candi(:,3)));
            errMat(vert_ind) = errMat(vert_ind)+abs(vert_err)';
            cntMat(vert_ind) = cntMat(vert_ind) + 1;
            vertLine{vertCnt,1} = vert_candi;
            vertCnt = vertCnt + 1;
        end
    end
    fitErr = [fitErr; vert_line_fit_err];
    kk = k;
    if use_van_plane
        try
            [~,~,hori_van] = svd([pextend(hori_paras(:,1:3)')'; [0 0 0 1]], 0);  % 扇形平面的法向应该共面
        catch
            continue;
            skags = 1;
        end
        
        try
            [~,~,vert_van] = svd([pextend(vert_paras(:,1:3)')'; [0 0 0 1]], 0);
        catch
            continue;
        end
        hori_van = hori_van(:,4);
        hori_van = hori_van./norm(hori_van(1:3));
        hori_van1 = hori_van./hori_van(4);
        hori_van_dir_err(kk,1) = sum(abs(dot(repmat(hori_van,1,size(hori_paras,1)), [pextend(hori_paras(:,1:3)')])));
        % hori_van_fir_err = ((dot(repmat(hori_van1,1,size(hori,1)), hori')))
        [~,~,hori_van_before] = svd([pextend(hori_paras(:,5:7)')'; [0 0 0 1]], 0);
        hori_van_before = hori_van_before(:,4);
        hori_van_before = hori_van_before./norm(hori_van_before(1:3));
        hori_van_before1 = hori_van_before./hori_van_before(4);
        hori_van_dir_err(kk,2) = sum(abs(dot(repmat(hori_van_before,1,size(hori_paras,1)), [pextend(hori_paras(:,5:7)')])));
        
        
        
        [~,~,vert_van] = svd([pextend(vert_paras(:,1:3)')'; [0 0 0 1]], 0);
        vert_van = vert_van(:,4);
        vert_van = vert_van./norm(vert_van(1:3));
        vert_van1 = vert_van./vert_van(4);
        vert_van_dir_err(kk,1) = sum(abs(dot(repmat(vert_van,1,size(vert_paras,1)),  [pextend(vert_paras(:,1:3)')])));
        % vert_van_fir_err = ((dot(repmat(vert_van1,1,size(vert,1)), vert')))
        [~,~,vert_van_before] = svd([pextend(vert_paras(:,5:7)')'; [0 0 0 1]], 0);
        vert_van_before = vert_van_before(:,4);
        vert_van_before = vert_van_before./norm(vert_van_before(1:3));
        vert_van_before1 = vert_van_before./vert_van_before(4);
        vert_van_dir_err(kk,2) = sum(abs(dot(repmat(vert_van_before,1,size(vert_paras,1)), [pextend(vert_paras(:,5:7)')])));
        
        
        hori_dir = inv(K)*hori_van(1:3);
        hori_dir = hori_dir./norm(hori_dir);
        vert_dir = inv(K)*vert_van(1:3);
        vert_dir = vert_dir./norm(vert_dir);
        
        hori_dir_before = inv(K_before)*hori_van_before(1:3);
        hori_dir_before = hori_dir_before./norm(hori_dir_before);
        vert_dir_before = inv(K_before)*vert_van_before(1:3);
        vert_dir_before = vert_dir_before./norm(vert_dir_before);
        
        err1(kk,:) = [dot(hori_van(1:3), vert_van(1:3)) dot(hori_van_before(1:3), vert_van_before(1:3))];
        err2(kk,:) = [CalcDegree(hori_van(1:3)', vert_van(1:3)') CalcDegree(hori_van_before(1:3)', vert_van_before(1:3)')];
        
    else
        
        try
            [~,~,hori_van] = svd(hori_paras(:,1:3));  % 求平行线的角点， 灭点
        catch
            continue;
            skags = 1;
        end
        
        try
            [~,~,vert_van] = svd(vert_paras(:,1:3));
        catch
            continue;
        end
        hori_van = hori_van(:,3);
        hori_van = hori_van./norm(hori_van(1:2));
        hori_van1 = hori_van./hori_van(3);
        hori_van_dir_err(kk,1) = sum(abs(dot(repmat(hori_van1,1,size(hori_paras,1)), hori_paras(:,1:3)')));
        % hori_van_fir_err = ((dot(repmat(hori_van1,1,size(hori,1)), hori')))
        [~,~,hori_van_before] = svd(hori_paras(:,4:6));
        hori_van_before = hori_van_before(:,3);
        hori_van_before = hori_van_before./norm(hori_van_before(1:2));
        hori_van_before1 = hori_van_before./hori_van_before(3);
        hori_van_dir_err(kk,2) = sum(abs(dot(repmat(hori_van_before1,1,size(hori_paras,1)), hori_paras(:,4:6)')));
        
        
        
        [~,~,vert_van] = svd(vert_paras(:,1:3));
        vert_van = vert_van(:,3);
        vert_van = vert_van./norm(vert_van(1:2));
        vert_van1 = vert_van./vert_van(3);
        vert_van_dir_err(kk,1) = sum(abs(dot(repmat(vert_van1,1,size(vert_paras,1)), vert_paras(:,1:3)')));
        % vert_van_fir_err = ((dot(repmat(vert_van1,1,size(vert,1)), vert')))
        [~,~,vert_van_before] = svd(vert_paras(:,4:6));
        vert_van_before = vert_van_before(:,3);
        vert_van_before = vert_van_before./norm(vert_van_before(1:2));
        vert_van_before1 = vert_van_before./vert_van_before(3);
        vert_van_dir_err(kk,2) = sum(abs(dot(repmat(vert_van_before1,1,size(vert_paras,1)), vert_paras(:,4:6)')));
        
        
        hori_dir = inv(K)*hori_van1;
        hori_dir = hori_dir./norm(hori_dir);
        vert_dir = inv(K)*vert_van1;
        vert_dir = vert_dir./norm(vert_dir);
        
        hori_dir_before = inv(K_before)*hori_van_before1;
        hori_dir_before = hori_dir_before./norm(hori_dir_before);
        vert_dir_before = inv(K_before)*vert_van_before1;
        vert_dir_before = vert_dir_before./norm(vert_dir_before);
        
        err1(kk,:) = [dot(hori_dir, vert_dir) dot(hori_dir_before, vert_dir_before)];
        err2(kk,:) = [CalcDegree(hori_dir', vert_dir') CalcDegree(hori_dir_before', vert_dir_before')];
        
    end
    
%     kk = kk+1;
if camId == 1
       eachCam{1,1} = eachCam{1,1} + 1;
       eachCam{1,2} = [eachCam{1,2}; k];
    end
    if camId == 2
        eachCam{2,1} = eachCam{2,1} + 1;
        eachCam{2,2} = [eachCam{2,2}; k];
    end
    if camId == 3
        eachCam{3,1} = eachCam{3,1} + 1;
        eachCam{3,2} = [eachCam{3,2}; k];
    end
    if camId == 4
        eachCam{4,1} = eachCam{4,1} + 1;
        eachCam{4,2} = [eachCam{4,2}; k];
    end
end

K_err = [(eachK{1,1} - eachK{1,2}) (eachK{2,1} - eachK{2,2}) (eachK{3,1} - eachK{3,2}) (eachK{4,1} - eachK{4,2})]

thr = 150;

goodness = find((hori_van_dir_err(:,1) < thr) & ( hori_van_dir_err(:,1) > 0) & (vert_van_dir_err(:,1) < thr) & (hori_van_dir_err(:,2) < thr) & (vert_van_dir_err(:,2) < thr));
% figure,hist(err2(goodness,1));hold on;hist(err2(goodness,2));
% figure,hist(err2(goodness,:), 300);legend('after','before');
avgErrMat = errMat./cntMat;

% figure,hist(fitErr,100);legend('after','before');
% figure,plot(sort(fitErr(:,1)));hold on;plot(sort(fitErr(:,2)));legend('after','before');

goodness1 = find((hori_van_dir_err(eachCam{1,2},1) < thr) & (hori_van_dir_err(eachCam{1,2},1) > 0) & (vert_van_dir_err(eachCam{1,2},1) < thr));
goodness2 = find((hori_van_dir_err(eachCam{2,2},1) < thr) & (hori_van_dir_err(eachCam{2,2},1) > 0) & (vert_van_dir_err(eachCam{2,2},1) < thr));
goodness3 = find((hori_van_dir_err(eachCam{3,2},1) < thr) & (hori_van_dir_err(eachCam{3,2},1) > 0) & (vert_van_dir_err(eachCam{3,2},1) < thr));
goodness4 = find((hori_van_dir_err(eachCam{4,2},1) < thr) & (hori_van_dir_err(eachCam{4,2},1) > 0) & (vert_van_dir_err(eachCam{4,2},1) < thr));


angErr1 = err2(eachCam{1,2},:);
angErr2 = err2(eachCam{2,2},:);
angErr3 = err2(eachCam{3,2},:);
angErr4 = err2(eachCam{4,2},:);


angThr = 88;

a1 = filterResult(angErr1(goodness1,:), angThr);
a2 = filterResult(angErr2(goodness2,:), angThr);
a3 = filterResult(angErr3(goodness3,:), angThr);
a4 = filterResult(angErr4(goodness4,:), angThr);


% figure,subplot(2,2,1);hist(angErr1(goodness1,:),10);title('camera 0');legend('after','before');
% subplot(2,2,2);hist(angErr2(goodness2,:),10);title('camera 1');legend('after','before');
% subplot(2,2,3);hist(angErr3(goodness3,:),10);title('camera 2');legend('after','before');
% % subplot(2,2,4);hist(angErr4(goodness4,:),10);title('camera 3');legend('after','before');

figure,subplot(2,2,1);hist(a1,100);title('camera 0');legend('after','before');
subplot(2,2,2);hist(a2,100);title('camera 1');legend('after','before');
subplot(2,2,3);hist(a3,100);title('camera 2');legend('after','before');
subplot(2,2,4);hist(a4,100);title('camera 3');legend('after','before');



end
function [hori_line_fit_err, line_para] = fitline(hori_line)

hori_line(:,3) = 1;
[~,~,V] = svd(hori_line);
line_para = V(:,3)';
line_para = sign(line_para(3)).*line_para;
line_para = line_para./norm(line_para(1:2));
hori_line_fit_err = dot(repmat(line_para', 1 ,size(hori_line,1)), hori_line');

end
function [hori_plane_fit_err, planeNorm] = fitline_plane(hori_line, K)


pt_bearing = (inv(K)*pextend(hori_line'))';

pt_bearing(:,4) = 1;

 [~,~,plane_norm] = svd([pt_bearing; [0 0 0 1]], 0);
planeNorm = plane_norm(:,4)';
planeNorm = planeNorm./norm(planeNorm(1:3));
planeNorm = sign(planeNorm(3)).*planeNorm;


hori_plane_fit_err = dot(repmat(planeNorm', 1 ,size(pt_bearing,1)), pt_bearing');


if 0
    hori_line(:,3) = 1;
    [~,~,V] = svd(hori_line);
    line_para = V(:,3)';
    line_para = sign(line_para(3)).*line_para;
    line_para = line_para./norm(line_para(1:2));
    hori_line_fit_err = dot(repmat(line_para', 1 ,size(hori_line,1)), hori_line');
end

end
function aa = filterResult(a, angThr)

% a = angErr3(goodness3,:);
b = find(a(:,1)>85);
aa = a(b,:);
% figure,hist(a(b,:),100);

end