function SimulateEpiSearchError()



close all



K = [400 0 320; 0 400 240; 0 0 1];
T = [rodrigues([0.1 -0.2 0.3]) [100 200 -300]'; 0 0 0 1];
xyz = [200 300 1200];
xyz = xyz./xyz(3);
depth = 1200;

Twc = {};

Tbc = [eye(3) [20 -20 100]'; 0 0 0 1];


angList = 10:10:60;
T0 = eye(4);
TwcMat =[];
transList = 1000.*rand(length(angList),3);
for i = 1 : length(angList)
    RRR = roty(angList(i));
    T_delt = [RRR transList(i,:)'; 0 0 0 1];
    k2c_comp{i,1} = inv(Tbc) * T_delt * Tbc;
    TwcMat = [TwcMat; [reshape(k2c_comp{i,1}(1:3,1:3),1,9) k2c_comp{i,1}(1:3,4)']];
    if i == 1
        deltaPoseMat = [reshape(eye(3), 1, 9) 0 0 0];
        T1 = inv(Tbc)*T_delt*T0;
        T1_0 = T1;
    else
        T1 = inv(Tbc)*T_delt*T0;
        deltaPose = inv(T1_0)*T1;
        deltaPoseMat = [deltaPoseMat; [reshape(deltaPose(1:3,1:3),1,9) deltaPose(1:3,4)']];
    end
    
    %     TwcMat = [TwcMat; [reshape(T1(1:3,1:3),1,9) T1(1:3,4)']];
end


host_pix = [50 100 1];
host_pix_unproj = pflat(inv(K)*host_pix')';
depth = 1500;
xyz = depth.*host_pix_unproj;

depth_err = -200 : 50 : 200;
depth_err = zeros(100,1);
for i = 1 : length(angList)
    
    [ptIcsGT(i,:), tgtPt3d] = TransformAndProject((depth).*host_pix_unproj, K, k2c_comp{i,1}(1:3, 1:3), k2c_comp{i,1}(1:3, 4));
    [ptIcs(i,:), tgtPt3d] = TransformAndProject((depth+depth_err(i)).*host_pix_unproj, K, k2c_comp{i,1}(1:3, 1:3), k2c_comp{i,1}(1:3, 4));
    F = inv(K')*SkewSymMat(k2c_comp{i,1}(1:3,4)) * k2c_comp{i,1}(1:3,1:3)*inv(K);
    epiline_target_in_host(:,i) = F'*[ptIcs(i,:) 1]';
    epiline_target_in_host(:,i) = epiline_target_in_host(:,i)./norm(epiline_target_in_host(1:2,i));
    
    
    
    
    target_unproject = inv(K)*[ptIcs(i,:) 1]';
    
    
    
end


[a,b,c] = svd(epiline_target_in_host*epiline_target_in_host');

intersection = c(:,3)./c(3,3);

err = host_pix - intersection';

Twc{1,1} = eye(4);
trace = [host_pix(1:2); ptIcs];
traceGT = [host_pix(1:2); ptIcsGT];



figure, imshow(ones(480, 640)); hold on; plot(host_pix(1), host_pix(2),'or');plot(ptIcs(:,1), ptIcs(:, 2),'x-g');
points = lineToBorderPoints(epiline_target_in_host', [480 640]);
line(points(:,[1,3])',points(:,[2,4])','Color', [0 0 1]);
plot(intersection(1),intersection(2),'xb');
plot(traceGT(1,1),traceGT(1,2),'sm');


for i = 1 : length(angList)
    Twc{1+i,1} = inv(k2c_comp{i,1});
    
    
end



for k = 1 : size(trace,1)
    base_id = k;
    
    cnt = 1;
    for i = 1 : length(Twc)
        
        if i == base_id
            continue;
        end
        
        T_th = inv(Twc{i})*Twc{base_id};
        F = inv(K')*SkewSymMat(T_th(1:3,4)) * T_th(1:3,1:3)*inv(K);
        E = SkewSymMat(T_th(1:3,4)) * T_th(1:3,1:3);
        epiplane_target_in_host2(:,cnt) = E'*(inv(K)*[trace(i,:) 1]');
        epiplane_target_in_host2(:,cnt) = epiplane_target_in_host2(:,cnt)./norm(epiplane_target_in_host2(:,cnt));
        epiline_target_in_host2(:,cnt) = F'*[trace(i,:) 1]';
        epiline_target_in_host2(:,cnt) = epiline_target_in_host2(:,cnt)./norm(epiline_target_in_host2(1:2,cnt));
        
        
        cnt = cnt + 1;
        
    end
    
    [a2,b2,c2] = svd(epiline_target_in_host2*epiline_target_in_host2');
%     [epi_norm, ~] = NormalizeVector(epiline_target_in_host2');
    [a3,b3,c3] = svd(epiplane_target_in_host2*epiplane_target_in_host2');
    intersection2 = c2(:,3)./c2(3,3);
    
    quality = b3(2,2)/sum(sum(b3));
    
    err2 = [k trace(base_id,:) - intersection2(1:2)'  traceGT(base_id,:) - intersection2(1:2)' quality]
%     trace - traceGT
    
    
    
    figure, imshow(ones(480, 640)); hold on; plot(trace(base_id,1), trace(base_id,2),'or');plot(trace(([1:base_id-1 base_id+1:end]),1), trace(([1:base_id-1 base_id+1:end]),2),'x-g');
    points2 = lineToBorderPoints(epiline_target_in_host2', [480 640]);
    line(points2(:,[1,3])',points2(:,[2,4])','Color', [0 0 1]);
    plot(intersection2(1),intersection2(2),'xb');
    plot(traceGT(base_id,1),traceGT(base_id,2),'sm');
    
end


end