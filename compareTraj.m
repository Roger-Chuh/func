function compareTraj()

close all

inputDirOld = 'G:\matlab\data\direct\sim\old';
inputDirNew = 'G:\matlab\data\direct\sim\new';

next_ = load('G:\matlab\data\direct\sim\output_next_5.txt');
old_ = load('G:\matlab\data\direct\sim\output_orca_5.txt');

[comm,a,b] = intersect(next_(:,1), old_(:,1));


next = next_(a,:);
old = old_(b,:);
poseNext = plotPose(next);
poseOld = plotPose(old);

figure,plotPath(poseNext);plotPath(poseOld);
figure,plot3(poseNext(:,10), poseNext(:,11), poseNext(:,12));
hold on;plot3(poseOld(:,10), poseOld(:,11), poseOld(:,12)); axis equal;legend('next','old');


poseCell = splitData(old,inputDirOld, 'old');
poseCell = splitData(next,inputDirNew, 'new');

end
function [poseMat, quatMat] = plotPose(Data)
poseMat = [];
quatMat = [];
for i = 1 : size(Data,1)
    data = Data(i,:);
    xyzw = data(5:8);
    trans = data(2:4);
    
%     R2 = quatern2rotMat(xyzw([4 1 2 3]));
    R = quat2rotm(xyzw([4 1 2 3]));
    if i == 1
        poseMat = [poseMat; [reshape(eye(3),1,9) [0 0 0]]];
        deltaT = inv([R trans';0 0 0 1]);
        quat = rotm2quat(eye(3));
        xyzw_new = quat([2 3 4 1]);
        quatMat = [quatMat; [data(1) [0 0 0] xyzw_new]];
    else
        T = deltaT * [R trans';0 0 0 1];
        poseMat = [poseMat; [reshape(T(1:3,1:3),1,9), T(1:3,4)']];
        quat = rotm2quat(T(1:3,1:3));
        xyzw_new = quat([2 3 4 1]);
        quatMat = [quatMat; [data(1) T(1:3,4)' xyzw_new]];
    end
%     rotm2quat(R)
%     rotMat2quatern(R2)
end
end
function poseCell = splitData(pose, inputDir, name)
data_ind = 1 : 100 : size(pose,1);

[~,quat] = plotPose(pose);

for i = 1 : length(data_ind)-1
    poseCell{i,1} = quat(data_ind(i):data_ind(i+1),:);
    write2file( poseCell{i,1}, inputDir, name, i);
end


end
function write2file(data, inputDir, name, id)
f_id=fopen(fullfile(inputDir,sprintf(strcat(name,'_%d.txt'),id)),'w');

for i = 1 : size(data,1)
    fprintf(f_id,'%f %f %f %f %f %f %f %f\n',data(i,1),data(i,2),data(i,3),data(i,4),data(i,5),data(i,6),data(i,7),data(i,8));
    
end
fclose(f_id);
end
