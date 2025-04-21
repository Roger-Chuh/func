function LoadYVRDataset()
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\SLAMRecord\controller_2024-12-13-01-15-45\IMU_QCOM\data.json';
% inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\SLAMRecord\controller_2024-10-08-05-05-48\IMU0\data.json';
inputDir2 = 'G:\matlab\data\direct\gt\D2_011\4\tbc\data_head.json';
inputDir = 'G:\matlab\data\direct\gt\D2_011\4\tbc\data_controller.json';
json_head = loadjson(fullfile(inputDir2));
json_controller = loadjson(fullfile(inputDir));

imudata_head = [];
for i = 1 : length(json_head.Sequence.Dataset.Data)
    gyro = [json_head.Sequence.Dataset.Data(i).g_x json_head.Sequence.Dataset.Data(i).g_y json_head.Sequence.Dataset.Data(i).g_z];
    acc = [json_head.Sequence.Dataset.Data(i).a_x json_head.Sequence.Dataset.Data(i).a_y json_head.Sequence.Dataset.Data(i).a_z];
    
%     gyro = gyro .* 2000 .* pi ./ 180 ./ 32768 .* 2.0;
%     acc = acc .* 9.807 ./ 2048 .* 2.0;
    imudata_head = [imudata_head;[json_head.Sequence.Dataset.Data(i).timestamp * 1e-9 gyro acc]];
end
imudata_controller = [];
for i = 1 : length(json_controller.Sequence.Dataset.Data)
    gyro = [json_controller.Sequence.Dataset.Data(i).g_x json_controller.Sequence.Dataset.Data(i).g_y json_controller.Sequence.Dataset.Data(i).g_z];
    acc = [json_controller.Sequence.Dataset.Data(i).a_x json_controller.Sequence.Dataset.Data(i).a_y json_controller.Sequence.Dataset.Data(i).a_z];
    
    gyro = gyro .* 2000 .* pi ./ 180 ./ 32768 .* 2.0;
    acc = acc .* 9.807 ./ 2048 .* 2.0;
    imudata_controller = [imudata_controller;[json_controller.Sequence.Dataset.Data(i).timestamp * 1e-9 gyro acc]];
end

figure,plot(imudata_head(:,2:4))
figure,plot(imudata_head(:,5:7))
figure,plot(1000 * diff(imudata_head(:,1)))

figure,plot(imudata_controller(:,2:4))
figure,plot(imudata_controller(:,5:7))
figure,plot(1000 * diff(imudata_controller(:,1)))







sensor_number = [];
for i = 1 : length(json_controller.Sequence.Dataset.Data)
sensor_number = [sensor_number; [json_controller.Sequence.Dataset.Data(i).sensor_number json_controller.Sequence.Dataset.Data(i).arrive_timestamp * 1e-6 json_controller.Sequence.Dataset.Data(i).timestamp * 1e-6]];
end



sensor_number = unique(sensor_number);

json1 = loadjson(fullfile(inputDir2));
imu = json.Sequence.Dataset.Data;

acc = [];
gyro = [];
acc_pred = [];
gyro_pred = [];
for i = 1 : length(imu)
    acc = [acc; [imu(i).a_x imu(i).a_y imu(i).a_z]];
    gyro = [gyro; [imu(i).g_x imu(i).g_y imu(i).g_z]];
    if (i >= 12 )
        acc_pred = [acc_pred; [SavitzkyGolayDerivative12(acc(i-11:i,:))]];
        gyro_pred = [gyro_pred; [SavitzkyGolayDerivative12(gyro(i-11:i,:))]];
    end
end

figure,subplot(2,1,1);plot(acc);title('acc');
subplot(2,1,2);plot(gyro);title('gyro');
figure,subplot(2,1,1);plot(acc_pred);title('acc pred');
subplot(2,1,2);plot(gyro_pred);title('gyro pred');

figure,subplot(2,1,1);hold on;plot(gyro);title('gyro');
subplot(2,1,2);hold on;plot(gyro_pred);title('gyro accel');

end
function ret = SavitzkyGolayDerivative12(PeekBack)


ret = 1000 * (PeekBack(1,:) * 0.03846 + PeekBack(2,:) * 0.03147 + PeekBack(3,:) * 0.02448 + PeekBack(4,:) * 0.01748 +...
    PeekBack(5,:) * 0.01049 + PeekBack(6,:) * 0.0035 - PeekBack(7,:) * 0.0035 - PeekBack(8,:) * 0.01049 -...
    PeekBack(9,:) * 0.01748 - PeekBack(10,:) * 0.02448 - PeekBack(11,:) * 0.03147 - PeekBack(12,:) * 0.03846);


end