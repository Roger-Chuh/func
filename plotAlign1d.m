function plotAlign1d(A)
% a = load('G:\matlab\data\direct\sim\align1d.txt');

% close all;

kf_num = zeros(size(A,1),1);

str = 'kf num: ';

figure,

for i = 1 : size(A,1)
    if(isempty(A{i,1}))
       continue; 
    end
   
    a = A{i,1};
    
    
    str = strcat(str, num2str(sum(a(:,2))));
    if(i < size(A,1))
        str = strcat(str, '-');
    end
    subplot(3,2,1);hold on; plot(a(:,3)); grid on;
    subplot(3,2,2);hold on; plot(a(:,4)); grid on; title('potential');
    subplot(3,2,3);hold on; plot(a(:,5)); grid on; title('stable');
    subplot(3,2,4);hold on; plot(a(:,6)); grid on; title('align1d success count');
    subplot(3,2,5);hold on; plot(10.*a(:,2));grid on;title('keyframe');
    
end

 subplot(3,2,1); title(str);
end