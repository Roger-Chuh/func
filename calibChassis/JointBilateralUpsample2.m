function upSampled = JointBilateralUpsample2(high, low, coarseMap, varargin)


if ndims(high) < 3
    high = double(cat(3, high, high, high));
    low = double(cat(3, low, low, low));
else
    high = double(high);
    low = double(low);
end

width = size(high, 2);
height = size(high, 1);
factor = 2;  % �������������ͼ��߶�֮��Ϊ2

if (nargin == 3)
    halfWindow = 5;
    sigma_d = 0.5;
    sigma_r = 0.1;
elseif (nargin == 4)
    halfWindow = varargin{1};
    sigma_d = 0.5;
    sigma_r = 0.1;
elseif (nargin == 5)
    halfWindow = varargin{1};
    sigma_d =  varargin{2};
    sigma_r = 0.1;
elseif (nargin == 6)
    halfWindow = varargin{1};
    sigma_d = varargin{2};
    sigma_r =  varargin{3};
else
    error('Too many input arguments');
end

type = 1;
if ndims(coarseMap) < 3
    upSampled = zeros(size(high(:,:,1)));
else
    type = 0;
    upSampled = zeros([size(high)]);
end

for i = 1:height
    for j = 1 : width
        p = high(i,j,:);
        low_i = i / factor;
        low_j = j / factor;
        iMax = floor(min(size(low,1) - 0, low_i + halfWindow));
        iMin = ceil(max(1, low_i - halfWindow));
        jMax = floor(min(size(low,2) - 0, low_j + halfWindow));
        jMin = ceil(max(1, low_j - halfWindow));
        
        lowWindow = coarseMap(iMin:iMax,jMin:jMax,:);
        if ndims(coarseMap) < 3
%             assert(min(lowWindow(:)) > 0);
        end
        
        iw = (iMin:iMax) - low_i;
        jw = (jMin:jMax) - low_j;
        [mx,my] = meshgrid(jw,iw);
        spatial = exp( -(mx.^2 + my.^2) ./ (2*sigma_d.^2) );
        
        %         highWindow������P�ڸ߷ֱ����µ�֧�Ŵ��ڣ�����ͷֱ���ʱ����Ϊ5*5
        %         �߷ֱ�����Ӧ�þ���10*10������ʽ��Ҫ�����һ�£���˸߷ֱ����´���Ҳ��5*5
        %         �����Ҫ�Դ��ڵ����ݽ���������һ�в���һ�У���һ�в���һ��
        highWindow = high(iMin * factor:factor:iMax * factor,jMin * factor:factor:jMax * factor,:);
        
        dR = highWindow(:,:,1) - p(1);
        dG = highWindow(:,:,2) - p(2);
        dB = highWindow(:,:,3) - p(3);
        range = exp( -(dR.^2 + dG.^2 + dB.^2) ./ (2*sigma_r.^2));
        
        
        spatial_range = spatial .* range;
        Kp = sum(spatial_range(:));
        
        if type
            depth = lowWindow.*spatial_range;
            depth = sum(depth(:))./Kp;
            upSampled(i,j) = depth;
            
        else
            
            Normal = lowWindow.*cat(3, spatial_range, spatial_range, spatial_range);
            normal = [sum(sum(Normal(:,:,1))); sum(sum(Normal(:,:,2))); sum(sum(Normal(:,:,3)))]./Kp;
            upSampled(i,j,:) = normal;
            
        end
        
    end
end

end
