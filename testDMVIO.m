function testDMVIO()

inputDir = 'G:\matlab\data\direct\gt\ke\3';


Gvec = [0:255];

min_=Gvec(1);
max_=Gvec(end);

  for i=1 : length(Gvec)
      Gvec(i) = 255.0 * (Gvec(i) - min_) / (max_-min_);		
  end


 for i = 0:length(Gvec)-1
 G1(i+1)=255.0*i/max_;
 end
 
veg = imread(fullfile(inputDir, 'vegnette.png'));
 
w = 512;
h = 512;
vm16 = zeros(1,w*h);

  maxV=0;
		for i=1:h
            for j = 1:w
            if(vm16->at(i) > maxV) maxV = vm16->at(i);

            end
		for(int i=0;i<w*h;i++)
			vignetteMap[i] = vm16->at(i) / maxV;
 
end