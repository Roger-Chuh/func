function [c,r,w,h] = FindLargestRectanglesPUS(I, crit, minSize)
% finds largest rectangle regions within all points set to 1.
% input: I       - B/W boolean matrix or output of FindLargestSquares
%        minSize - [height width] - minimum width and height on regions of 
%                  interest (used to restrict final choise)
%        crit    - Optimazation Criteria parameters to optimize:
%                   crit(1)*height + crit(2)*width + crit(3)*height*width
% output (changed by P.Seibold): 
%         c    -  left column of largest all-white rectangle of the image
%               
%         r    - top row of largest all-white rectangle of the image
%                 
%         w    - width of largest all-white rectangle of the image
%                 
%         h    - height of largest all-white rectangle of the image
%Other modifications by P.Seibold:
%   FindLargestSquares(I) from Jaroslaw Tuszynski is included in this file

%This a script from Jaroslaw Tuszynski with small modifications
%see: http://www.mathworks.com/matlabcentral/fileexchange/28155-inscribed-rectangle


% Copyright (c) 2010, Jaroslaw Tuszynski
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.


 
if (nargin<2)
  crit = [1 1 0];
end
if (nargin<3)
  minSize = [1 1];
end
p = crit;
[nR nC] = size(I);
if (minSize(1)<1), minSize(1)= floor(minSize(1)*nR); end
if (minSize(2)<1), minSize(2)= floor(minSize(2)*nC); end
if (max(I(:)) - min(I(:))==1),
  S = FindLargestSquares(I);
else
  S = I;
end
n = max(S(:));
W = S; % make a carbon copy of the matrix data
H = S;
C = ((p(1)+p(2)) + p(3)*S) .* S; % p(1)*width + p(2)*height + p(3)*height*width for height=width=S;
d = round((3*n)/4);
minH = max(min(minSize(1), d),1);
minW = max(min(minSize(2), d),1);

%% look for rectangles with width>height
hight2width = zeros(n+1,1);  % Store array with largest widths aviable for a given height
for r = 1 : nR               % each row is processed independently
  hight2width(:) = 0;        % reset the List
  for c = nC: -1 : 1         % go through all pixels in a row right to left
    s = S(r,c);              % s is a size of a square with its corner at (r,c)
    if (s>0)                 % if pixel I(r,c) is true
      MaxCrit = C(r,c);      % initialize the Max Criteria using square
      for hight = s:-1:1     % go through all possible width&hight combinations. Start with more likely to be the best
        width = hight2width(hight); % look up width for a given hight
        width = max(width+1,s);
        hight2width(hight) = width;
        Crit = p(1)*hight + p(2)*width + p(3)*width*hight;
        if (Crit>MaxCrit),   % check if it produces larger Criteria
          MaxCrit = Crit;    % if it does than save the results
          W(r,c)  = width;
          H(r,c)  = hight;
        end % if Crit
      end % for hight
      C(r,c)  = MaxCrit;
    end % if s
    hight2width((s+1):end) = 0;    % hights>s will not be aviable for the next pixel
  end % for c
end
clear hight2width

%% look for rectangles with width<height
width2hight = zeros(n+1,1);  % Store array with largest widths aviable for a given height
for c = 1 : nC               % each column is processed independently
  width2hight(:) = 0;        % reset the List
  for r = nR: -1 : 1         % go through all pixels in a column bottom to top
    s = S(r,c);              % s is a size of a square with its corner at (r,c)
    if (s>0)                 % if pixel I(r,c) is true
      MaxCrit = C(r,c);      % initialize the Max Criteria using square
      for width = s:-1:1     % go through all possible width&hight combinations. Start with more likely to be the best
        hight = width2hight(width); % look up hight for a given width
        hight = max(hight+1,s);
        width2hight(width) = hight;
        Crit = p(1)*hight + p(2)*width + p(3)*width*hight;
        if (Crit>MaxCrit),   % check if it produces larger Criteria
          MaxCrit = Crit;    % if it does than save the results
          W(r,c)  = width;
          H(r,c)  = hight;
        end % if Crit
      end % for width
      C(r,c)  = MaxCrit;
    end % if s
    width2hight((s+1):end) = 0;    % hights>s will not be aviable for the next pixel
  end % for r
end

%% Create container mask
CC = C;
CC( H<minH | W<minW ) = 0; % first try to obey size restrictions
[~, pos] = max(CC(:));
if (isempty(pos)), [~, pos] = max(C(:)); end % but when it fails than drop them
[r c] = ind2sub(size(C), pos);
M = false(size(C));
M( r:(r+H(r,c)-1), c:(c+W(r,c)-1) ) = 1;
%%change by P.Seibold:
w=W(r,c)-1;%width
h=H(r,c)-1;%height


%%Find largest square
function S = FindLargestSquares(I)
%FindLargestSquares - finds largest sqare regions with all points set to 1.
%input:  I - B/W boolean matrix
%output: S - for each pixel I(r,c) return size of the largest all-white square with its upper -left corner at I(r,c)  
[nr nc] = size(I);
S = double(I>0);
for r=(nr-1):-1:1
  for c=(nc-1):-1:1
    if (S(r,c))
      a = S(r  ,c+1);
      b = S(r+1,c  );
      d = S(r+1,c+1);
      S(r,c) = min([a b d]) + 1;
    end
  end  
end  