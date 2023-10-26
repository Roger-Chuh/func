%ROBUST_GM

% Copyright Snap Inc. 2020
% This sample code is made available by Snap Inc. for informational
% purposes only.  It is provided as-is, without warranty of any kind,
% express or implied, including any warranties of merchantability, fitness
% for a particular purpose, or non-infringement.  In no event will Snap
% Inc. be liable for any damages arising from the sample code or your use
% thereof.

function [s, W] = robust_gm(s, width)
tau_sq = width * width;
a = tau_sq ./ (s + tau_sq);
s = s .* a;
W = a .* a;
end