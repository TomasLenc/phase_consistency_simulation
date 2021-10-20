function [wvfOut] = myReson(wvfIn,resonatorFreq,fs,Q)

if nargin==3
    Q = 13; % number of periods for decay 
end

bw = resonatorFreq./Q;

R = 1 - bw./fs*pi;

psi = resonatorFreq./fs*(2*pi);    

gainFactor = 1; 

% approximating actual peak magnitude response from theta
% (Steiglitz, 1994)
%  'reson_z'
B_unscaled = [1 0 -1];
theta = acos((1+R^2)./(2*R).*cos(psi));
%     %  'reson'
%     B_unscaled = [1];
%     theta =  acos((2*R)./(1+R^2).*cos(psi));
%     %  'reson_r' %we have no peak response correction equation for reson_r
%     theta = psi;
%     B_unscaled = [1 0 -R];

realFreq = psi./(2*pi).*fs;

B = B_unscaled.*gainFactor;
A = [1 -2*R*cos(theta) R^2];

wvfOut = filter(B,A,wvfIn);

