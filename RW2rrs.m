function [U, rrs, nm] = RW2rrs(RW, nm, L)

% MW algorithm and for SPM and associated uncertainty estimates
%
% Juliana Tavora, University of Maine, 2020
%
% See the following publication for details on the method:
% Tavora, J, et al., An algorithm to estimate Suspended Particulate Matter
% concentrations and associated uncertainties from Remote Sensing Reflectance
% in Coastal Environments
%
% INPUTS:
%
% nm       -  wavelengths associated with measured Remote sensing reflectance
% RW       -  measured water-leaving reflectance
% L        -  constants to calculate U by Gordon et al (1988)
%
% OUTPUTS:
%
% U        -  function of IOPs by Gordon et al (1988) at 630 nm and longer
% rrs      -  measured below water rrs(nm) at 630 nm and and longer
% nm       -  wavelengths from 630 nm and longer
%
%-------------------------------------------------------------------------%

% water-leaving reflectance to above water remote sensing reflectance
Rrs = RW ./ pi;

% Lee et al., 2002 approach to transfer above water remote sensing
% reflectance to below water remote sensing reflectance
rrs = Rrs ./ (0.52+ (1.7 .* Rrs));

% Gordon et al (1988) approach
raiz_delta = sqrt((L(1)).^2 - (4*L(2).*(-rrs)));

if imag(raiz_delta) ~= 0
    raiz_delta = NaN;
end

U = (-L(1) + raiz_delta)./(2*L(2));

rrs          = rrs(:,(nm>=630));
U            = U(:,(nm>=630));
nm           = nm(nm>=630);

end

