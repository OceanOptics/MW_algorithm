function [Weight] = weight_asses(std, rrs, SPM_model, IOP_matrix, L, U, nm, Q, Q_filter)

% Wheight assessment function for SPM and associated uncertainty estimates
%
% Juliana Tavora, University of Maine, 2020
%
% See the following publication for details on the method:
% Tavora, J, et al., An algorithm to estimate Suspended Particulate Matter 
% concentrations and associated uncertaintiesfrom Remote Sensing Reflectance
%
% INPUTS:
%
% nm         -  wavelengths associated with measured Remote sensing reflectance
% std        -  standard deviation of mesaured below water rrs(nm) 
% rrs        -  mesaured below water rrs(nm)
% U          -  function of IOPs by Gordon et al (1988) 
% L          -  constants to calculate U by Gordon et al (1988)
% IOP_matrix -  matrix of all IOPs combinations
% SPM_model  -  all SPM estimates
% Q          -  saturation estimates 
% Q_filter   -  saturation threshold 
%
% OUTPUTS:
%
% Weight     - Maximum uncertainty taking into account absolute and relative uncertainties propagated to SPM 
%
%----------------------------------------------------------------------------------------------------------------------------------------------------%

SPM_median = (squeeze(nanmedian(SPM_model,2)));

for i=1:length(rrs(:,1))
    
%----------------------------------------------------------------------------------------------------------------------------------------------------%
% uncertainty
    
    relative = (sqrt(2).* 0.05 .* rrs(i,:));
    
    if isnan(std(i,:))
        absolute = repmat(nanstd(abs((smooth(rrs(i,:),10))' - rrs(i,:))),1,length(nm));
    else
        absolute = std(i,:);
    end
    
    max_noise(i,:) = (max([absolute;relative],[],1));
    
end

deltaU =  max_noise ./ (L(1) + 2.*L(2).* U);
iop =  repmat((squeeze(IOP_matrix(1,:,:))+squeeze(IOP_matrix(2,:,:))./squeeze(IOP_matrix(2,:,:))),1,1,length(rrs(:,1)));
iop(Q < 0 | Q > Q_filter) = NaN;
iop = squeeze(nanmedian((iop),2));

% final weight
Weight = 1./  ((deltaU' .* SPM_median) ./ ( U' - (U.^2)'  .* iop));

end

