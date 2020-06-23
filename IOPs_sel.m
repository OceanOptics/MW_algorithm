function [max_values_IOP] = IOPs_sel(SPM_model, SPM_mw, coefs_matrix, shape_matrix, nm, dim, Nmax)

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
% nm              -  wavelengths associated with measured Remote sensing reflectance
% dim             -  total number combination of IOPs
% SPM_model       -  all SPM estimates
% SPM_mw          -  final SPM estimate (SPM weighted median)
% coefs_matrix    -  matrix of mass-specific coefficients (ap443, ap750, bbp700)
% shape_matrix    -  matrix of shape parameters (S and Y)
% Nmax            -  number of most common entries
%
% OUTPUTS:
%
% max_values_IOP  -  most common 20 IOP combinations and percentage of occurrence
%
%-------------------------------------------------------------------------%

% Selects SPM solutions within 10%

mask_SPM = (squeeze(SPM_model)  >= (0.9.*SPM_mw) & squeeze(SPM_model) <= (1.1.*SPM_mw));
mask_SPM = double(mask_SPM);  mask_SPM(mask_SPM == 0) = NaN;
valid_nm = (nansum(mask_SPM,2)); valid_nm(valid_nm == 0) = NaN;
index = find(~isnan(valid_nm));
valid_nm = nm(index)';

bbp700_mask = squeeze(coefs_matrix(3, :, :)).*mask_SPM; bbp700_mask = bbp700_mask(index,:);
ap443_mask = squeeze(coefs_matrix(1, :, :)).*mask_SPM; ap443_mask = ap443_mask(index,:);
ap750_mask = squeeze(coefs_matrix(2, :, :)).*mask_SPM; ap750_mask = ap750_mask(index,:);
Y_mask = squeeze(shape_matrix(2, :, :)).*mask_SPM; Y_mask = Y_mask(index,:);
S_mask = squeeze(shape_matrix(1, :, :)).*mask_SPM; S_mask = S_mask(index,:);


working_IOPs = [];
for kk=1:dim
    
    if any(~isnan(bbp700_mask(:,kk)))
        sel_nm = find(~isnan(bbp700_mask(:,kk)));
        
        working_IOPs = [working_IOPs; valid_nm(sel_nm), bbp700_mask(sel_nm,kk), ap443_mask(sel_nm,kk), ap750_mask(sel_nm,kk), Y_mask(sel_nm,kk), S_mask(sel_nm,kk)];
        
    end
end

if isempty(working_IOPs)
    
    max_values = NaN(20,6);
    
else
    
    c = unique(working_IOPs(:,1)); % finds the unique IOP combinations
    for ii = 1:length(c)
        counts(ii,1) = sum(working_IOPs(:,1) ==c(ii)); % number of times each unique combination is repeated
    end
    
    [B, ~, iA] = unique(working_IOPs(:,2:end), 'rows');
    numoccurences = accumarray(iA, 1);
    indices = accumarray(iA, find(iA), [], @(rows){rows});  %the find(iA) simply generates (1:size(a,1))'
    
    [Avec, Ind ] = sort(numoccurences(:),1,'descend');
    max_values = Avec(1:Nmax);
    [ind_row, ~] = ind2sub(size(numoccurences),Ind(1:Nmax)); % fetch indices
    
    perc = (numoccurences(ind_row).*100) ./length(valid_nm);
    
    max_values_IOP = [perc, B(ind_row,:)];
    
end


end

