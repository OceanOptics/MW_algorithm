function [SPM_mw, err_mw, IOP_table] = MW_algorithm(nm, std, rrs, U, L, S, Y, a_nap443, bbp700, a_nap750, temp, Q_filter)

% MW algorithm and for SPM and associated uncertainty estimates
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
% S          -  exponente of exponential spectral shape for a_NAP :[0.006:0.001:0.014]
% Y          -  exponente of power-law spectral shape for bbp : [0:0.15:1.8]
% a_nap443   -  mass specific absorption at the reference at 443 nm : [0.01:0.01:0.06]
% a_nap750   -  mass specific absorption at the reference at 750 nm :[0.013:0.001:0.015]
% bbp700     -  mass specific backscattering at 700 nm : [0.002:0.001:0.021]
% temp       -  temperature of the water in Celsius which remote sensing
%               reflectance was measured
% Q_filter   -  saturation threshold 
%
% This function calls other 3 functions:  
%   1) 'asw_correction' Water absorption for temperature correction 
%   2) 'weight_assess' for estimates of maximum uncertainties and
%       estabilishing the weight used for final SPM and uncertainty estimates
%   3) 'IOP_sel' constrains combinations of IOPs and returns a table with
%   the 20 most common IOPs for each measurement
%
% OUTPUTS:
%
% SPM_mw    -  SPM estimates
% err_mw    -  uncertainties associated with the method in percentage          
% IOP_table -  table of the 20 most common combinations of IOPs. The
%              following code line can be added to estimate median and ranges
%              of IOP
%
% The SPM retrivals and uncertainties can also be applied to ocean color 
% remote sensing data such as Landsat8/OLI, MODIS/aqua and others. 
% Adjustments must be done however to the 'weight_assess' function 
% regarding atmospheric correction uncertanties (see section uncertainties,
% reference above).

%-------------------------------------------------------------------------%

% number of combination of IOPs
dim = (length(S)*length(Y)*length(a_nap443)*length(a_nap750)*length(bbp700)); 

%generate vectors for IOP's -> water absorption:

[absorp_cor] = asw_correction(temp,nm);
a_sw = repmat(absorp_cor',1,1,dim);
a_sw = permute(a_sw, [1 3 2]);
clear absorp_cor

%-------------------------------------------------------------------------%

%generate vectors for IOP's -> absorption:

for i = 1:length(S)
    for ii = 1:length(nm)
        for iii = 1:length(a_nap443)
            for iiii = 1:length(a_nap750)
                a_p(ii,iii,iiii,i) =  a_nap443(iii) .* (exp(-S(i)*(nm(ii)-443))  - exp(-S(i)*(750-443)) )  + a_nap750(iiii);
                S_slope(ii,iii,iiii,i) = S(i);
                a_nap443_matrix(ii,iii,iiii,i) = a_nap443(iii);
                a_nap750_matrix(ii,iii,iiii,i) =a_nap750(iiii);
            end
        end
    end
end

ap = reshape(a_p,[length(nm),(length(S)*length(a_nap443)*length(a_nap750))])';
S_res = reshape(S_slope,[length(nm),(length(S)*length(a_nap443)*length(a_nap750))])';
a_nap443_res = reshape(a_nap443_matrix,[length(nm),(length(S)*length(a_nap443)*length(a_nap750))])';
a_nap750_res = reshape(a_nap750_matrix,[length(nm),(length(S)*length(a_nap443)*length(a_nap750))])';
clear i ii iii a_p

%-------------------------------------------------------------------------%
%generate vectors for IOP's -> backscattering:

for n = 1:length(Y)
    for nn = 1:length(nm)
        for nnn = 1:length(bbp700)
            bb_p(nn,nnn,n) = bbp700(nnn).*((700./nm(nn)).^(Y(n)));
            Y_slope(nn,nnn,n) = Y(n);
            bbp_matrix(nn,nnn,n) = bbp700(nnn);
        end
    end
end

bbp = reshape(bb_p,[length(nm),(length(Y)*length(bbp700))])';
Y_res = reshape(Y_slope,[length(nm),(length(Y)*length(bbp700))])';
bbp_res = reshape(bbp_matrix,[length(nm),(length(Y)*length(bbp700))])';
clear bb_p n nn nnn

%-------------------------------------------------------------------------%
%generate all possible combination of eigenvectors

k = 0;
for i = 1:length(S)*length(a_nap443)*length(a_nap750)
    for n = 1:length(Y)*length(bbp700)
        k = k+1;
        IOP_matrix(:, :, k) = [ap(i, :); bbp(n, :)];
        shape_matrix(:, :, k) = [S_res(i, :); Y_res(n, :)];
        coefs_matrix(:, :, k) = [a_nap443_res(i, :); a_nap750_res(i, :); bbp_res(n, :)];
    end
end

%--------------------------------------------------------------------------%
%% SPM calculation

Q = NaN(length(nm),(dim),length(rrs(:,1)));


U_model = (squeeze(IOP_matrix(2,:,:)) ./ (squeeze(IOP_matrix(1,:,:)) + squeeze(IOP_matrix(2,:,:))));


ii=1;
for i=1:length(rrs(:,1))
    
    %----------------------------------------------------------------------------------------------------------------------------------------------------%
    % SPM estimates
    
    if length(temp) ==1
        SPM_model(:,:,ii) = (1./squeeze(IOP_matrix(2,:,:)) .* (a_sw ./  ( (((1-U(i,:))./U(i,:))' - (squeeze(IOP_matrix(1,:,:))./squeeze(IOP_matrix(2,:,:)))))));
    else
        SPM_model(:,:,ii) = (1./squeeze(IOP_matrix(2,:,:)) .* (squeeze(a_sw(:,:,i)) ./  ( (((1-U(i,:))./U(i,:))' - (squeeze(IOP_matrix(1,:,:))./squeeze(IOP_matrix(2,:,:)))))));
    end
    
    %----------------------------------------------------------------------------------------------------------------------------------------------------%
    % saturation checking
    Q(:,:,ii) = U(i,:)' ./ U_model;
    
    %----------------------------------------------------------------------------------------------------------------------------------------------------%
    fprintf('%0.0f/%0.0f\n',i,length(rrs(:,1)));
    
    ii=ii+1;
    
end

SPM_model(SPM_model < 0 | Q < 0 | Q > Q_filter) = NaN;

%---------------------------------------------------------------------------------------------------------------------------------------------------------%
%% Solutions

% ----------------------------------------------------------------------- %
% Wheighting estimates 

[Weight] = weight_asses(std, rrs, SPM_model, IOP_matrix, L, U, nm, Q, Q_filter);
 W = repmat(Weight,1,1,(dim));
 W = permute(W, [1 3 2]);

% filtering results for saturation 
W(SPM_model < 0 | Q < 0 | Q > Q_filter)=NaN;
Q(SPM_model < 0 | Q < 0 | Q > Q_filter)=NaN;
coefs_matrix(:, SPM_model < 0 | Q < 0 | Q > Q_filter) = NaN;
shape_matrix(:, SPM_model < 0 | Q < 0 | Q > Q_filter) = NaN;

% ----------------------------------------------------------------------- %
% SPM retrieval

SPM_mw = (squeeze(nansum(nansum(SPM_model.*W),2)) ./ squeeze(nansum(nansum(W,2))))';

SPM_84 = squeeze(prctile(squeeze(SPM_model),84,2));
SPM_16 = squeeze(prctile(squeeze(SPM_model),16,2));

valid_nm = sum(~isnan(SPM_84),1);

SPM_max = ((nansum(SPM_84.*Weight)) ./ (nansum(Weight))).*(1./sqrt(valid_nm));
SPM_min = ((nansum(SPM_16.*Weight)) ./ (nansum(Weight))).*(1./sqrt(valid_nm));

err_mw = (((SPM_max - SPM_min)./2)*100)./SPM_mw;

%% constraining IOPs

disp('Constraining IOPs')

for i= 1:(length(rrs(:,1)))
    
    [max_values_IOP] = IOPs_sel(squeeze(SPM_model(:,:,i)), SPM_mw(i), coefs_matrix, shape_matrix, nm, dim, 20);
    
    
    IOP_table(i,:,:) = max_values_IOP;
    
    %----------------------------------------------------------------------------------------------------------------------------------------------------%
    fprintf('%0.0f/%0.0f\n',i,length(rrs(:,1)));
end


end
