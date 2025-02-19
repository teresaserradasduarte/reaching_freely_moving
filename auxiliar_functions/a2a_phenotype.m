% IDENTITY OF A2A AND WT
function phenotype_group = a2a_phenotype(mouse)

idx_a2_temp = [32, 33, 34, 35, 36];
idx_wt_temp = [31, 37, 38, 39, 40];

idthis_mouse = str2num(mouse(7:8));

if ismember(idthis_mouse,idx_wt_temp)
    phenotype_group = 'ctr';
elseif ismember(idthis_mouse,idx_a2_temp)
    phenotype_group = 'a2a';
else
    phenotype_group = 'unknown';
end


end

% Str31 GFP
% Str32 CASP
% Str33 CASP
% Str34 CASP
% Str35 CASP
% Str36 CASP
% Str37 GFP
% Str38 GFP
% Str39 GFP
% Str40 GFP
