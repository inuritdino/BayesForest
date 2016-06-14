function scat_out = extract_scatter(scat_in,scat_names,order,trunk)
% USAGE: SCAT_OUT = EXTRACT_SCATTER(SCAT_IN,SCAT_NAMES,ORDER)
%
% SCAT_IN: the input struct of the scatters obtained from GEN_SCATTER().
%
% SCAT_NAMES: the names of scatters to extract from SCAT_IN. Available
% names are: 'taper' (tapering function), 'bra' (branching angles),'curv' 
% (3D curvature in space), 'lchi_lapar' (length of children branches as a 
% function of length along the parent), 'lchi_bra_lapar' (the same as 
% previous but augmented with the branching angles). See READ_SCATTER_DAT()
%
% ORDER: array of topological order values of branches which are used in the
% scatter extraction (Gravelius order with trunk order = 0).
%
% SCAT_OUT: cell array of the output scatter in a specified order:
% according to the SCAT_NAMES, if numel(ORDER) > 1, then in the following
% order: ORDER(1)SCAT_NAMES(1), ORDER(1)SCAT_NAMES(2) ...
% ORDER(end)SCAT_NAMES(1), ORDER(end)SCAT_NAMES(2) ...
% ORDER(end)SCAT_NAMES(end)

if(nargin < 4)
    trunk = [];
end
if(~isa(scat_names,'cell'))
    scat_names = {scat_names};
end
scat_out = cell(1,length(order)*length(scat_names));
for ii = 1:length(order)
    for jj = 1:length(scat_names)
        if(strcmp(scat_names{jj},'taper'))
            scat_out{(ii-1)*length(scat_names)+jj} = (scat_in(order(ii)).tap)';
        elseif(strcmp(scat_names{jj},'bra'))
            scat_out{(ii-1)*length(scat_names)+jj} = (scat_in(order(ii)).bra)';
        elseif(strcmp(scat_names{jj},'curv'))
            scat_out{(ii-1)*length(scat_names)+jj} = (scat_in(order(ii)).curv)';
        elseif(strcmp(scat_names{jj},'lchi_lapar'))
            scat_out{(ii-1)*length(scat_names)+jj} = (scat_in(order(ii)).lchi_lapar)';
        elseif(strcmp(scat_names{jj},'lchi_bra_lapar'))
            scat_out{(ii-1)*length(scat_names)+jj} = (scat_in(order(ii)).lchi_bra_lapar)';
        elseif(strcmp(scat_names{jj},'ltot_rini'))
            scat_out{(ii-1)*length(scat_names)+jj} = (scat_in(order(ii)).ltot_rini)';
        elseif(strcmp(scat_names{jj},'az'))
            scat_out{(ii-1)*length(scat_names)+jj} = (scat_in(order(ii)).az)';
        end
    end
end
if(~isempty(trunk))
    m = length(order)*length(scat_names);
    scat_out{m+1} = trunk.scatter.tap';
    scat_out{m+2} = trunk.scatter.curv';
end

end