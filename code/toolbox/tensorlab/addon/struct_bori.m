function [x,state] = struct_bori(z, task)
%STRUCT_BORI short discription
%   long description
%   
    
% Author(s): Nico Vervliet       (Nico.Vervliet@esat.kuleuven.be)
%
% Version History:
% - 2015/11/10   NV      Initial version
    
    state = [];
    
    if isempty(task)
        x = nan(size(z{2},1), size(z{1},1));
        for r = 1:size(z{1},1)
            x(:,r) = polyval(z{1}(r,:), z{2}(:,r));
        end
    elseif ~isempty(task.r)
        x = nan(size(z{2},1), size(z{1},1));
        for r = 1:size(z{1},1)
            x(:,r) = polyval(task.r{1}(r,:), z{2}(:,r)) + ...
                     polyval(polyder(z{1}(r,:)),z{2}(:,r)).*task.r{2}(:,r);
        end
    elseif ~isempty(task.l)
        c = nan(size(z{1}));
        t = nan(size(z{2}));
        for r = 1:size(z{1},1)
            c(r,:) = task.l(:,r).'*conj(bsxfun(@power, z{2}(:,r), size(z{1},2)-1:-1:0));
            t(:,r) = conj(polyval(polyder(z{1}(r,:)), z{2}(:,r))).*task.l(:,r);
        end
        x = {c,t};
    end
    
end
