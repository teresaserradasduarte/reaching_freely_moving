
function [dist_trav, displacement, tortuosity] = tortuosity_calc_2D(xyz_vec,range)
% returns the power sspectrum of a vector
% input: 3D signal (across time)
%        range = [start stop], optional, default = 1:end
% output: dist_trav,displacement,tortuosity
% teresa, 23/03/2023 (last update: 5/8/2023)

% Default inputs
if nargin<2
    range = [1 size(xyz_vec,1)];
end

% Convert to mm - distance is not he same in all directions
% distance travelled
start = range(1);
stop = range(2);
diff_r = diff(xyz_vec(start:stop,:));
dist_points = squeeze(sqrt(...
    diff_r(:,1).^2 + ...
    diff_r(:,2).^2 ...
    ));
dist_trav = sum(dist_points,'omitnan');
% displacement
displacement = sqrt(...
    (xyz_vec(stop,1)-xyz_vec(start,1)).^2+...
    (xyz_vec(stop,2)-xyz_vec(start,2)).^2);
% tortuosity
tortuosity = dist_trav/displacement;

end