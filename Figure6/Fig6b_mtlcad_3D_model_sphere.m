figure;         clf;            base_dir = cd;
%%
% cd('mtlcad_3D_spheres-FF11-final');     data = readtable('mtlcad_3D_spheres_FF11_final.csv');
% cd('mtlcad_3D_spheres-FF37-final');     data = readtable('mtlcad_3D_spheres_FF37_final.csv');
cd('mtlcad_3D_spheres-FF50-final');       data = readtable('mtlcad_3D_spheres_FF50_final.csv');
positions = data{:,1:3};
radii = data{:,4};
num_spheres = length(radii);

hold on;                                            % Allow multiple spheres in the same plot
for i = 1:num_spheres
    % Generate sphere coordinates
    [x, y, z] = sphere(100);    
    x = x * radii(i) + positions(i, 1);
    y = y * radii(i) + positions(i, 2);
    z = z * radii(i) + positions(i, 3);
    surface = surf(x, y, z);
    % Choose a specific range for grey color (closer to white, e.g., 0.7 to 0.9)
    random_gray = 0.7 + 0.1 * rand();  % Random gray color in range [0.7, 0.9]
    % Set the color of the surface
    surface.FaceColor = [random_gray, random_gray, random_gray];  % Set grayscale color (RGB)
    surface.EdgeColor = 'none';  % Remove mesh edges for smooth appearance
%     shading flat;              % Smooth interpolation of colors
end

colormap(gray);                                        % Grayscale color map
light('Position', [-1, -1, 10], 'Style', 'infinite');  % Add a light source
lighting gouraud;               % Use Gouraud shading for smooth lighting
material([0.8, 0.6, 0.1]);      % Adjust material properties [ambient, diffuse, specular]

% % Set figure and axes background to transparent
% set(gcf, 'Color', 'none');  % Transparent figure background
% set(gca, 'Color', 'none');  % Transparent axes background
axis off;       % Remove axis, grid, and labels
axis equal;     % Equal aspect ratio
view(3);        % 3D view
rotate3d on;    % Enable rotation of the plot
cd(base_dir);
% exportgraphics(gcf,'Fig1c_mtlcad_3D_model_sphere.pdf','BackgroundColor','none');
exportgraphics(gcf, 'Fig1c_mtlcad_3D_model_sphere.png', 'Resolution', 400);
