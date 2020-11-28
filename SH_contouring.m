function SH_contouring(order)

% Developer: Alvin Pei-Yan, Li
% E-mail: d05548014@ntu.edu.tw / a4624393@gmail.com

fprintf('This is for the %dth SH coils\n', order);

R = 14;
N = 1107;
hsize = 0.2;

N = N * (1e6 / 9998) * R ^ 2;
R = R * 10;
a = 4 * pi * R ^ 2 / N;
d = sqrt(a);
M_th = round(pi / d);
d_th = pi / M_th;
d_phi = a / d_th;

Loci_Matrix = zeros(1, 3); % [x1, y1, z1; x2, y2, z2; ....]
Sph_matrix = zeros(1, 3); % [phi, theta, b_value]

% x = R * sin(theta) * cos(phi)
% y = R * sin(theta) * sin(phi)
% z = R * cos(theta)

N_count = 0;
for m = 0 : round(M_th - 1)
    theta = pi * (m + 0.5) / M_th;
    if theta <= pi * (1 - hsize)
        M_phi = round(2 * pi * sin(theta) / d_phi);

        for n = 0 : M_phi - 1
            phi = 2 * pi * n / M_phi;

            x = R * sin(theta) * cos(phi);
            y = R * sin(theta) * sin(phi);
            z = R * cos(theta);

            Loci_Matrix(N_count + 1, :) = [x, y, z];
            Sph_matrix(N_count+1,:) = [phi, theta, 0];

            N_count=N_count+1;
        end
    end
end

filename = 'CurrentPatternOnLoops.mat';
load(filename,'current_patterns');

fprintf('contouring\n');

% mode 1
X = Sph_matrix(:, 1);
Y = Sph_matrix(:, 2);

for n = 1 : (order + 1) ^ 2 - 1
%for n = order
    fprintf('This is the %d SH coil\n', n);
    V = current_patterns(:, n);
    Xnew = X - pi * 0.5;
    Xnew(X > pi * 1.5) = Xnew(X > pi * 1.5) - pi * 2;

    g = 1;
    for k = 1 : size(X, 1)
         if (X(k) == 0)
            indicator(g) = k;
            g = g + 1;
         end
    end

    Vnew = V;
 %   indicator(length(indicator) + 1) = length(X) + 1;
 %   for k = 1 : length(indicator) - 1
 %       start_ = indicator(k) + 1;
 %       stop_ = indicator(k + 1) - 1;
 %       seg_X = Xnew(start_ : stop_);
 %       for m = start_ : stop_
 %           angle = Xnew(m);
 %           if nnz(seg_X == -angle)
 %              seg_V = V(start_ : stop_);
 %              Vnew(m) = 0.5 * (seg_V(seg_X == angle) + seg_V(seg_X == -angle));
 %           end
 %       end
 %   end
    
    V = Vnew / max([abs(max(Vnew)), abs(min(Vnew))]);

    fprintf('Generate the heavy mesh tasks\n');

    X2 = [X; X + 2 * pi]; Y2 = [Y; Y]; V2 = [V; V];
    X3 = [X2; X2(X2 >= pi) - pi; X2(X2 < pi) + 3 * pi]; 
    Y3 = [Y2; -Y2(X2 >= pi); -Y2(X2 < pi)]; 
    V3 = [V2; V2(X2 >= pi); V2(X2 < pi)];
    [XI, YI] = meshgrid(-pi/2 : 0.005 : pi * 1.5, 0 : 0.005 : max(Y));
    [XI2, YI2] = meshgrid(0 : 0.005 : 4 * pi, -max(Y) : 0.005 : max(Y));
    G3 = griddata(X3, Y3, V3, XI2, YI2, 'v4');
    Gsymm = G3(1 : size(XI, 1), round(0.75 * size(XI, 2)) + 1 : round(1.75 * size(XI, 2)));
    figure;
    values = reshape(V, [numel(V), 1]);
    bin_width = 2*iqr(values(:))*numel(values)^(-1/3);
    bin_counts = round((max(values) - min(values)) / bin_width);
    [C, h] = contourf(XI, YI, Gsymm, bin_counts, 'LineWidth', 2); colormap('jet'); colorbar    
    axis equal
    xlabel('azimuth angle')
    ylabel('elevation angle')
    c=colorbar;
    c.FontSize=60;
    c.Label.String = 'current density (A/mm^2)';
    box off;
    ax=gca;
    ax.FontSize=60;
    ax.TickDir='out';
    caxis([-1, 1]);
    c.Ticks=[-1, 0, 1];
    switch n
        case 1
            m = -1; l = 1;
        case 2
            m = 0; l = 1;
        case 3
            m = 1; l = 1;
        case 4
            m = -2; l = 2;
        case 5
            m = -1; l = 2;
        case 6
            m = 0; l = 2;
        case 7
            m = 1; l = 2;
        case 8
            m = 2; l = 2;
        case 9
            m = -3; l = 3;
        case 10
            m = -2; l = 3;
        case 11
            m = -1; l = 3;
        case 12
            m = 0 ; l = 3;
        case 13
            m = 1; l = 3;
        case 14
            m = 2; l = 3;
        case 15
            m = 3; l = 3;
    end
    title(sprintf('(%d, %d)', m, l));
                                 
    figure;
    fprintf('Translate the coordinates system\n');

    XI_1d = reshape(XI, [numel(XI), 1]);
    YI_1d = reshape(YI, [numel(YI), 1]);
    Gsymm_1d = reshape(Gsymm, [numel(Gsymm), 1]);
    Loci_Matrix_sph = cat(2, XI_1d, YI_1d);
    Loci_Matrix_xyzv = zeros(length(XI_1d), 4);

    el = -Loci_Matrix_sph(:, 2) + max(Y);
    az = Loci_Matrix_sph(:, 1) + pi;
    Loci_Matrix_xyzv(:, 1) = R * sin(el) .* cos(az);
    Loci_Matrix_xyzv(:, 2) = R * sin(el) .* sin(az);
    Loci_Matrix_xyzv(:, 3) = R * cos(el);
    Loci_Matrix_xyzv(:, 4) = Gsymm_1d;

    [B, I] = sort(YI_1d);
    Loci_Matrix_(:, 1) = Loci_Matrix_xyzv(I, 1);
    Loci_Matrix_(:, 2) = Loci_Matrix_xyzv(I, 2);
    Loci_Matrix_(:, 3) = Loci_Matrix_xyzv(I, 3);
    sorted_G = Loci_Matrix_xyzv(I, 4);
    
    border_theta = max(B);
    DT = delaunayTriangulation(Loci_Matrix_);
    [K,q] = convexHull(DT);
    
    emc=trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),sorted_G);
    axis equal;axis on
    xlabel('X-axis (mm)');
    ylabel('Y-axis (mm)');
    zlabel('Z-axis (mm)');
    view(-44,17);
    colormap(jet)
    c=colorbar;
    c.FontSize=60;
    ax=gca;
    ax.FontSize=60;
    c.AxisLocation='in';
    emc.FaceColor='interp';
    emc.EdgeColor='interp';
    box off
    caxis([-1, 1]);
    c.Ticks=[-1, 0, 1]; % include data of 99.8% 
    c.Label.String = 'current density (A/mm^2)';
    title(sprintf('(%d, %d)', m, l));
    hold on

    k = 1;
    cnum = 1;
    while k < size(C, 2)
       kl = C(2, k);
       v = k + 1 : k + kl;
       el = -C(2, v) + max(Y);
       az = C(1, v) + pi;
       xv = R * sin(el) .* cos(az) * 1.005;
       yv = R * sin(el) .* sin(az) * 1.005;
       zv = R * cos(el) * 1.005;
       plot3(xv,yv,zv,'LineWidth',2, 'Color','k'), hold on;
       cnum = cnum + 1;
       k = k + kl + 1;
    end

    fprintf('Done!\n');
    fprintf('Saving the data\n');
    fprintf('at contour_points_%d.mat C\n', n);
    %eval(sprintf('save contour_points_%d.mat C bin_counts', n));
end

end
