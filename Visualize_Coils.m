function Visualize_Coils(order)

% Developer: Pei-Yan, Li
% E-mail: d05548014@ntu.edu.tw

for n = 1 : (order + 1) ^ 2 - 1
    fprintf('This is the %d SH coil\n', n);
    contours_matrix = sprintf('contour_points_%d.mat', n)
    load(contours_matrix, 'C'); % C is denoted as the matrix of contours' coordinates in Cardesian system.
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
    axis equal;axis on
    xlabel('X-axis (mm)');
    ylabel('Y-axis (mm)');
    zlabel('Z-axis (mm)');
    view(-44,17);
    colormap(jet)
    c=colorbar;
    c.FontSize=15;
    ax=gca;
    ax.FontSize=15;
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
end
end
