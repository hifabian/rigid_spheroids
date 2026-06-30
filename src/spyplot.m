[L2,Gxz,iLy,Gyz,iLx] = build_matrix(16, "store", false);


figure('Color','w','Position',[100 100 900 700]);
tiledlayout(1,5,'TileSpacing','compact','Padding','compact');

M = {L2,Gxz,iLy,Gyz,iLx};
titles = {'$\mathbf{L^2}$','$\mathbf{\Gamma_{xz}}$','$\mathbf{L_y}$', ...
    '$\mathbf{\Gamma_{yz}}$','$\mathbf{L_x}$'};

for i = 1:5
    ax = nexttile;
    spy(M{i},'k',8);
    title(titles{i},'Interpreter','latex','FontSize',24);
    %set(gca,'FontSize',11,'Box','on');
    ax.XTick = [];
    ax.YTick = [];
end

exportgraphics(gcf,'outputs/sparse_matrices.pdf','ContentType','vector');