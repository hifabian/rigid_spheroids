[L2, G, iLy, W] = build_matrix(16, "store", false);


figure('Color','w','Position',[100 100 900 700]);
tiledlayout(1,4,'TileSpacing','compact','Padding','compact');

M = {L2,G,iLy,W};
titles = {'$\mathbf{L^2}$','$\mathbf{\Gamma}$','$\mathbf{L_y}$','$\mathbf{\Omega_p}$'};

for i = 1:4
    ax = nexttile;
    spy(M{i},'k',8);
    title(titles{i},'Interpreter','latex','FontSize',24);
    %set(gca,'FontSize',11,'Box','on');
    ax.XTick = [];
    ax.YTick = [];
end

exportgraphics(gcf,'outputs/sparse_matrices.pdf','ContentType','vector');