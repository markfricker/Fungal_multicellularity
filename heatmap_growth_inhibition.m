map = colorcet('L3');

immobile = ones(res1);
autolytic = ones(res1);
motile = ones(res1);
fungi = ones(res1);

% if the fastest growing colony has grown by a factor of
% of n, and a slower growing colony has grown by a factor 
% of m, the channels for the colonies are coloured according 
% with intensity 1 and m/n respectively. Higher n increases the
% apparent diffence in colour, making a lighter CYM merge.
n = 1000;

for i = 1:res1
    for j = 1:res1
        
        immobile(i,j) = Mu_immobile(i,j)/Mu_cell(i,j);
        
        autolytic(i,j) = Mu_autolytic(i,j)/Mu_cell(i,j);
        
        motile(i,j) = Mu_motile(i,j)/Mu_cell(i,j);
        
        fungi(i,j) = Mu_fungal(i,j)/Mu_cell(i,j);

    end
end

immobile_im = applycolourmap(immobile, map, [0, 1]);
autolytic_im = applycolourmap(autolytic, map, [0, 1]);
motile_im = applycolourmap(motile, map, [0, 1]);
fungi_im = applycolourmap(fungi, map, [0, 1]);

figure(1)
imagesc(log(Recalcitrance)/log(10), log(C_to_N_vector)/log(10), immobile_im)
set(gca,'YDir','normal')
xlabel('log_{10}(Recalcitrance), hours')
ylabel('log_{10}(C:N ratio)')
title('Immobile Growth Rate', 'Fontsize', 14)
axis square
axis off

figure(2)
imagesc(log(Recalcitrance)/log(10), log(C_to_N_vector)/log(10), autolytic_im)
set(gca,'YDir','normal')
xlabel('log_{10}(Recalcitrance), hours')
ylabel('log_{10}(C:N ratio)')
title('Autolytic Growth Rate', 'Fontsize', 14)
axis square
axis off

figure(3)
imagesc(log(Recalcitrance)/log(10), log(C_to_N_vector)/log(10), motile_im)
set(gca,'YDir','normal')
xlabel('log_{10}(Recalcitrance), hours')
ylabel('log_{10}(C:N ratio)')
title('Motile Growth Rate', 'Fontsize', 14)
axis square
axis off

figure(4)
imagesc(log(Recalcitrance)/log(10), log(C_to_N_vector)/log(10), fungi_im)
set(gca,'YDir','normal')
xlabel('log_{10}(Recalcitrance), hours')
ylabel('log_{10}(C:N ratio)')
title('Fungal Growth Rate', 'Fontsize', 14)
axis square
axis off