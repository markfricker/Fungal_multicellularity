res1 = length(Mu_fungal);
MG = ones(res1, res1, 3);
CYM = ones(res1, res1, 3);

% if the fastest growing colony has grown by a factor of
% of n, and a slower growing colony has grown by a factor 
% of m, the channels for the colonies are coloured according 
% with intensity 1 and m/n respectively. Higher n increases the
% apparent diffence in colour, making a lighter CYM merge.
n = 1000;

for i = 1:res1
    for j = 1:res1
        
        max_g =  max([Mu_motile(i,j), Mu_fungal(i,j), Mu_autolytic(i,j)]);
       
        max_f_or_i = max([Mu_immobile(i,j), Mu_fungal(i,j)]);

        if max_g > 0
            
            
            CYM(i,j,2) = ...
                1 - (1/n)*exp(Mu_fungal(i,j)*log(n)/max_g);
            
            CYM(i,j,3) = ...
                1 - (1/n)*exp(Mu_autolytic(i,j)*log(n)/max_g);
            
            CYM(i,j,1) = ...
                1 - (1/n)*exp(Mu_motile(i,j)*log(n)/max_g);
        end
        
        if max_f_or_i > 0
            
            MG(i,j,1) = ...
                1 - (1/n)*exp(Mu_immobile(i,j)*log(n)/max_f_or_i);
            
            MG(i,j,3) = MG(i,j,1);
            
            MG(i,j,2) = ...
                1 - (1/n)*exp(Mu_fungal(i,j)*log(n)/max_f_or_i);
        end
    end
end

% figure(2)
% imagesc(log(Recalcitrance)/log(10), log(C_to_N_vector)/log(10), MG)
% set(gca,'YDir','normal')
% xlabel('log_{10}(Recalcitrance), hours')
% ylabel('log_{10}(C:N ratio)')
% title('Fungi vs Persistent', 'Fontsize', 14)
% axis square
% axis off

figure
imagesc(log(Recalcitrance)/log(10), log(C_to_N_vector)/log(10), CYM)
set(gca,'YDir','normal')
xlabel('log_{10}(Recalcitrance), hours')
ylabel('log_{10}(C:N ratio)')
title('Fungi vs Mobile vs Degradable', 'Fontsize', 14)
axis square
axis off
