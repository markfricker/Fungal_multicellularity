mu1 = Mu_fungal;
mu2 = Mu_autolytic;

% max_ratio is the ratio between growth rate results in a colour that is
% indistinguishable from the colour given when the ratio between growth
% rates is infinite 
max_ratio = 2;

colour_growth_ratio = zeros(res1, res1, 3);
 
map = colorcet('R3');

growth_ratio = zeros(res1);

for i = 1:res1
    for j = 1:res1
        
        log_ratio = log(mu1(i,j)) - log(mu2(i,j));
        
        if log_ratio > 0
        
            d = 2 + 2*(log_ratio/log(max_ratio))^0.8;
            growth_ratio(i,j) = d;
        else
            d = 2 - 2*(-log_ratio/log(max_ratio))^0.8;
            growth_ratio(i,j) = d;
        end

        if (mu1(i,j) == 0) && ...
           (mu2(i,j) == 0)
            
            c1 = 1;
            c2 = 1;
            c3 = 1;
            
        elseif d >= 4
            
            c1 = 1;
            c2 = 0;
            c3 = 0;
            
        elseif d > 3
            
            c1 = 1;
            c2 = 4 - d;
            c3 = 0;
            
        elseif d > 2
            
            c1 = d - 2;
            c2 = 1;
            c3 = 0;
            
        elseif d > 1
            
            c1 = 0;
            c2 = 1;
            c3 = 2 - d;
            
        elseif d > 0
            
            c1 = 0;
            c2 = d;
            c3 = 1;
            
        else
            c1 = 0;
            c2 = 0;
            c3 = 1;
        end
        
        colour_growth_ratio(i,j,:) = [c1,c2,c3];
        
    end
end
 
rgbim = applycolourmap(growth_ratio, map, [0, 4]);

for i = 1:res1
    for j = 1:res1
        
        if (rgbim(i,j,1) == 0) && ...
           (rgbim(i,j,2) == 0) && ...
           (rgbim(i,j,3) == 0) 
       
            rgbim(i,j,:) = 1;
        end
    end
end

figure
imagesc(log(Recalcitrance), log(C_to_N_vector)/log(10), colour_growth_ratio)
set(gca,'YDir','normal')
xlabel('log_{10}(Recalcitrance), hours')
ylabel('log_{10}(C:N ratio)')
title('Log ratio growth rates', 'Fontsize', 14)
axis square
% axis off

figure
imagesc(log(Recalcitrance), log(C_to_N_vector)/log(10), rgbim)
set(gca,'YDir','normal')
xlabel('log_{10}(Recalcitrance), hours')
ylabel('log_{10}(C:N ratio)')
title('Log ratio growth rates', 'Fontsize', 14)
axis square
% axis off
