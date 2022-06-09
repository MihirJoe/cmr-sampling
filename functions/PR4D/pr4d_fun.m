function [PEInd] = pr4d_fun(param)

% Output
% ky(i,e) = PEInd(i,1,e); kz(i,e) = PEInd(i,2,e) where 'e' is encoding and
% 'i' is the order of acquisition


N   = param.PE; %Size of the square grid
FR  = param.FR; % Nominal number of frames
n   = param.n; % Nominal numbers of samples per frame per encoding
M   = param.M; % Number of samples per encoding
s   = param.s/2; % controls sampling density
ar  = param.ar; % as>=0; controls the aspect ratio/orientation of the high density ellipsoid region in the middle
E   = param.E; % Number of encodings 1<E<5
cg  = param.cg; % radial shift; cg=n implies central n pixels will not be sampled
gr  = 2 - (sqrt(5) + 1)/2; % golden ratio for radial increment
gs  = rem(param.gs, 1); % ratio for angular increment; other options: 35, 16, 13
dsp = param.dsp; % plot figures or not
p   = 0; % A larger value avoids a large number of samples at the center. The sampling density is not sensitity to this parameter


PEInd = zeros(M,2,E);
% avgSamp = zeros(N(1), N(2), E); % time average sampling
for e=1:E
%     kk=e;
    R  = zeros(M,1);
    %     RS = zeros(M,1);
    T  = zeros(M,1);
    x  = zeros(M,1);%rep(0,M);
    y  = zeros(M,1);%rep(0,M);
    Ri = floor((max(N)-1)/2);
    
    for j=1:M
        theta=2*pi*gs*(j+(e-1)/E);
        %     theta = theta+fa*e;%(e-1)/E;
        if j==1, r = R(j);
        else,    r  = rem((R(j-1) + Ri*gr)-0, Ri) + 0;
        end
        R(j) = r;
        theta_new = atan2((N(2)/N(1))^ar * sin(theta), cos(theta));
        x(j) = (r+p)^s   *  cos(theta_new);
        y(j) = (r+p)^s   *  sin(theta_new);
        T(j) = theta_new;
    end
    PEInd(:, 1, e) = round((N(1)/2-(cg+1))*x/max(abs(x)) + cg*cos(T)) + floor(N(1)/2)+1;
    PEInd(:, 2, e) = round((N(2)/2-(cg+1))*y/max(abs(y)) + cg*sin(T)) + floor(N(2)/2)+1;
end

% Plotting
if dsp == 1
%     E = 1; % comment out
    tiFont = 20; % title font
    axFont = 14; % axis font
    laFont = 18; % label font
    figure(1);
    col = 5;
    tiledlayout(E,col,'TileSpacing','compact', 'Padding', 'compact')
    for e=1:E % number of rows
        rang = [1, n]; % the number of samples plot in each frame
        for i=1:col % position in subplot space
            nexttile;
            if i < 4
                %Plot samples
                tSamp = zeros(N);
                for j=rang(1):rang(2)
                    tSamp(PEInd(j,1,e), PEInd(j,2,e)) = tSamp(PEInd(j,1,e), PEInd(j,2,e)) + 1;
                end
                imagesc(tSamp.^0.25); axis('image'); colormap(gray);
                set(gca, 'FontSize', axFont, 'FontName','times');
                if i == 1, ylabel(['$k_y$ (Encoding ' num2str(e),')'], 'FontSize', laFont,'Interpreter','latex');
                else i>1 && i <4; set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]); end
                if e==E, xlabel('$k_z$', 'FontSize', laFont,'Interpreter','latex');
                else, set(gca, 'xtick',[]); end
                if e==1, title([num2str(rang(1)) '-' num2str(rang(2)) ' Samples'], 'FontSize', tiFont,'Interpreter','latex'); end
                
                rang(1) = rang(2) + 1;
                rang(2) = rang(2) + n;
            
            elseif i == 4
                tSamp = zeros(N);
                for j=1:M
                    tSamp(PEInd(j,1,e), PEInd(j,2,e)) = tSamp(PEInd(j,1,e), PEInd(j,2,e)) + 1;
                end
                %Plot accumulation
                imagesc(tSamp.^0.25), axis('image'); colormap(gray);
                set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
                if e==E, xlabel('$k_z$', 'FontSize', laFont,'Interpreter','latex');
                else, set(gca,'xtick',[]); end
                if e==1, title([num2str(1) '-' num2str(M) ' Samples'], 'FontSize', tiFont,'Interpreter','latex'); end
              
            elseif i == 5
                %Plot binary mask
                imagesc(logical(tSamp)); axis('image'); colormap(gray);
                set(gca, 'FontSize', axFont, 'FontName','times','ytick',[]);
                if e==E, xlabel('$k_z$', 'FontSize', laFont,'Interpreter','latex');
                else, set(gca,'xtick',[]); end 
                if e==1, title(['Binary Mask'], 'FontSize', tiFont,'Interpreter','latex'); end
            end
        end
    end
    set(gcf,'color','w','units','points','position',[10,10,900,1070]); %export_fig('pr4d',gcf,'-m4','-png');
end