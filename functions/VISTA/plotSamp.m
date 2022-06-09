function [] = plotSamp(samp, param, ind)
% Author: Rizwan Ahmad (ahmad.46@osu.edu)

p = param.PE; % PE Index
t = param.FR; % Number of Frames
n = param.n; % Number of samples in each frame
N = length(ind); % Number of Samples

% Plot the sampling pattern marginals
pp = find(samp==1);
[x,~]=ind2sub([p,t],pp);
x = diff(x(:));
x(x<0)=[];

figure; 
        subplot(211); plot((sum(samp,2))); xlabel('k_y'); ylabel('Number of samples');
        subplot(212); plot((sum(samp,1))); xlabel('frames'); ylabel('Number of samples');
%         subplot(313); plot(1:numel(x),x,1:numel(x), sort(x),'r'); xlabel('samples'); ylabel('k-space jumps');



% 2D Binary map
figure; imagesc(squeeze(samp(:,:))); colormap(gray); axis('image');
title('VISTA Sampling', 'FontSize', 28); ylabel('k_y', 'FontSize', 28); xlabel('t', 'FontSize', 28);
ax = gca;
ax.FontSize = 18;

% Plot PE Index vs Acquisition Order for VISTA


PEind = zeros(1,n);

for i=1:size(samp, 2)
    if mod(i, 2) == 0
        % Sort in descending order (even frames)
        PEind((i-1)*n+1 : i*n) = sort(find(samp(:,i) == 1), 'descend');
    else
        PEind((i-1)*n+1 : i*n) = sort(find(samp(:,i) == 1), 'ascend');
    end
end

%acqOrderSamp

figure;
plot(1:N, PEind, '-');
hold on
plot(1:N, PEind, '.', 'MarkerEdgeColor', [0 0.4470 0.7410], 'MarkerFaceColor', [0 0.4470 0.7410], 'MarkerSize', 15);
title('PE Index vs Acquisition Order', 'FontSize', 20);
xlabel('Acquisition Order', 'FontSize', 20);
ylabel('PE Index', 'FontSize', 20);
ax = gca;
ax.FontSize = 18;
end

% 2D plot
% figure; 
% for i=1:param.t
%     k = (i-1)*param.tr + 1 : i*param.tr;
%     col = [rand,rand,rand];
%     col = col - 0.5*min(col);
%     col = (2/3 + 1/3*rand)*col/max(col);
%     plot(T(k),P(k),'o','markersize',param.sz,'color',[0 0 0],'markerfacecolor',col); axis('image'); hold on;
% end
% title([param.typ, ', Re = ' num2str(p*t/sum(samp(:)))]);
% axis([-floor(param.t/2), ceil(param.t/2)-1, -floor(param.p/2), floor(param.p/2)-1]);
