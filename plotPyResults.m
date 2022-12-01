function [] = plotPyResults(space, many, newFirst)
close all
clearvars -except space many newFirst
clc

if nargin < 3
    newFirst = 0;
    if nargin < 2
        many = 'all';
        if nargin < 1
            space = 'S2';
        end
    end
end

[sols, dist] = readpy(space,many,newFirst);

for i = 1:max(size(sols))
    n = size(sols{i});
    m = size(dist{i});

    if n(2) ~= 2
        print('something is wrong here!')
    end

    x = cell(n(1)*n(2),1);

    for k = 1:n(1)
        x{2*k-1} = reshape(sols{i}(k,1,:,:),n(3),n(4));
        x{2*k  } = reshape(sols{i}(k,2,:,:),n(3),n(4));
        if strcmp(space,'TS2')
            dist{i}(k,:) = riemannianDistance(x(2*k-1:2*k));
        end
    end

    % print exact flow
    plotGradFlow(x{1}, x{2})
    % explicit
%     plotGradFlow(x{3}, x{4})
%     hold on
%     plot3(0,0,-1,'ok',MarkerSize=5,MarkerFaceColor='k')
    % implicit
%     plotGradFlow(x{5}, x{6})
%     hold on
%     plot3(0,0,-1,'ok',MarkerSize=5,MarkerFaceColor='k')

    % print distance
    color = linspecer(min(m)+1);
    h = linspace(0.1,15,max(m));
    figure()
    hold on
    for k = 2:min(m)
        plot(h,dist{i}(k,:),'LineWidth',2.5,'Color',color(k,:))
    end
    plot(h,dist{i}(1,:),'LineWidth',2.5,'Color',color(1,:))
    plot(h,ones(max(m),1)*dist{1}(1,1),'LineWidth',2.5,'Color',color(end,:))
    legend('explicit Lie-Euler','implicit Lie-Euler','exact flow','FontSize',15)
    grid on
    xlabel('Time step size', 'FontSize',18)
    ylabel('Riemannian distance', 'FontSize',18)
    hold off

% %     checking contractivity condition
%     param.D = [3.91079856 0 0; 0 3.91079856 0; 0 0 1.07126021];
%     param.D = [1.87722313 0 0; 0 1.87722313 0; 0 0 1.07020357];
%     x10 = x{1}(:,1);
%     x20 = x{2}(:,1);
%     j = 0;
%     for tpar = 0:0.001:1
%         j = j+1;
%         y0 = cart2sph(ProjS2(x10,x20,tpar));
%         [contr(j,:), lam(j,:)] = contrCond(space,y0(2:3),param);
%     end
%     
%     figure()
%     plot(lam(:,1))
%     hold on
%     plot(lam(:,2))
%     hold off
%     if contr
%         disp('contractive')
%     else
%         disp('not so much')
%     end
    
end

function rslt = ProjS2(p,q,t)
    rslt = ((1-t)*p+t*q)/norm((1-t)*p+t*q);
end

end
