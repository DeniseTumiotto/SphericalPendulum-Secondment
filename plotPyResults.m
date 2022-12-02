function [] = plotPyResults(space, many, newFirst)
% plots the solutions of a gradient flow on S2 and TS2
% 
% :param space: which manifold (S2 or TS2)
% :param many: how many couple of solutions read
% :param newFirst: reading fisrt the newest solutions
%
% :returns: plots of interest
%

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
if strcmp(space,'S2')
    D = readpy(space,many,newFirst,'D');
elseif strcmp(space,'TS2')
    [D, ~] = readpy(space,many,newFirst,'D','M');
end

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

%     checking contractivity condition
    if strcmp(space,'S2')
        param.D = D{i};
        nParam = 100;
        contr = zeros(nParam,1);
        lam = zeros(nParam,2);
        tPar = linspace(0,1,nParam);
        j = 0;
        for tpar = tPar
            j = j+1;
            y0 = cart2sph(ProjS2(x{1}(:,1),x{2}(:,1),tpar));
            [contr(j,:), lam(j,:)] = contrCond(space,y0(2:3),param);
        end
        if contr
            disp('Logarithmic norm less or equal to 0!')
        else
            disp('Logarithmic norm positive in at least one point!')
        end
    end
    plotEig(tPar,lam)
end

    function rslt = ProjS2(p,q,t)
        rslt = ((1-t)*p+t*q)/norm((1-t)*p+t*q);
    end

    function [] = plotEig(space, l)
        ll = max(l,[],2);
        figure()
        plot(space, ll,'b','LineWidth',3)
%         hold on
%         plot(space, 0*space,'r','LineWidth',3)
%         plot(space, l(:,2),'r','LineWidth',3)
%         title('Eigenvalues of the symmetric part of the Jacobian', 'FontSize',20)
        title('Max eigenvalue of the symmetric part of the Jacobian', 'FontSize',20)
        xlabel('Parametrer','FontSize',18)
        ylabel('Eigenvalue','FontSize',18)
        grid on
        hold off
    end
end
