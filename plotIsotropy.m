function plotIsotropy(s, D)

n = size(s,1);
m = size(D{1},2);
h = linspace(0.1,10,m);

color = linspecer(n+2);

figure()
plot(h(1:floor(m/2)),ones(floor(m/2),1)*D{1}(1,1),'LineWidth',2.5,'Color',color(n+1,:))
hold on
for i = 1:n
    plot(h(1:floor(m/2)),D{i}(3,1:floor(m/2)),'LineWidth',2.5,'Color',color(i,:))
end
legend('initial distance','c=0','c=0.1','c=0.5','c=1','FontSize',15)
grid on
xlabel('Time step size', 'FontSize',18)
ylabel('Riemannian distance', 'FontSize',18)
title('Isotropy study -- Implicit Lie-Euler','FontSize',24)
hold off


end