x = linspace(-9,9,181);
g = exp(-x.^2/10);
c = cos(3*x);
figure(1)
plot(x,g,'LineWidth',7);
hold on
plot(x,g.*c*.986,'LineWidth',5)
hold off
axis([-9,9,-1.1, 1.1])
set(gca,'visible','off')

