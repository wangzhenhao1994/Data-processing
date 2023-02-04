% Using the exact colors from the SourceForge web site
x = linspace(-9,9,181);
g = exp(-x.^2/10);
c = cos(3*x);
figure(1)
area(x,g,'EdgeColor',[0 170 221]/255,'LineWidth',6,'FaceColor',[255 102 0]/255);
hold on
plot(x,g.*c*.986,'Color',[51 51 51]/255,'LineWidth',6)
hold off
axis([-9,9,-1.1, 1.1])
set(gca,'visible','off')

