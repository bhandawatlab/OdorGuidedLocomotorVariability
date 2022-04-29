% set up the axis and hunger-valence function
x = (0:0.01:1)';
y = @(x) tanh(x.*3);

%% plot the hunger-valence curve
figure;set(gcf,'Position',[2 42 487 924])
subplot(3,1,1);
plot(x,y(x),'k','linewidth',1)
xlabel('hunger');ylabel('valence');
title('Hunger-Valence curve (y = tanh(x.*3))')

%% plot hunger and valence distribution for higher hunger based on the curve
subplot(3,1,2);
pd = makedist('normal','mu',0.8,'sigma',0.05);
xx = (0:0.01:1)';
yy = pdf(pd,xx);
plot(xx,yy,'linewidth',1);hold on;

xx = normrnd(0.8,0.05,1,10000);
parmhat = lognfit(y(xx));
pd = makedist('Lognormal','mu',parmhat(1),'sigma',parmhat(2));
xx = (0:0.01:1)';
yy = pdf(pd,x);
plot(xx,yy,'linewidth',1)
xlabel('value');ylabel('pdf')
legend({'hunger','valence'})
title('Higher hunger')

%% plot hunger and valence distribution for lower hunger  based on the curve
subplot(3,1,3);
pd = makedist('normal','mu',0.3,'sigma',0.05);
xx = (0:0.01:1)';
yy = pdf(pd,x);
plot(xx,yy,'linewidth',1);hold on;

xx = normrnd(0.3,0.05,1,10000);
parmhat = lognfit(y(xx));
pd = makedist('Lognormal','mu',parmhat(1),'sigma',parmhat(2));
xx = (0:0.01:1)';
yy = pdf(pd,x);
plot(xx,yy,'linewidth',1)
xlabel('value');ylabel('pdf')
legend({'hunger','valence'})
title('Lower hunger')

%% printing figures
print('-painters','-dpdf',['HungerValence.pdf']);
