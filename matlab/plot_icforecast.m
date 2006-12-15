function plot_icforecast(Variable)
% stephane.adjemian@ens.fr
    
load conditional_forecasts;


eval(['ci1 = forecasts.cond.ci.' Variable ';'])
eval(['m1 = forecasts.cond.mean.' Variable ';'])
eval(['ci2 = forecasts.uncond.ci.' Variable ';'])
eval(['m2 = forecasts.uncond.mean.' Variable ';'])



H = length(m1);

% area(1:H,ci1(2,:),'FaceColor',[.9 .9 .9],'BaseValue',min([min(ci1(1,:)),min(ci2(1,:))]))

h1 = area(1:H,ci1(2,1:H))
set(h1,'BaseValue',min([min(ci1(1,:)),min(ci2(1,:))]))
set(h1,'FaceColor',[.9 .9 .9])

hold on
% area(1:H,ci1(1,:),'FaceColor',[1 1 1],'BaseValue',min([min(ci1(1,:)),min(ci2(1,:))]))
h2 = area(1:H,ci1(1,1:H));
set(h2,'BaseValue',min([min(ci1(1,:)),min(ci2(1,:))]))
set(h2,'FaceColor',[1 1 1])
plot(1:H,m1,'-k','linewidth',3)
plot(1:H,m2,'--k','linewidth',3)
plot(1:H,ci2(1,:),'--k','linewidth',1)
plot(1:H,ci2(2,:),'--k','linewidth',1)
axis tight
hold off