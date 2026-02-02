x=0:5:80;
% y is the monitoring statistics
y = [-1.830760,  -1.704123,  -1.602857,  -1.517437,  -1.442527,  -1.031871,  -1.148359,   3.070736,   6.793267,  14.367242, 29.501673,  44.450258,  71.990199, 124.102520, 162.030532];
plot(x(1:16),[1,y],'k');
hold on
h=9.20*ones(16,1);
light_blue = [0 0.447058823529412 0.741176470588235];
plot(x(1:16),h,'Color', light_blue);
xlim([0 75]);
ylim([-20 170]);
xlabel('Time t');
ylabel('Monitoring Statistics');