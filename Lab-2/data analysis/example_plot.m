% plot
plot(t,fit1,'-g',t2,fit2,'-k')
xlabel('time (s)','FontSize',13,'FontWeight','bold');
ylabel('natural log [T-Tinf/Ti-Tinf]','FontSize',13,'FontWeight','bold');
title('semi-log profile force convection');
format short
txt1 = ['slope Run 1 ' num2str(p(1,1)) ' s  h = ' sprintf('%0.3g',h1_fc),'W/m^2K']; 
text(400, -0.65, txt1);
txt2 = ['slope Run 2 ' num2str(p2(1,1))  ' s  h = ' sprintf('%0.3g',h2_fc),'W/m^2K']; 
text(500, -1, txt2);

legend('Run 1','Run 2','fit 1','fit 2')

hold