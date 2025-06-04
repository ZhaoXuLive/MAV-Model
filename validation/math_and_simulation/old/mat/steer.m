% time_step = time;
% STRa_out_step = STRa_out;

sty={'r:','r-','r--','b:','b-','b--'};
lines = ["r-"; ":"; "-."; "--"; "--."; ":."];
axle_str={'Axle 1','Axle 3','Axle 4','Axle 5','Axle 6','Axle 8'};
unit_str={'Unit 1 (math model)','Unit 2 (math model)','Unit 3 (math model)', 'Unit 1 (truckmaker)', 'Unit 2 (truckmaker)', 'Unit 3 (truckmaker)'};

% axle steer angles
subplot(1,2,1)
j = 1;
for i=1:6
    j = i;
    if j == 2 || j == 7
        j = j+1;
    end
    plot(time_step(1:2091),STRa_out_step(1:2091,j)*180/pi,lines(j),'LineWidth', 1);hold on
end
grid on
xlabel('Time [s]','FontName','Times New Roman','FontSize',12)
ylabel('Axle steer angle [deg]','FontName','Times New Roman','FontSize',12)
axis([0 20 -0.5 1.5])
legend(axle_str,'FontName','Times New Roman','FontSize',9);
title('Steering Angle of Step Response','FontName','Times New Roman','FontSize',12)

subplot(1,2,2)
j = 1;
for i=1:6
    j = i;
    if j == 2 || j == 7
        j = j+1;
    end
    plot(time(1:3800),STRa_out(1:3800,j)*180/pi,lines(j),'LineWidth', 1);hold on
end
grid on
xlabel('Time [s]','FontName','Times New Roman','FontSize',12)
ylabel('Axle steer angle [deg]','FontName','Times New Roman','FontSize',12)
axis([0 20 -0.3 0.4])
legend(axle_str,'FontName','Times New Roman','FontSize',9);
title('Steering Angle of Sinusoidal Response','FontName','Times New Roman','FontSize',12)