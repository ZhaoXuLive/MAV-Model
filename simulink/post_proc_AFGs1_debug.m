%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%
% function []=post_proc(savepath)
%post processing for CRRC Qingdao stage 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

close all
time = out.tdig_out;
ranbegin=1;%range for output trim 
ranend=1;%range for output trim 
dt=time(2)-time(1);

%%
figure('Name','vel_lat')
plot(time(ranbegin:end-ranend),vehicle_pva9(ranbegin:end-ranend,5)*3.6,'-','LineWidth', 2);hold on
plot(time(ranbegin:end-ranend),trailer1_pva9(ranbegin:end-ranend,5)*3.6,'-.','LineWidth', 2);
plot(time(ranbegin:end-ranend),trailer2_pva9(ranbegin:end-ranend,5)*3.6,'--','LineWidth', 2);
grid on
title('Lateral Velocity of units')
xlabel('Time [s]')
xlim([ranbegin/1000 inf])
ylabel('Velocity [kph]')
legend('1st unit','2nd unit','3rd unit','Location','Best');

ref_speed=squeeze(ref_speed_out);
figure('Name','vel_long')
plot(time(ranbegin:end-ranend),vehicle_pva9(ranbegin:end-ranend,4)*3.6,'-','LineWidth', 2);hold on
plot(time(ranbegin:end-ranend),trailer1_pva9(ranbegin:end-ranend,4)*3.6,'-.','LineWidth', 2);
plot(time(ranbegin:end-ranend),trailer2_pva9(ranbegin:end-ranend,4)*3.6,'--','LineWidth', 2);
plot(time(ranbegin:end-ranend),ref_speed(ranbegin:end-ranend)*3.6,'r','LineWidth', 1);
grid on
% grid minor
title('Longitudinal Velocity of Units')
xlabel('Time [s]')
xlim([ranbegin/1000 inf])
ylabel('Velocity [kph]')
legend('1st unit','2nd unit','3rd unit','vel_{ref}','Location','Best');

%% 
artc1=q_out*[0 0 1 -1 0]'*180/pi;
artc2=q_out*[0 0 0 1 -1]'*180/pi;
artc1_v=qd_out*[0 0 1 -1 0]'*180/pi;
artc2_v=qd_out*[0 0 0 1 -1]'*180/pi;

figure('Name','articulation');

subplot(2,1,1),
yyaxis left,plot(time(ranbegin:end-ranend),artc1(ranbegin:end-ranend),'b-','LineWidth', 2),ylabel('Angle[deg]');
hold on
yyaxis right,plot(time(ranbegin:end-ranend),artc1_v(ranbegin:end-ranend),'r-.','LineWidth', 2);ylabel('Velocity[dps]');
title('Articulation Angle and Velocity of Joint 2');
xlim([ranbegin/1000 inf])

grid on

subplot(2,1,2);
yyaxis left,plot(time(ranbegin:end-ranend),artc2(ranbegin:end-ranend),'b-','LineWidth', 2),ylabel('Angle [deg]');
hold on
yyaxis right,plot(time(ranbegin:end-ranend),artc2_v(ranbegin:end-ranend),'r-.','LineWidth', 2);ylabel('Velocity [dps]');
title('Articulation Angle and Velocity of Joint 3');
xlabel('Time [s]')
xlim([ranbegin/1000 inf])
grid on

%% 
%call offtrack of each axel and joint
% N=size(r12_out_sync,3);
% s12_out=zeros(12,2,N);
% 
% for i=1:12
%     for j=1:N
%         r=r12_out_sync(i,:,j);r=squeeze(r);
%         s12_out(i,:,j)=x2ss(r,TRACK,sxy_ref(i,:),20);
%     end
% end
% 
% 
% off_1st_JOINT=squeeze(s12_out(1,1:2,:));
% off_2nd_JOINT=squeeze(s12_out(5,1:2,:));
% off_3rd_JOINT=squeeze(s12_out(8,1:2,:));
% off_4th_JOINT=squeeze(s12_out(12,1:2,:));
% 
% off_fx_1st=squeeze(s12_out(2,1:2,:));
% off_mx_1st=squeeze(s12_out(3,1:2,:));
% off_rx_1st=squeeze(s12_out(4,1:2,:));
% off_fx_2nd=squeeze(s12_out(6,1:2,:));
% off_rx_2nd=squeeze(s12_out(7,1:2,:));
% off_fx_3rd=squeeze(s12_out(9,1:2,:));
% off_mx_3rd=squeeze(s12_out(10,1:2,:));
% off_rx_3rd=squeeze(s12_out(11,1:2,:));
% 
% dis_sequence=linspace(0,off_fx_1st(1,end),length(time));
% curvature_est=interp1(TRACK(:,1),TRACK(:,8),dis_sequence');

off_1st_JOINT=squeeze(offtrack_out(1:2,1,:));
off_2nd_JOINT=squeeze(offtrack_out(1:2,5,:));
off_3rd_JOINT=squeeze(offtrack_out(1:2,8,:));
off_4th_JOINT=squeeze(offtrack_out(1:2,12,:));

off_fx_1st=squeeze(offtrack_out(1:2,2,:));
off_mx_1st=squeeze(offtrack_out(1:2,3,:));
off_rx_1st=squeeze(offtrack_out(1:2,4,:));
off_fx_2nd=squeeze(offtrack_out(1:2,6,:));
off_rx_2nd=squeeze(offtrack_out(1:2,7,:));
off_fx_3rd=squeeze(offtrack_out(1:2,9,:));
off_mx_3rd=squeeze(offtrack_out(1:2,10,:));
off_rx_3rd=squeeze(offtrack_out(1:2,11,:));

figure('Name','offtrack')
ranbegin2=1;
plot(off_fx_1st(1,ranbegin2:end-ranend),off_fx_1st(2,ranbegin2:end-ranend),'-','LineWidth', 2);hold on
plot(off_mx_1st(1,ranbegin2:end-ranend),off_mx_1st(2,ranbegin2:end-ranend),'-','LineWidth', 2);hold on
plot(off_rx_1st(1,ranbegin2:end-ranend),off_rx_1st(2,ranbegin2:end-ranend),'-.','LineWidth', 2);
plot(off_fx_2nd(1,ranbegin2:end-ranend),off_fx_2nd(2,ranbegin2:end-ranend),':','LineWidth', 2);
plot(off_rx_2nd(1,ranbegin2:end-ranend),off_rx_2nd(2,ranbegin2:end-ranend),'--','LineWidth', 2);
plot(off_fx_3rd(1,ranbegin2:end-ranend),off_fx_3rd(2,ranbegin2:end-ranend),'-','LineWidth', 2);
plot(off_mx_3rd(1,ranbegin2:end-ranend),off_mx_3rd(2,ranbegin2:end-ranend),'-','LineWidth', 2);
plot(off_rx_3rd(1,ranbegin2:end-ranend),off_rx_3rd(2,ranbegin2:end-ranend),'-d','LineWidth',1,'MarkerSize',2);
ylabel('Offtrack [m]');

grid on
legend('a1',...
    'a2',...
    'a3',...
    'a4',...
    'a5',...
    'a6',...
    'a7',...
    'a8',...
    'Location',...
    'Best');

title('Offtrack of axles');
xlabel('Distance [m]')

%%
%side slip
figure('Name','sideslip')

plot(time(ranbegin:end-ranend),-atan2(vehicle_pva9(ranbegin:end-ranend,5),vehicle_pva9(ranbegin:end-ranend,4))*180/pi,'-','LineWidth', 2),hold on
plot(time(ranbegin:end-ranend),-atan2(trailer1_pva9(ranbegin:end-ranend,5),trailer1_pva9(ranbegin:end-ranend,4))*180/pi,'-.','LineWidth', 2)
plot(time(ranbegin:end-ranend),-atan2(trailer2_pva9(ranbegin:end-ranend,5),trailer2_pva9(ranbegin:end-ranend,4))*180/pi,'--','LineWidth', 2)

legend('1st unit','2nd unit','3rd unit','Location','Best')
xlabel('Time[s]')
xlim([ranbegin/1000 inf])

ylabel('Sideslip angle[deg]')
grid on
title('Sideslip Angle')
%%
%steering angle
figure('Name','steer_angle')
% time_seq=steer6_out.Time;
steer_angle=steer8_out;
sty={'r','g--','b-.','c-.','m:','k--','y--','b:'};

for i=1:8
    plot(time(ranbegin:end-ranend),steer_angle(ranbegin:end-ranend,i)*180/pi,sty{i},'LineWidth', 2);hold on
end
grid on
xlabel('Time [s]');
xlim([ranbegin/1000 inf])
ylabel('Steering angle [deg]');
title('Steering Angle of Axles')
legend('a1',...
    'a2',...
    'a3',...
    'a4',...
    'a5',...
    'a6',...
    'a7',...
    'a8',...
    'Location',...
    'Best');

%%
% %%
% figure('Name','brake')
% plot(time(ranbegin:end-ranend),brake6(ranbegin:end-ranend,:)','LineWidth', 2);
% grid on
% xlabel('Time [s]');
% xlim([ranbegin/1000 inf])
% 
% ylabel('Brake Torque [Nm]');
% title('Brake Torque of Axles')
% legend('front axle of 1st unit',...
%     'rear axle of 1st unit',...
%     'front axle of 2nd unit',...
%     'rear axle of 2nd unit',...
%     'front axle of 3rd unit',...
%     'rear axle of 3rd unit',...
%     'Location',...
%     'Best');


% %%
% figure('Name','drive');
% plot(time(ranbegin:end-ranend),drive6(ranbegin:end-ranend,1),'r-','LineWidth', 2);
% hold on
% plot(time(ranbegin:end-ranend),drive6(ranbegin:end-ranend,6),'b','LineWidth', 2);
% % plot(time,drive_torque_actual(:,1),'b-.');hold on
% % plot(time,drive_torque_actual(:,3),'m-.');
% title('Drive Torque');
% grid on
% % grid minor
% legend('1st','6th','Location','best');
% 
% % legend('front','rear','front(actual)','rear(actual)','Location','best');
% xlabel('Time [s]')
% xlim([ranbegin/1000 inf])
% 
% ylabel('Drive Torque [Nm]')
%%
figure('Name','yaw_rate');
plot(time(ranbegin:end-ranend),vehicle_att9_rpy(ranbegin:end-ranend,8)*180/pi,'r-','LineWidth', 2);
hold on
plot(time(ranbegin:end-ranend),trailer1_att9_rpy(ranbegin:end-ranend,8)*180/pi,'g-.','LineWidth', 2);
hold on
plot(time(ranbegin:end-ranend),trailer2_att9_rpy(ranbegin:end-ranend,8)*180/pi,'b--','LineWidth', 2);

% plot(time,drive_torque_actual(:,1),'b-.');hold on
% plot(time,drive_torque_actual(:,3),'m-.');
title('Yaw Rate');
grid on
% grid minor
legend('1st unit','2nd unit','3rd unit','Location','best');

% legend('front','rear','front(actual)','rear(actual)','Location','best');
xlabel('Time [s]')
xlim([ranbegin/1000 inf])

ylabel('Yaw Rate [deg/s]')

%%
figure('Name','yaw_angle');
plot(time(ranbegin:end-ranend),vehicle_att9_rpy(ranbegin:end-ranend,7)*180/pi,'r-','LineWidth', 2);
hold on
plot(time(ranbegin:end-ranend),trailer1_att9_rpy(ranbegin:end-ranend,7)*180/pi,'g-.','LineWidth', 2);
hold on
plot(time(ranbegin:end-ranend),trailer2_att9_rpy(ranbegin:end-ranend,7)*180/pi,'b--','LineWidth', 2);

% plot(time,drive_torque_actual(:,1),'b-.');hold on
% plot(time,drive_torque_actual(:,3),'m-.');
title('Yaw Angle');
grid on
% grid minor
legend('1st unit','2nd unit','3rd unit','Location','best');

% legend('front','rear','front(actual)','rear(actual)','Location','best');
xlabel('Time [s]')
xlim([ranbegin/1000 inf])

ylabel('Yaw Angle [deg]')

%%
% %%
% %acc command
% figure('Name','acc_command')
% % time_seq=steer6_out.Time;
% acc_command=acc_command_out.Data;
% % plot(time_seq,steer_angle,'LineWidth', 2);
% plot(time,acc_command,'b','LineWidth', 2);
% grid on
% xlabel('Time [s]');
% ylabel('Acceleration command');
% title('Acceleration command')
% legend('Acceleration command',...
%     'Location',...
%     'Best');

% % 
%%
figure('Name','angerf')
% subplot(2,1,1)
% yyaxis left;
plot(time(ranbegin:end-ranend),ANGref(ranbegin:end-ranend,:),'-','LineWidth', 2);
ylabel('ANGref [rad]');

grid on
legend('front axle of 1st unit',...
    'rear axle of 1st unit',...
    'front axle of 2nd unit',...
    'rear axle of 2nd unit',...
    'front axle of 3rd unit',...
    'rear axle of 3rd unit',...
    'Location',...
    'Best');

title('ANGref');
xlabel('time [s]')
xlim([ranbegin/1000 inf])

%%
figure('Name','fb_signal')
% subplot(2,1,1)
% yyaxis left;
plot(time(ranbegin:end-ranend),fb_error6(ranbegin:end-ranend,:),'-','LineWidth', 2);
ylabel('Feedback signal [rad]');

grid on
legend('front axle of 1st unit',...
    'rear axle of 1st unit',...
    'front axle of 2nd unit',...
    'rear axle of 2nd unit',...
    'front axle of 3rd unit',...
    'rear axle of 3rd unit',...
    'Location',...
    'Best');

title('Feedback Signal');
xlabel('time [s]')
xlim([ranbegin/1000 inf])

%%
figure('Name','ff_signal')
% subplot(2,1,1)
% yyaxis left;
plot(time(ranbegin:end-ranend),ff_error6(ranbegin:end-ranend,:),'-','LineWidth', 2);
ylabel('Feedforward signal [rad]');

grid on
legend('front axle of 1st unit',...
    'rear axle of 1st unit',...
    'front axle of 2nd unit',...
    'rear axle of 2nd unit',...
    'front axle of 3rd unit',...
    'rear axle of 3rd unit',...
    'Location',...
    'Best');

title('Feedforward Signal');
xlabel('time [s]')
xlim([ranbegin/1000 inf])

%% 
figure('Name','Trajectory')
for i = 1:12
plot(squeeze(r12_out(i,1,:)),squeeze(r12_out(i,2,:)),'-','LineWidth', 2);
hold on
end
ylabel('Y  [m]');

grid on
legend('Joint 1',...
    'Axle 1',...
    'Axle 2',...
    'Axle 3',...
    'Joint 2',...
    'Axle 4',...
    'Axle 5',...
    'Joint 3',...
    'Axle 6',...
    'Axle 7',...
    'Axle 8',...
    'Joint 4');

title('Trajectory of each key points');
xlabel('X [m]')
xlim([ranbegin/1000 inf])

%%
% figsav('figures');

%%
if 0
file_max = fopen([savepath '/MATRICS.txt'], 'w');
fprintf(file_max,'Max acc lat 1st:   ');
fprintf(file_max,'%4.2f ', max(cg1st_a_lat));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Max acc lon 1st:   ');
fprintf(file_max,'%4.2f ', max(cg1st_a_long));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Max acc lat 2nd:   ');
fprintf(file_max,'%4.2f ', max(cg2nd_a_lat));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Max acc lon 2nd:   ');
fprintf(file_max,'%4.2f ', max(cg2nd_a_long));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Max acc lat 3rd:   ');
fprintf(file_max,'%4.2f ', max(cg3rd_a_lat));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Max acc lon 3rd:   ');
fprintf(file_max,'%4.2f ', max(cg3rd_a_long));
fprintf(file_max,' \n');
fprintf(file_max,' \n');


fprintf(file_max,'Max acc tot 1st:   ');
fprintf(file_max,'%4.2f ', max(hypot(cg1st_a_lat,cg1st_a_long)));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Max acc tot 2nd:   ');
fprintf(file_max,'%4.2f ', max(hypot(cg2nd_a_lat,cg2nd_a_long)));
fprintf(file_max,' \n');
fprintf(file_max,' \n');


fprintf(file_max,'Max acc tot 3rd:   ');
fprintf(file_max,'%4.2f ', max(hypot(cg3rd_a_lat,cg3rd_a_long)));
fprintf(file_max,' \n');
fprintf(file_max,' \n');


fprintf(file_max,'Max vel tot:       ');
fprintf(file_max,'%4.2f ', max(velocity/kph));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Yaw rate 1st:      ');
fprintf(file_max,'%4.2f ', max(cg1st_a_lat*9.81)/max(velocity));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Yaw rate 2nd:      ');
fprintf(file_max,'%4.2f ', max(cg2nd_a_lat*9.81)/max(velocity));
fprintf(file_max,' \n');
fprintf(file_max,' \n');

fprintf(file_max,'Yaw rate 3rd:      ');
fprintf(file_max,'%4.2f ', max(cg3rd_a_lat*9.81)/max(velocity));
fprintf(file_max,' \n');
fprintf(file_max,' \n');


fclose(file_max);


end

%%
% save([savepath '/sim_data.mat']);

%close all
%clear all



