% STability Evaluation Tool (STET.m)
% Calls articABCD2 to build state space (linear) model of ART
% (includes joint damping parameters)
% Reads data from excel file and restructures it for analysis
% updte to include cases with different number of units

clear
close all
g=9.81;deg=pi/180;kph=1/3.6;tiny=1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ DATA FROM FILE
infile='STET input data2.xlsx';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A0=readtable(infile, 'format', 'auto', 'VariableNamingRule', 'preserve');
AA=table2cell(A0);
%find locations of UNIT etc
temp = find(strcmp( AA, 'CASE' ));[rc,cc] = ind2sub(size(AA), temp); 
temp = find(strcmp( AA, 'SPEED' ));[rs,cs] = ind2sub(size(AA), temp); 
temp = find(strcmp( AA, 'UNIT' ));[ru,cu] = ind2sub(size(AA), temp); 
temp = find(strcmp( AA, 'AXLE' ));[ra,ca] = ind2sub(size(AA), temp); 
Ncases=length(ru); %how many boxes in the spreadsheet

%%%%%%% START LOOP OVER CASES IN THE SPREADSHEET %%%%%%%%%%%%%%%%%%%%%%%
Nu=zeros(Ncases,1);
% Na=zeros()
for k=1:Ncases
    kase=AA(rc(k),cc(k)+1);
    PARSu=AA(ru(k)+1:ru(k)+5,cu(k)+1:cu(k)+5);
    PARSa=AA(ra(k)+1:ra(k)+10,ca(k)+1:ca(k)+3);
    speed=AA(rs(k),cs(k)+1);
    
    %convert to numeric
    kase=str2double(kase);
    temp = cellfun(@str2num,PARSu,'UniformOutput',false); PARSu=cell2mat(temp);
    nu=size(PARSu,1);Nu(kase)=nu; % !!!! correction - Nu array now set up according to the case number
    temp = cellfun(@str2num,PARSa,'UniformOutput',false); PARSa=cell2mat(temp);Na=size(PARSa,1);
    U=str2double(speed)*kph;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CREATE STATE SPACE MATRICES INCLUDING THE EFFECT OF DAMPING
    [A,B,BM,C,D]=articABCD2(PARSu,PARSa,U);
    
    Bj=PARSu(:,5); %joint damping
    KM=zeros(nu,2*nu);
    S=[-1,1;1,-1]; %structure of the individual damping sub-matrix
    for i=1:nu-1
        KM(i:i+1,i+1:i+2)=KM(i:i+1,i+2:i+2)+Bj(i)*S;
    end
    disp(KM);
    disp(BM);
    disp(BM*KM);
    A=A+BM*KM %close the loop via joint damping
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %analysis ...
    [eigvec,eigval]=eig(A);
    eigval=diag(eigval);
    
    %collect results into structure ART, ordering according to kase number
    ART(kase).case=kase;
    ART(kase).U=U;
    ART(kase).PARSu=PARSu;
    ART(kase).PARSa=PARSa;
    ART(kase).eigvec=eigvec;
    ART(kase).eigval=eigval;
    
end
%%%%%%% END LOOP OVER CASES IN THE SPREADSHEET %%%%%%%%%%%%%%%%%%%%%%%

%% tabulate results in terms of frequency and settling time

for k=1:Ncases % 5 percent settling time, frequency in Hz, damping ratio
    val=ART(k).eigval;
    ts=-3*ones(size(val))./real(val);ts=round(100*ts)/100;
    fn=abs(imag(val))/2/pi;fn=round(100*fn)/100;
    zeta=-real(val)./abs(val);zeta=round(1000*zeta)/1000;
%     ts=log(0.05.*sqrt(1-zeta.^2))./real(val);ts=round(100*ts)/100;
    ART(k).ts=ts;
    ART(k).fn= fn;
    ART(k).zeta=zeta;
end

%array I picks out non-negative frequency terms
for k=1:Ncases
    vals=ART(k).eigval;
    I=find(imag(ART(k).eigval)>-tiny);
    ART(k).I=I;
end

%% reconstruction
%%% reconstruct the ART structure, assigning mode number from 1 to ~
Nmodes=length(ART(1).I);

for i=1:Ncases
    current_I=ART(i).I;
    for j=1:length(current_I)
        ART_rec(i).case=ART(i).case;
        ART_rec(i).ts(j)=ART(i).ts(current_I(j));
        ART_rec(i).eigvec(:,j)=ART(i).eigvec(:,current_I(j));
        ART_rec(i).eigval(j)=ART(i).eigval(current_I(j));
        ART_rec(i).I(j)=j;
        ART_rec(i).fn(j)=ART(i).fn(current_I(j));
        ART_rec(i).zeta(j)=ART(i).zeta(current_I(j));
    end
end


%% put into Matlab table format
%%assign missing modes for cases with less units with NaNs
%%%
    ART_table(1).case=ART_rec(1).case;
    ART_table(1).ts=ART_rec(1).ts;
    ART_table(1).fn=ART_rec(1).fn;
    ART_table(1).eigvec=ART_rec(1).eigvec;
    ART_table(1).eigval=ART_rec(1).eigval;
    ART_table(1).zeta=ART_rec(1).zeta;
    ART_table(1).I=ART_rec(1).I;

% 
for l=2:Ncases
    
    dist=zeros(Nmodes,length(ART_rec(l).I));

    for i=1:Nmodes
        val=ART_rec(1).eigval(ART_rec(1).I(i));
        for j=1:length(ART_rec(l).I)
            dist(i,j)=norm(val-ART_rec(l).eigval(ART_rec(l).I(j)));%buid the distance matrix
        end    
    end

    [~,min_i] =min(dist);%%find the smallet distance in each row, and row correponds to latter cases
    diffs=setdiff(ART_rec(1).I,min_i);

    for i=1:length(diffs)
        ART_table(l).ts(diffs(i))=NaN;
        ART_table(l).fn(diffs(i))=NaN;
        ART_table(l).zeta(diffs(i))=NaN;
        ART_table(l).eigval(diffs(i))=NaN;
        [r,c]=size(ART_rec(l).eigvec);
        ART_table(l).eigvec(:,diffs(i))=zeros(r,1);
        ART_table(l).I(diffs(i))=diffs(i);
    end

    for k=1:length(ART_rec(l).I)
        ART_table(l).case=ART_rec(l).case;
        ART_table(l).ts(min_i(k))=ART_rec(l).ts(k);
        ART_table(l).eigvec(:,min_i(k))=ART_rec(l).eigvec(:,k);
        ART_table(l).eigval(min_i(k))=ART_rec(l).eigval(k);
        ART_table(l).I(min_i(k))=min_i(k);
        ART_table(l).fn(min_i(k))=ART_rec(l).fn(k);
        ART_table(l).zeta(min_i(k))=ART_rec(l).zeta(k);      

    end

end

headings=cell(1,Ncases+1);
headings{1}='Mode';
for i=1:Ncases
    headings{i+1}=['Case ',num2str(i)];
end


I1=ART_table(1).I; %use to denote the mode numbers in the table
TS=zeros(Nmodes,Ncases);ZETA=zeros(Nmodes,Ncases);FN=zeros(Nmodes,Ncases);
for k=1:Ncases
    I=ART_table(k).I;
    TS(:,k)=ART_table(k).ts(I);
    ZETA(:,k)=ART_table(k).zeta(I);
    FN(:,k)=ART_table(k).fn(I);
end



TS_table=array2table([I1',TS],'VariableNames',headings);
disp('Settling times (s)'),disp(TS_table)
ZETA_table=array2table([I1',ZETA],'VariableNames',headings);
disp('Damping ratios'),disp(ZETA_table)
FN_table=array2table([I1',FN],'VariableNames',headings);
disp('Natural frequencies (Hz)'),disp(FN_table)


%%
%%%assign missing modes for cases with less units with zeros, cause NaNs
%%%will cut the line


    ART_new(1).case=ART_rec(1).case;
    ART_new(1).ts=ART_rec(1).ts;
    ART_new(1).fn=ART_rec(1).fn;
    ART_new(1).eigvec=ART_rec(1).eigvec;
    ART_new(1).eigval=ART_rec(1).eigval;
    ART_new(1).zeta=ART_rec(1).zeta;
    ART_new(1).I=ART_rec(1).I;
    ART_new(1).PARSu=ART(1).PARSu;

% 
for l=2:Ncases    
    dist=zeros(Nmodes,length(ART_rec(l).I));
    for i=1:Nmodes
        val=ART_rec(1).eigval(ART_rec(1).I(i));
        for j=1:length(ART_rec(l).I)
            dist(i,j)=norm(val-ART_rec(l).eigval(ART_rec(l).I(j)));
        end    
    end

    [~,min_i] =min(dist);
    diffs=setdiff(ART_rec(1).I,min_i);

    for i=1:length(diffs)
        ART_new(l).ts(diffs(i))=0;
        ART_new(l).fn(diffs(i))=0;
        ART_new(l).zeta(diffs(i))=0;
        ART_new(l).eigval(diffs(i))=0;
        [r,c]=size(ART_rec(l).eigvec);
        ART_new(l).eigvec(:,diffs(i))=zeros(r,1);
        ART_new(l).I(diffs(i))=diffs(i);
%            ART_new(l).ts(diffs(i))=NaN;
%             ART_new(l).fn(diffs(i))=NaN;
%             ART_new(l).zeta(diffs(i))=NaN;
%             ART_new(l).eigval(diffs(i))=NaN;
%             [r,c]=size(ART_rec(l).eigvec);
%             ART_new(l).eigvec(:,diffs(i))=zeros(r,1);
%             ART_new(l).I(diffs(i))=diffs(i);
    end

    for k=1:length(ART_rec(l).I)
        ART_new(l).case=ART_rec(l).case;
        ART_new(l).ts(min_i(k))=ART_rec(l).ts(k);
        ART_new(l).eigvec(:,min_i(k))=ART_rec(l).eigvec(:,k);
        ART_new(l).eigval(min_i(k))=ART_rec(l).eigval(k);
        ART_new(l).I(min_i(k))=min_i(k);
        ART_new(l).fn(min_i(k))=ART_rec(l).fn(k);
        ART_new(l).zeta(min_i(k))=ART_rec(l).zeta(k);
        ART_new(l).PARSa=ART(l).PARSa;
        ART_new(l).PARSu=ART(l).PARSu;
    end

end

%%

%% plot the mode shapes using artic_shape_01
% % arrange in order from least stable to most stable
% % plot initial mode shape (green) plus ART shape after T seconds
% T=2;
% 
% Nmodes_i=zeros(Ncases,1);%%no of modes for each case
% 
% for k=1:Ncases
%     I=find(abs(ART_new(k).eigval)>tiny);%%find modes whose eigen value not being zero
%     Nmodes_i(k)=length(I);
% %     figure('NumberTitle', 'off', 'Name', ['Mode Shapes- Case ',num2str(k)],'Position',[0   0   843   794]);
%     figure('NumberTitle', 'off', 'Name', ['Mode Shapes: Case ',num2str(k)],'Position',[680   184   843   794]);
% %     figure('NumberTitle', 'off', 'Name', ['Mode Shapes  Case ',num2str(k)]);
%     for j=1:Nmodes_i(k)
%         m=I(j); %mode number ref. original array
%         % tiny=1e-5;
%         Q=complex(zeros(1,Nu(k)+2));
%         
%         ev=ART_new(k).eigval(m);
%         fn=ART_new(k).fn(m);
%         zeta=ART_new(k).zeta(m);
%         ts=ART_new(k).ts(m);
% 
%         X=ART_new(k).eigvec(:,m)';%complex mode shape as row vector
%         %flip sign if necessary so the max real part is positive
%         [~,p]=max(abs(real(X)));
%         if X(p)<0, X=-X; end
%         
%         %convert eigenvector X to standard q vector (Q = complex amplitude)
%         PARSu=ART_new(k).PARSu;
%         Xr=-PARSu(:,4)';
%         Xf=PARSu(:,3)'+Xr;
%         a1=Xf(1); %distance from CG1 to front joint
%         
%         Q(2)=(U*X(1)+a1*X(2))/ev; %integration of lateral velocity at joint
%         Q(3)=X(2)/ev; %integration of front yaw rate.
%         for i=2:Nu(k) %use articulation angles to find other yaw angles
%             Q(2+i)=Q(1+i)+X(Nu(k)+i); %e.g. i=2: Q(4)=Q(3)+X(5)
%         end
%         Q1=real(Q);
%         Q2=real(Q*exp(ev*T));
%         % plot results
%         subplot(ceil(Nmodes_i(k)/2),2,j)
%         artic_shape_01(Q2,PARSu,'k')
%         axis('equal'),hold on
%         artic_shape_01(Q1,PARSu,'g')
%         
%         yl =  ylim;%get current axis property and set text position
%         xl =  xlim;
%         text(mean(xl)-20,mean(yl)+25,['mode: ',num2str(m),'   f= ',num2str(fn),' Hz' ])
%         text(mean(xl)-20,mean(yl)+18,['\zeta = ',num2str(zeta),'      t_s = ',num2str(ts) ])
%         axis([xl(1)-5,xl(2)+5,yl(1)-5,yl(2)+15])
%     end
% end
% 
% 
% 
% 
% 
% 
% %%
% figure('NumberTitle', 'off', 'Name', 'Trends for key stability parameters')
% tls={'settling time (5%)','damping ratio (\zeta)','natural frequency (Hz)'};
% % set up for plotting
% xplot=(1:Ncases)';
% 
% for p=1:3 %plot number
%     subplot(3,1,p),title(tls{p})
%     hold on
% end
% subplot(3,1,3),xlabel('Case number')
% 
% %fix the number of modes based on case 1
% % n=sum(imag(ART_new(1).eigval)>-tiny); %number of modes in case 1
% 
% 
% 
% TS=zeros(Ncases,Nmodes);
% ZETA=zeros(Ncases,Nmodes);
% FN=zeros(Ncases,Nmodes);
% modelabels=cell(1,Nmodes);
% for k=1:Ncases
%     for j=1:Nmodes
%         TS(k,j)=ART_new(k).ts(j);
%         ZETA(k,j)=ART_new(k).zeta(j);
%         FN(k,j)=ART_new(k).fn(j);
%         modelabels{j}=['Mode ', num2str(j)];
%     end
% end
% 
% subplot(3,1,1),set(gca,'XTickMode','manual','Xtick',(1:Ncases))
% plot(xplot,TS)
% subplot(3,1,2),set(gca,'XTickMode','manual','Xtick',(1:Ncases))
% plot(xplot,ZETA)
% subplot(3,1,3),set(gca,'XTickMode','manual','Xtick',(1:Ncases))
% plot(xplot,FN)
% 
% for p=1:3 %plot number
%     subplot(3,1,p)
%     ax = gca;
%     ax.ColorOrder = [0 0 0;1 0 0;0 0 1;0 1 0];
%     ax.LineStyleOrder = {'-','--',':'};
% end
% legend(modelabels)

