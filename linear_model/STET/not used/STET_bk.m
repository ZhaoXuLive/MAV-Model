% STability Evaluation Tool (STET.m)
% Calls articABCD2 to build state space (linear) model of ART
% (includes joint damping parameters)
% Reads data from excel file and restructures it for analysis
clear
close all
g=9.81;deg=pi/180;kph=1/3.6;tiny=1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%READ DATA FROM FILE
infile='STET input data.xlsx';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A0=readtable(infile);AA=table2cell(A0);
%find locations of UNIT etc
temp = find(strcmp( AA, 'CASE' ));[rc,cc] = ind2sub(size(AA), temp ) ;
temp = find(strcmp( AA, 'SPEED' ));[rs,cs] = ind2sub(size(AA), temp ) ;
temp = find(strcmp( AA, 'UNIT' ));[ru,cu] = ind2sub(size(AA), temp ) ;
temp = find(strcmp( AA, 'AXLE' ));[ra,ca] = ind2sub(size(AA), temp ) ;
Ncases=length(ru); %how many boxes in the spreadsheet


%%%%%%% START LOOP OVER CASES IN THE SPREADSHEET %%%%%%%%%%%%%%%%%%%%%%%
Nu=zeros(Ncases,1);
for k=1:Ncases
    kase=AA(rc(k),cc(k)+1);
    PARSu=AA(ru(k)+1:ru(k)+5,cu(k)+1:cu(k)+5);
    PARSa=AA(ra(k)+1:ra(k)+10,ca(k)+1:ca(k)+3);
    speed=AA(rs(k),cs(k)+1);
    
    %convert to numeric
    kase=str2double(kase);
    temp = cellfun(@str2num,PARSu,'UniformOutput',false); PARSu=cell2mat(temp);
    nu=size(PARSu,1);Nu(kase)=nu;
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
    A=A+BM*KM; %close the loop via joint damping
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


%% put into Matlab table format
headings=cell(1,Ncases+1);
headings{1}='Mode';
for i=1:Ncases
    headings{i+1}=['Case ',num2str(i)];
end


I1=ART(1).I; %use to denote the mode numbers in the table
n=length(I1);n2=ceil(n/2); %for plotting later
TS=zeros(n,Ncases);ZETA=zeros(n,Ncases);FN=zeros(n,Ncases);
for k=1:Ncases
    I=ART(k).I;
    TS(:,k)=ART(k).ts(I);
    ZETA(:,k)=ART(k).zeta(I);
    FN(:,k)=ART(k).fn(I);
end



TS_table=array2table([I1,TS],'VariableNames',headings);
disp('Settling times (s)'),disp(TS_table)
ZETA_table=array2table([I1,ZETA],'VariableNames',headings);
disp('Damping ratios'),disp(ZETA_table)
FN_table=array2table([I1,FN],'VariableNames',headings);
disp('Natural frequencies (Hz)'),disp(FN_table)


%% plot the mode shapes using artic_shape_01
% arrange in order from least stable to most stable
% plot initial mode shape (green) plus ART shape after T seconds
T=2;


for k=1:Ncases
    
    I=ART(k).I;

    %figure('NumberTitle', 'off', 'Name', ['Mode Shapes: Case ',num2str(k)],'Position',[680   184   843   794]);
    figure('NumberTitle', 'off', 'Name', ['Mode Shapes: Case ',num2str(k)],'Position',[680   184   843   794]);

    for j=1:n
        m=I(j); %mode number ref. original array
        % tiny=1e-5;
        Q=complex(zeros(1,Nu+2));
        
        ev=ART(k).eigval(m);
        fn=ART(k).fn(m);
        zeta=ART(k).zeta(m);
        ts=ART(k).ts(m);

        X=ART(k).eigvec(:,m)';%complex mode shape as row vector
        %flip sign if necessary so the max real part is positive
        [~,p]=max(abs(real(X)));
        if X(p)<0, X=-X; end
        
        %convert eigenvector X to standard q vector (Q = complex amplitude)
        PARSu=ART(k).PARSu;
        Xr=-PARSu(:,4)';
        Xf=PARSu(:,3)'+Xr;
        a1=Xf(1); %distance from CG1 to front joint
        
        Q(2)=(U*X(1)+a1*X(2))/ev; %integration of lateral velocity at joint
        Q(3)=X(2)/ev; %integration of front yaw rate.
        for i=2:Nu %use articulation angles to find other yaw angles
            Q(2+i)=Q(1+i)+X(Nu+i); %e.g. i=2: Q(4)=Q(3)+X(5)
        end
        Q1=real(Q);
        Q2=real(Q*exp(ev*T));
        
        % plot results
        subplot(n2,2,j)
        artic_shape_01(Q2,PARSu,'k')
        axis('equal'),hold on
        artic_shape_01(Q1,PARSu,'g')
        text(-28,15,['mode: ',num2str(m),' : ',num2str(fn),' Hz' ])
        text(-28,10,['\zeta = ',num2str(zeta),' :  t_s = ',num2str(ts) ])
        
    end
end

%% PLOT TRENDS FOR FN, TS, ZETA

figure('NumberTitle', 'off', 'Name', 'Trends for key stability parameters')
tls={'settling time (5%)','damping ratio (\zeta)','natural frequency (Hz)'};
% set up for plotting
xplot=(1:Ncases)';

for p=1:3 %plot number
    subplot(3,1,p),title(tls{p})
    hold on
end
subplot(3,1,3),xlabel('Case number')


%fix the number of modes based on case 1
n=sum(imag(ART(1).eigval)>-tiny); %number of modes in case 1
TS=zeros(Ncases,n);
ZETA=zeros(Ncases,n);
FN=zeros(Ncases,n);
modelabels=cell(1,n);
for k=1:Ncases
    for j=1:n
        m=I(j); %mode number ref. original array
        TS(k,j)=ART(k).ts(m);
        ZETA(k,j)=ART(k).zeta(m);
        FN(k,j)=ART(k).fn(m);
        modelabels{j}=['Mode ', num2str(m)];
    end
end

subplot(3,1,1),set(gca,'XTickMode','manual','Xtick',(1:Ncases))
plot(xplot,TS)
subplot(3,1,2),set(gca,'XTickMode','manual','Xtick',(1:Ncases))
plot(xplot,ZETA)
subplot(3,1,3),set(gca,'XTickMode','manual','Xtick',(1:Ncases))
plot(xplot,FN)

for p=1:3 %plot number
    subplot(3,1,p)
    ax = gca;
    ax.ColorOrder = [0 0 0;1 0 0;0 0 1;0 1 0];
    ax.LineStyleOrder = {'-','--',':'};
end
legend(modelabels)


