function track=trackgen4(i,R)


if nargin ==1
    R=50;
end

tiny=1e-6;
p0=[0 0];t0=[1 0]; %default start position and direction
switch i %set track definition data

    case 1 %closed big track

        %length of straight, radius of curve
        L=[700 100 100 200 50 50 450 50 150 100];
        %swept angles for curved segments (+ = left);
        Q=[0 pi/2 pi/2 0 -pi/2 pi/2 0 pi/2 0 pi/2]; 
        track=build5(L,Q,p0,t0);


    case 2 %closed circle (define as 4 segments)

        L2=[1 1 1 1 ]*R;
        Q2=[pi/2 pi/2 pi/2 pi/2];
        L=[250 L2 L2 100];
        Q=[0 Q2 Q2 0];

        track=build5(L,Q,p0,t0);

        
    case 3 %open straight

        L=5000;
        Q=0;
        track=build5(L,Q,p0,t0);
        
    case 4 %open wiggle
        L=[250 100 300 300 200 600 200 100 50]*R/50;
        q1=3*20*pi/180;q2=3*35*pi/180;
        
        Q=[0 q2 -0.2*q2 -0.8*q2 -q1 q1 -q2 q2 0];
        track=build5(L,Q,p0,t0);
        
    case 5 %Hockenheim (ISO)
        R=[0 68 0 45 150 150 0 160 160 35 0 95 0 40 0 90 100 0 50 0 69.5 0]; %radii
        Q=-[0 69 0 102 35 -35 0 -47 47 89 0 84 0 -155 0 -52 27 0 77 0 119 0]*pi/180; %turn angles
        L=[200 0 220 0 0 0 200 0 0 0 60 0 220 0 60 0 0 45 0 45 0 297]; %lengths of straights
        L=L+R; %all lengths in one array
        
        track=build5c(L,Q,p0,t0); %forces closure
                  
        
    case 7 %simple repeating 90 deg arcs
        K=2; %number of arc-groups

        Larcs=ones(1,4*K)*R;
        Qarcs=kron(ones(1,K),[1 -1 -1 1]*pi/2);
        L= [250 ,Larcs, 100];
        Q=[0  Qarcs 0];
        
        track=build5(L,Q,p0,t0);
		
		
    case 8 %figure 8
        L0=[R R R R R R R R ];
        Q0=[pi/2 -pi/2 -pi/2 -pi/2 -pi/2 pi/2 pi/2 pi/2 ];
        Q1=[pi/2 -pi/2];
        L1=[R R];
        L=[250 L0 L0  L1 250];
        Q=[0 Q0 Q0  Q1 0];
        track=build5(L,Q,p0,t0);
        
    case 9 %lane change
        W=3; %metres lateral, R is curve radius
        q=acos(1-W/(2*R));
        L=[250,R,R,200];
        Q=[0,q,-q,0];
        track=build5(L,Q,p0,t0);

                    
end

end






