function track=pts2trk_v2(points)

    points = downsample(points,100);
    ptsNum = length(points);
    track = zeros(ptsNum,8);
    
    delta = diff(points);
    delta = [delta;delta(ptsNum-1,:)];
    N_delta = zeros(ptsNum,1);
    for i=1:ptsNum
        N_delta(i)=norm(delta(i,:));
    end
    track(:,1) = [0;cumsum(N_delta(1:ptsNum-1))];
    track(:,2:3) = points;
    track(:,4) = delta(:,1) ./ N_delta;
    track(:,5) = delta(:,2) ./ N_delta;
    track(:,6:7) = [-track(:,5),track(:,4)];
    for i=1:ptsNum-1
        ang = atan2(track(i+1,5),track(i+1,4))-atan2(track(i,5),track(i,4));
        if ang > pi
            ang = ang-2*pi;
        elseif ang < -pi
            ang = ang+2*pi;
        end
        track(i,8)=ang/N_delta(i);
    end
end