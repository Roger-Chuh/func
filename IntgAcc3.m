function [vel1, rot1, t1, angle,ang,dis1,aWorld] = IntgAcc3(a0,a1,w0,w1,vel0,dis0,dt,RTworld)
% % vel(1)=vel0;
% % dis(1)=dis0;
% %  T = [rot1 t1;0 0 0 1]; pt3d1 = T*pt3d0

% vel0: world velocity
% RTworld£ºlocal to world rotation 






if size(vel0,2) == 3;
    vel0 = vel0';
end
% % % for j = 1:3
% % % 
% % % % % for i=2:length(a)
% % % %  numerically integrate acceleration to find velocity
% % %     vel1(j)=vel0(j)+1*(a0(j))*dt;
% % % %  numerically integrate angular rate to find roll pitch yaw( rotation around the x, y, z axis)
% % % %     RollPitchYaw(j)=w0(j)+0.5*(w0(j)+w1(j))*dt;
% % % RollPitchYaw(j)= 1*(w0(j))*dt;
% % % %  numerically integrate velocity to find position
% % % %     dis1(j)=dis0(j)+0.5*(vel0(j)+vel1(j))*dt;
% % % dis1(j)=1*(vel0(j))*dt;
% % % % % end
% % % 
% % % 
% % % 
% % % end





% RTworld : rotation from local-1 to world

RollPitchYaw = w0*dt; %RollPitchYaw = RollPitchYaw';  
rot1 = rodrigues(RollPitchYaw); % rotation from local to local - 1 ;

vel1 = vel0+ RTworld*a0.*dt; vel1 = vel1';     % local velocity to world velocity


dis1 = dis0 + vel0*dt + RTworld*0.5*a0.*dt*dt; %dis1 = dis1';
% % % % % % dis1=dis0 + 0.5*(vel0+vel1')*dt;

aWorld = RTworld*a0;
% % % rot1 = rotz(-180*RollPitchYaw(3)/pi)*roty(-180*RollPitchYaw(2)/pi)*rotx(-180*RollPitchYaw(1)/pi);

angle = [180*RollPitchYaw(1)/pi 180*RollPitchYaw(2)/pi 180*RollPitchYaw(3)/pi];
t1 = 1;   %rot1*(-dis1');
% % % % t1 = RTworld'*(dis1);


q = v2q(rodrigues(RTworld*rot1));
q2 = v2q(RollPitchYaw);
vec = [dis1;q;vel1']; % world frame
vec2 = [0.5*a0.*dt*dt;q2;a0.*dt]; % local frame

vec3 = [vec2 vec];

% T = [rot1 t1; 0 0 0 1] transforms pt3d in t0 cs to t1 cs;
% pt2d2 = T*pt3d1;

ang = rad2deg(norm(RollPitchYaw));


% % % % % rotz(45)*roty(45)
% % % % % ans*[sqrt(2);0;0]
% % % % % roty(45)*rotz(45)
% % % % % ans*[sqrt(2);0;0]
% % % % % rotz(45)*roty(-45)
% % % % % ans*[sqrt(2);0;0]




end