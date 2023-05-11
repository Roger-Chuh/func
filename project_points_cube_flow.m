function [xp,dxpdom,dxpdT,dxpdf,dxpdc,dxpdk,dxpdalpha] = project_points_cube_flow(iImg, cbGrid, r,t,foc,cen,k,alpha,cubeIdx,config)

gridHomo1 = [cbGrid{iImg}{1,2}; zeros(1, length(cubeIdx{1}))];
gridHomo2 = [cbGrid{iImg}{2,2}; zeros(1, length(cubeIdx{2}))];
gridHomo3 = [cbGrid{iImg}{3,2}; zeros(1, length(cubeIdx{3}))];
T1 = [rodrigues(r) t; 0 0 0 1];
T2 = T1*config.rt2;
T3 = T1*config.rt3;
[xp1,dxp1dom,dxp1dT,dxp1df,dxp1dc,dxp1dk,dxp1dalpha] = project_points2(gridHomo1,r,t,foc,cen,k,alpha);



%  [rotVecRight,transVecRight,drrdrl,drrdtl,drrdrref,drrdtref,dtrdrl,dtrdtl,dtrdrref,dtrdtref] = compose_motion(rodrigues(config.rt2(1:3,1:3)),config.rt2(1:3,4),r,t);
%  [rotVecRight,transVecRight,drrdrl, drrdtl, drrdrref,drrdtref,dtrdrl, dtrdtl, dtrdrref,dtrdtref] = compose_motion(rotVecLeft,                    transVecLeft,     rotVecRef,   transVecRef );
[rotVec2,    transVec2,    dr2dr12,dr2dt12,dr2dr1,  dr2dt1,  dt2dr12,dt2dt12,dt2dr1,  dt2dt1]   = compose_motion(rodrigues(config.rt2(1:3,1:3)),config.rt2(1:3,4),r, t);
[xp2,dxp2dom,dxp2dT,dxp2df,dxp2dc,dxp2dk,dxp2dalpha] = project_points2(gridHomo2, rotVec2, transVec2, foc, cen, k, alpha);

dxp2dr1 = dxp2dom * dr2dr1 + dxp2dT * dt2dr1;
dxp2dt1 = dxp2dom * dr2dt1 + dxp2dT * dt2dt1;
if 1
    dxp2dr12 = dxp2dom * dr2dr12 + dxp2dT * dt2dr12;
    dxp2dt12 = dxp2dom * dr2dt12 + dxp2dT * dt2dt12;
end


[rotVec3,    transVec3,    dr3dr13,dr3dt13,dr3dr1,  dr3dt1,  dt3dr13,dt3dt13,dt3dr1,  dt3dt1]   = compose_motion(rodrigues(config.rt3(1:3,1:3)),config.rt3(1:3,4),r, t);
[xp3,dxp3dom,dxp3dT,dxp3df,dxp3dc,dxp3dk,dxp3dalpha] = project_points2(gridHomo3, rotVec3, transVec3, foc, cen, k, alpha);

dxp3dr1 = dxp3dom * dr3dr1 + dxp3dT * dt3dr1;
dxp3dt1 = dxp3dom * dr3dt1 + dxp3dT * dt3dt1;
if 1
    dxp3dr13 = dxp3dom * dr3dr13 + dxp3dT * dt3dr13;
    dxp3dt13 = dxp3dom * dr3dt13 + dxp3dT * dt3dt13;
end




xp = [xp1 xp2 xp3];
dxpdom = [dxp1dom; dxp2dr1; dxp3dr1];
dxpdT = [dxp1dT; dxp2dt1; dxp3dt1];
dxpdf = [dxp1df; dxp2df; dxp3df];
dxpdc = [dxp1dc; dxp2dc; dxp3dc];
dxpdk = [dxp1dk; dxp2dk; dxp3dk];
dxpdalpha = [dxp1dalpha; dxp2dalpha; dxp3dalpha];

if 0
    figure(1),subplot(2,2,3);plot(xp1(1,:),xp1(2,:),'or');plot(xp2(1,:),xp2(2,:),'og');plot(xp3(1,:),xp3(2,:),'ob');
end

end