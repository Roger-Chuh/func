function rMatNew = CorrectRot(axis, gVec)
ang = CalcDegree(axis',gVec');
xVec = cross(axis,gVec);
xVec = xVec./norm(xVec);
zVec = cross(xVec, gVec);
zVec = zVec./norm(zVec);

rMatNew = [xVec'; gVec'; zVec']';
end