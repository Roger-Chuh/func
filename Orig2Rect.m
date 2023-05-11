function pixRect = Orig2Rect(pix, KOrig, KRect, R,kc)

[pixUndist] = normalize_pixel(pix',[KOrig(1,1);KOrig(2,2)],[KOrig(1,3);KOrig(2,3)],kc,0);
pixUndistR = R*pextend(pixUndist);
pixRect = pflat(KRect*pixUndistR);
pixRect = pixRect(1:2,:)';



end