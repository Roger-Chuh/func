function PcShow(xyz)


        pcshow(xyz, [0 0 1], 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
        xlabel('X (mm)');
        ylabel('Y (mm)');
        zlabel('Z (mm)');
        set(gca,  'CameraUpVector',[0 -1 0],'DataAspectRatio',[1 1 1]);
        title('Reconstructed 3-D Scene');



end