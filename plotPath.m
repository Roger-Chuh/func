function plotPath(data, clr)
hold on;
scale = 0.5;
if 0
    
    % %     pcshow(UU(:,UU(1,:)<XX & UU(3,:)<ZZ)', [1 0 0], 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
    % %     legend('direction','proposed method','sfm reference');
    % %     xlabel('X (mm)');
    % %     ylabel('Y (mm)');
    % %     zlabel('Z (mm)');
    % %     set(gca,  'CameraUpVector',[0 -1 0],'DataAspectRatio',[1 1 1]);
    % %     title('Reconstructed 3-D Scene');
    quiver3(data(:,10), data(:,11), data(:,12), scale.*data(:,1), scale.*data(:,2), scale.*data(:,3),'Color', [1 0 0]);
    quiver3(data(:,10), data(:,11), data(:,12), scale.*data(:,4), scale.*data(:,5), scale.*data(:,6),'Color', [0 1 0]);
    quiver3(data(:,10), data(:,11), data(:,12), scale.*data(:,7), scale.*data(:,8), scale.*data(:,9),'Color', [0 0 1]);
    legend('x','y','z');
    xlabel('x(mm)');ylabel('y(mm)');zlabel('z(mm)')
    % hold on;
    % for i = 1 : size(data,1)
    %     plot3(data(i,1), data(i,2), data(i,3), 'o', 'Color', color, 'LineWidth', lineWidth);
    % %     text(data(i,1), data(i,2), data(i,3), num2str(i));
    % end
    axis equal;
else
    if size(data,2) == 12
        
        pcshow([data(1,10), data(1,11), data(1,12)], [0 0 1], 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
        xlabel('X (mm)');
        ylabel('Y (mm)');
        zlabel('Z (mm)');
        if 0
            plot3(data(:,10), data(:,11),data(:,12),'-r');
        end
        set(gca,  'CameraUpVector',[0 -1 0],'DataAspectRatio',[1 1 1]);
        title('Reconstructed 3-D Scene');
        hold on
    else
%         if flag
            pcshow([data(:,1), data(:,2), data(:,3)], clr, 'VerticalAxis', 'Y', 'VerticalAxisDir', 'Down');
%         end
        plot3(data(1,1), data(1,2), data(1,3),'sb' ,'MarkerSize',8,'LineWidth',8);
        plot3(data(:,1), data(:,2), data(:,3), '-o', 'Color',clr);
        xlabel('X (mm)');
        ylabel('Y (mm)');
        zlabel('Z (mm)');
%         if flag
            set(gca,  'CameraUpVector',[0 -1 0],'DataAspectRatio',[1 1 1]);
%         end
        title('Reconstructed 3-D Scene');
        hold on
        return;
    end
    % %    scatter3(XYZtemp(:,1),XYZtemp(:,2),XYZtemp(:,3),'filled'); axis equal;
    
    
    %     quiver3(poseVec11(1:kl,1),poseVec11(1:kl,2),poseVec11(1:kl,3),poseVec11(1:kl,4),poseVec11(1:kl,5),poseVec11(1:kl,6),0.1);
    %     quiver3(poseVec11(1:kl,1),poseVec11(1:kl,2),poseVec11(1:kl,3),poseVec11(1:kl,7),poseVec11(1:kl,8),poseVec11(1:kl,9),0.1);
    %     quiver3(poseVec11(1:kl,1),poseVec11(1:kl,2),poseVec11(1:kl,3),poseVec11(1:kl,10),poseVec11(1:kl,11),poseVec11(1:kl,12),0.1);
    quiver3(data(:,10), data(:,11), data(:,12), scale.*data(:,1), scale.*data(:,2), scale.*data(:,3),scale);
    quiver3(data(:,10), data(:,11), data(:,12), scale.*data(:,4), scale.*data(:,5), scale.*data(:,6),scale);
    quiver3(data(:,10), data(:,11), data(:,12), scale.*data(:,7), scale.*data(:,8), scale.*data(:,9),scale);
    % %    quiver3(poseVec11Robot(1:kl,1),poseVec11Robot(1:kl,2),poseVec11Robot(1:kl,3),poseVec11Robot(1:kl,4),poseVec11Robot(1:kl,5),poseVec11Robot(1:kl,6),0.1);
    % %    quiver3(poseVec11Robot(1:kl,1),poseVec11Robot(1:kl,2),poseVec11Robot(1:kl,3),poseVec11Robot(1:kl,7),poseVec11Robot(1:kl,8),poseVec11Robot(1:kl,9),0.1);
    % %    quiver3(poseVec11Robot(1:kl,1),poseVec11Robot(1:kl,2),poseVec11Robot(1:kl,3),poseVec11Robot(1:kl,10),poseVec11Robot(1:kl,11),poseVec11Robot(1:kl,12),0.1);
    
    legend('world origin','pnp x','pnp y', 'pnp z');
    plot3(data(1,10), data(1,11), data(1,12),'sb','MarkerSize',15,'LineWidth',15);
end

end