function movefiles(flag)

if flag == 1
    system('D:\ForDump\apriltag\cmake-build-debug2\camtest2.exe  -W 960 -H 540 -w');
else
    system('D:\ForDump\apriltag\cmake-build-debug2\camtest2.exe  -W 960 -H 540 -i 3');
    pause(2);
%     movefile('D:\ForDump\apriltag\cmake-build-debug2\test_002.avi', 'D:\ForDump\apriltag\cmake-build-debug2\test_022.avi');
    copyfile('D:\ForDump\apriltag\cmake-build-debug2\test_002.avi','\\192.168.50.172\nextvpu\3.CTO办公室\Roger\temp\log');
    copyfile('D:\ForDump\apriltag\cmake-build-debug2\test_002.log','\\192.168.50.172\nextvpu\3.CTO办公室\Roger\temp\log');
%     movefile('D:\ForDump\apriltag\cmake-build-debug2\test_022.avi', 'D:\ForDump\apriltag\cmake-build-debug2\test_002.avi');
%     copyfile('D:\ForDump\apriltag\cmake-build-debug2\test_002.avi','\\192.168.50.172\nextvpu\3.CTO办公室\Roger\temp\log');
end

end