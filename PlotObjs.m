function PlotObjs(pop)


Objs=[pop.Objectives];

    plot3(Objs(1,:),-1*Objs(2,:), Objs(3,:),'s','MarkerSize',6, 'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]);
    
    xlabel('1st Objective');
    ylabel('2nd Objective');
    zlabel('3rd Objective');
    
    grid on;
    

end