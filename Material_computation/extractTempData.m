    %Usage of 4340 AISI data in MATLAB Figure
    fig = open('AnalyseThermoMagnetique_4340ferrite.fig');
    axObjs = fig.Children;
    
    %Wh
    dataObjs =  axObjs(1).Children;
    Temp1 = dataObjs(2).XData;
    Whdata = dataObjs(2).YData;
    
    %Hc
    dataObjs =  axObjs(2).Children;
    Temp2 = dataObjs(2).XData;
    Hcdata = dataObjs(2).YData;
    
    %Br
    dataObjs =  axObjs(3).Children;
    Temp3 = dataObjs(2).XData;
    Brdata = dataObjs(2).YData;
    
    %Bsat
    dataObjs =  axObjs(4).Children;
    Temp4 = dataObjs(2).XData;
    Bsatdata = dataObjs(2).YData;
    
    %murmax
    dataObjs =  axObjs(5).Children;
    Temp5 = dataObjs(2).XData;
    murmaxdata = dataObjs(2).YData;
    
    minlen = min([length(Temp1),length(Temp2),length(Temp3),length(Temp4),length(Temp5)]);

    
    Br_list = Brdata(1:minlen);    %T
    Bsat_list = Bsatdata(1:minlen);  %T
    Hc_list = Hcdata(1:minlen); %A/m
    Wh_list = Whdata(1:minlen);
    mur_max_list = murmaxdata(1:minlen);
    a_list = zeros(size(param_phy.Br_list));
    s_list = zeros(size(param_phy.Br_list));
    b_list = zeros(size(param_phy.Br_list));