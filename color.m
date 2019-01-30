function [co] = color()
% COLOR Custom color maps for plotting
% http://jdherman.github.io/colormap/
% RGB color codes
    % 8 for alkanes increasing with C
    % 7 for alkenes increasing with C 
    % 4 for BTEX compounds 
    % 5 for additive compounds 
    % 1 for oxygen 
co = [ % alkene - teals 
%         20 166 164
%         18 153 151
%         17 140 139
%         15 128 126 
%         14 117 116 
%         12 99 98
%         9 79 78
%         7 59 58 
        132,213,241;
        26,188,236;
        0,170,229;
        0,155,222;
        1,143,213;
        3,135,202;
        3,121,188;
        2,83,167;
            
        % alkene - purples 
%         230 105 239
%         206 94 154
%         179 81 186
%         154 70 160
%         132 60 137
%         105 48 109
%         81 37 47 
        220,121,234;
        202,4,228;
        183,2,211;
        163,0,194;
        134,0,168;
        106,0,143;
        89,0,121;    
        % BTEX - bright color - summer (4)
%         255,255,7;
%         93,198,112;
%         46,129,168;
%         0,0,179;
         
        232,248,0;
        101,235,138;
        41,198,128;
        30,149,53
        % Additives - brights - jet(6)
%         34 117 165
%         247 92 2
%         242 196 14
%         216 2 102
%         0 216 108
      	255,189,0;
        255,34,237;
        175,0,248;
        0,188,234;
        0,250,217

        % Oxygen - grey ?
        216 3 0];
    
co = co./255;    
    
end

