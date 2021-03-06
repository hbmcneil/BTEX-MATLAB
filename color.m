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
        132,213,241;
        26,188,236;
        0,170,229;
        0,155,222;
        1,143,213;
        3,135,202;
        3,121,188;
        2,83,167;
            
        % alkene - purples 
        220,121,234;
        202,4,228;
        183,2,211;
        163,0,194;
        134,0,168;
        106,0,143;
        89,0,121;    
        
        % BTEX - greens
        232,248,0;
        101,235,138;
        41,198,128;
        30,149,53
        
        % Additives - brights 
      	255,189,0;
        255,34,237;
        175,0,248;
        0,188,234;
        0,250,217
        
        % Oxygen - red
        216 3 0];
    
co = co./255;    
    
end

