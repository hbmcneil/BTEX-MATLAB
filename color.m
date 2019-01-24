function [co] = color()
% COLOR Custom color maps for plotting
% RGB color codes
    % 8 for alkanes increasing with C
    % 7 for alkenes increasing with C 
    % 4 for BTEX compounds 
    % 5 for additive compounds 
    % 1 for oxygen 
co = [ % alkene - teals 
        20 166 164
        18 153 151
        17 140 139
        15 128 126 
        14 117 116 
        12 99 98
        9 79 78
        7 59 58 
        % alkene - purples 
        230 105 239
        206 94 114
        179 81 186
        154 70 160
        132 60 137
        105 48 109
        81 37 47 
        % BTEX - bright color
        178 76 99 
        255 140 66 
        83 255 69 
        30 46 222
        % Additives - brights?
        34 117 165
        247 92 2
        242 196 14
        216 2 102
        0 216 108
        % Oxygen - grey ?
        216 3 0];
    
co = co./255;    
    
end

