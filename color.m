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
        % alkene - dusky purples 
        218 214 214
        189 187 195 
        175 159 165
        166 144 164 
        150 127 153 
        120 81 122
        95 68 108 
        % BTEX - bright color
        178 76 99 
        255 140 66 
        83 255 69 
        30 46 222
        % Additives - reds?
        157 68 181
        181 68 110
        84 8 4 
        70 99 102
        171 52 40
        % Oxygen - grey ?
        57 57 58];
    
co = co./255;    
    
end

