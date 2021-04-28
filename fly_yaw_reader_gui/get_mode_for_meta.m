function [ mode ] = get_mode_for_meta( mode_meta )

if mode_meta > 0 && mode_meta < 3.5
    mode = 'V';
elseif mode_meta >= 3.5 && mode_meta < 7
    mode = 'I';
end
        
end

