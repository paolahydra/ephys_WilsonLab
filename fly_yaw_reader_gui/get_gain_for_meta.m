function [ gain ] = get_gain_for_meta( gain_meta )

if gain_meta > 0 && gain_meta < 2.34
    gain = 0.5;
elseif gain_meta >= 2.34 && gain_meta < 2.85
    gain = 1;
elseif gain_meta >= 2.85 && gain_meta < 3.34
    gain = 2;
elseif gain_meta >= 3.34 && gain_meta < 3.85
    gain = 5;
elseif gain_meta >= 3.85 && gain_meta < 4.37
    gain = 10;
elseif gain_meta >= 4.37 && gain_meta < 4.85
    gain = 20;
elseif gain_meta >= 4.85 && gain_meta < 5.34
    gain = 50;
elseif gain_meta >= 5.34 && gain_meta < 5.85
    gain = 100;
elseif gain_meta >= 5.85 && gain_meta < 6.37
    gain = 200;
elseif gain_meta >= 6.37 && gain_meta < 6.85
    gain = 500;
end

end

